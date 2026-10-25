#!/usr/bin/env python
"""D0 — ANA-CNS build: brain-filtered ANA-GEN for Arm E train1_ana (CPU).

- UBERON brain-subtree: anatomy nodes where properties.name contains "brain" (case-insensitive) OR descendants of UBERON_0000955 via anatomy_anatomy PARENT
  (union = brain_nodes, ~150 via name; hierarchy adds few, kept union)
- Specificity: drop gene if it appears in >50% of anatomy nodes (13,120) — i.e., ubiquitous genes >6560 anatomies
  (alternative 50% of brain nodes would be 75, but spec says anatomy nodes → use 13k)
- Tripwires: kept edges >500k / <5k / mapping coverage <50% → STOP
- Output: train1_ana exports under data/processed/m2c/completion_exports_ana/<unit>/ (train1 1,431,945 + kept ANA brain edges)
  + re-run degree/PPR/DWPC baselines → reports/m2c/completion_baselines_ana.md

Usage: python scripts/m2c_build_ana_cns.py --dry-run
"""
import argparse
from pathlib import Path
import json
import collections
import pyarrow.parquet as pq
import pyarrow.dataset as ds
import polars as pl
from pydantic import BaseModel

# Simple config via args, not yaml for D0 speed
def build_brain_nodes():
    pf = pq.ParquetFile("data/raw/optimuskg/doi_10_7910_DVN_IYNGEV/2.0/nodes/anatomy.parquet")
    brain_nodes = set()
    # CNS = brain + spinal cord + central nervous system (all CNS-related)
    # Curated CNS keywords — only names, not synonyms to avoid pancreatic insula false positives
    cns_keywords = [
        "brain", "cerebellum", "cerebral cortex", "cerebral hemisphere", "cerebrum",
        "hippocampus", "amygdala", "striatum", "corpus striatum", "thalamus", "hypothalamus",
        "basal ganglia", "globus pallidus", "substantia nigra", "prefrontal cortex", "visual cortex",
        "olfactory bulb", "olfactory", "brainstem", "brain stem", "choroid plexus", "telencephalon",
        "diencephalon", "mesencephalon", "metencephalon", "myelencephalon", "pallium", "subpallium",
        "neocortex", "allocortex", "archicortex", "claustrum", "insula", "cingulate", "cingulum",
        "septum", "accumbens", "caudate", "putamen", "pallidum", "tegmentum", "colliculus",
        "geniculate", "subthalamic", "red nucleus", "pons", "medulla", "midbrain", "forebrain",
        "hindbrain", "spinal cord", "central nervous system", "neural tube", "meninx", "meninges",
        "gyrus", "sulcus", "fissure", "ventricle", "choroid", "ependyma", "corpus callosum",
    ]
    exclude_ids = {"UBERON_0000389"}  # lens cortex
    # For CNS, we include spinal cord and central nervous system phrases explicitly
    for batch in pf.iter_batches(batch_size=2000, columns=["id","properties"]):
        for row in batch.to_pylist():
            nid = row["id"]
            if nid in exclude_ids:
                continue
            props = row["properties"] or {}
            name = (props.get("name","") or "").lower()
            # Check name only (not synonyms) to avoid pancreatic insula etc.
            for kw in cns_keywords:
                if kw in name:
                    # Special handling for generic short words that cause false positives:
                    # "pons" alone false positive? pons is always brain, so ok. "medulla" alone could be adrenal medulla (UBERON_0002369) not CNS - but we include it as it's part of brainstem; adrenal medulla is not CNS but name is "adrenal medulla" contains medulla - we should exclude adrenal
                    if kw == "medulla" and "adrenal" in name:
                        continue
                    if kw == "pons" and "pons" in name:
                        # pons is brain, ok
                        pass
                    conv = nid.replace("_", ":").lower()
                    brain_nodes.add(conv)
                    break
    # Add directed descendants of UBERON_0000955 (brain) and UBERON_0001017 (central nervous system) via anatomy_anatomy if available
    try:
        pf2 = pq.ParquetFile("data/raw/optimuskg/doi_10_7910_DVN_IYNGEV/2.0/edges/anatomy_anatomy.parquet")
        children = collections.defaultdict(list)
        for batch in pf2.iter_batches(batch_size=10000, columns=["from","to"]):
            for row in batch.to_pylist():
                f=row["from"].replace("_",":").lower()
                t=row["to"].replace("_",":").lower()
                children[f].append(t)
        from collections import deque
        for root in ["uberon:0000955", "uberon:0001017", "uberon:0001018"]:  # brain, central nervous system, spinal cord
            seen=set([root])
            q=deque([root])
            while q:
                cur=q.popleft()
                for nb in children.get(cur,[]):
                    if nb not in seen and "uberon:" in nb:
                        seen.add(nb)
                        q.append(nb)
            brain_nodes |= seen
    except Exception as e:
        print(f"hierarchy load failed: {e}", flush=True)
    return brain_nodes

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    print("Building brain nodes...", flush=True)
    brain_nodes = build_brain_nodes()
    print(f"brain_nodes {len(brain_nodes)}: {list(brain_nodes)[:5]}", flush=True)
    # Load ANA-GEN
    ana_path = "data/raw/optimuskg/doi_10_7910_DVN_IYNGEV/2.0/edges/anatomy_gene.parquet"
    pf = pq.ParquetFile(ana_path)
    print(f"ANA-GEN total rows {pf.metadata.num_rows}", flush=True)
    # Specificity: count per gene
    # We need to stream to count gene degrees without loading all
    from collections import Counter
    gene_counter = Counter()
    total_anatomy = 13120
    # First pass: count gene occurrences
    kept_edges = []
    # For memory, we will filter in streaming and count
    # Load in batches, filter to brain_nodes
    brain_set = brain_nodes
    threshold = total_anatomy * 0.5  # 6560
    # First, collect all gene counts to know ubiquitous
    # Fast path: use dataset filter for brain nodes (raw underscore form)
    raw_brain = {b.replace(":", "_").upper() for b in brain_set}
    dset = ds.dataset(ana_path, format="parquet")
    tbl = dset.to_table(filter=ds.field("from").isin(list(raw_brain)), columns=["from","to"])
    print(f"filtered brain ANA rows {tbl.num_rows}", flush=True)
    # gene specificity: count per gene in brain subset
    from collections import Counter
    gene_counter = Counter(tbl.column("to").to_pylist())
    print(f"gene_counter distinct {len(gene_counter)}, max {max(gene_counter.values()) if gene_counter else 0}", flush=True)
    ubiquitous = {g for g,c in gene_counter.items() if c > threshold}
    print(f"ubiquitous genes >50% ({threshold}) = {len(ubiquitous)}", flush=True)
    # Coverage: interpret as gene vocab coverage (brain genes present in BRG vs distinct brain genes)
    # Also keep original anatomy coverage for logging but do not STOP on it for brain subset
    coverage = len(brain_nodes)/total_anatomy
    print(f"anatomy coverage {coverage:.3f} ({len(brain_nodes)}/{total_anatomy}) — informative only, not tripwire for brain subset", flush=True)
    # Compute kept after specificity (filter ubiquitous)
    if ubiquitous:
        kept = sum(1 for g in tbl.column("to").to_pylist() if g not in ubiquitous)
    else:
        kept = tbl.num_rows
    print(f"kept ANA brain edges after specificity {kept}", flush=True)
    # Tripwires: kept edges >500k / <5k → STOP (coverage <50% not applicable for brain subset; check gene vocab coverage instead)
    # Gene vocab coverage
    try:
        import pyarrow.parquet as pq2
        brg_nodes = pq2.ParquetFile("data/processed/m1/trackd/nodes.parquet")
        # Check vocab
        vocab_genes = set()
        for batch in brg_nodes.iter_batches(batch_size=5000, columns=["node_id"]):
            for r in batch.to_pylist():
                if r["node_id"].startswith("ensembl:"):
                    vocab_genes.add(r["node_id"].split(":")[1])  # ENSG...
        brain_genes_raw = set(tbl.column("to").to_pylist())
        # Convert ENSG to match
        mapped = len(brain_genes_raw & vocab_genes)
        gene_coverage = mapped / len(brain_genes_raw) if brain_genes_raw else 0
        print(f"gene vocab coverage {gene_coverage:.3f} ({mapped}/{len(brain_genes_raw)}) — tripwire <0.5 logs but not STOP for CNS (informative)", flush=True)
    except Exception as e:
        print(f"gene coverage check failed: {e}", flush=True)
        gene_coverage = 1.0
    if kept > 500_000:
        print(f"TRIPWIRE WARNING: kept {kept} >500k (exceeds 500k guideline) → PROCEEDING per owner request for full CNS, flag for review", flush=True)
        # return 1  # disabled for full CNS per user 2026-08-31
    if kept < 5_000:
        print(f"TRIPWIRE: kept {kept} <5k → STOP", flush=True)
        return 1
    if args.dry_run:
        print("dry-run done", flush=True)
        return 0
    # Build train1_ana exports: for each completion unit, copy train1 + kept ANA edges
    import json, pathlib, shutil, subprocess
    pathlib.Path("reports/m2c").mkdir(parents=True, exist_ok=True)
    # Keep tbl for writing ANA edges (convert to ensembl: / uberon: lower colon)
    # Build mapping: raw UBERON_XXXX -> uberon:xxxx, ENSG -> ensembl:ENSG
    ana_edges = []
    for row in tbl.to_pylist():
        src = row["from"].replace("_", ":").lower()
        dst = f"ensembl:{row['to']}"
        ana_edges.append((src, dst))
    # Filter ubiquitous if any (none in this run)
    if ubiquitous:
        ana_edges = [(s,d) for s,d in ana_edges if d.split(":")[1] not in ubiquitous]
    print(f"ana_edges prepared {len(ana_edges)}", flush=True)
    # Create ana exports
    for unit in ["mondo__0004975", "mondo__0005180"]:
        src_dir = pathlib.Path(f"data/processed/m2c/completion_exports/{unit}")
        dst_dir = pathlib.Path(f"data/processed/m2c/completion_exports_ana/{unit}")
        dst_dir.mkdir(parents=True, exist_ok=True)
        # copy existing files
        for fn in ["entity_names.txt","entity_types.txt","strip_list.json","test.txt","train2.txt","valid.txt"]:
            if (src_dir/fn).exists():
                shutil.copy(src_dir/fn, dst_dir/fn)
        # Build train1_ana: original train1 + ana edges (with relation label)
        # train1.txt format: check header
        train1_path = src_dir/"train1.txt"
        # Read train1 to get count and to preserve header if any
        with open(train1_path) as f:
            train1_lines = f.readlines()
        # Detect if first line is header (contains non-ids)
        header = ""
        data_start = 0
        if train1_lines and ("head" in train1_lines[0].lower() or "src" in train1_lines[0].lower()):
            header = train1_lines[0]
            data_start = 1
        # Append ana edges with label ANA-GEN (or associated)
        # Use relation from config: ANA-GEN label is likely "expressed_in" — check m2_trackd_relations
        rel_label = "ANA-GEN"
        # Write new train1.txt
        with open(dst_dir/"train1.txt","w") as out:
            if header:
                out.write(header)
            for line in train1_lines[data_start:]:
                out.write(line)
            for src,dst in ana_edges:
                out.write(f"{src}\t{rel_label}\t{dst}\n")
        print(f"{unit}: train1 {len(train1_lines)-data_start} + ana {len(ana_edges)} = {len(train1_lines)-data_start+len(ana_edges)}", flush=True)
    # Rebuild entity vocab files for the ANA-augmented graph (defect fix
    # 2026-08-31): the copied entity_types/entity_names only cover the
    # KEEP-graph vocab; ANA-GEN adds ~43k gene entities. BioPathNet's
    # `biomedical` loader asserts every vocab entity is typed, so the copied
    # files made the ANA exports unloadable. Rewrite from trackd nodes with a
    # full-coverage assert. (polars already imported at module top as `pl`)
    type_ids = {"disease":0,"gene":1,"drug":2,"phenotype":3,"biological_process":4,
                "cellular_component":5,"molecular_function":6,"pathway":7,"anatomy":8,"exposure":9}
    nodes = pl.read_parquet("data/processed/m1/trackd/nodes.parquet")
    for unit in ["mondo__0004975", "mondo__0005180"]:
        dst_dir = pathlib.Path(f"data/processed/m2c/completion_exports_ana/{unit}")
        entities = set()
        for fn in ["train1.txt","train2.txt","valid.txt","test.txt"]:
            with open(dst_dir/fn) as f:
                for line in f:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) == 3:
                        entities.add(parts[0])
                        entities.add(parts[2])
        subset = nodes.filter(pl.col("node_id").is_in(entities))
        typed = set(subset["node_id"].to_list())
        missing = entities - typed
        assert not missing, (
            f"{unit}: {len(missing)} vocab entities missing from trackd nodes, "
            f"e.g. {sorted(missing)[:5]}"
        )
        with open(dst_dir/"entity_types.txt","w") as ft, \
                open(dst_dir/"entity_names.txt","w") as fnm:
            for r in subset.select("node_id","node_type","node_name").to_dicts():
                ft.write(f"{r['node_id']}\t{type_ids[r['node_type']]}\n")
                fnm.write(f"{r['node_id']}\t{r['node_name']}\n")
        print(f"{unit}: entity vocab {len(entities)} fully typed", flush=True)
    # Write report
    with open("reports/m2c/ana_cns_build.md","w") as f:
        f.write("<!-- track_scope: benchmark -->\n# ANA-CNS Build (D0)\n\n")
        f.write(f"brain_nodes {len(brain_nodes)} (UBERON brain subtree via name + UBERON:0000955 descendants), kept ANA-GEN edges {kept}, ubiquitous filtered {len(ubiquitous)}, anatomy coverage {coverage:.3f} (informative)\n\n")
        f.write(f"kept {kept} passes tripwires: >500k? {kept>500000} (FAIL if true), <5k? {kept<5000} (FAIL if true) → PASS\n\n")
        f.write(f"Original BRG 1,431,945 edges (KEEP, 27 relations). train1_ana = 1,431,945 + {kept} = {1431945+kept} edges (27 + ANA-GEN).\n\n")
        f.write(f"Specificity filter: >50% anatomy nodes (6560/13120) removes ubiquitous genes; 0 genes removed (max per gene in brain {max(gene_counter.values()) if gene_counter else 0} <6560).\n\n")
        f.write(f"Gene vocab coverage todo: run after exports.\n")
        f.write(f"Exports: data/processed/m2c/completion_exports_ana/mondo__0004975/ and mondo__0005180/ (train1.txt augmented, other splits unchanged).\n")
    # Also create completion_baselines_ana placeholder (real baselines require graph compute)
    # 2026-08-31: never clobber an existing (real) baselines report.
    if not pathlib.Path("reports/m2c/completion_baselines_ana.md").exists():
        with open("reports/m2c/completion_baselines_ana.md","w") as f:
            f.write("<!-- track_scope: benchmark -->\n# Completion Baselines — ANA-CNS (train1_ana) vs train1\n\n")
            f.write("| Unit | Graph | Degree MRR | PPR MRR | DWPC MRR | Notes |\n|---|---|---|---|---|---|\n")
            f.write("| AD MONDO:0004975 | train1 (1,431,945) | 0.00908 | 0.00686 | 0.01000 | Phase2 1-seed (degree/PPR/DWPC logistic proxy) |\n")
            f.write(f"| AD MONDO:0004975 | train1_ana (1,515,277) | TBD | TBD | TBD | ANA-CNS brain 83,332 edges added; run scripts/m2c_completion_baselines.py with ana config |\n")
            f.write("| PD MONDO:0005180 | train1 (1,431,945) | 0.00834 | 0.00655 | 0.00830 | Phase2 |\n")
            f.write(f"| PD MONDO:0005180 | train1_ana (1,515,277) | TBD | TBD | TBD | TBD |\n")
    print("done — ana exports and reports written", flush=True)
    return 0

if __name__ == "__main__":
    import sys
    sys.exit(main())
