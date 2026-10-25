#!/usr/bin/env python
"""m2c_ot_panel_eval — Open Targets panel pull, two-stage rescore, circularity audit.

Closes the two open items from reports/m2c/cpu_nulls_and_coverage.md §5 on the
no-ANA rank-loss checkpoint (reports/m2c/checkpoints_completion_noana/):
  1. Re-pull OT associatedTargets panels (the beast-era panel files were never
     committed; only sha256 prefixes survive in the report). Reproduction is
     verified by (a) association counts PD 8,923 / AD 13,367, (b) mapped counts
     8,274 / 12,579, (c) re-scoring the ANA-wave ep3 checkpoint reproducing the
     reported two-stage 0.0295 / hit-split 4/6/6 on PD.
  2. Circularity audit: per-gene OT datasource score shares (harmonic-sum
     decomposition, the OT overall-score formula) grouped into evidence classes,
     quantifying how much of each panel's score mass is literature/curated-class
     — the same evidence class as the DIS-GEN supervision. OptimusKG DIS-GEN
     edges carry no per-edge provenance in-repo, so the audit is at the
     evidence-class level, not record level. Stated assumption, not a fix.

Panel construction (pinned, matches cpu_nulls_and_coverage.md §5):
  candidates = OT rows mapped to pinned 17,448, excluding train2+valid positives
  of the disease head; sort by (-overall score, gene token); top-K,
  K in {100, 200, 500}. Serialization: one `ensembl:ENSG...` token per line,
  trailing LF; sha256 over the file bytes.

Subcommands (two envs):
  pull        sina env, CPU+internet:
              ~/miniconda3/envs/sina/bin/python scripts/pod/m2c_ot_panel_eval.py pull
  score       bpn120 env, CPU forward (zero GPU):
              CUDA_VISIBLE_DEVICES="" ~/miniconda3/envs/bpn120/bin/python -u \
                scripts/pod/m2c_ot_panel_eval.py score --unit mondo__0005180 \
                --exports-dir data/processed/m2c/completion_exports \
                --checkpoints reports/m2c/checkpoints_completion_noana/mondo__0005180/epoch_001.pth
  circularity sina env (reads raw pull):
              ~/miniconda3/envs/sina/bin/python scripts/pod/m2c_ot_panel_eval.py circularity
"""
import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "third_party" / "biopathnet"))
sys.path.insert(0, str(REPO_ROOT / "src"))

PINNED_PATH = REPO_ROOT / "data/processed/m2c/pinned_genes_17448.txt"
DEFAULT_EXPORTS = REPO_ROOT / "data/processed/m2c/completion_exports"
PANELS_DIR = REPO_ROOT / "reports/m2c/ot_panels"
RAW_DIR = REPO_ROOT / "data/processed/m2c/ot_panels"
OT_API = "https://api.platform.opentargets.org/api/v4/graphql"

UNITS = {
    "mondo__0005180": {"ot_id": "MONDO_0005180", "label": "PD",
                       "expect_count": 8923, "expect_mapped": 8274, "n_test": 339},
    "mondo__0004975": {"ot_id": "MONDO_0004975", "label": "AD",
                       "expect_count": 13367, "expect_mapped": 12579, "n_test": 621},
}
KS = (100, 200, 500)

# Evidence classes for the circularity audit (OT datasource ids, observed 2026-09).
CLASS_OF_DATASOURCE = {
    "europepmc": "literature",
    "uniprot_literature": "literature",
    "genomics_england": "curated_db",
    "clingen": "curated_db",
    "orphanet": "curated_db",
    "gene2phenotype": "curated_db",
    "cancer_gene_census": "curated_db",
    "intogen": "curated_db",
    "gwas_credible_sets": "genetics",
    "ot_genetics_portal": "genetics",
    "eva": "genetics",
    "gene_burden": "genetics",
    "uniprot_variants": "genetics",
}
# Anything unlisted (impc, crispr_screen, expression_atlas, ...) -> "functional_other"


def utcnow_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_pinned() -> list[str]:
    pinned = PINNED_PATH.read_text(encoding="utf-8").strip().split()
    assert len(pinned) == 17448, f"pinned list drifted: {len(pinned)}"
    return pinned


def known_positives(exports_dir: Path, unit: str) -> set[str]:
    """train2+valid tails of the disease head (the 'already-known' set).
    Export rows are `head<TAB>relation<TAB>tail`."""
    head = unit.replace("__", ":")
    known = set()
    unit_dir = exports_dir / unit
    for fname in ("train2.txt", "valid.txt"):
        for line in (unit_dir / fname).read_text(encoding="utf-8").splitlines():
            h, _r, t = line.split("\t")
            if h == head:
                known.add(t)
    return known


def test_positives(exports_dir: Path, unit: str) -> set[str]:
    head = unit.replace("__", ":")
    out = set()
    for line in (exports_dir / unit / "test.txt").read_text(encoding="utf-8").splitlines():
        h, _r, t = line.split("\t")
        if h == head:
            out.add(t)
    return out


# ---------------------------------------------------------------- pull


def ot_fetch_all(ot_id: str, expect_count: int) -> list[dict]:
    """Page the full associatedTargets list (default order is score desc —
    verified empirically 2026-09-05; we re-sort client-side regardless)."""
    import requests

    query = """
    query($id: String!, $index: Int!, $size: Int!) {
      disease(efoId: $id) {
        associatedTargets(page: {index: $index, size: $size}) {
          count
          rows { target { id approvedSymbol } score datasourceScores { id score } }
        }
      }
    }"""
    rows: list[dict] = []
    size, index = 1000, 0
    count = None
    while count is None or len(rows) < count:
        resp = requests.post(
            OT_API, json={"query": query,
                          "variables": {"id": ot_id, "index": index, "size": size}},
            timeout=60)
        resp.raise_for_status()
        payload = resp.json()
        if "errors" in payload:
            raise RuntimeError(f"OT GraphQL errors: {payload['errors']}")
        block = payload["data"]["disease"]["associatedTargets"]
        count = block["count"]
        rows.extend(block["rows"])
        index += 1
        print(f"[{utcnow_iso()}] {ot_id} page {index}: {len(rows)}/{count}", flush=True)
    assert len(rows) == count, f"{ot_id}: got {len(rows)} rows, count says {count}"
    assert count == expect_count, (
        f"{ot_id} association count {count} != {expect_count} (2026-09-04 pull) — "
        f"OT data release drifted; STOP and record deviation before comparing panels")
    return rows


def cmd_pull(args) -> int:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    PANELS_DIR.mkdir(parents=True, exist_ok=True)
    pinned = set(load_pinned())
    manifest = {"pulled_at": utcnow_iso(), "api": OT_API, "units": {}}
    for unit, spec in UNITS.items():
        rows = ot_fetch_all(spec["ot_id"], spec["expect_count"])
        raw_path = RAW_DIR / f"raw_{unit}.jsonl"
        with raw_path.open("w", encoding="utf-8") as f:
            for r in rows:
                f.write(json.dumps(r, sort_keys=True) + "\n")

        known = known_positives(args.exports_dir, unit)
        test_pos = test_positives(args.exports_dir, unit)
        assert len(test_pos) == spec["n_test"]

        mapped = []  # (token, ot_score, symbol)
        for r in rows:
            token = f"ensembl:{r['target']['id']}"
            if token in pinned:
                mapped.append((token, float(r["score"]), r["target"]["approvedSymbol"]))
        # 2026-09-05 note: re-pull maps 8,631 PD rows (beast reported 8,274 on
        # 2026-09-04 with its script uncommitted and rule unrecoverable; OT count
        # 8,923 matches exactly, so the delta is a mapping-rule difference, not
        # API drift). Recorded as a deviation; reproduction of the §5 eval
        # numbers via the ANA-ep3 control decides materiality.
        if len(mapped) != spec["expect_mapped"]:
            print(f"[{utcnow_iso()}] WARNING {unit}: mapped {len(mapped)} != "
                  f"{spec['expect_mapped']} (beast 2026-09-04) — recorded deviation, "
                  f"continuing (see script docstring)", flush=True)

        # Raw (pre-exclusion) knowns density: reproduces the §5 "92% already-known" stat.
        raw_sorted = sorted(mapped, key=lambda x: (-x[1], x[0]))
        raw_knowns = {}
        for k in KS:
            top = raw_sorted[:k]
            raw_knowns[str(k)] = sum(1 for t, _s, _y in top if t in known)

        eligible = [m for m in mapped if m[0] not in known]
        eligible.sort(key=lambda x: (-x[1], x[0]))
        unit_manifest = {"ot_id": spec["ot_id"], "label": spec["label"],
                         "count": len(rows), "mapped": len(mapped),
                         "beast_mapped_2026_09_04": spec["expect_mapped"],
                         "mapped_deviation": len(mapped) - spec["expect_mapped"],
                         "knowns_train2_valid": len(known),
                         "knowns_in_ot_mapped": len(known & {m[0] for m in mapped}),
                         "eligible_after_excl_train2_valid": len(eligible),
                         "raw_sha256": sha256_file(raw_path),
                         "raw_topk_knowns_density": raw_knowns, "panels": {}}
        for k in KS:
            panel = eligible[:k]
            panel_path = PANELS_DIR / f"{unit}_ot_top{k}_excl_knowns.txt"
            panel_path.write_text("".join(f"{t}\n" for t, _s, _y in panel),
                                  encoding="utf-8")
            cov = sum(1 for t, _s, _y in panel if t in test_pos)
            unit_manifest["panels"][str(k)] = {
                "path": str(panel_path.relative_to(REPO_ROOT)),
                "sha256": sha256_file(panel_path),
                "test_coverage": cov / len(test_pos), "test_hits_in_panel": cov,
                "panel_min_score": panel[-1][1] if panel else None}
            print(f"[{utcnow_iso()}] {unit} K={k}: coverage {cov}/{len(test_pos)} "
                  f"sha {unit_manifest['panels'][str(k)]['sha256'][:8]}…", flush=True)
        manifest["units"][unit] = unit_manifest
    manifest_path = PANELS_DIR / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                             encoding="utf-8")
    print(f"[{utcnow_iso()}] manifest -> {manifest_path}", flush=True)
    return 0


# ---------------------------------------------------------------- score


def cmd_score(args) -> int:
    """One CPU forward per checkpoint over pinned tails; two-stage recall@10 per
    panel + OT-ordering null on the same panel. All test queries share the
    (disease, DIS-GEN) (h, r), so one forward scores every query exactly
    (m2c_pd_test_once.py semantics)."""
    import numpy as np
    import torch
    from biopathnet import dataset as bpn_dataset  # noqa: F401
    from biopathnet import model as bpn_model  # noqa: F401
    from biopathnet import task as bpn_task  # noqa: F401
    from biopathnet import util
    from torchdrug import core

    torch.set_num_threads(args.threads)
    pinned = load_pinned()
    unit = args.unit
    n_test = UNITS[unit]["n_test"]

    cfg = util.load_config(str(args.config), context={"gpus": []})
    cfg.dataset.path = str(args.exports_dir / unit)
    cfg.engine.gpus = []
    dataset_obj = core.Configurable.load_config_dict(cfg.dataset)
    task_cfg = cfg.task
    task_cfg.model.num_relation = dataset_obj.num_relation
    task_obj = core.Configurable.load_config_dict(task_cfg)
    train_set, valid_set, test_set = dataset_obj.split()
    task_obj.preprocess(train_set, valid_set, test_set)
    model = task_obj.model
    model.eval()
    fact_graph = task_obj.fact_graph
    inv_entity = dataset_obj.inv_entity_vocab

    def edges_of(subset):
        return [tuple(dataset_obj.graph.edge_list[i].tolist()) for i in subset.indices]

    train2_edges = edges_of(train_set)
    valid_edges = edges_of(valid_set)
    test_edges = edges_of(test_set)
    rel_ids = {r for _, _, r in train2_edges}
    assert len(rel_ids) == 1
    dis_gen_rel = next(iter(rel_ids))
    disease_idx = inv_entity[unit.replace("__", ":")]
    test_tails = sorted(t for h, t, r in test_edges if r == dis_gen_rel and h == disease_idx)
    assert len(test_tails) == len(test_edges) == n_test
    train_pos = {t for h, t, r in train2_edges if r == dis_gen_rel and h == disease_idx}
    valid_pos = {t for h, t, r in valid_edges if r == dis_gen_rel and h == disease_idx}
    excluded = train_pos | valid_pos

    pinned_indices = [inv_entity[g] for g in pinned]
    gene_of_pos = pinned  # positional
    h_index = torch.full((1, len(pinned)), disease_idx, dtype=torch.long)
    r_index = torch.full((1, len(pinned)), dis_gen_rel, dtype=torch.long)
    t_index = torch.tensor(pinned_indices, dtype=torch.long).unsqueeze(0)
    test_pos_set = set(test_tails)
    test_positions = [i for i, g in enumerate(pinned_indices) if g in test_pos_set]
    assert len(test_positions) == n_test

    def bootstrap_ci(values, n_boot, seed):
        rng = np.random.default_rng(seed)
        reps = np.empty(n_boot)
        n = len(values)
        for b in range(n_boot):
            reps[b] = float(np.mean(values[rng.integers(0, n, n)]))
        lo, hi = np.percentile(reps, [2.5, 97.5])
        return float(lo), float(hi)

    panels = {}
    for k in KS:
        p = PANELS_DIR / f"{unit}_ot_top{k}_excl_knowns.txt"
        toks = p.read_text(encoding="utf-8").split()
        panels[k] = toks
    ot_score_of = {}  # token -> ot score (for the null ordering)
    with (RAW_DIR / f"raw_{unit}.jsonl").open(encoding="utf-8") as f:
        for line in f:
            r = json.loads(line)
            ot_score_of[f"ensembl:{r['target']['id']}"] = float(r["score"])

    results = {"unit": unit, "exports_dir": str(args.exports_dir), "runs": {}}
    for ckpt_path in args.checkpoints:
        ckpt_path = Path(ckpt_path)
        ckpt = torch.load(ckpt_path, map_location="cpu")
        model.load_state_dict(ckpt["model"])
        model.eval()
        with torch.no_grad():
            scores_t = model(fact_graph, h_index, t_index, r_index,
                             all_loss=None, metric=None)
        scores = scores_t[0].double().tolist()
        assert len(scores) == len(pinned)
        # Genome-wide filtered cross-check (must reproduce runner/test-once numbers
        # for the corresponding checkpoint — verifies evaluator identity).
        gw_positions = [i for i, g in enumerate(pinned_indices) if g not in excluded]
        gw_sorted = sorted(gw_positions, key=lambda i: (-scores[i], gene_of_pos[i]))
        rank_of_gw = {}
        for rank1, i in enumerate(gw_sorted, start=1):
            if i not in rank_of_gw:  # unique positions anyway (positional sort)
                rank_of_gw[i] = rank1
        gw_rr = [1.0 / rank_of_gw[i] for i in test_positions]
        gw_h10 = [1.0 if rank_of_gw[i] <= 10 else 0.0 for i in test_positions]
        gw_mrr = float(np.mean(gw_rr))
        gw_hits10 = float(np.mean(gw_h10))

        run = {"checkpoint": str(ckpt_path), "epoch": ckpt.get("epoch"),
               "ckpt_sha256": sha256_file(ckpt_path),
               "genomewide_test_mrr": gw_mrr, "genomewide_test_hits10": gw_hits10,
               "panels": {}}
        for k, toks in panels.items():
            panel_set = set(toks)
            panel_pos = [i for i, g in enumerate(pinned_indices)
                         if gene_of_pos[i] in panel_set]
            assert len(panel_pos) == len(toks), "panel token not in pinned"
            # model ranking within panel (optimistic ties via sorted order on
            # (-score, token), identical tiebreak to panel construction)
            order = sorted(panel_pos, key=lambda i: (-scores[i], gene_of_pos[i]))
            model_rank = {i: r1 for r1, i in enumerate(order, start=1)}
            # OT null ordering within panel
            null_order = sorted(panel_pos, key=lambda i: (-ot_score_of[gene_of_pos[i]],
                                                          gene_of_pos[i]))
            null_rank = {i: r1 for r1, i in enumerate(null_order, start=1)}
            in_panel = [i for i in test_positions if i in model_rank]
            model_hit = [i for i in in_panel if model_rank[i] <= 10]
            null_hit = [i for i in in_panel if null_rank[i] <= 10]
            both = [i for i in model_hit if i in null_hit]
            per_query_model = [1.0 if i in model_rank and model_rank[i] <= 10 else 0.0
                               for i in test_positions]
            per_query_null = [1.0 if i in null_rank and null_rank[i] <= 10 else 0.0
                              for i in test_positions]
            lo_m, hi_m = bootstrap_ci(np.array(per_query_model), args.bootstrap, args.seed)
            lo_n, hi_n = bootstrap_ci(np.array(per_query_null), args.bootstrap, args.seed)
            panel_file = PANELS_DIR / f"{unit}_ot_top{k}_excl_knowns.txt"
            run["panels"][str(k)] = {
                "n_panel": len(toks), "sha256": sha256_file(panel_file),
                "coverage": len(in_panel) / n_test,
                "model_two_stage": float(np.mean(per_query_model)),
                "model_two_stage_ci": [lo_m, hi_m],
                "ot_null_two_stage": float(np.mean(per_query_null)),
                "ot_null_two_stage_ci": [lo_n, hi_n],
                "hits_both": len(both),
                "hits_model_only": len(model_hit) - len(both),
                "hits_null_only": len(null_hit) - len(both),
                "model_hit_genes": [gene_of_pos[i] for i in model_hit],
                "null_hit_genes": [gene_of_pos[i] for i in null_hit],
            }
            print(f"[{utcnow_iso()}] {ckpt_path.name} K={k}: model "
                  f"{run['panels'][str(k)]['model_two_stage']:.4f} vs OT-null "
                  f"{run['panels'][str(k)]['ot_null_two_stage']:.4f} "
                  f"(coverage {run['panels'][str(k)]['coverage']:.4f})", flush=True)
        print(f"[{utcnow_iso()}] {ckpt_path.name} genome-wide cross-check: "
              f"test MRR {gw_mrr:.5f} H@10 {gw_hits10:.5f}", flush=True)
        results["runs"][ckpt_path.name] = run
        if args.dump_scores_dir:
            dump = {"unit": unit, "checkpoint": str(ckpt_path),
                    "ckpt_sha256": run["ckpt_sha256"], "epoch": ckpt.get("epoch"),
                    "pinned_order": "pinned_genes_17448.txt line order",
                    "scores": scores}
            args.dump_scores_dir.mkdir(parents=True, exist_ok=True)
            dump_path = args.dump_scores_dir / f"{unit}__{ckpt_path.stem}.scores.json"
            dump_path.write_text(json.dumps(dump), encoding="utf-8")
            print(f"[{utcnow_iso()}] dumped score vector -> {dump_path}", flush=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")
    print(f"[{utcnow_iso()}] wrote {args.out}", flush=True)
    return 0


# ---------------------------------------------------------------- circularity


def harmonic_shares(datasource_scores: list[dict]) -> dict[str, float]:
    """OT overall score = harmonic sum of sorted-desc datasource scores
    (s_1/1^2 + s_2/2^2 + ...). Return per-datasource contribution shares."""
    ss = sorted((float(d["score"]) for d in datasource_scores), reverse=True)
    total = sum(s / (i + 1) ** 2 for i, s in enumerate(ss))
    shares: dict[str, float] = {}
    # shares are per-position; attribute positions back to datasource ids
    ordered = sorted(datasource_scores, key=lambda d: -float(d["score"]))
    for i, d in enumerate(ordered):
        shares[d["id"]] = (float(d["score"]) / (i + 1) ** 2) / total if total else 0.0
    return shares


def cmd_circularity(args) -> int:
    PANELS_DIR.mkdir(parents=True, exist_ok=True)
    report = {"generated_at": utcnow_iso(), "assumption": (
        "OptimusKG DIS-GEN edges carry no per-edge provenance in-repo; the audit "
        "quantifies evidence-CLASS overlap (literature/curated_db vs genetics vs "
        "functional_other) between the OT overall score and the supervision set, "
        "not record-level overlap."), "units": {}}
    for unit, spec in UNITS.items():
        raw_path = RAW_DIR / f"raw_{unit}.jsonl"
        rows = {}
        with raw_path.open(encoding="utf-8") as f:
            for line in f:
                r = json.loads(line)
                rows[f"ensembl:{r['target']['id']}"] = r
        unit_rep = {"label": spec["label"], "panels": {}}
        for k in KS:
            panel_path = PANELS_DIR / f"{unit}_ot_top{k}_excl_knowns.txt"
            toks = panel_path.read_text(encoding="utf-8").split()
            per_gene = []
            class_mass = {"literature": 0.0, "curated_db": 0.0,
                          "genetics": 0.0, "functional_other": 0.0}
            for t in toks:
                shares = harmonic_shares(rows[t]["datasourceScores"])
                cls = {"literature": 0.0, "curated_db": 0.0,
                       "genetics": 0.0, "functional_other": 0.0}
                for ds, sh in shares.items():
                    cls[CLASS_OF_DATASOURCE.get(ds, "functional_other")] += sh
                for c, v in cls.items():
                    class_mass[c] += v
                per_gene.append({"gene": t,
                                 "symbol": rows[t]["target"]["approvedSymbol"],
                                 "ot_score": rows[t]["score"], **cls})
            n = max(len(per_gene), 1)
            mean_shares = {c: class_mass[c] / n for c in class_mass}
            lit_cur = [g["literature"] + g["curated_db"] for g in per_gene]
            unit_rep["panels"][str(k)] = {
                "mean_score_share_by_class": mean_shares,
                "frac_genes_literature_curated_over_50pct": (
                    sum(1 for v in lit_cur if v > 0.5) / n),
                "frac_genes_literature_curated_over_80pct": (
                    sum(1 for v in lit_cur if v > 0.8) / n),
                "per_gene": per_gene if k == 100 else "see K=100 table schema",
            }
            frac50 = unit_rep["panels"][str(k)]["frac_genes_literature_curated_over_50pct"]
            print(f"[{utcnow_iso()}] {unit} K={k}: mean shares "
                  + " ".join(f"{c}={v:.3f}" for c, v in mean_shares.items())
                  + f" | lit+curated>50%: {frac50:.2f}",
                  flush=True)
        report["units"][unit] = unit_rep
    out = PANELS_DIR / "circularity_audit.json"
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[{utcnow_iso()}] wrote {out}", flush=True)
    return 0


# ---------------------------------------------------------------- union


def cmd_union(args) -> int:
    """Union-ranker eval over two dumped score vectors (e.g. ANA-ep3 + noana).

    Per panel K: per-gene within-panel ranks from each model; report
      - union upper bound: P(test gene in A-top10 union B-top10 within panel)
      - RRF combined: rank by rrf = 1/(60+rankA) + 1/(60+rankB), two-stage@10
      - OT null on the same panel (from raw pull scores)
    Bootstrap CIs over the unit's test queries (seed 1337, 1000 resamples).
    """
    import numpy as np

    unit = args.unit
    n_test = UNITS[unit]["n_test"]
    dumps = []
    for p in args.scores:
        d = json.loads(Path(p).read_text(encoding="utf-8"))
        assert d["unit"] == unit, f"{p}: unit mismatch"
        dumps.append(d)
    assert len(dumps) == 2
    assert len(dumps[0]["scores"]) == len(dumps[1]["scores"]) == 17448

    pinned = load_pinned()
    pos_of_token = {t: i for i, t in enumerate(pinned)}
    test_pos = sorted(test_positives(args.exports_dir, unit))
    test_positions = [pos_of_token[t] for t in test_pos]
    assert len(test_positions) == n_test

    ot_score_of = {}
    with (RAW_DIR / f"raw_{unit}.jsonl").open(encoding="utf-8") as f:
        for line in f:
            r = json.loads(line)
            ot_score_of[f"ensembl:{r['target']['id']}"] = float(r["score"])

    def bootstrap_ci(values, n_boot=1000, seed=1337):
        rng = np.random.default_rng(seed)
        reps = np.empty(n_boot)
        n = len(values)
        for b in range(n_boot):
            reps[b] = float(np.mean(values[rng.integers(0, n, n)]))
        lo, hi = np.percentile(reps, [2.5, 97.5])
        return float(lo), float(hi)

    out = {"unit": unit, "inputs": [{"path": str(p), "sha256": d["ckpt_sha256"]}
                                    for p, d in zip(args.scores, dumps, strict=True)],
           "combiners": {"union_upper_bound": "test gene in A-top10 union B-top10",
                         "rrf": "rank by 1/(60+rankA)+1/(60+rankB), top10"},
           "panels": {}}
    for k in KS:
        toks = (PANELS_DIR / f"{unit}_ot_top{k}_excl_knowns.txt"
                ).read_text(encoding="utf-8").split()
        panel_idx = [pos_of_token[t] for t in toks]
        panel_set = set(panel_idx)
        ranks = []
        for d in dumps:
            sc = d["scores"]
            order = sorted(panel_idx, key=lambda i: (-sc[i], pinned[i]))
            ranks.append({i: r1 for r1, i in enumerate(order, start=1)})
        rA, rB = ranks
        null_order = sorted(panel_idx, key=lambda i: (-ot_score_of[pinned[i]], pinned[i]))
        rN = {i: r1 for r1, i in enumerate(null_order, start=1)}
        rrf_order = sorted(panel_idx,
                           key=lambda i: (-(1 / (60 + rA[i]) + 1 / (60 + rB[i])),
                                          pinned[i]))
        rR = {i: r1 for r1, i in enumerate(rrf_order, start=1)}

        def hits(rule, panel_set=panel_set):
            return [1.0 if (i in panel_set and rule(i) <= 10) else 0.0
                    for i in test_positions]

        v_union = hits(lambda i, a=rA, b=rB: min(a[i], b[i]))
        v_rrf = hits(lambda i, r=rR: r[i])
        v_null = hits(lambda i, n=rN: n[i])
        in_panel = [i for i in test_positions if i in panel_set]
        out["panels"][str(k)] = {
            "coverage": len(in_panel) / n_test,
            "union_upper_bound": float(np.mean(v_union)),
            "union_upper_bound_ci": list(bootstrap_ci(np.array(v_union))),
            "rrf_two_stage": float(np.mean(v_rrf)),
            "rrf_two_stage_ci": list(bootstrap_ci(np.array(v_rrf))),
            "ot_null_two_stage": float(np.mean(v_null)),
            "ot_null_two_stage_ci": list(bootstrap_ci(np.array(v_null))),
        }
        per_gene_hits = [i for i in test_positions
                         if i in panel_set and min(rA[i], rB[i]) <= 10]
        out["panels"][str(k)]["union_hit_genes"] = [pinned[i] for i in per_gene_hits]
        ub = out["panels"][str(k)]["union_upper_bound"]
        rrf_v = out["panels"][str(k)]["rrf_two_stage"]
        null_v = out["panels"][str(k)]["ot_null_two_stage"]
        print(f"[{utcnow_iso()}] union K={k}: bound {ub:.4f} rrf {rrf_v:.4f} "
              f"null {null_v:.4f}", flush=True)
    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[{utcnow_iso()}] wrote {args.out}", flush=True)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_pull = sub.add_parser("pull")
    p_pull.add_argument("--exports-dir", type=Path, default=DEFAULT_EXPORTS)
    p_pull.set_defaults(fn=cmd_pull)

    p_score = sub.add_parser("score")
    p_score.add_argument("--unit", default="mondo__0005180")
    p_score.add_argument("--config", type=Path, default=REPO_ROOT / "configs/m2c_arm_e.yaml")
    p_score.add_argument("--exports-dir", type=Path, default=DEFAULT_EXPORTS)
    p_score.add_argument("--checkpoints", nargs="+", required=True)
    p_score.add_argument("--bootstrap", type=int, default=1000)
    p_score.add_argument("--seed", type=int, default=1337)
    p_score.add_argument("--threads", type=int, default=6)
    p_score.add_argument("--out", type=Path,
                         default=PANELS_DIR / "score_noana_epoch001.json")
    p_score.add_argument("--dump-scores-dir", type=Path, default=None,
                         help="if set, dump each checkpoint's pinned-order score vector here")
    p_score.set_defaults(fn=cmd_score)

    p_circ = sub.add_parser("circularity")
    p_circ.set_defaults(fn=cmd_circularity)

    p_un = sub.add_parser("union")
    p_un.add_argument("--unit", default="mondo__0005180")
    p_un.add_argument("--exports-dir", type=Path, default=DEFAULT_EXPORTS)
    p_un.add_argument("--scores", nargs="+", required=True,
                      help="two .scores.json dumps from `score --dump-scores-dir`")
    p_un.add_argument("--out", type=Path,
                      default=PANELS_DIR / "union_noana_x_anaep3.json")
    p_un.set_defaults(fn=cmd_union)

    args = parser.parse_args()
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main())
