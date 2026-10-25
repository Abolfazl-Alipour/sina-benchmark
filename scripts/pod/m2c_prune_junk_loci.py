#!/usr/bin/env python
"""m2c_prune_junk_loci — build KEEP exports with non-pinned GEN nodes' edges removed.

Owner directive 2026-09-04: "kill them if you can prove that it is not helping (like ANA)".
This builds the pruned graph; the proof is the ±junk floor comparison
(reports/m2c/junk_loci_null.md), NOT the pruning itself.

Junk loci (C5): GEN-typed entities NOT in the pinned 17,448 candidate list —
database loci (pseudogenes, lncRNA, small-RNA, IG/TCR segments) that are never
scored, never negatives, never targets, but sit in the message graph. This script
removes every train1 edge touching them. train2/valid/test are copied
byte-identical (supervision is never pruned; C1 asserts their tails ⊆ pinned).

Usage:
    ~/miniconda3/envs/sina/bin/python scripts/pod/m2c_prune_junk_loci.py
"""
import argparse
import hashlib
import json
import shutil
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
PINNED_PATH = REPO_ROOT / "data/processed/m2c/pinned_genes_17448.txt"


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--exports-root", type=Path,
                    default=REPO_ROOT / "data/processed/m2c/completion_exports")
    ap.add_argument("--out-root", type=Path,
                    default=REPO_ROOT / "data/processed/m2c/completion_exports_nojunk")
    ap.add_argument("--units", nargs="+", default=["mondo__0005180", "mondo__0004975"])
    args = ap.parse_args()

    pinned = set(PINNED_PATH.read_text(encoding="utf-8").split())
    assert len(pinned) == 17448, f"pinned len {len(pinned)} != 17448"

    for unit in args.units:
        src = args.exports_root / unit
        dst = args.out_root / unit
        assert src.exists(), f"missing {src}"
        dst.mkdir(parents=True, exist_ok=True)

        # GEN type id for THIS export = type of the pinned genes (first-appearance map)
        type_of = {}
        with open(src / "entity_types.txt", encoding="utf-8") as f:
            for line in f:
                eid, tid = line.split()
                type_of[eid] = int(tid)
        pinned_types = {type_of[g] for g in pinned if g in type_of}
        assert len(pinned_types) == 1, f"pinned genes span multiple type ids {pinned_types}"
        gene_type = pinned_types.pop()
        missing = [g for g in pinned if g not in type_of]
        assert not missing, f"{len(missing)} pinned genes missing e.g. {missing[:3]}"

        junk = {e for e, t in type_of.items() if t == gene_type and e not in pinned}

        kept = removed = 0
        with open(src / "train1.txt", encoding="utf-8") as fin, \
                open(dst / "train1.txt", "w", encoding="utf-8") as fout:
            for line in fin:
                parts = line.split()
                if parts[0] in junk or parts[1] in junk:
                    removed += 1
                    continue
                fout.write(line)
                kept += 1

        # Supervision splits copied byte-identical — never pruned
        copies = {}
        for name in ("train2.txt", "valid.txt", "test.txt",
                     "entity_types.txt", "entity_names.txt", "strip_list.json"):
            shutil.copyfile(src / name, dst / name)
            copies[name] = sha256_file(dst / name)
            assert copies[name] == sha256_file(src / name), f"copy mismatch {name}"

        manifest = {
            "unit": unit,
            "source": str(src),
            "gene_type_id": gene_type,
            "pinned_genes": len(pinned),
            "gen_nodes_total": sum(1 for t in type_of.values() if t == gene_type),
            "junk_loci": len(junk),
            "train1_edges_in": kept + removed,
            "train1_edges_kept": kept,
            "train1_edges_removed": removed,
            "pct_removed": round(100.0 * removed / (kept + removed), 2),
            "copies_sha256": copies,
            "train1_sha256": sha256_file(dst / "train1.txt"),
        }
        (dst / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        print(f"{unit}: GEN nodes {manifest['gen_nodes_total']}, junk {len(junk)}, "
              f"train1 {kept + removed} -> {kept} "
              f"(removed {removed}, {manifest['pct_removed']}%)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
