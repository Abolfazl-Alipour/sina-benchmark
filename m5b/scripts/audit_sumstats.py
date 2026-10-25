"""Phase 0a: disease-side summary-stats availability audit.

For each M4B panel locus, decide whether full regional disease stats are
retrievable (GWAS Catalog harmonised sumstats keyed by study accession).
Outputs a coverage fraction + per-locus table; no modeling, no downloads
beyond the harmonised file list.
"""

import json
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
M4B_DIR = REPO_ROOT / "data" / "processed" / "m4b"
M5B_DIR = REPO_ROOT / "m5b"

GCST_RE = re.compile(r"GCST\d+")


def load_harmonised_gcsts(path: Path) -> set:
    gcsts = set()
    with path.open(encoding="utf-8", errors="replace") as f:
        for line in f:
            found = GCST_RE.findall(line)
            if found:
                gcsts.add(found[-1])
    return gcsts


def main() -> int:
    panels = json.loads((M4B_DIR / "panels.json").read_text(
        encoding="utf-8"))["loci"]
    ot_scores = json.loads((M4B_DIR / "ot_scores.json").read_text(
        encoding="utf-8"))["scores"]
    harm = load_harmonised_gcsts(M5B_DIR / "cache" / "harmonised_list.txt")
    by_locus: dict = {}
    for r in json.loads((M4B_DIR / "positives.json").read_text(
            encoding="utf-8"))["rows"]:
        if r.get("study_acc"):
            by_locus.setdefault((r["chr"], r["pos"]), set()).add(
                r["study_acc"])

    rows = []
    for loc in panels:
        accs = sorted(by_locus.get((loc["chr"], loc["pos"]), set()))
        have = sorted({a for acc in accs for a in GCST_RE.findall(acc)
                       if a in harm} | (set(accs) & harm))
        rows.append({"locus_id": loc["locus_id"], "chr": loc["chr"],
                     "pos": loc["pos"], "n_study_accs": len(accs),
                     "has_sumstats": bool(have), "gcsts": have,
                     "ot_covered": loc["locus_id"] in ot_scores})
    out = {"n_loci": len(rows),
           "frac_with_sumstats": round(
               sum(r["has_sumstats"] for r in rows) / max(len(rows), 1), 4),
           "covered": _sub(rows, True), "uncovered": _sub(rows, False),
           "rows": rows}
    (M5B_DIR / "reports" / "phase0a_sumstats_audit.json").write_text(
        json.dumps(out, indent=1), encoding="utf-8")
    print(f"loci {out['n_loci']}: with-sumstats "
          f"{out['frac_with_sumstats']}", flush=True)
    print(f"covered {out['covered']['n']}: "
          f"{out['covered']['frac_with_sumstats']}", flush=True)
    print(f"uncovered {out['uncovered']['n']}: "
          f"{out['uncovered']['frac_with_sumstats']}", flush=True)
    return 0


def _sub(rows: list, covered: bool) -> dict:
    s = [r for r in rows if r["ot_covered"] == covered]
    return {"n": len(s), "frac_with_sumstats": round(
        sum(r["has_sumstats"] for r in s) / max(len(s), 1), 4)}


if __name__ == "__main__":
    raise SystemExit(main())
