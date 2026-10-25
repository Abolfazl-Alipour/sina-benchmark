#!/usr/bin/env python
"""m2c_loci_headtohead — SINA model vs Open Targets Genetics L2G at locus ranking.

Owner task brief 2026-09-11 (benchmark track, CPU-only, zero GPU spend): on GWAS
loci for PD/AD, restricted to the pinned 17,448 universe, which ranks the known
disease gene higher within its locus — the saved model score vector or the OT
Genetics L2G score? Pre-registered decision rule keys on per-locus top-1 rate
(model minus OT >= 5 points with bootstrap CI excluding 0 = WIN), PD and AD
separately, no pooling.

Pipeline (subcommands):
  gate   reproduce the published per-locus model-only top-1 (PD 0.904 / AD 0.954,
         `reports/m2c/prospective_per_locus.md` Probe B) from the SAME score
         dumps + loci + definitions. Mismatch beyond +/-2 points -> exit 2, STOP.
  fetch  pull OT Genetics L2G via the Platform credibleSets API
         (api.genetics.opentargets.org has no DNS record; L2G now lives at
         api.platform.opentargets.org/api/v4/graphql). Per study: enumerate
         credible sets; per GWAS member (study_acc, chr, pos): match the
         credible set with the same study and lead-variant chr/pos; cache each
         matched set's l2GPredictions rows. Raw JSON cached under
         reports/m2c/ot_l2g/<unit>/ with SHA-256 manifest; re-pull is
         byte-identical (cache-first + --verify-fetch spot check).
  run    gate first (abort on discrepancy), then per-locus head-to-head over
         model / OT-L2G / degree / chance arms with pessimistic ties in every
         arm, identical known-gene definition (train2 head positives) in every
         arm, bootstrap CIs, pre/post-2020 strata, OT-missing loci counted as
         OT misses in the primary view + covered-subset view.

Usage (sina env, CPU-only, no GPU job is authorized):
  ~/miniconda3/envs/sina/bin/python scripts/pod/m2c_loci_headtohead.py gate --unit mondo__0005180
  ~/miniconda3/envs/sina/bin/python scripts/pod/m2c_loci_headtohead.py fetch --unit mondo__0005180
  ~/miniconda3/envs/sina/bin/python scripts/pod/m2c_loci_headtohead.py run --unit mondo__0005180 \\
      --out reports/m2c/loci_headtohead_pd.json --report reports/m2c/loci_headtohead.md
"""

import argparse
import hashlib
import json
import time
from collections import defaultdict
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
PINNED_PATH = REPO_ROOT / "data" / "processed" / "m2c" / "pinned_genes_17448.txt"
DEFAULT_EXPORTS = REPO_ROOT / "data" / "processed" / "m2c" / "completion_exports"
DEFAULT_OT_CACHE = REPO_ROOT / "reports" / "m2c" / "ot_l2g"
OT_API = "https://api.platform.opentargets.org/api/v4/graphql"
GATE_EXPECTED = {"mondo__0005180": 0.904, "mondo__0004975": 0.954}
GATE_TOL = 0.02


def utcnow_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_text(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()


def pessimistic_ranks(scores: Mapping[str, float],
                      genes: Sequence[str]) -> dict[str, int]:
    """Worst-rank for every gene: 1 + strictly-greater + tied-others.

    Deterministic (no RNG); exact ties share the worst rank of the tie block.
    """
    order = sorted(genes, key=lambda g: (-scores[g], g))
    ranks: dict[str, int] = {}
    i = 0
    n = len(order)
    while i < n:
        j = i
        while j + 1 < n and scores[order[j + 1]] == scores[order[i]]:
            j += 1
        for g in order[i:j + 1]:
            ranks[g] = j + 1
        i = j + 1
    return ranks


def best_known_rank(ranks: Mapping[str, int], known: set[str]) -> int | None:
    """Minimum rank among known genes present in the ranking; None if absent."""
    vals = [r for g, r in ranks.items() if g in known]
    return min(vals) if vals else None


def hit_at_k(best_rank: int | None, k: int) -> float:
    return 1.0 if best_rank is not None and best_rank <= k else 0.0


def bootstrap_ci(values: np.ndarray, n_boot: int, seed: int) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    n = len(values)
    means = np.array([float(np.mean(values[rng.integers(0, n, n)]))
                      for _ in range(n_boot)])
    lo, hi = np.percentile(means, [2.5, 97.5])
    return float(lo), float(hi)


def chance_draw_means(probs: Sequence[float], n_draws: int,
                      seed: int) -> np.ndarray:
    """One seeded top-1 pick per locus per draw; returns mean hit rate per draw.

    A random top-1 pick among n_genes with n_known knowns hits with
    probability n_known/n_genes per locus.
    """
    rng = np.random.default_rng(seed)
    probs = np.asarray(probs, dtype=float)
    draws = rng.random((n_draws, len(probs))) < probs[None, :]
    return draws.mean(axis=1)


def classify_locus(n_pinned_genes: int, n_known: int) -> str:
    """gradeable | ungradeable-few-genes | no-known (listed, never pooled)."""
    if n_pinned_genes < 2:
        return "ungradeable-few-genes"
    if n_known < 1:
        return "no-known"
    return "gradeable"


def load_pinned(pinned_arg: Path | None, exports_dir: Path,
                unit: str) -> tuple[list[str], str]:
    """Canonical pinned list; fallback derives token-1 from entity_types.

    Deviation (used on hosts missing the canonical artifact): derive from
    {exports_dir}/{unit}/entity_types.txt type token 1, assert 17,448.
    Proven equivalent 2026-09-04 (exact Probe 1.5 reproduction); the
    regression gate re-proves it on every run. Returns (genes, source).
    """
    if pinned_arg is not None:
        genes = pinned_arg.read_text(encoding="utf-8").strip().split()
        assert len(genes) == 17448, len(genes)
        return genes, f"explicit:{pinned_arg}"
    if PINNED_PATH.exists():
        genes = PINNED_PATH.read_text(encoding="utf-8").strip().split()
        assert len(genes) == 17448, len(genes)
        return genes, "canonical:pinned_genes_17448.txt"
    print("WARNING: canonical pinned file absent; deriving token-1 from "
          "entity_types.txt (deviation, see results JSON pinned_source)",
          flush=True)
    genes = [e for e, t in (ln.split() for ln in
             (exports_dir / unit / "entity_types.txt").read_text(
                 encoding="utf-8").splitlines()) if t == "1"]
    assert len(genes) == 17448, len(genes)
    return genes, f"derived:entity_types-token-1:{unit}"


def load_splits(exports_dir: Path, unit: str, disease: str):
    def tails(fname):
        out = []
        with (exports_dir / unit / fname).open(encoding="utf-8") as f:
            for line in f:
                h, _, t = line.rstrip("\n").split("\t")
                if h == disease:
                    out.append(t)
        return out
    return set(tails("train2.txt")), set(tails("valid.txt")), tails("test.txt")


def load_scores(scores_path: Path):
    dump = json.loads(scores_path.read_text(encoding="utf-8"))
    for key in ("scores", "pinned_scores", "score_vector"):
        if key in dump:
            skey = key
            break
    else:
        raise SystemExit(f"unknown score-dump keys: {sorted(dump.keys())}")
    assert dump.get("pinned_order") == "pinned_genes_17448.txt line order", \
        dump.get("pinned_order")
    scores = np.array(dump[skey], dtype=float)
    assert len(scores) == 17448, len(scores)
    return scores, dump.get("checkpoint"), skey


def load_loci(loci_path: Path):
    doc = json.loads(loci_path.read_text(encoding="utf-8"))
    return doc["loci"], doc


def train1_degree(exports_dir: Path, unit: str, universe: set[str]):
    """Undirected degree over train1.txt restricted to universe genes."""
    deg: dict[str, int] = defaultdict(int)
    with (exports_dir / unit / "train1.txt").open(encoding="utf-8") as f:
        for line in f:
            h, _, t = line.rstrip("\n").split("\t")
            if h in universe:
                deg[h] += 1
            if t in universe:
                deg[t] += 1
    return deg


def locus_gene_dist(loc: dict) -> dict:
    """Nearest Catalog SNP distance per gene (verbatim Probe B definition)."""
    dist = {}
    for a in loc["members"]:
        for g, roles in a["genes"].items():
            if "snp_gene" in roles and (g not in dist or dist[g] > 0):
                dist[g] = 0
        for col, dcol in (("upstream_id", "upstream_dist"),
                          ("downstream_id", "downstream_dist")):
            gid, dval = a.get(col), a.get(dcol)
            if (gid and gid.startswith("ensembl:") and dval is not None
                    and (gid not in dist or dval < dist[gid])):
                dist[gid] = dval
    return dist


def probe_b_rows(loci: list, scores: np.ndarray, pos_of: dict,
                 train_pos: set, deg: dict) -> list[dict]:
    """Verbatim Probe B row construction (token tiebreak, NOT pessimistic).

    Mirrors scripts/pod/m2c_prospective_model.py Probe B exactly so the
    regression gate compares like with like. The head-to-head uses
    pessimistic_ranks instead (see locus_arm_metrics).
    """
    b_rows = []
    for loc in loci:
        genes = loc["genes_pinned"]
        if len(genes) < 2:
            continue
        known = [g for g in genes if g in train_pos]
        if not known:
            continue
        by_model = sorted(genes, key=lambda g: (-scores[pos_of[g]], g))
        dist = locus_gene_dist(loc)
        with_d = [g for g in genes if g in dist]
        by_near = sorted(with_d, key=lambda g: (dist[g], g)) if with_d else []
        ranks_in = {}
        for i, g in enumerate(by_model):
            ranks_in[g] = i + 1
        top1 = by_model[0]
        by_deg = sorted(genes, key=lambda g: (-deg.get(g, 0), g))
        b_rows.append({
            "locus": loc["id"], "n_genes": len(genes), "n_known": len(known),
            "model_top1": top1, "model_top1_known": top1 in train_pos,
            "degree_top1": by_deg[0],
            "degree_top1_known": by_deg[0] in train_pos,
            "nearest_top1": by_near[0] if by_near else None,
            "nearest_top1_known": (by_near[0] in train_pos) if by_near else None,
            "known_ranks_model": sorted(ranks_in[g] for g in known),
            "chance_top1": len(known) / len(genes),
        })
    return b_rows


def cmd_gate(args) -> int:
    """Reproduce published per-locus model-only top-1 within GATE_TOL."""
    disease = args.unit.replace("__", ":")
    pinned, pinned_source = load_pinned(args.pinned, args.exports_dir, args.unit)
    pos_of = {g: i for i, g in enumerate(pinned)}
    scores, ckpt, skey = load_scores(args.scores)
    print(f"score vector: {skey}, ckpt={ckpt}", flush=True)
    loci, _doc = load_loci(args.loci_json)
    train_pos, _valid_pos, _test = load_splits(args.exports_dir, args.unit,
                                               disease)
    deg = train1_degree(args.exports_dir, args.unit, set(pinned))
    rows = probe_b_rows(loci, scores, pos_of, train_pos, deg)
    rate = float(np.mean([r["model_top1_known"] for r in rows])) if rows else 0.0
    expected = GATE_EXPECTED[args.unit]
    print(f"GATE {args.unit}: {len(rows)} loci, model top1-known {rate:.3f} "
          f"(expected {expected:.3f} ± {GATE_TOL}, pinned: {pinned_source})",
          flush=True)
    if abs(rate - expected) > GATE_TOL:
        print(f"GATE FAIL: discrepancy {rate:.3f} vs {expected:.3f} — STOP, "
              f"do not proceed (definitions diverged)", flush=True)
        return 2
    print("GATE PASS", flush=True)
    return 0


def ot_post(query: str, variables: dict, timeout_s: int = 90) -> dict:
    """Single GraphQL POST with retries (requests imported lazily)."""
    import requests
    last: Exception | None = None
    for attempt in range(4):
        try:
            resp = requests.post(OT_API, json={"query": query,
                                               "variables": variables},
                                 timeout=timeout_s)
            resp.raise_for_status()
            return resp.json()
        except Exception as exc:
            last = exc
            time.sleep(2 + attempt * 2)
    raise RuntimeError(f"OT GraphQL failed after retries: {last}")


CS_LIST_Q = ("query($sids:[String!]!,$idx:Int!){credibleSets(studyIds: $sids,"
             "page: {index: $idx, size: 100})"
             "{count rows{studyLocusId studyId studyType "
             "variant{id chromosome position}}}}")
CS_ONE_Q = ("query($slid:String!){credibleSet(studyLocusId: $slid)"
            "{studyLocusId studyId "
            "variant{id chromosome position} "
            "l2GPredictions{count rows{score target{id approvedSymbol}}}}}")


def cmd_fetch(args) -> int:
    """Pull OT L2G per (study, lead variant), cache raw JSON + manifest."""
    cache = Path(args.ot_cache_dir) / args.unit
    studies_dir = cache / "studies"
    sets_dir = cache / "credible_sets"
    studies_dir.mkdir(parents=True, exist_ok=True)
    sets_dir.mkdir(parents=True, exist_ok=True)
    loci, _doc = load_loci(args.loci_json)
    members = [(loc["id"], m) for loc in loci for m in loc["members"]]
    study_accs = sorted({m["study_acc"] for _, m in members if m.get("study_acc")})
    print(f"fetch {args.unit}: {len(members)} members, "
          f"{len(study_accs)} distinct studies", flush=True)
    if args.dry_run:
        print(f"DRY-RUN: would query {len(study_accs)} studies + matched "
              f"credible sets only; write under {cache}", flush=True)
        return 0
    # 1. enumerate credible sets per study (cached per study).
    study_sets: dict[str, list] = {}
    for acc in study_accs:
        fpath = studies_dir / f"{acc}.json"
        if fpath.exists():
            study_sets[acc] = json.loads(fpath.read_text(encoding="utf-8"))
            continue
        rows, idx, total = [], 0, None
        while True:
            d = ot_post(CS_LIST_Q, {"sids": [acc], "idx": idx})
            try:
                cs = d["data"]["credibleSets"]
            except (KeyError, TypeError) as err:
                raise RuntimeError(f"bad credibleSets response for {acc}: "
                                   f"{json.dumps(d)[:300]}") from err
            total = cs["count"]
            rows.extend(cs["rows"])
            if len(rows) >= total or not cs["rows"]:
                break
            idx += 1
        fpath.write_text(json.dumps(rows, indent=1, sort_keys=True) + "\n",
                         encoding="utf-8")
        study_sets[acc] = rows
        print(f"  study {acc}: {len(rows)} credible sets", flush=True)
        time.sleep(0.2)
    # 2. match members to credible sets by (study, chr, lead pos).
    index: dict[tuple, list] = defaultdict(list)
    for acc, rows in study_sets.items():
        for r in rows:
            v = r.get("variant") or {}
            index[(acc, str(v.get("chromosome")),
                   v.get("position"))].append(r["studyLocusId"])
    matched: dict[str, list] = defaultdict(list)
    unmatched: list[dict] = []
    for loc_id, m in members:
        key = (m.get("study_acc"), str(m.get("chr")), m.get("pos"))
        hits = index.get(key, [])
        if hits:
            matched[loc_id].extend(hits)
        else:
            unmatched.append({"locus": loc_id, "study": m.get("study_acc"),
                              "chr": m.get("chr"), "pos": m.get("pos"),
                              "reason": "no-positional-match"})
    # 3. pull l2G rows per matched studyLocusId (cache-first, resumable).
    wanted = sorted({s for v in matched.values() for s in v})
    have = sum((sets_dir / f"{s}.json").exists() for s in wanted)
    print(f"matched {len(wanted)} credible sets ({have} cached), "
          f"{len(unmatched)} unmatched members", flush=True)
    for s in wanted:
        fpath = sets_dir / f"{s}.json"
        if fpath.exists():
            continue
        d = ot_post(CS_ONE_Q, {"slid": s})
        try:
            cs = d["data"]["credibleSet"]
        except (KeyError, TypeError) as err:
            raise RuntimeError(f"bad credibleSet response for {s}: "
                               f"{json.dumps(d)[:300]}") from err
        fpath.write_text(json.dumps(cs, indent=1, sort_keys=True) + "\n",
                         encoding="utf-8")
        time.sleep(0.2)
    if args.verify_fetch:
        import random
        rng = random.Random(1337)
        for s in rng.sample(wanted, min(3, len(wanted))):
            d = ot_post(CS_ONE_Q, {"slid": s})
            fresh = json.dumps(d["data"]["credibleSet"], indent=1,
                               sort_keys=True) + "\n"
            if sha256_text(fresh) != sha256_file(sets_dir / f"{s}.json"):
                print(f"VERIFY FAIL: {s} differs on re-pull", flush=True)
                return 2
        print("VERIFY PASS: re-pull byte-identical on sample", flush=True)
    manifest = {"pulled_at": utcnow_iso(), "api": OT_API, "unit": args.unit,
                "studies": len(study_accs),
                "matched_loci": len(matched),
                "matched_sets": len(wanted),
                "unmatched_members": len(unmatched),
                "unmatched": unmatched,
                "matched": {k: sorted(set(v)) for k, v in matched.items()},
                "cache_sha256": {s: sha256_file(sets_dir / f"{s}.json")
                                 for s in wanted}}
    (cache / "manifest.json").write_text(json.dumps(manifest, indent=2,
                                                    sort_keys=True) + "\n",
                                         encoding="utf-8")
    print(f"fetch done: {len(wanted)} sets cached under {cache}", flush=True)
    return 0


def load_ot_gene_scores(cache: Path, fetch_manifest: dict) -> dict[str, dict]:
    """Per-locus OT gene scores = max L2G over matched credible sets."""
    sets_dir = cache / "credible_sets"
    per_locus: dict[str, dict] = {}
    for loc_id, sids in fetch_manifest["matched"].items():
        agg: dict[str, float] = {}
        for s in sids:
            cs = json.loads((sets_dir / f"{s}.json").read_text(encoding="utf-8"))
            for row in ((cs.get("l2GPredictions") or {}).get("rows") or []):
                tok = "ensembl:" + row["target"]["id"]
                if tok not in agg or row["score"] > agg[tok]:
                    agg[tok] = float(row["score"])
        per_locus[loc_id] = agg
    return per_locus


def grade_known(genes: list, grade_pos: set) -> set:
    """Known genes for grading = locus genes ∩ the selected split's
    positives (train2 for recovery, test for unseen-gene grading)."""
    return {g for g in genes if g in grade_pos}


def locus_arm_metrics(score_of: dict, genes: list, known: set) -> dict:
    """Pessimistic ranks + best-known top-1/top-5 for one arm, one locus."""
    sub = {g: score_of.get(g, float("-inf")) for g in genes}
    ranks = pessimistic_ranks(sub, genes)
    best = best_known_rank(ranks, known)
    return {"best_rank": best, "top1": hit_at_k(best, 1),
            "top5": hit_at_k(best, 5)}


def ot_arm_metrics(ot_sub: dict, genes: list, known: set) -> dict:
    """OT arm with the missing-data rule: no OT predictions for in-locus
    pinned genes -> miss (top-1/top-5 = 0.0, rank None, excluded from the
    mean-rank denominator downstream)."""
    if not ot_sub:
        return {"best_rank": None, "top1": 0.0, "top5": 0.0}
    return locus_arm_metrics(ot_sub, genes, known)


def summarize(rows: list[dict], key: str, n_boot: int, seed: int) -> dict:
    vals = np.array([r[key] for r in rows if r[key] is not None])
    return {"n": len(vals), "mean": round(float(np.mean(vals)), 4),
            "ci": [round(float(x), 4) for x in bootstrap_ci(vals, n_boot, seed)]} \
        if len(vals) else {"n": 0, "mean": None, "ci": [None, None]}


def verdict_for_diff(diff: float, lo: float, hi: float) -> str:
    if diff >= 0.05 and lo > 0:
        return "WIN"
    if diff <= -0.05 and hi < 0:
        return "LOSS"
    return "TIE"


def cmd_run(args) -> int:
    """Gate, then head-to-head; writes results JSON (+ report md)."""
    from sina.manifests import Run
    unit, disease = args.unit, args.unit.replace("__", ":")
    manifest_path = Path(args.manifest or
                         (REPO_ROOT / "reports" / "m2c" /
                          f"loci_headtohead_{args.unit}.manifest.json"))
    try:
        import subprocess
        git_ref = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO_ROOT,
                                 capture_output=True, text=True,
                                 timeout=30).stdout.strip()
    except Exception:
        git_ref = ""
    cfg = {}
    if args.config:
        cfg = json.loads(Path(args.config).read_text(encoding="utf-8"))
    run = Run(script="scripts/pod/m2c_loci_headtohead.py",
              manifest_path=manifest_path, config=cfg, git_ref=git_ref)
    pinned, pinned_source = load_pinned(args.pinned, args.exports_dir,
                                        args.unit)
    pos_of = {g: i for i, g in enumerate(pinned)}
    scores, ckpt, _skey = load_scores(args.scores)
    loci, _doc = load_loci(args.loci_json)
    train_pos, _valid_pos, _test_list = load_splits(args.exports_dir,
                                                     args.unit, disease)
    test_pos = set(_test_list)
    grade_pos = test_pos if args.known_split == "test" else train_pos
    # Stage 1: regression gate (structural — always runs first).
    if "gate" not in run.done_stages:
        with run.stage("gate"):
            deg_gate = train1_degree(args.exports_dir, args.unit, set(pinned))
            rows_gate = probe_b_rows(loci, scores, pos_of, train_pos, deg_gate)
            rate = float(np.mean([r["model_top1_known"] for r in rows_gate]))
            expected = GATE_EXPECTED[args.unit]
            print(f"GATE {args.unit}: {len(rows_gate)} loci, {rate:.3f} vs "
                  f"{expected:.3f}", flush=True)
            if abs(rate - expected) > GATE_TOL:
                raise SystemExit(f"GATE FAIL {rate:.3f} vs {expected:.3f}")
    # Stage 2: metrics (OT cache must exist — fetch is a separate subcommand).
    results_path = Path(args.out or (REPO_ROOT / "reports" / "m2c" /
                                     f"loci_headtohead_{args.unit}.json"))
    if "metrics" not in run.done_stages:
        with run.stage("metrics", outputs=[results_path]):
            cache = Path(args.ot_cache_dir) / args.unit
            fmanifest_path = cache / "manifest.json"
            if not fmanifest_path.exists():
                raise SystemExit(f"no OT cache at {cache} — run fetch first")
            fmanifest = json.loads(fmanifest_path.read_text(encoding="utf-8"))
            ot_gene = load_ot_gene_scores(cache, fmanifest)
            deg = train1_degree(args.exports_dir, args.unit, set(pinned))
            deg_score = {g: float(deg.get(g, 0)) for g in pinned}
            model_score = {g: float(scores[pos_of[g]]) for g in pinned}
            graded, ungradeable, noknown = [], [], []
            for loc in loci:
                genes = [g for g in loc["genes_pinned"] if g in pos_of]
                known = grade_known(genes, grade_pos)
                cls = classify_locus(len(genes), len(known))
                if cls != "gradeable":
                    (ungradeable if cls.startswith("ungradeable") else noknown
                     ).append({"locus": loc["id"], "n_genes": len(genes),
                               "n_known": len(known), "class": cls})
                    continue
                ot_sub = {g: v for g, v in ot_gene.get(loc["id"], {}).items()
                          if g in set(genes)}
                row = {"locus": loc["id"], "chr": loc.get("chr"),
                       "n_genes": len(genes), "n_known": len(known),
                       "year_min": loc.get("year_min"),
                       "year_max": loc.get("year_max"),
                       "ot_covered": bool(ot_sub),
                       "model": locus_arm_metrics(model_score, genes, known),
                       "ot": ot_arm_metrics(ot_sub, genes, known),
                       "degree": locus_arm_metrics(deg_score, genes, known)}
                row["chance_p"] = len(known) / len(genes)
                graded.append(row)
            if args.dry_run:
                print(f"DRY-RUN: {len(graded)} gradeable, "
                      f"{len(ungradeable)} ungradeable, {len(noknown)} no-known; "
                      f"no writes", flush=True)
                return 0
            chance_means = chance_draw_means([r["chance_p"] for r in graded],
                                             100, args.seed)
            out = {
                "unit": unit, "disease": disease, "generated_at": utcnow_iso(),
                "seed": args.seed, "bootstrap": args.bootstrap,
                "checkpoint": ckpt, "pinned_source": pinned_source,
                "git_ref": git_ref, "known_split": args.known_split,
                "n_loci_total": len(loci), "n_gradeable": len(graded),
                "ungradeable": ungradeable, "no_known": noknown,
                "chance_draws": {"n_draws": 100,
                                 "mean": round(float(np.mean(chance_means)), 4),
                                 "ci": [round(float(x), 4) for x in
                                        np.percentile(chance_means, [2.5, 97.5])]},
                "strata": {},
            }
            for sname, srows in (
                    ("all", graded),
                    ("new_loci", [r for r in graded
                                  if (r["year_min"] or 0) >= args.cutoff]),
                    ("established", [r for r in graded
                                     if (r["year_min"] or 0) < args.cutoff])):
                if not srows:
                    out["strata"][sname] = {"n": 0}
                    continue
                m_top1 = np.array([r["model"]["top1"] for r in srows])
                o_top1 = np.array([r["ot"]["top1"] for r in srows])
                diffs = m_top1 - o_top1
                dmean = float(np.mean(diffs))
                dlo, dhi = bootstrap_ci(diffs, args.bootstrap, args.seed)
                cov = [r for r in srows if r["ot_covered"]]
                entry = {"n": len(srows), "n_ot_covered": len(cov)}
                entry["model"] = {
                    "top1": summarize(
                        [{"x": r["model"]["top1"]} for r in srows], "x",
                        args.bootstrap, args.seed),
                    "top5": summarize(
                        [{"x": r["model"]["top5"]} for r in srows], "x",
                        args.bootstrap, args.seed + 1),
                    "mean_best_rank": summarize(
                        [{"x": float(r["model"]["best_rank"])} for r in srows],
                        "x", args.bootstrap, args.seed + 2),
                }
                entry["ot_primary"] = {
                    "top1": summarize(
                        [{"x": r["ot"]["top1"]} for r in srows], "x",
                        args.bootstrap, args.seed),
                    "top5": summarize(
                        [{"x": r["ot"]["top5"]} for r in srows], "x",
                        args.bootstrap, args.seed + 1),
                    "mean_best_rank": summarize(
                        [{"x": float(r["ot"]["best_rank"])}
                         for r in cov], "x", args.bootstrap, args.seed + 2),
                }
                entry["ot_covered_subset"] = {
                    "n": len(cov),
                    "top1": summarize(
                        [{"x": r["ot"]["top1"]} for r in cov], "x",
                        args.bootstrap, args.seed),
                    "top5": summarize(
                        [{"x": r["ot"]["top5"]} for r in cov], "x",
                        args.bootstrap, args.seed + 1),
                }
                entry["degree"] = {
                    "top1": summarize(
                        [{"x": r["degree"]["top1"]} for r in srows], "x",
                        args.bootstrap, args.seed),
                    "top5": summarize(
                        [{"x": r["degree"]["top5"]} for r in srows], "x",
                        args.bootstrap, args.seed + 1),
                }
                entry["diff_model_minus_ot_top1"] = {
                    "mean": round(dmean, 4), "ci": [round(dlo, 4), round(dhi, 4)],
                    "verdict": verdict_for_diff(dmean, dlo, dhi),
                }
                out["strata"][sname] = entry
            results_path.parent.mkdir(parents=True, exist_ok=True)
            results_path.write_text(json.dumps(out, indent=2, sort_keys=True)
                                    + "\n", encoding="utf-8")
    # Stage 3: report md (renders every unit results file matching THIS run's
    # naming — a test-knowns run must not render the full-panel numbers).
    report_path = Path(args.report or (REPO_ROOT / "reports" / "m2c" /
                                       "loci_headtohead.md"))
    sibling_results = [
        results_path.with_name(results_path.name.replace(args.unit, u))
        for u in ("mondo__0005180", "mondo__0004975")
    ]
    if "report" not in run.done_stages:
        with run.stage("report", outputs=[report_path]):
            render_report(report_path, sibling_results)
    v = json.loads(results_path.read_text(encoding="utf-8"))
    s = v["strata"]["all"]
    print(f"{args.unit}: model top1 {s['model']['top1']['mean']} vs OT "
          f"{s['ot_primary']['top1']['mean']} "
          f"(diff {s['diff_model_minus_ot_top1']['mean']} "
          f"{s['diff_model_minus_ot_top1']['ci']} "
          f"{s['diff_model_minus_ot_top1']['verdict']})", flush=True)
    return 0


def _fmt(m: dict, key: str = "mean", ci_key: str = "ci") -> str:
    if m.get("n", 1) == 0 or m.get("mean") is None:
        return "n/a (n=0)"
    return f"{m[key]:.4f} [{m[ci_key][0]:.4f},{m[ci_key][1]:.4f}] (n={m['n']})"


def render_report(report_path: Path,
                  results_paths: list[Path] | None = None) -> None:
    """Render the report from the given unit results files.

    Defaults to the full-panel files; callers running a non-default split
    (e.g. known_split=test) must pass their own results paths so the report
    can never silently render the wrong panel's numbers.
    """
    if results_paths is None:
        results_paths = [
            REPO_ROOT / "reports" / "m2c" / f"loci_headtohead_{u}.json"
            for u in ("mondo__0005180", "mondo__0004975")
        ]
    units = {}
    for p in results_paths:
        if p.exists():
            v = json.loads(p.read_text(encoding="utf-8"))
            units[v["unit"]] = v
    if not units:
        raise SystemExit("no unit results files to render")
    L = ["<!-- track_scope: benchmark -->",
         "# Locus-ranking head-to-head — SINA model vs Open Targets Genetics L2G",
         "",
         f"Generated {utcnow_iso()}. Units present: "
         + ", ".join(sorted(units)) + ".",
         "",
         "## Question and pre-registered decision rule",
         "",
         "On GWAS loci, restricted to the pinned 17,448 universe, which ranks "
         "the known disease gene higher within its locus — the saved model "
         "score vector or the OT Genetics L2G score? Locus and known-gene "
         "(train2 head positives) definitions are verbatim Probe B "
         "(`m2c_prospective_model.py`); ties pessimistic in every arm; "
         "identical known definition in every arm.",
         "",
         "WIN = model top-1 minus OT top-1 >= 0.05 with bootstrap 95% CI of "
         "the per-locus difference excluding 0. Within ±0.05 points or CI "
         "overlapping 0 = TIE. OT ahead with CI clear = LOSS. Per disease, "
         "no pooling. Primary stratum: all gradeable loci.",
         ""]
    for unit, v in sorted(units.items()):
        L += [f"## {unit} ({v['disease']})", "",
              f"Checkpoint: `{v.get('checkpoint')}`; pinned: "
              f"{v.get('pinned_source')}; known_split: "
              f"**{v.get('known_split', 'train')}**; loci total "
              f"{v['n_loci_total']}, "
              f"gradeable {v['n_gradeable']}; ungradeable "
              f"{len(v['ungradeable'])}; no-known {len(v['no_known'])}; "
              f"chance draws mean {v['chance_draws']['mean']} "
              f"{v['chance_draws']['ci']}.", ""]
        L += ["| stratum | n | model top-1 | OT top-1 (primary) | "
              "OT covered-subset top-1 | degree top-1 | diff (verdict) |",
              "|---|---|---|---|---|---|---|"]
        for sname in ("all", "new_loci", "established"):
            s = v["strata"].get(sname, {"n": 0})
            if not s.get("n"):
                L.append(f"| {sname} | 0 | n/a | n/a | n/a | n/a | n/a |")
                continue
            d = s["diff_model_minus_ot_top1"]
            L.append(f"| {sname} | {s['n']} | {_fmt(s['model']['top1'])} | "
                     f"{_fmt(s['ot_primary']['top1'])} | "
                     f"{_fmt(s['ot_covered_subset']['top1'])} "
                     f"(n={s['ot_covered_subset']['n']}) | "
                     f"{_fmt(s['degree']['top1'])} | {d['mean']:.4f} "
                     f"{d['ci']} {d['verdict']} |")
        L += ["",
              "| stratum | model top-5 | OT top-5 | degree top-5 | "
              "model mean-best-rank | OT mean-best-rank (covered n) |",
              "|---|---|---|---|---|---|"]
        for sname in ("all", "new_loci", "established"):
            s = v["strata"].get(sname, {"n": 0})
            if not s.get("n"):
                L.append(f"| {sname} | n/a | n/a | n/a | n/a | n/a |")
                continue
            L.append(f"| {sname} | {_fmt(s['model']['top5'])} | "
                     f"{_fmt(s['ot_primary']['top5'])} | "
                     f"{_fmt(s['degree']['top5'])} | "
                     f"{_fmt(s['model']['mean_best_rank'])} | "
                     f"{_fmt(s['ot_primary']['mean_best_rank'])} |")
        L += ["",
              f"Ungradeable loci ({len(v['ungradeable'])}): "
              + (", ".join(r["locus"] for r in v["ungradeable"]) or "none"),
              f"No-known loci ({len(v['no_known'])}): "
              + (", ".join(r["locus"] for r in v["no_known"]) or "none"), ""]
    L += ["## Limitations (mandatory)",
          "",
          "1. Known genes are train positives — recovery/ranking quality, "
          "not novel discovery.",
          "2. Model weights include a 0.64% DisGeNET-encumbered slice — "
          "internal go/no-go only.",
          "3. OT L2G and the model train on overlapping literature "
          "ecosystems — product-relevance comparison, not a clean "
          "scientific one.",
          "4. One seed (1337), frozen checkpoints, no retraining — CIs are "
          "over loci, not seeds.",
          ""]
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(L), encoding="utf-8")
    print(f"report -> {report_path}", flush=True)


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="mode", required=True)
    for name in ("gate", "fetch", "run"):
        p = sub.add_parser(name)
        p.add_argument("--unit", required=True,
                       choices=("mondo__0005180", "mondo__0004975",
                                "mondo__0005148", "mondo__0005090",
                                "mondo__0007254", "mondo__0005575",
                                "mondo__0008903", "mondo__0008315"))
        p.add_argument("--scores", type=Path, default=None)
        p.add_argument("--loci-json", type=Path, default=None)
        p.add_argument("--exports-dir", type=Path, default=DEFAULT_EXPORTS)
        p.add_argument("--pinned", type=Path, default=None)
        p.add_argument("--ot-cache-dir", type=Path, default=DEFAULT_OT_CACHE)
        p.add_argument("--cutoff", type=int, default=2020)
        p.add_argument("--known-split", choices=["train2", "test"],
                       default="train2")
        p.add_argument("--bootstrap", type=int, default=1000)
        p.add_argument("--seed", type=int, default=1337)
        p.add_argument("--out", type=Path, default=None)
        p.add_argument("--report", type=Path, default=None)
        p.add_argument("--manifest", type=Path, default=None)
        p.add_argument("--config", type=Path, default=None)
        p.add_argument("--dry-run", action="store_true")
        p.add_argument("--verify-fetch", action="store_true")
    return ap


def main() -> int:
    args = build_parser().parse_args()
    if args.scores is None:
        args.scores = (REPO_ROOT / "reports" / "m2c" / "ot_panels" /
                       "score_vectors" /
                       f"{args.unit}__epoch_001.scores.json")
    if args.loci_json is None:
        args.loci_json = (REPO_ROOT / "reports" / "m2c" /
                          f"gwas_loci_{args.unit}.json")
    if args.mode == "gate":
        return cmd_gate(args)
    if args.mode == "fetch":
        return cmd_fetch(args)
    return cmd_run(args)


if __name__ == "__main__":
    raise SystemExit(main())
