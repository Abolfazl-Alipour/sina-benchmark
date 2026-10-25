#!/usr/bin/env python
"""compute_ppr — personalized PageRank from EFO trait nodes (CPU).

V3 network feature for the M4B control arm (owner-ordered amendment:
disease-anchored propagation, not anonymous degree). Personalization =
mapped EFO trait node; random walk on the undirected cleaned full graph
(full_exports train1+train2, DISGENET-free); power iteration to tol.
Output: per (trait, panel-gene) PPR scores.

Usage (sina env, CPU, ~1h for 120 traits, resumable by trait):
  ~/miniconda3/envs/sina/bin/python scripts/m4b/compute_ppr.py
"""

import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
FULL = REPO_ROOT / "data" / "processed" / "m2c" / "full_exports"
M4B_DIR = REPO_ROOT / "data" / "processed" / "m4b"
SRC_TSV = Path("/tmp/m4b_src/genetics-gold-standards/gold_standards"
               "/processed/gwas_gold_standards.191108.tsv")
OUT = M4B_DIR / "ppr_scores.parquet"
ALPHA = 0.15
TOL = 1e-6
MAX_ITERS = 100


def utcnow_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def load_graph():
    names = [ln.split("\t")[0] for ln in
             (FULL / "entity_names.txt").read_text(
                 encoding="utf-8").splitlines()]
    idx = {e: i for i, e in enumerate(names)}
    rows, cols = [], []
    for name in ("train1.txt", "train2.txt"):
        for chunk in pd.read_csv(FULL / name, sep="\t", header=None,
                                 chunksize=500_000, dtype=str):
            for _, r in chunk.iterrows():
                a, b = idx.get(r[0]), idx.get(r[2])
                if a is not None and b is not None and a != b:
                    rows.append(a)
                    cols.append(b)
    n = len(names)
    adj = sparse.coo_matrix((np.ones(len(rows)), (rows, cols)),
                            shape=(n, n)).tocsr()
    adj = ((adj + adj.T) > 0).astype(np.float64)
    deg = np.asarray(adj.sum(axis=1)).ravel()
    deg[deg == 0] = 1.0
    walk = sparse.diags(1.0 / deg).dot(adj).tocsr()
    return names, idx, walk


def mapped_traits(idx: dict) -> dict:
    g = pd.read_csv(SRC_TSV, sep="\t", low_memory=False,
                    usecols=["trait_info.ontology"])
    efos = set()
    for x in g["trait_info.ontology"].dropna().unique():
        parts = str(x).split("_")
        if len(parts) == 2 and parts[1]:
            efos.add("efo:" + parts[1].lower())
    return {e: idx[e] for e in sorted(efos) if e in idx}


def panel_genes(idx: dict) -> set:
    feat = pd.read_parquet(M4B_DIR / "features" / "features.parquet",
                           columns=["ensg"])
    return {"ensembl:" + e for e in feat["ensg"].unique()
            if "ensembl:" + e in idx}


def ppr(walk, seed: int, n: int) -> np.ndarray:
    e = np.zeros(n)
    e[seed] = 1.0
    p = e.copy()
    for _ in range(MAX_ITERS):
        new = ALPHA * e + (1.0 - ALPHA) * (walk.T.dot(p))
        if np.abs(new - p).sum() < TOL:
            p = new
            break
        p = new
    return p


def main() -> int:
    from sina.manifests import Run
    run = Run(script="scripts/m4b/compute_ppr.py",
              manifest_path=M4B_DIR / "ppr.manifest.json",
              config={}, git_ref="")
    with run.stage("ppr"):
        names, idx, walk = load_graph()
        print(f"graph: {walk.shape[0]} nodes, {walk.nnz} directed edges",
              flush=True)
        traits = mapped_traits(idx)
        print(f"mapped traits: {len(traits)}", flush=True)
        genes = panel_genes(idx)
        gidx = sorted(idx[e] for e in genes)
        gnames = [names[i].split(":")[1] for i in gidx]
        have: set[str] = set()
        if OUT.exists():
            old = pd.read_parquet(OUT, columns=["trait"])
            have = set(old["trait"].unique().tolist())
            print(f"resuming: {len(have)} traits done", flush=True)
        todo = [t for t in traits if t not in have]
        frames = []
        if have:
            frames.append(pd.read_parquet(OUT))
        for k, t in enumerate(todo):
            vec = ppr(walk, traits[t], walk.shape[0])
            frames.append(pd.DataFrame({"trait": t, "ensg": gnames,
                                        "ppr": vec[gidx]}))
            print(f"  [{k + 1}/{len(todo)}] {t}", flush=True)
        pd.concat(frames, ignore_index=True).to_parquet(OUT, index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
