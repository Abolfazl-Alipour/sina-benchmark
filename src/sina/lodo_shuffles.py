"""Degree-preserving supervision shuffles — B3 dual negative controls.

Directives v2 Phase B3 + CHANGES/0022 Light Gate 0 permutation control: rewrite
the DIS-GEN supervision (train2.txt) of a LODO export unit under one of two
degree-preserving nulls, keeping every other export file byte-identical:

- ``across_disease`` (popularity null): each GENE keeps its exact supervision
  degree; its edges are reassigned to uniformly sampled diseases. Gene
  popularity is kept, disease pairing destroyed.
- ``within_disease`` (specificity null): each DISEASE keeps its exact
  supervision degree; its edges are reassigned to uniformly sampled genes of
  the same unit's train2. Disease degree is kept, gene identity destroyed.

Leakage contract (hard guards, mirroring the A2 tuple-disjointness guard of
scripts/m2_build_lodo_exports.py): no duplicate (disease, relation, gene)
tuples in the shuffled output, and the shuffled tuples stay disjoint from the
unit's valid.txt/test.txt tuples. Collisions are resampled; after
``max_attempts`` failures for one edge the shuffle raises — it never silently
keeps a colliding tuple.

All functions are deterministic: same triples + same seed → same output, and
output rows are returned sorted for byte-stable file writes.
"""

import random
from collections import Counter
from collections.abc import Sequence
from pathlib import Path

Triple = tuple[str, str, str]  # (src_id, relation, dst_id) = (disease, relation, gene)

MAX_ATTEMPTS = 1000

MODES = ("across_disease", "within_disease")


class ShuffleError(RuntimeError):
    """Raised when no collision-free assignment exists within max_attempts."""


def curie_safe(curie: str) -> str:
    """Filesystem-safe unit directory name ('mondo:0010196' -> 'mondo__0010196')."""
    return curie.replace(":", "__")


def read_triples(path: Path) -> list[Triple]:
    """Read a headerless TSV triple file (src_id, relation, dst_id)."""
    rows: list[Triple] = []
    with Path(path).open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line:
                src, rel, dst = line.split("\t")
                rows.append((src, rel, dst))
    return rows


def write_triples(triples: Sequence[Triple], path: Path) -> None:
    """Write triples as sorted, headerless TSV lines (byte-stable output)."""
    with Path(path).open("w", encoding="utf-8") as fh:
        for t in sorted(triples):
            fh.write("\t".join(t) + "\n")


def gene_marginal(triples: Sequence[Triple]) -> dict[str, int]:
    """Supervision degree per gene (dst_id)."""
    return dict(Counter(t[2] for t in triples))


def disease_marginal(triples: Sequence[Triple]) -> dict[str, int]:
    """Supervision degree per disease (src_id)."""
    return dict(Counter(t[0] for t in triples))


def across_disease_shuffle(
    triples: Sequence[Triple],
    disease_pool: Sequence[str],
    forbidden: set[Triple],
    seed: int | str,
    max_attempts: int = MAX_ATTEMPTS,
) -> list[Triple]:
    """Popularity null: preserve each gene's exact degree, resample its diseases.

    Every edge ``(d, rel, g)`` keeps its relation and gene; the disease is drawn
    uniformly from ``disease_pool``. Draws colliding with ``forbidden``
    (valid/test tuples) or with an already-emitted tuple are resampled, up to
    ``max_attempts`` per edge. Genes and their edges are visited in sorted order
    so RNG consumption is input-order independent.
    """
    pool = sorted(set(disease_pool))
    if not pool:
        raise ShuffleError("across_disease: empty disease pool")
    rng = random.Random(seed)
    by_gene: dict[str, list[Triple]] = {}
    for t in triples:
        by_gene.setdefault(t[2], []).append(t)
    used: set[Triple] = set()
    out: list[Triple] = []
    for gene in sorted(by_gene):
        for _, rel, _ in sorted(by_gene[gene]):
            for _ in range(max_attempts):
                cand = (rng.choice(pool), rel, gene)
                if cand in forbidden or cand in used:
                    continue
                used.add(cand)
                out.append(cand)
                break
            else:
                raise ShuffleError(
                    f"across_disease: no collision-free disease for gene {gene!r} "
                    f"(rel {rel!r}) after {max_attempts} attempts"
                )
    return sorted(out)


def within_disease_shuffle(
    triples: Sequence[Triple],
    forbidden: set[Triple],
    seed: int | str,
    max_attempts: int = MAX_ATTEMPTS,
) -> list[Triple]:
    """Specificity null: preserve each disease's exact degree, resample its genes.

    Every edge ``(d, rel, g)`` keeps its disease and relation; the gene is drawn
    uniformly from the genes present in ``triples``. Same collision/resample
    contract as :func:`across_disease_shuffle`.
    """
    genes = sorted({t[2] for t in triples})
    if not genes:
        raise ShuffleError("within_disease: no genes in train2")
    rng = random.Random(seed)
    by_disease: dict[str, list[Triple]] = {}
    for t in triples:
        by_disease.setdefault(t[0], []).append(t)
    used: set[Triple] = set()
    out: list[Triple] = []
    for disease in sorted(by_disease):
        for _, rel, _ in sorted(by_disease[disease]):
            for _ in range(max_attempts):
                cand = (disease, rel, rng.choice(genes))
                if cand in forbidden or cand in used:
                    continue
                used.add(cand)
                out.append(cand)
                break
            else:
                raise ShuffleError(
                    f"within_disease: no collision-free gene for disease {disease!r} "
                    f"(rel {rel!r}) after {max_attempts} attempts"
                )
    return sorted(out)


def verify_shuffle(
    before: Sequence[Triple],
    after: Sequence[Triple],
    mode: str,
    forbidden: set[Triple],
) -> dict:
    """Post-shuffle verification row; raises AssertionError on any failed check.

    The preserved marginal depends on the mode: gene-degree for
    ``across_disease``, disease-degree for ``within_disease``. Edge count,
    tuple uniqueness and valid/test disjointness are always checked.
    """
    if mode not in MODES:
        raise ValueError(f"unknown shuffle mode {mode!r}")
    if mode == "across_disease":
        preserved = gene_marginal(before) == gene_marginal(after)
        marginal = "gene_degree"
    else:
        preserved = disease_marginal(before) == disease_marginal(after)
        marginal = "disease_degree"
    n_duplicates = len(after) - len(set(after))
    disjoint = set(after).isdisjoint(forbidden)
    row = {
        "mode": mode,
        "marginal_checked": marginal,
        "edges_before": len(before),
        "edges_after": len(after),
        "genes": len(gene_marginal(after)),
        "diseases": len(disease_marginal(after)),
        "degree_marginal_preserved": preserved,
        "n_duplicate_tuples": n_duplicates,
        "valid_test_disjoint": disjoint,
    }
    assert len(after) == len(before), f"{mode}: edge count {len(after)} != {len(before)}"
    assert preserved, f"{mode}: {marginal} marginal not preserved"
    assert n_duplicates == 0, f"{mode}: {n_duplicates} duplicate tuples"
    assert disjoint, f"{mode}: shuffled train2 collides with valid/test tuples"
    return row
