"""Minimal OBO 1.2 parser for MONDO (crosswalk + is_a closure).

Only what M2A needs: term id, name, is_a parents, xrefs, obsoletion flag.
Streaming — mondo.obo is ~53 MB / ~50k terms.
"""

from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path

import polars as pl


@dataclass
class OboTerm:
    id: str
    name: str = ""
    parents: list[str] = field(default_factory=list)
    xrefs: list[str] = field(default_factory=list)
    is_obsolete: bool = False


def parse_obo(path: Path | str) -> Iterator[OboTerm]:
    """Yield OboTerm for every [Term] stanza with an id."""
    path = Path(path)
    term: OboTerm | None = None
    in_term = False
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line == "[Term]":
                if term is not None:
                    yield term
                term = None
                in_term = True
                continue
            if line.startswith("["):
                if term is not None:
                    yield term
                term = None
                in_term = False
                continue
            if not in_term:
                continue
            if line.startswith("id: "):
                term = OboTerm(id=line[4:].strip())
            elif term is None:
                continue
            elif line.startswith("name: "):
                term.name = line[6:].strip()
            elif line.startswith("is_a: "):
                # "is_a: MONDO:0000001 {source="..."} ! disease" → id only
                term.parents.append(line[6:].split(" ")[0].strip())
            elif line.startswith("xref: "):
                # "xref: DOID:10652 {source="..."}"
                term.xrefs.append(line[6:].split(" ")[0].split(" {")[0].strip())
            elif line.startswith("is_obsolete: "):
                term.is_obsolete = line.split(": ", 1)[1].strip() == "true"
    if term is not None:
        yield term


def build_xref_map(terms: list[OboTerm], prefixes: set[str]) -> pl.DataFrame:
    """external CURIE → MONDO id for xrefs whose prefix is in `prefixes`.

    Obsolete terms are excluded. One external id can map to several MONDO ids;
    all rows are kept (downstream decides).
    """
    rows = []
    for t in terms:
        if t.is_obsolete:
            continue
        for x in t.xrefs:
            prefix = x.split(":", 1)[0] if ":" in x else ""
            if prefix in prefixes:
                rows.append((prefix, x, t.id))
    return pl.DataFrame(
        rows,
        schema=["external_prefix", "external_id", "mondo_id"],
        orient="row",
    )


def descendant_closure(terms: list[OboTerm], root: str) -> set[str]:
    """root + all transitive descendants via is_a edges (child → parent)."""
    children: dict[str, list[str]] = {}
    for t in terms:
        for p in t.parents:
            children.setdefault(p, []).append(t.id)
    seen: set[str] = set()
    stack = [root]
    while stack:
        node = stack.pop()
        if node in seen:
            continue
        seen.add(node)
        stack.extend(children.get(node, []))
    return seen
