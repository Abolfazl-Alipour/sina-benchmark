"""Unit tests for sina.obo (fixture OBO text only)."""

import polars as pl

from sina.obo import build_xref_map, descendant_closure, parse_obo

FIXTURE = """format-version: 1.2
data-version: releases/test

[Term]
id: MONDO:0000001
name: disease

[Term]
id: MONDO:0004975
name: Alzheimer disease
is_a: MONDO:0000001 ! disease
xref: DOID:10652 {source="MONDO:equivalentTo"}
xref: EFO:0000249 {source="MONDO:equivalentTo"}

[Term]
id: MONDO:1000001
name: familial Alzheimer disease
is_a: MONDO:0004975 ! Alzheimer disease
xref: Orphanet:1020

[Term]
id: MONDO:0005180
name: Parkinson disease
is_a: MONDO:0000001 ! disease
xref: DOID:14330

[Term]
id: MONDO:9999999
name: obsolete thing
is_obsolete: true
xref: DOID:99999
"""


def _terms(tmp_path):
    p = tmp_path / "mini.obo"
    p.write_text(FIXTURE, encoding="utf-8")
    return list(parse_obo(p))


def test_parse_obo_fields(tmp_path):
    terms = {t.id: t for t in _terms(tmp_path)}
    assert terms["MONDO:0004975"].name == "Alzheimer disease"
    assert terms["MONDO:0004975"].parents == ["MONDO:0000001"]
    assert "DOID:10652" in terms["MONDO:0004975"].xrefs
    assert terms["MONDO:9999999"].is_obsolete


def test_xref_map_respects_prefixes_and_obsoletes(tmp_path):
    xw = build_xref_map(_terms(tmp_path), prefixes={"DOID", "EFO"})
    assert isinstance(xw, pl.DataFrame)
    pairs = set(zip(xw["external_id"].to_list(), xw["mondo_id"].to_list(), strict=True))
    assert ("DOID:10652", "MONDO:0004975") in pairs
    assert ("EFO:0000249", "MONDO:0004975") in pairs
    assert ("DOID:99999", "MONDO:9999999") not in pairs  # obsolete excluded
    assert not any(x.startswith("Orphanet") for x in xw["external_id"])  # prefix filtered


def test_descendant_closure(tmp_path):
    terms = _terms(tmp_path)
    ad = descendant_closure(terms, "MONDO:0004975")
    assert ad == {"MONDO:0004975", "MONDO:1000001"}
    pd = descendant_closure(terms, "MONDO:0005180")
    assert pd == {"MONDO:0005180"}
