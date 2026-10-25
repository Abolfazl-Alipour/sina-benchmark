# V1 — Licenses (recorded verbatim from primary sources, 2026-07-24)

**Status: DRAFT — awaiting owner signature**

## PrimeKG (Track-R) — Harvard Dataverse doi:10.7910/DVN/IXA7BM

License field of the dataset record, verbatim from the Dataverse API
(`GET /api/datasets/:persistentId/?persistentId=doi:10.7910/DVN/IXA7BM`):

> `"license":{"name":"CC0 1.0","uri":"http://creativecommons.org/publicdomain/zero/1.0","rightsIdentifier":"CC0-1.0","rightsIdentifierScheme":"SPDX","schemeUri":"https://spdx.org/licenses/","languageCode":"en"}`

No `termsOfUse` / `termsOfAccess` set. Dataset v2.1, released 2022-05-02.

⚠️ **Discrepancy with the spec:** §3.3-V1 and Appendix A state PrimeKG is
"CC BY-NC-SA 4.0". The primary source (Dataverse record) says **CC0 1.0**.
Recorded verbatim as found — see `CHANGES/0006_primekg_license_verbatim.md`.
**Working assumption until owner sign-off:** treat PrimeKG as
research/benchmark-only anyway, because integrated sources (notably
DrugBank-derived content) carry their own restrictions regardless of the
Dataverse-level license. The three-track architecture is unaffected.

## OptimusKG (Track-D) — Dataverse doi:10.7910/DVN/IYNGEV + GitHub mims-harvard/optimuskg

Dataset terms of use, verbatim from the Dataverse record (v2.0, RELEASED,
published 2026-04-14):

> "The OptimusKG codebase is released under the MIT License. OptimusKG integrates
> multiple primary data resources, each of which is subject to its own license and
> terms of use. These terms may impose restrictions on redistribution, commercial
> use, or downstream applications of the resulting knowledge graph or its subsets.
> Some resources provide data under academic or noncommercial licenses, while
> others may impose attribution or usage requirements. As a result, use of
> OptimusKG may be partially restricted depending on the specific data components
> included in a given instantiation. Users are responsible for reviewing and
> complying with the license and terms of use of each primary dataset, as
> specified by the original data providers. OptimusKG does not alter or override
> these source-specific licensing conditions."

Code license: **MIT** (GitHub repo LICENSE). PyPI client `optimuskg` 1.0.0: MIT.

## Consequence for architecture (unchanged from spec §0.0-B.15)

- Track-R (PrimeKG) and Track-D (OptimusKG): research/benchmark use only.
- Track-C (SINA-KG-Open): permissive-source rebuild via the MIT OptimusKG
  pipeline — required for any commercial use. See `configs/m1_trackc_sources.yaml`
  (M3) and the narrative-separation lint (spec §13.13).

Owner signature: ______________________  Date: __________
