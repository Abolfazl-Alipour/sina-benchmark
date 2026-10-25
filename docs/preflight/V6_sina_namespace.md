# V6 — `sina` namespace disambiguation

**Status: DRAFT — awaiting owner signature**

The PyPI/import name **`sina`** is used by this project for its own package
(`src/sina/`, `pip install -e .`). It is a private, non-published distribution name
installed only inside the project conda env and CI; it is not uploaded to PyPI.

Disambiguation notes:

- `sina` is also the name of unrelated public projects (e.g. the Sina Weibo
  ecosystem, an old Python signals library). This project has no dependency on, and
  no code derived from, any of them.
- KG write-back edges (spec §11.4) use the `sina:` CURIE-style namespace prefix for
  SINA-generated evidence edges; that prefix is project-internal and disjoint from
  all upstream source prefixes (`ncbigene:`, `mondo:`, `drugbank:`, …).
- If the package is ever published, the distribution name will change to a
  collision-free name (e.g. `sina-platform`) before publication; the import package
  name may remain `sina`.

Owner signature: ______________________  Date: __________
