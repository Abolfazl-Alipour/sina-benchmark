# Dependency compatibility matrix (spec §3 item 1)

**Status: DRAFT — living document; updated whenever `environment.yml` changes**

Target platform: Python 3.12.x, Linux x86-64. Machine: 12 threads / 49 GB
(see `CHANGES/0001_threads_12.md`). All facts below verified against live
publisher APIs on 2026-07-24; sources are the URLs in each row.

## 1. Installed now (conda env `sina`, Phase 0/1/2)

Pinned in `environment.yml`; full transitive export in `environment.lock.yml`.

| Package | Version | Source |
|---|---|---|
| python | 3.12.13 | conda-forge |
| pydantic | 2.13.4 | PyPI |
| pyyaml | 6.0.3 | PyPI |
| structlog | 26.1.0 | PyPI |
| polars | 1.43.0 | PyPI (verified latest, 2026-07-24) |
| duckdb | 1.5.5 | PyPI (verified latest) |
| pyarrow | 25.0.0 | PyPI (verified latest) |
| httpx | 0.28.1 | PyPI |
| tenacity | 9.1.4 | PyPI |
| psutil | 7.2.2 | PyPI |
| pytest | 9.1.1 | PyPI |
| ruff | 0.16.0 | PyPI |
| optimuskg | 1.0.0 | PyPI (M2) |
| torch | 2.13.0+cpu | pytorch CPU index (M2A-0; import-verified with PyKEEN) |
| pykeen | 1.11.1 | PyPI (M2A-0) |
| scipy | 1.18.0 | PyPI (M2A-0) |

## 2. Deferred — enter the env when their phase starts (CHANGES/0003)

Hashes are publisher-provided SHA-256 for the named artifact. Where the publisher
ships **no checksum**, the hash is computed at first download and recorded in
`data/raw/REGISTRY.json` / this matrix at that time.

| Package / tool | Pinned version | Phase | URL | SHA-256 (artifact) |
|---|---|---|---|---|
| ~~torch (CPU build)~~ → installed M2A-0 | 2.13.0+cpu | 2 | `https://download.pytorch.org/whl/cpu` | `4ca4a9394b0c771238a4f73590fdbbc4debad85ed0fa63d026ae1b085da7d6e2` (`torch-2.13.0+cpu-cp312-cp312-manylinux_2_28_x86_64.whl`) |
| ~~pykeen~~ → installed M2A-0 | 1.11.1 | 2 | `https://pypi.org/project/pykeen/` | `b6cc0ab012fc6ba4f37b58fb7081a409a5c1f04cbbb0e0e9c2256e911757d4ba` (wheel) |
| rdkit | 2026.3.4 | 4–5 | `https://pypi.org/project/rdkit/` | `bec605934ad781e91aa10d434ac986ee18e191ae9dfa65fcabb1784aadf33edb` (`rdkit-2026.3.4-cp312-cp312-manylinux_2_28_x86_64.whl`) |
| dimorphite_dl | 2.0.2 | 5 | `https://pypi.org/project/dimorphite_dl/` | `44125dfe2a5e1a37239f61914428d7f51d8fa61ec701a1337e71c293514133a8` (wheel) |
| openjdk | conda-forge latest LTS | 2 (AnyBURL) | conda-forge | conda-forge pin at install time |
| ~~optimuskg~~ → installed M2 | 1.0.0 | 1 | `https://pypi.org/project/optimuskg/` | `47f108371830c8d1f052b389df964eadec8fd2c9c9e7b3b74ef65464af4ff078` (wheel) |
| AnyBURL | 23-1x (bugfix 2025-11-12) | 2 | `https://web.informatik.uni-mannheim.de/AnyBURL/AnyBURL-23-1x.jar` | none published — record on download |
| AutoDock Vina | 1.2.7 | 5 | `https://github.com/ccsb-scripps/AutoDock-Vina/releases/tag/v1.2.7` | none published — record on download |
| smina | static build | 5 | `https://sourceforge.net/projects/smina/files/smina.static/download` | none published — record on download |
| Open Babel | 3.2.1 | 5 | conda-forge `openbabel` / `https://github.com/openbabel/openbabel/releases/tag/openbabel-3-2-1` | none published — conda-forge pin |
| fpocket | 4.2.3 | 5 | conda-forge `fpocket` / `https://github.com/Discngine/fpocket/releases/tag/4.2.3` | none published — conda-forge pin |

## 3. Notes

- **torch**: install the CPU build locally via
  `pip install torch==2.13.0+cpu --index-url https://download.pytorch.org/whl/cpu`;
  the plain PyPI wheel is the CUDA-bundled build. Cloud GPU images pin separately
  per §12 job YAML.
- **AnyBURL** real upstream is `github.com/OpenBioLink/AnyBURL` (no GitHub
  releases; jar from the Mannheim site).
- **smina** SourceForge static binary has no versioned URL or checksum; treat as a
  pinned-by-hash asset after first download.
- Docker image builds happen in CI (`docker` job), not locally
  (`CHANGES/0002_no_local_docker.md`).
