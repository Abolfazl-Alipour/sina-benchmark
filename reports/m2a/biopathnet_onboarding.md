<!-- track_scope: benchmark -->
# M2B BioPathNet onboarding (spec §7.1)

**Status:** DONE — official repo pinned, environment built, mock run green on CPU.

## Pinned artifacts

- **Repo:** `github.com/emyyue/BioPathNet`, commit `b051a9003a1363b5655cdb7dfdad60609c0e4fe3` (main HEAD, 2026-07-25), cloned to `third_party/biopathnet` (git-ignored).
- **License:** MIT (verified, `third_party/biopathnet/LICENSE`, "Copyright (c) 2021 MilaGraph").
- **Paper:** *Enhancing link prediction in biomedical knowledge graphs with BioPathNet*, Nature Biomedical Engineering, DOI 10.1038/s41551-025-01598-z (peer-reviewed version of record, published Jan 2026).
- **Framework:** NBFNet-family GNN over TorchDrug; learns path representations; supports a background regulatory graph (BRG) for message passing — matches spec §7.4's BRG requirements.
- **TorchDrug dep:** `github.com/DeepGraphLearning/torchdrug`, commit `6066fbd82360abb5f270cba1eca560af01b8cc90`, cloned to `third_party/torchdrug`, installed editable.

## Environment (separate conda env `biopathnet`)

- Python 3.10; torch 2.0.1+cpu; torchvision 0.15.2; torchaudio 2.0.2.
- PyG CPU wheels for torch 2.0.1: `pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv`.
- `ogb easydict pyyaml numpy==1.26.4 ninja`.
- **Gotchas fixed during bring-up:**
  1. `setuptools>=81` removed `pkg_resources` (torchdrug needs it) → pinned `setuptools<81` (80.10.2).
  2. torchdrug imports RDKit → X11 libs needed → conda `xorg-libxrender xorg-libxext xorg-libx11`.
  3. `create_working_directory` in `biopathnet/util.py` does `os.chdir` into the experiment dir **before** loading the dataset → relative dataset paths break. Fix per repo README: absolute paths in config (`configs/m2b_mock_run.yaml`). No code patch needed.
  4. torchdrug JIT-compiles a C++ extension → needs `ninja` on PATH (run via `conda run -n biopathnet`, not the bare env python).

## Mock run (smallest real config, CPU)

Command:

```bash
cd third_party/biopathnet
conda run -n biopathnet --no-capture-output python script/run.py \
  -s 1024 -c /home/reza/projects/SinaProject/configs/m2b_mock_run.yaml --gpus null
```

Result (mock data, 5 epochs, NBFNet 6×32 dims): **test MRR 0.4358, hits@1 0.125, hits@3 0.708, AUROC 0.486** — training, evaluation, and checkpointing all work on CPU.

## What this unlocks / what it does not

- The full BioPathNet pipeline (BRG + supervision files, NBFNet training, eval) is proven runnable locally on small graphs.
- PrimeKG reproduction configs exist upstream (`config/primekg/*.yaml`, e.g. `alzheimers.yaml` with `num_epoch: 4`, `batch_size: 16`, `num_negative: 64`) — these become the Task-D reproduction baseline (§7.2).
- Real-graph runs (PrimeKG 8.1M edges, OptimusKG curated Task-G 11.7M edges) need cloud GPU (Phase 0.5); no VRAM assumptions made here (§7.1).
