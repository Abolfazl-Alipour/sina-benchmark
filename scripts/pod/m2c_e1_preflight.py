"""m2c_e1_preflight — Arm E E1 pre-flight / batch-ladder probe on the ANA-CNS completion exports.

Runs N training batches + 1 validation pass of the Arm E base config
(configs/m2c_arm_e.yaml, CHANGES/0026 frozen hyperparams) on one completion
unit, recording peak VRAM (torch.cuda.max_memory_allocated), wall-clock and
ms/batch. Used for (a) the sm_120 bring-up pre-flight (50 batches + 1 val,
v6 E1) and (b) the batch ladder (v6 §4, ceiling 91.2 GB on the 96 GB
RTX PRO 6000).

NOTE: this drives the stock 64-negative BioPathNet path. Message passing
(the dominant VRAM term) is identical to the full-vocab/multi-task Arm E
variant; tail scoring adds <1 GB. Final ladder numbers are re-confirmed with
the full-stack runner before E2.

Usage:
    python scripts/pod/m2c_e1_preflight.py --unit mondo__0005180 \
        --batch-size 32 --batches 50 [--skip-eval] [--output reports/m2c/e1_preflight_g7e.json]
"""

import argparse
import hashlib
import json
import random
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import torch
from torchdrug import core
from torchdrug.utils import comm

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.append(str(REPO_ROOT / "third_party" / "biopathnet"))
from biopathnet import dataset, layer, model, task, util  # noqa: E402, F401

EXPORTS_DIR = REPO_ROOT / "data/processed/m2c/completion_exports_ana"


def utcnow_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")  # noqa: UP017


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-c", "--config", type=Path, default=REPO_ROOT / "configs/m2c_arm_e.yaml")
    parser.add_argument("--unit", default="mondo__0005180",
                        choices=["mondo__0005180", "mondo__0004975"])
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--batches", type=int, default=50)
    parser.add_argument("--seed", type=int, default=1337)
    parser.add_argument("--eval-batch-size", type=int, default=4)
    parser.add_argument("--skip-eval", action="store_true")
    parser.add_argument("--output", type=Path,
                        default=REPO_ROOT / "reports/m2c/e1_preflight_g7e.json")
    args = parser.parse_args()

    if torch.cuda.is_available():
        torch.cuda.init()

    gpus = [0]
    device = gpus[0]
    started = utcnow_iso()

    random.seed(args.seed)
    np.random.seed(args.seed)
    torch.manual_seed(args.seed + int(comm.get_rank()))

    cfg = util.load_config(str(args.config), context={"gpus": gpus})
    cfg.dataset.path = str(EXPORTS_DIR / args.unit)
    cfg.engine.gpus = gpus
    cfg.engine.batch_size = args.batch_size
    cfg.task.sample_weight = False
    cfg.train.num_epoch = 1
    cfg.train.batch_per_epoch = args.batches

    torch.cuda.empty_cache()
    torch.cuda.reset_peak_memory_stats(device)

    t_load = time.perf_counter()
    dataset_obj = core.Configurable.load_config_dict(cfg.dataset)
    load_sec = round(time.perf_counter() - t_load, 1)
    print(f"dataset loaded: {dataset_obj.num_entity} entities, "
          f"{dataset_obj.num_relation} relations ({load_sec}s)", flush=True)

    solver = util.build_solver(cfg, dataset_obj)

    solver.model.split = "train"
    t0 = time.perf_counter()
    solver.train(num_epoch=1, batch_per_epoch=args.batches)
    train_sec = time.perf_counter() - t0
    train_peak_gb = round(torch.cuda.max_memory_allocated(device) / 1e9, 2)
    ms_per_batch = round(train_sec / args.batches * 1000, 1)
    print(f"train: {args.batches} batches in {train_sec:.1f}s "
          f"({ms_per_batch} ms/batch), peak VRAM {train_peak_gb} GB", flush=True)

    valid_metrics = {}
    eval_sec = 0.0
    if not args.skip_eval:
        solver.model.split = "valid"
        orig_bs = solver.batch_size
        solver.batch_size = args.eval_batch_size
        torch.cuda.empty_cache()
        t0 = time.perf_counter()
        metric = solver.evaluate("valid")
        eval_sec = round(time.perf_counter() - t0, 1)
        solver.batch_size = orig_bs
        for k, v in metric.items():
            if isinstance(v, torch.Tensor):
                valid_metrics[k] = (
                    float(v.item()) if v.numel() == 1 else [float(x) for x in v.flatten()]
                )
            elif isinstance(v, (int, float)):
                valid_metrics[k] = float(v)
        print(f"valid: mrr={valid_metrics.get('mrr', 0):.5f} "
              f"hits@10={valid_metrics.get('hits@10', 0):.4f} ({eval_sec}s)", flush=True)

    peak_gb = round(torch.cuda.max_memory_allocated(device) / 1e9, 2)
    result = {
        "script": "m2c_e1_preflight",
        "unit": args.unit,
        "batch_size": args.batch_size,
        "batches": args.batches,
        "seed": args.seed,
        "node": "g7e.2xlarge RTX PRO 6000 Blackwell 96GB sm_120",
        "env": "bpn120 torch 2.7.1+cu128 (owner-authorized divergence, CHANGES/0026)",
        "dataset_load_sec": load_sec,
        "train_seconds": round(train_sec, 1),
        "ms_per_batch": ms_per_batch,
        "train_peak_vram_gb": train_peak_gb,
        "valid_metrics": valid_metrics,
        "valid_seconds": eval_sec,
        "peak_vram_gb": peak_gb,
        "started_at": started,
        "finished_at": utcnow_iso(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(result, indent=2) + "\n"
    args.output.write_text(payload, encoding="utf-8")
    sha = hashlib.sha256(payload.encode("utf-8")).hexdigest()
    manifest = {"script": "m2c_e1_preflight", "status": "done",
                "outputs": {str(args.output): sha}}
    args.output.with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=2) + "\n",
                                                         encoding="utf-8")
    print(f"saved {args.output} (peak VRAM {peak_gb} GB)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
