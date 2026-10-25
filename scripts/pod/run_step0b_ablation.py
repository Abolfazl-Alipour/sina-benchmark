"""Phase C Step 0b GPU Ablation Runner (SINA v5.1 §7, Directives v2 Phase C).

Runs the combined Step 0b ablation matrix across 3 probe units:
  1. Werner syndrome (mondo:0010196)
  2. HGPS (mondo:0008310)
  3. Renal cysts and diabetes syndrome (mondo:0007669)

Evaluates:
  - C1: DIS-PHE {drop (1,420,095) vs keep (1,431,945)} x pass-1 x seed 1337
  - C2: Weighting {as-is vs capped} x {1 vs 2 passes} x seed 1337

Features:
  - Resumable: skips training if checkpoint already exists.
  - Evaluation micro-batching (batch_size=4 during eval) to prevent OOM
    during full entity candidate ranking on test sets.
  - Combined 2-epoch evaluation (scoring pass 1 & pass 2 in same run).
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

BIOPATHNET_DIR = Path("/opt/sina/third_party/biopathnet")
if not BIOPATHNET_DIR.is_dir():
    BIOPATHNET_DIR = REPO_ROOT / "third_party" / "biopathnet"
sys.path.append(str(BIOPATHNET_DIR))
from biopathnet import dataset, layer, model, task, util  # noqa: E402, F401

PROBE_UNITS = [
    {"unit": "mondo:0010196", "dir": "mondo__0010196", "name": "Werner syndrome"},
    {"unit": "mondo:0008310", "dir": "mondo__0008310", "name": "HGPS"},
    {
        "unit": "mondo:0007669",
        "dir": "mondo__0007669",
        "name": "Renal cysts and diabetes syndrome",
    },
]


def utcnow_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")  # noqa: UP017


def evaluate_solver(solver, split: str, eval_batch_size: int = 4) -> dict:
    solver.model.split = split
    orig_batch_size = solver.batch_size
    solver.batch_size = eval_batch_size
    torch.cuda.empty_cache()
    metric = solver.evaluate(split)
    solver.batch_size = orig_batch_size
    torch.cuda.empty_cache()

    cleaned = {}
    for k, v in metric.items():
        if isinstance(v, torch.Tensor):
            cleaned[k] = float(v.item()) if v.numel() == 1 else [float(x) for x in v.flatten()]
        elif isinstance(v, (int, float)):
            cleaned[k] = float(v)
    return cleaned


def run_arm_combined_2epochs(
    config_template_path: Path,
    dataset_dir: Path,
    batch_size: int,
    seed: int,
    gpus: list[int],
    checkpoints_dir: Path,
    unit_curie: str,
) -> dict:
    if torch.cuda.is_available():
        torch.cuda.init()

    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed + int(comm.get_rank()))

    cfg = util.load_config(str(config_template_path), context={"gpus": gpus})
    cfg.dataset.path = str(dataset_dir)
    cfg.engine.batch_size = batch_size
    cfg.engine.gpus = gpus
    cfg.task.sample_weight = False
    cfg.train.num_epoch = 1

    device = gpus[0]
    torch.cuda.empty_cache()
    torch.cuda.reset_peak_memory_stats(device)

    dataset_obj = core.Configurable.load_config_dict(cfg.dataset)
    solver = util.build_solver(cfg, dataset_obj)

    ckpt1 = checkpoints_dir / f"{unit_curie.replace(':', '__')}_drop_pass1.pth"
    ckpt2 = checkpoints_dir / f"{unit_curie.replace(':', '__')}_drop_pass2.pth"

    # --- Epoch 1 ---
    if ckpt1.exists():
        print(f"    [Resume] Loading existing checkpoint {ckpt1.name}...")
        solver.load(str(ckpt1))
        ep1_sec = 0.0
        peak1 = round(torch.cuda.max_memory_allocated(device) / 1e9, 2)
    else:
        t0 = time.perf_counter()
        solver.model.split = "train"
        solver.train(num_epoch=1)
        ep1_sec = round(time.perf_counter() - t0, 1)
        peak1 = round(torch.cuda.max_memory_allocated(device) / 1e9, 2)
        solver.save(str(ckpt1))

    val1 = evaluate_solver(solver, "valid")
    test1 = evaluate_solver(solver, "test")

    # --- Epoch 2 ---
    if ckpt2.exists():
        print(f"    [Resume] Loading existing checkpoint {ckpt2.name}...")
        solver.load(str(ckpt2))
        ep2_sec = 0.0
        peak2 = round(torch.cuda.max_memory_allocated(device) / 1e9, 2)
    else:
        t1 = time.perf_counter()
        solver.model.split = "train"
        solver.train(num_epoch=1)
        ep2_sec = round(time.perf_counter() - t1, 1)
        peak2 = round(torch.cuda.max_memory_allocated(device) / 1e9, 2)
        solver.save(str(ckpt2))

    val2 = evaluate_solver(solver, "valid")
    test2 = evaluate_solver(solver, "test")

    del solver, dataset_obj
    torch.cuda.empty_cache()

    return {
        "pass1": {
            "train_seconds": ep1_sec,
            "peak_vram_gb": peak1,
            "valid_mrr": val1.get("mrr"),
            "valid_hits10": val1.get("hits@10"),
            "valid_metrics": val1,
            "test_mrr": test1.get("mrr"),
            "test_hits10": test1.get("hits@10"),
            "test_metrics": test1,
        },
        "pass2": {
            "train_seconds": ep2_sec,
            "peak_vram_gb": peak2,
            "valid_mrr": val2.get("mrr"),
            "valid_hits10": val2.get("hits@10"),
            "valid_metrics": val2,
            "test_mrr": test2.get("mrr"),
            "test_hits10": test2.get("hits@10"),
            "test_metrics": test2,
        },
    }


def run_arm_1epoch(
    config_template_path: Path,
    dataset_dir: Path,
    sample_weight: bool,
    batch_size: int,
    seed: int,
    gpus: list[int],
    checkpoints_dir: Path,
    ckpt_name: str,
) -> dict:
    if torch.cuda.is_available():
        torch.cuda.init()

    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed + int(comm.get_rank()))

    cfg = util.load_config(str(config_template_path), context={"gpus": gpus})
    cfg.dataset.path = str(dataset_dir)
    cfg.engine.batch_size = batch_size
    cfg.engine.gpus = gpus
    cfg.task.sample_weight = sample_weight
    cfg.train.num_epoch = 1

    device = gpus[0]
    torch.cuda.empty_cache()
    torch.cuda.reset_peak_memory_stats(device)

    dataset_obj = core.Configurable.load_config_dict(cfg.dataset)
    solver = util.build_solver(cfg, dataset_obj)

    ckpt_path = checkpoints_dir / ckpt_name
    if ckpt_path.exists():
        print(f"    [Resume] Loading existing checkpoint {ckpt_name}...")
        solver.load(str(ckpt_path))
        train_sec = 0.0
        peak_vram = round(torch.cuda.max_memory_allocated(device) / 1e9, 2)
    else:
        t0 = time.perf_counter()
        solver.model.split = "train"
        solver.train(num_epoch=1)
        train_sec = round(time.perf_counter() - t0, 1)
        peak_vram = round(torch.cuda.max_memory_allocated(device) / 1e9, 2)
        solver.save(str(ckpt_path))

    val = evaluate_solver(solver, "valid")
    test = evaluate_solver(solver, "test")

    del solver, dataset_obj
    torch.cuda.empty_cache()

    return {
        "train_seconds": train_sec,
        "peak_vram_gb": peak_vram,
        "valid_mrr": val.get("mrr"),
        "valid_hits10": val.get("hits@10"),
        "valid_metrics": val,
        "test_mrr": test.get("mrr"),
        "test_hits10": test.get("hits@10"),
        "test_metrics": test,
    }


def main() -> int:
    if torch.cuda.is_available():
        torch.cuda.init()

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-c", "--config", default="configs/m2b_trimmed_ad_long.yaml", type=Path)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--seed", type=int, default=1337)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--output", type=Path, default="reports/m2b/step0b_ablation_results.json")
    args = parser.parse_args()

    started = utcnow_iso()
    output_path = REPO_ROOT / args.output
    output_path.parent.mkdir(parents=True, exist_ok=True)
    checkpoints_dir = output_path.parent / "checkpoints_step0b"
    checkpoints_dir.mkdir(parents=True, exist_ok=True)

    cfg_template = REPO_ROOT / args.config
    exports_dir = Path("/opt/sina/data/processed/m2/lodo_exports")

    print(
        f"Phase C Step 0b (Combined & Resumable): {len(PROBE_UNITS)} probe units x 3 runs = "
        f"{len(PROBE_UNITS) * 3} runs (yielding all 4 arms)"
    )

    if args.dry_run:
        for u in PROBE_UNITS:
            u_unit = u["unit"]
            u_dir = u["dir"]
            d_drop = exports_dir / u_dir
            d_keep = exports_dir / f"{u_dir}_keep_disphe"
            print(f"  [DRY-RUN] {u_unit}:")
            print(f"    - drop (2 epochs): {d_drop} (exists: {d_drop.exists()})")
            print(f"    - keep (1 epoch):  {d_keep} (exists: {d_keep.exists()})")
            print(f"    - capped (1 epoch): {d_drop} (exists: {d_drop.exists()})")
        return 0

    all_results = []
    for u in PROBE_UNITS:
        u_unit = u["unit"]
        u_name = u["name"]
        u_dir = u["dir"]
        curie_prefix = u_unit.replace(":", "__")
        d_drop = exports_dir / u_dir
        d_keep = exports_dir / f"{u_dir}_keep_disphe"

        print("\n==================================================")
        print(f"=== Probe Unit: {u_unit} ({u_name}) ===")
        print("==================================================")
        unit_res = {"unit": u_unit, "name": u_name, "arms": {}}

        # 1. Combined DROP 2 epochs (scores both pass 1 and pass 2)
        print("\n[1/3] Running DROP unweighted 2 epochs (evaluating pass 1 & pass 2)...")
        res_drop_2ep = run_arm_combined_2epochs(
            config_template_path=cfg_template,
            dataset_dir=d_drop,
            batch_size=args.batch_size,
            seed=args.seed,
            gpus=[0],
            checkpoints_dir=checkpoints_dir,
            unit_curie=u_unit,
        )
        p1 = res_drop_2ep["pass1"]
        p2 = res_drop_2ep["pass2"]
        p1_val_mrr = f"{p1['valid_mrr']:.4f}" if p1["valid_mrr"] is not None else "N/A"
        p1_test_mrr = f"{p1['test_mrr']:.4f}" if p1["test_mrr"] is not None else "N/A"
        p2_val_mrr = f"{p2['valid_mrr']:.4f}" if p2["valid_mrr"] is not None else "N/A"
        p2_test_mrr = f"{p2['test_mrr']:.4f}" if p2["test_mrr"] is not None else "N/A"
        print(
            f"  -> Pass 1: {p1['train_seconds']}s | peak {p1['peak_vram_gb']} GB | "
            f"valid MRR: {p1_val_mrr} | test MRR: {p1_test_mrr}"
        )
        print(
            f"  -> Pass 2: {p2['train_seconds']}s | peak {p2['peak_vram_gb']} GB | "
            f"valid MRR: {p2_val_mrr} | test MRR: {p2_test_mrr}"
        )
        unit_res["arms"]["drop_pass1_unweighted"] = {
            "disphe": "drop",
            "epochs": 1,
            "sample_weight": False,
            **p1,
        }
        unit_res["arms"]["drop_pass2_unweighted"] = {
            "disphe": "drop",
            "epochs": 2,
            "sample_weight": False,
            **p2,
        }

        # 2. KEEP DIS-PHE 1 epoch
        print("\n[2/3] Running KEEP DIS-PHE 1 epoch (unweighted)...")
        res_keep = run_arm_1epoch(
            config_template_path=cfg_template,
            dataset_dir=d_keep,
            sample_weight=False,
            batch_size=args.batch_size,
            seed=args.seed,
            gpus=[0],
            checkpoints_dir=checkpoints_dir,
            ckpt_name=f"{curie_prefix}_keep_pass1.pth",
        )
        k_val_mrr = f"{res_keep['valid_mrr']:.4f}" if res_keep["valid_mrr"] is not None else "N/A"
        k_test_mrr = f"{res_keep['test_mrr']:.4f}" if res_keep["test_mrr"] is not None else "N/A"
        print(
            f"  -> Keep Pass 1: {res_keep['train_seconds']}s | "
            f"peak {res_keep['peak_vram_gb']} GB | "
            f"valid MRR: {k_val_mrr} | "
            f"test MRR: {k_test_mrr}"
        )
        unit_res["arms"]["keep_pass1_unweighted"] = {
            "disphe": "keep",
            "epochs": 1,
            "sample_weight": False,
            **res_keep,
        }

        # 3. DROP CAPPED 1 epoch
        print("\n[3/3] Running DROP CAPPED 1 epoch (sample_weight=True)...")
        res_capped = run_arm_1epoch(
            config_template_path=cfg_template,
            dataset_dir=d_drop,
            sample_weight=True,
            batch_size=args.batch_size,
            seed=args.seed,
            gpus=[0],
            checkpoints_dir=checkpoints_dir,
            ckpt_name=f"{curie_prefix}_drop_capped_pass1.pth",
        )
        c_val_mrr = (
            f"{res_capped['valid_mrr']:.4f}" if res_capped["valid_mrr"] is not None else "N/A"
        )
        c_test_mrr = (
            f"{res_capped['test_mrr']:.4f}" if res_capped["test_mrr"] is not None else "N/A"
        )
        print(
            f"  -> Capped Pass 1: {res_capped['train_seconds']}s | "
            f"peak {res_capped['peak_vram_gb']} GB | "
            f"valid MRR: {c_val_mrr} | "
            f"test MRR: {c_test_mrr}"
        )
        unit_res["arms"]["drop_pass1_capped"] = {
            "disphe": "drop",
            "epochs": 1,
            "sample_weight": True,
            **res_capped,
        }

        all_results.append(unit_res)

    finished = utcnow_iso()
    payload_data = {
        "experiment": "phase_c_step0b_ablation_combined",
        "batch_size": args.batch_size,
        "seed": args.seed,
        "started_at": started,
        "finished_at": finished,
        "results": all_results,
    }

    payload = json.dumps(payload_data, indent=2) + "\n"
    output_path.write_text(payload, encoding="utf-8")
    sha256 = hashlib.sha256(payload.encode("utf-8")).hexdigest()

    manifest_path = output_path.with_suffix(".manifest.json")
    manifest = {
        "script": "run_step0b_ablation.py",
        "started_at": started,
        "finished_at": finished,
        "status": "done",
        "outputs": {str(output_path): sha256},
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print(f"\nStep 0b ablation complete. Results saved to {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
