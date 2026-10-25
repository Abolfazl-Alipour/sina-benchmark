"""Phase D Wave 1 Main Benchmark Runner with Live Heartbeat (13 Primary Disease Units).

Features:
- Live Heartbeat file (/tmp/sina_heartbeat.json) updated at every unit and step.
- Resumable: skips already computed checkpoints in reports/m2b/checkpoints_wave1/.
- Unbuffered output and explicit error trap.
- Evaluates full valid and test set metrics with batch micro-batching.
"""

import argparse
import hashlib
import json
import os
import random
import sys
import time
import traceback
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

HEARTBEAT_PATH = Path("/tmp/sina_heartbeat.json")
ERROR_LOG_PATH = Path("/tmp/sina_last_error.log")

PRIMARY_UNITS = [
    {"unit": "mondo:0010196", "dir": "mondo__0010196_keep_disphe", "name": "Werner syndrome"},
    {"unit": "mondo:0008310", "dir": "mondo__0008310_keep_disphe", "name": "HGPS"},
    {
        "unit": "mondo:0011913",
        "dir": "mondo__0011913_keep_disphe",
        "name": "Mandibuloacral dysplasia type A",
    },
    {
        "unit": "mondo:0012853",
        "dir": "mondo__0012853_keep_disphe",
        "name": "Mandibuloacral dysplasia type B",
    },
    {
        "unit": "mondo:0014157",
        "dir": "mondo__0014157_keep_disphe",
        "name": "Mandibuloacral dysplasia",
    },
    {
        "unit": "mondo:0957576",
        "dir": "mondo__0957576_keep_disphe",
        "name": "Wiedemann-Rautenstrauch syndrome",
    },
    {"unit": "mondo:0008199", "dir": "mondo__0008199_keep_disphe", "name": "Parkinson disease 1"},
    {
        "unit": "mondo:0007669",
        "dir": "mondo__0007669_keep_disphe",
        "name": "Renal cysts and diabetes syndrome",
    },
    {
        "unit": "mondo:0010785",
        "dir": "mondo__0010785_keep_disphe",
        "name": "Maternally-inherited diabetes and deafness",
    },
    {
        "unit": "mondo:0012192",
        "dir": "mondo__0012192_keep_disphe",
        "name": ("Permanent neonatal diabetes mellitus with pancreatic and cerebellar agenesis"),
    },
    {
        "unit": "mondo:0012522",
        "dir": "mondo__0012522_keep_disphe",
        "name": "Transient neonatal diabetes mellitus 3",
    },
    {
        "unit": "mondo:0100164",
        "dir": "mondo__0100164_keep_disphe",
        "name": "Permanent neonatal diabetes mellitus",
    },
    {
        "unit": "mondo:0100165",
        "dir": "mondo__0100165_keep_disphe",
        "name": "Permanent neonatal diabetes mellitus 1",
    },
]

# Wave 2 secondary panel (16 units, CHANGES/0024 wave schedule): 5
# secondary-evidence + 11 secondary-low-context, 1 seed each, descriptive.
SECONDARY_UNITS = [
    {
        "unit": "mondo:0007088",
        "dir": "mondo__0007088_keep_disphe",
        "name": "Alzheimer disease type 1",
    },
    {
        "unit": "mondo:0011150",
        "dir": "mondo__0011150_keep_disphe",
        "name": "acroosteolysis-keloid-like lesions-premature aging syndrome",
    },
    {
        "unit": "mondo:0013523",
        "dir": "mondo__0013523_keep_disphe",
        "name": "Nestor-Guillermo progeria syndrome",
    },
    {
        "unit": "mondo:0014527",
        "dir": "mondo__0014527_keep_disphe",
        "name": "progeroid features-hepatocellular carcinoma predisposition",
    },
    {
        "unit": "orphanet:75496",
        "dir": "orphanet__75496_keep_disphe",
        "name": "Ehlers-Danlos syndrome progeroid type",
    },
    {
        "unit": "mondo:0000828",
        "dir": "mondo__0000828_keep_disphe",
        "name": "juvenile-onset Parkinson disease",
    },
    {
        "unit": "mondo:0017279",
        "dir": "mondo__0017279_keep_disphe",
        "name": "young-onset Parkinson disease",
    },
    {
        "unit": "mondo:0100087",
        "dir": "mondo__0100087_keep_disphe",
        "name": "familial Alzheimer disease",
    },
    {
        "unit": "mondo:0015140",
        "dir": "mondo__0015140_keep_disphe",
        "name": "early-onset autosomal dominant Alzheimer disease",
    },
    {
        "unit": "mondo:0018911",
        "dir": "mondo__0018911_keep_disphe",
        "name": "maturity-onset diabetes of the young",
    },
    {
        "unit": "mondo:0020525",
        "dir": "mondo__0020525_keep_disphe",
        "name": "transient neonatal diabetes mellitus",
    },
    {
        "unit": "mondo:0012480",
        "dir": "mondo__0012480_keep_disphe",
        "name": "diabetes mellitus transient neonatal 2",
    },
    {
        "unit": "mondo:0030089",
        "dir": "mondo__0030089_keep_disphe",
        "name": "diabetes mellitus permanent neonatal 4",
    },
    {"unit": "efo:0005422", "dir": "efo__0005422_keep_disphe", "name": "skin aging"},
    {
        "unit": "mondo:0018203",
        "dir": "mondo__0018203_keep_disphe",
        "name": "LMNA-related cardiocutaneous progeria syndrome",
    },
    {
        "unit": "mondo:0019321",
        "dir": "mondo__0019321_keep_disphe",
        "name": "atypical Werner syndrome",
    },
]

# Dual-shuffle negative controls (CHANGES/0022 B3, keep-graph variants):
# 3 probe units x {across_disease (popularity null), within_disease
# (specificity null)}. Same pinned config + seed 1337; expectation is MRR
# collapse to the degree null (~0.01). Dirs generated by m2_lodo_shuffles.py
# (configs/m2_lodo_shuffles_keep.yaml).
SHUFFLE_UNITS = [
    {
        "unit": "mondo:0010196",
        "dir": "mondo__0010196_keep_shuffle_across_disease",
        "name": "Werner syndrome (across-disease shuffle)",
        "mode": "across_disease",
    },
    {
        "unit": "mondo:0010196",
        "dir": "mondo__0010196_keep_shuffle_within_disease",
        "name": "Werner syndrome (within-disease shuffle)",
        "mode": "within_disease",
    },
    {
        "unit": "mondo:0008310",
        "dir": "mondo__0008310_keep_shuffle_across_disease",
        "name": "HGPS (across-disease shuffle)",
        "mode": "across_disease",
    },
    {
        "unit": "mondo:0008310",
        "dir": "mondo__0008310_keep_shuffle_within_disease",
        "name": "HGPS (within-disease shuffle)",
        "mode": "within_disease",
    },
    {
        "unit": "mondo:0007669",
        "dir": "mondo__0007669_keep_shuffle_across_disease",
        "name": "Renal cysts and diabetes syndrome (across-disease shuffle)",
        "mode": "across_disease",
    },
    {
        "unit": "mondo:0007669",
        "dir": "mondo__0007669_keep_shuffle_within_disease",
        "name": "Renal cysts and diabetes syndrome (within-disease shuffle)",
        "mode": "within_disease",
    },
]


def utcnow_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")  # noqa: UP017


def write_heartbeat(
    unit_idx: int,
    total_units: int,
    unit_curie: str,
    unit_name: str,
    status: str,
    extra: dict | None = None,
):
    try:
        vram_mb = 0
        if torch.cuda.is_available():
            vram_mb = int(torch.cuda.memory_allocated() / 1e6)
        data = {
            "pid": os.getpid(),
            "timestamp_utc": utcnow_iso(),
            "epoch_time": time.time(),
            "unit_idx": unit_idx,
            "total_units": total_units,
            "unit_curie": unit_curie,
            "unit_name": unit_name,
            "status": status,
            "gpu_vram_mb": vram_mb,
            "extra": extra or {},
        }
        HEARTBEAT_PATH.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    except Exception as e:
        print(f"[Heartbeat Error] Failed to write heartbeat: {e}", file=sys.stderr)


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


def main():
    if torch.cuda.is_available():
        torch.cuda.init()

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-c", "--config", default="configs/m2b_trimmed_ad_long.yaml", type=Path)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--seed", type=int, default=1337)
    parser.add_argument(
        "--unit-set",
        choices=["primary", "shuffle", "secondary"],
        default="primary",
        help=(
            "primary = Wave 1/Wave 2 multi-seed (13 primaries); shuffle = 6 dual-shuffle "
            "negative controls; secondary = Wave 2 secondary panel (16 units)"
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="results JSON (default depends on --unit-set and --seed)",
    )
    args = parser.parse_args()

    units = {
        "primary": PRIMARY_UNITS,
        "shuffle": SHUFFLE_UNITS,
        "secondary": SECONDARY_UNITS,
    }[args.unit_set]
    # Seed-scoped outputs/checkpoints: Wave 2 reruns primaries on seeds
    # 1339/5027 and must not collide with the seed-1337 Wave 1 artifacts.
    # The seed-1337 primary naming is legacy — sentinel resume depends on it.
    seed_tag = "" if args.seed == 1337 else f"_seed{args.seed}"
    default_output = {
        "primary": f"reports/m2b/wave1_benchmark_results{seed_tag}.json",
        "shuffle": f"reports/m2b/shuffle_control_results{seed_tag}.json",
        "secondary": f"reports/m2b/wave2_secondary_results{seed_tag}.json",
    }[args.unit_set]
    if args.output is None:
        args.output = Path(default_output)
    experiment = {
        "primary": "phase_d_wave1_primary_benchmark",
        "shuffle": "phase_d_shuffle_collapse_controls",
        "secondary": "phase_d_wave2_secondary_benchmark",
    }[args.unit_set] + seed_tag

    started = utcnow_iso()
    output_path = REPO_ROOT / args.output
    output_path.parent.mkdir(parents=True, exist_ok=True)
    checkpoints_dir = (
        output_path.parent
        / (
            {
                "primary": "checkpoints_wave1",
                "shuffle": "checkpoints_shuffle_controls",
                "secondary": "checkpoints_wave2_secondary",
            }[args.unit_set]
        )
    )
    checkpoints_dir.mkdir(parents=True, exist_ok=True)

    cfg_template = REPO_ROOT / args.config
    exports_dir = Path("/opt/sina/data/processed/m2/lodo_exports")

    print("===============================================================")
    print(f"=== {experiment} ({len(units)} units) =====")
    print(f"=== PID: {os.getpid()} | Heartbeat: {HEARTBEAT_PATH} ===")
    print("===============================================================\n", flush=True)

    # Load existing results if resuming
    results = []
    if output_path.exists():
        try:
            existing = json.loads(output_path.read_text(encoding="utf-8"))
            results = existing.get("results", [])
            print(f"[Resume] Loaded {len(results)} evaluated units from {output_path.name}")
        except Exception:
            results = []

    if args.unit_set == "primary":
        # Legacy Wave 1 results are keyed by curie only (no "dir" field)
        evaluated_keys = {r["curie"] for r in results}
    else:
        evaluated_keys = {r.get("dir") or r["curie"] for r in results}
    gpus = [0]
    device = gpus[0]

    try:
        for idx, u in enumerate(units, 1):
            curie = u["unit"]
            name = u["name"]
            d_dir = exports_dir / u["dir"]
            if args.unit_set == "primary":
                # seed-1337 naming is legacy — existing checkpoints + sentinel resume
                # depend on it; later seeds get seed-scoped names (Wave 2)
                ckpt_file = checkpoints_dir / f"{curie.replace(':', '__')}_keep{seed_tag}_pass1.pth"
            else:
                ckpt_file = checkpoints_dir / f"{u['dir']}{seed_tag}_pass1.pth"

            unit_key = curie if args.unit_set == "primary" else u["dir"]
            if unit_key in evaluated_keys:
                print(f"[{idx}/{len(units)}] {curie} ({name}) already evaluated, skipping.")
                write_heartbeat(idx, len(units), curie, name, "already_evaluated")
                continue

            print(f"\n[{idx}/{len(units)}] Starting {curie} ({name})...", flush=True)
            write_heartbeat(idx, len(units), curie, name, "starting_unit")

            random.seed(args.seed)
            np.random.seed(args.seed)
            torch.manual_seed(args.seed + int(comm.get_rank()))

            cfg = util.load_config(str(cfg_template), context={"gpus": gpus})
            cfg.dataset.path = str(d_dir)
            cfg.engine.gpus = gpus
            cfg.engine.batch_size = args.batch_size
            cfg.task.sample_weight = False
            cfg.train.num_epoch = 1

            torch.cuda.empty_cache()
            torch.cuda.reset_peak_memory_stats(device)

            dataset_obj = core.Configurable.load_config_dict(cfg.dataset)
            solver = util.build_solver(cfg, dataset_obj)

            if ckpt_file.exists() and ckpt_file.stat().st_size > 0:
                print(f"    [Resume] Loading existing checkpoint {ckpt_file.name}...", flush=True)
                write_heartbeat(idx, len(PRIMARY_UNITS), curie, name, "loading_checkpoint")
                solver.load(str(ckpt_file))
                train_sec = 0.0
                peak_vram = round(torch.cuda.max_memory_allocated(device) / 1e9, 2)
            else:
                write_heartbeat(idx, len(PRIMARY_UNITS), curie, name, "training_epoch")
                t0 = time.perf_counter()
                solver.model.split = "train"
                solver.train(num_epoch=1)
                train_sec = round(time.perf_counter() - t0, 1)
                peak_vram = round(torch.cuda.max_memory_allocated(device) / 1e9, 2)
                solver.save(str(ckpt_file))
                print(
                    f"    Saved checkpoint to {ckpt_file.name}"
                    f" (train time: {train_sec}s | peak: {peak_vram} GB)",
                    flush=True,
                )

            write_heartbeat(idx, len(PRIMARY_UNITS), curie, name, "evaluating_valid")
            val = evaluate_solver(solver, "valid")

            write_heartbeat(idx, len(PRIMARY_UNITS), curie, name, "evaluating_test")
            test = evaluate_solver(solver, "test")

            del solver, dataset_obj
            torch.cuda.empty_cache()

            val_mrr = float(val.get("mrr", 0.0))
            test_mrr = float(test.get("mrr", 0.0))
            test_h10 = float(test.get("hits@10", 0.0))
            test_h100 = float(test.get("hits@100", 0.0))
            test_auroc = float(test.get("auroc", 0.0))
            test_ap = float(test.get("ap", 0.0))

            print(
                f"    -> Val MRR: {val_mrr:.5f} | Test MRR: {test_mrr:.5f} | "
                f"Test Hits@10: {test_h10 * 100:.2f}% | Test AUROC: {test_auroc:.3f}"
                f" | Test AP: {test_ap:.3f}",
                flush=True,
            )

            unit_result = {
                "curie": curie,
                "name": name,
                "dir": u["dir"],
                "mode": u.get("mode", "primary"),
                "train_seconds": train_sec,
                "peak_vram_gb": peak_vram,
                "valid_mrr": val_mrr,
                "test_mrr": test_mrr,
                "test_hits10": test_h10,
                "test_hits100": test_h100,
                "test_auroc": test_auroc,
                "test_ap": test_ap,
                "valid_metrics": val,
                "test_metrics": test,
            }
            results.append(unit_result)
            evaluated_keys.add(unit_key)

            # Incremental save
            payload_data = {
                "experiment": experiment,
                "batch_size": args.batch_size,
                "seed": args.seed,
                "started_at": started,
                "updated_at": utcnow_iso(),
                "results": results,
            }
            payload = json.dumps(payload_data, indent=2) + "\n"
            output_path.write_text(payload, encoding="utf-8")
            write_heartbeat(idx, len(units), curie, name, "unit_completed", {"test_mrr": test_mrr})

        finished = utcnow_iso()
        payload_data["finished_at"] = finished
        payload = json.dumps(payload_data, indent=2) + "\n"
        output_path.write_text(payload, encoding="utf-8")
        sha256 = hashlib.sha256(payload.encode("utf-8")).hexdigest()

        manifest_path = output_path.with_suffix(".manifest.json")
        manifest = {
            "script": "run_wave1_benchmark.py",
            "unit_set": args.unit_set,
            "started_at": started,
            "finished_at": finished,
            "status": "done",
            "outputs": {str(output_path): sha256},
        }
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        write_heartbeat(len(units), len(units), "ALL", "Completed", "benchmark_finished")
        print(
            f"\n{experiment} complete! All {len(units)} units evaluated. Saved to {output_path}",
            flush=True,
        )
        return 0

    except Exception as e:
        err_msg = (
            f"Fatal exception in Wave 1 benchmark at {utcnow_iso()}:\n{traceback.format_exc()}"
        )
        print(err_msg, file=sys.stderr, flush=True)
        ERROR_LOG_PATH.write_text(err_msg, encoding="utf-8")
        write_heartbeat(-1, len(PRIMARY_UNITS), "ERROR", "Crash", "fatal_error", {"error": str(e)})
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
