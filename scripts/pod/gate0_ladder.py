"""Light Gate 0 VRAM/throughput ladder probe (CHANGES/0022).

Single-batch forward+backward at the pre-registered batch rungs
300→280→250→200→150→100 on one LODO unit export of the 1,431,945-edge BRG
(default: Werner `mondo_0010196`). Per rung records peak VRAM
(torch.cuda.max_memory_allocated, reset per rung), mean ms/batch over timed
steps after warmup, and extrapolated pass time
  pass_minutes = ms/batch / 1000 / 60 * (train2 rows / batch_size).
torch.cuda.OutOfMemoryError (or RuntimeError matching "out of memory") at a
rung is caught, the cache emptied, and the ladder continues at the next lower
rung; the below-100 fallback (gradient accumulation) is out of scope here.

Engine setup mirrors scripts/pod/run_trimmed_long.py (same third_party/
biopathnet imports, util.load_config with the `{{ gpus }}` template context,
core.Configurable.load_config_dict for the dataset, util.build_solver per
rung), except the repo root is derived from __file__ instead of a hardcoded
sys.path, so the script works from any checkout location and any cwd. Results
are written as JSON (+ a small sha256-manifest) to the config's
ladder.output_json, resolved relative to the repo root. Logging is JSON-lines
on stdout, mirroring sina.logging conventions without importing sina (the
biopathnet env does not have the sina package installed).

Usage (from any cwd; biopathnet conda env, py3.10 / torch 2.0.1+cu118):
  cd <repo>/third_party/biopathnet && conda run -n biopathnet --no-capture-output \
    python <repo>/scripts/pod/gate0_ladder.py \
    --config <repo>/configs/m2b_gate0_ladder.yaml [--dry-run] [--rungs 300,280,250]
"""

import argparse
import hashlib
import json
import re
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
BIOPATHNET_DIR = REPO_ROOT / "third_party" / "biopathnet"

DEFAULT_RUNGS = [300, 280, 250, 200, 150, 100]
DEFAULT_OUTPUT_JSON = "reports/m2b/gate0_ladder_results.json"
SEED = 1337  # plan 2026-08-10: single seed for Gate 0
REQUIRED_SECTIONS = ("dataset", "task", "optimizer", "engine", "ladder")


def log(event, **fields):
    """JSON-lines stdout logging, mirroring sina.logging's structlog shape."""
    record = {
        "event": event,
        "level": "info",
        "timestamp": datetime.now(timezone.utc).isoformat(),  # noqa: UP017
    }
    record.update(fields)
    print(json.dumps(record, default=str), flush=True)


def utcnow_iso():
    return datetime.now(timezone.utc).isoformat()  # noqa: UP017


def git_ref():
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
        return out.stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def resolve_repo_path(value):
    """Resolve a config path; relative paths anchor at the repo root."""
    path = Path(value).expanduser()
    return path if path.is_absolute() else (REPO_ROOT / path).resolve()


def parse_rungs(text):
    if not text:
        return None
    return [int(part.strip()) for part in text.split(",") if part.strip()]


def count_lines(path):
    with path.open("r", encoding="utf-8") as fin:
        return sum(1 for _ in fin)


def load_probe_config(config_path):
    """Parse the config torch-free (dry-run and pre-flight).

    The engine section templates `gpus: {{ gpus }}` from the ladder section,
    so first render with an empty jinja context (placeholder -> null) to read
    `ladder.gpus`; the real run re-renders via biopathnet util.load_config
    with that value, exactly like run_trimmed_long.py.
    """
    import yaml  # deferred: dry-run stays torch-free

    raw = config_path.read_text(encoding="utf-8")
    try:
        import jinja2

        rendered = jinja2.Template(raw).render()
    except ImportError:
        rendered = re.sub(r"\{\{.*?\}\}", "", raw)
    return yaml.safe_load(rendered)


def validate_config(cfg):
    missing = [section for section in REQUIRED_SECTIONS if section not in cfg]
    if missing:
        raise ValueError(f"config missing required sections: {', '.join(missing)}")
    if "path" not in cfg["dataset"]:
        raise ValueError("config missing dataset.path")


def is_oom(exc):
    import torch

    return isinstance(exc, torch.cuda.OutOfMemoryError) or (
        isinstance(exc, RuntimeError) and "out of memory" in str(exc).lower()
    )


def print_table(rows):
    print()
    print(f"{'batch':>6}  {'status':<7}  {'peak_vram_gb':>12}  {'ms/batch':>10}  {'pass_min':>9}")
    for row in rows:
        ms = f"{row['ms_per_batch']:.1f}" if "ms_per_batch" in row else "-"
        pm = f"{row['pass_minutes']:.1f}" if row.get("pass_minutes") is not None else "-"
        print(
            f"{row['batch_size']:>6}  {row['status']:<7}  {row['peak_vram_gb']:>12.2f}"
            f"  {ms:>10}  {pm:>9}"
        )


def run_ladder(config_path, cfg, rungs, warmup_batches, timed_batches, gpus, output_json):
    """GPU path: build the BioPathNet solver per rung and time train steps."""
    import random

    import numpy as np
    import torch
    from torch.utils import data as torch_data
    from torchdrug import core
    from torchdrug import data as torchdrug_data
    from torchdrug import utils as torchdrug_utils

    if not BIOPATHNET_DIR.is_dir():
        raise RuntimeError(f"BioPathNet checkout not found at {BIOPATHNET_DIR}")
    sys.path.append(str(BIOPATHNET_DIR))
    from biopathnet import dataset, layer, model, task, util  # noqa: F401 (class registration)

    if not torch.cuda.is_available():
        raise RuntimeError("CUDA device required for the ladder probe")
    torch.cuda.init()  # memory-stats calls require an initialized context (torch 2.0.1)

    random.seed(SEED)
    np.random.seed(SEED)
    torch.manual_seed(SEED)

    dataset_path = resolve_repo_path(cfg["dataset"]["path"])
    train2_path = dataset_path / "train2.txt"
    if not train2_path.is_file():
        raise RuntimeError(f"train2.txt not found at {train2_path} (LODO unit export missing?)")
    num_train2 = count_lines(train2_path)

    # Same construction as run_trimmed_long.py: render `{{ gpus }}` then load.
    cfg_ez = util.load_config(str(config_path), context={"gpus": gpus})
    cfg_ez.dataset.path = str(dataset_path)
    log("dataset_load", path=str(dataset_path), num_train2_triples=num_train2)
    dataset_obj = core.Configurable.load_config_dict(cfg_ez.dataset)

    device = gpus[0]  # int index for torch.cuda memory/synchronize calls
    rows = []
    for rung in rungs:
        log("rung_start", batch_size=rung)
        row = {"batch_size": rung, "status": "ok"}
        solver = None
        try:
            torch.cuda.empty_cache()
            torch.cuda.reset_peak_memory_stats(device)
            cfg_ez.engine.batch_size = rung
            solver = util.build_solver(cfg_ez, dataset_obj)
            # Replicates torchdrug Engine.train's inner loop (world_size 1).
            sampler = torch_data.DistributedSampler(
                solver.train_set, solver.world_size, solver.rank
            )
            dataloader = torchdrug_data.DataLoader(
                solver.train_set, solver.batch_size, sampler=sampler, num_workers=solver.num_worker
            )
            model = solver.model
            model.split = "train"
            model.train()
            timings = []
            for batch_id, batch in enumerate(dataloader):
                if batch_id >= warmup_batches + timed_batches:
                    break
                batch = torchdrug_utils.cuda(batch, device=solver.device)
                torch.cuda.synchronize(device)
                start = time.perf_counter()
                loss, _metric = model(batch)
                loss.backward()
                solver.optimizer.step()
                solver.optimizer.zero_grad()
                torch.cuda.synchronize(device)
                if batch_id >= warmup_batches:
                    timings.append((time.perf_counter() - start) * 1000.0)
            mean_ms = sum(timings) / len(timings)
            row["ms_per_batch"] = round(mean_ms, 1)
            row["per_batch_ms"] = [round(t, 1) for t in timings]
            row["pass_minutes"] = round(mean_ms / 1000.0 / 60.0 * (num_train2 / rung), 1)
        except Exception as exc:
            row["status"] = "oom" if is_oom(exc) else "error"
            if row["status"] == "error":
                row["error"] = f"{type(exc).__name__}: {exc}"
            log("rung_exception", batch_size=rung, status=row["status"], error=row.get("error", ""))
        finally:
            row["peak_vram_gb"] = round(torch.cuda.max_memory_allocated(device) / 1e9, 2)
            del solver
            torch.cuda.empty_cache()
            torch.cuda.reset_peak_memory_stats(device)
        log("rung_done", **{k: v for k, v in row.items() if k != "per_batch_ms"})
        rows.append(row)
    return rows, num_train2, dataset_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-c", "--config", required=True, type=Path)
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="validate config + print planned rungs, exit 0 without touching GPU",
    )
    parser.add_argument("--rungs", default=None, help="comma-separated batch-size override")
    args = parser.parse_args()

    try:
        cfg = load_probe_config(args.config)
        validate_config(cfg)
    except (ValueError, OSError) as exc:
        print(f"config error: {exc}", file=sys.stderr)
        return 1

    ladder = cfg["ladder"]
    rungs = parse_rungs(args.rungs) or [int(r) for r in ladder.get("rungs", DEFAULT_RUNGS)]
    warmup_batches = int(ladder.get("warmup_batches", 1))
    timed_batches = int(ladder.get("timed_batches", 3))
    gpus = [int(g) for g in ladder.get("gpus", [0])]
    output_json = resolve_repo_path(ladder.get("output_json", DEFAULT_OUTPUT_JSON))
    dataset_path = resolve_repo_path(cfg["dataset"]["path"])

    log(
        "gate0_ladder_config",
        config=str(args.config),
        rungs=rungs,
        warmup_batches=warmup_batches,
        timed_batches=timed_batches,
        gpus=gpus,
        output_json=str(output_json),
        dataset_path=str(dataset_path),
    )

    if args.dry_run:
        train2_path = dataset_path / "train2.txt"
        if train2_path.is_file():
            num_train2 = count_lines(train2_path)
            log("dry_run_train2", num_train2_triples=num_train2)
        else:
            log("dry_run_warning", warning=f"train2.txt not found at {train2_path}")
        print(f"planned rungs: {rungs} (warmup {warmup_batches}, timed {timed_batches})")
        print(f"output: {output_json}")
        return 0

    started = utcnow_iso()
    try:
        rows, num_train2, dataset_path = run_ladder(
            args.config, cfg, rungs, warmup_batches, timed_batches, gpus, output_json
        )
    except RuntimeError as exc:
        print(f"gate0_ladder: {exc}", file=sys.stderr)
        return 1

    import torch  # for the metadata below; already loaded by run_ladder

    results = {
        "script": "scripts/pod/gate0_ladder.py",
        "change": "CHANGES/0022",
        "git_ref": git_ref(),
        "timestamp": utcnow_iso(),
        "gpu_name": torch.cuda.get_device_name(gpus[0]),
        "torch_version": torch.__version__,
        "config": str(args.config),
        "dataset_path": str(dataset_path),
        "num_train2_triples": num_train2,
        "warmup_batches": warmup_batches,
        "timed_batches": timed_batches,
        "rungs": rows,
    }
    output_json.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(results, indent=2) + "\n"
    output_json.write_text(payload, encoding="utf-8")
    sha256 = hashlib.sha256(payload.encode("utf-8")).hexdigest()

    manifest_path = output_json.with_suffix(".manifest.json")
    manifest = {
        "run_id": datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"),  # noqa: UP017
        "script": "gate0_ladder",
        "config": {"config": str(args.config), "rungs": rungs, "gpus": gpus},
        "git_ref": results["git_ref"],
        "started_at": started,
        "finished_at": utcnow_iso(),
        "status": "done",
        "outputs": {str(output_json): sha256},
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print_table(rows)
    print(f"\nresults: {output_json}\nsha256:  {sha256}\nmanifest: {manifest_path}")
    log("gate0_ladder_done", output_json=str(output_json), sha256=sha256)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
