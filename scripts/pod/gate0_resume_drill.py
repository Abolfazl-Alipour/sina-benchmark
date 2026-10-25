"""Phase B Gate 0 Resume Drill on 1.43M BRG (CHANGES/0022, Directives v2 Phase B4).

Verifies checkpoint save, solver state reloading, weight fidelity, and
post-resume gradient step execution on GPU.
"""

import argparse
import hashlib
import json
import sys
from datetime import UTC, datetime
from pathlib import Path

import torch
from torch.utils import data as torch_data
from torchdrug import core
from torchdrug import data as torchdrug_data
from torchdrug import utils as torchdrug_utils

REPO_ROOT = Path(__file__).resolve().parent.parent.parent

BIOPATHNET_DIR = Path("/opt/sina/third_party/biopathnet")
if not BIOPATHNET_DIR.is_dir():
    BIOPATHNET_DIR = REPO_ROOT / "third_party" / "biopathnet"
sys.path.append(str(BIOPATHNET_DIR))
from biopathnet import dataset, layer, model, task, util  # noqa: E402, F401


def utcnow_iso() -> str:
    return datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def solver_load_state(solver, checkpoint_path: str) -> None:
    state = torch.load(checkpoint_path, map_location=solver.device)
    for key in [
        "fact_graph",
        "fact_graph_supervision",
        "graph",
        "train_graph",
        "valid_graph",
        "test_graph",
        "full_valid_graph",
        "full_test_graph",
    ]:
        state["model"].pop(key, 0)
    solver.model.load_state_dict(state["model"], strict=False)
    if "optimizer" in state:
        solver.optimizer.load_state_dict(state["optimizer"])
        for s in solver.optimizer.state.values():
            for k, v in s.items():
                if isinstance(v, torch.Tensor):
                    s[k] = v.to(solver.device)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-c", "--config", default="configs/m2b_gate0_ladder.yaml", type=Path)
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument(
        "--output", type=Path, default="reports/m2b/gate0_resume_drill_results.json"
    )
    args = parser.parse_args()

    started = utcnow_iso()
    cfg_path = REPO_ROOT / args.config
    output_path = REPO_ROOT / args.output
    output_path.parent.mkdir(parents=True, exist_ok=True)

    gpus = [0]
    cfg = util.load_config(str(cfg_path), context={"gpus": gpus})
    cfg.engine.batch_size = args.batch_size
    cfg.engine.gpus = gpus

    dataset_obj = core.Configurable.load_config_dict(cfg.dataset)
    solver1 = util.build_solver(cfg, dataset_obj)

    # 1. Step solver 1 for 1 batch
    sampler1 = torch_data.DistributedSampler(solver1.train_set, solver1.world_size, solver1.rank)
    loader1 = torchdrug_data.DataLoader(solver1.train_set, solver1.batch_size, sampler=sampler1)
    solver1.model.split = "train"
    solver1.model.train()

    batch = next(iter(loader1))
    batch = torchdrug_utils.cuda(batch, device=solver1.device)
    loss1, _ = solver1.model(batch)
    loss1.backward()
    solver1.optimizer.step()
    solver1.optimizer.zero_grad()

    # Save checkpoint
    ckpt_path = output_path.parent / "drill_checkpoint.pth"
    solver1.save(str(ckpt_path))
    ckpt_bytes = ckpt_path.stat().st_size

    # 2. Build fresh solver 2 and load checkpoint
    solver2 = util.build_solver(cfg, dataset_obj)
    solver_load_state(solver2, str(ckpt_path))
    solver2.model.split = "train"
    solver2.model.train()

    # 3. Verify weight fidelity
    params_match = True
    max_diff = 0.0
    for (_n1, p1), (_n2, p2) in zip(
        solver1.model.named_parameters(), solver2.model.named_parameters(), strict=True
    ):
        diff = (p1 - p2).abs().max().item()
        if diff > max_diff:
            max_diff = diff
        if diff > 1e-6:
            params_match = False

    # 4. Post-resume forward/backward step
    batch2 = next(iter(loader1))
    batch2 = torchdrug_utils.cuda(batch2, device=solver2.device)
    loss2, _ = solver2.model(batch2)
    loss2.backward()
    solver2.optimizer.step()
    solver2.optimizer.zero_grad()

    # Clean up checkpoint
    if ckpt_path.exists():
        ckpt_path.unlink()

    finished = utcnow_iso()
    results = {
        "drill": "gate0_resume_drill",
        "batch_size": args.batch_size,
        "checkpoint_bytes": ckpt_bytes,
        "params_match": params_match,
        "max_param_diff": max_diff,
        "pre_save_loss": float(loss1.item()),
        "post_resume_loss": float(loss2.item()),
        "status": "PASS" if params_match else "FAIL",
        "started_at": started,
        "finished_at": finished,
    }

    payload = json.dumps(results, indent=2) + "\n"
    output_path.write_text(payload, encoding="utf-8")
    sha256 = hashlib.sha256(payload.encode("utf-8")).hexdigest()

    manifest_path = output_path.with_suffix(".manifest.json")
    manifest = {
        "script": "gate0_resume_drill.py",
        "started_at": started,
        "finished_at": finished,
        "status": results["status"],
        "outputs": {str(output_path): sha256},
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print(f"Resume drill: {results['status']} (max diff {max_diff:.2e})")
    print(f"Results written to: {output_path}")
    return 0 if params_match else 1


if __name__ == "__main__":
    raise SystemExit(main())
