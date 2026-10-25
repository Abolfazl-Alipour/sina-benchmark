"""Long-run driver for the trimmed-graph learning run (CHANGES/0015/0016).

Same engine setup as third_party/biopathnet/script/run.py, but:
- validation scored EVERY epoch (run.py validates every num_epoch/10),
- checkpoint saved every epoch,
- early-stop after N consecutive epochs without validation MRR improvement
  (N = cfg.early_stop_patience, default 20 i.e. effectively off per
  CHANGES/0017 owner directive 2026-08-06; the run trains all 20 epochs),
- optional resume: set top-level `checkpoint:` in the config to continue
  training from a saved epoch,
- final: best-valid checkpoint reloaded and scored on full valid + test.

Usage (from third_party/biopathnet):
  conda run -n biopathnet --no-capture-output \
    python /opt/sina/scripts/pod/run_trimmed_long.py -s 1339 \
    -c /opt/sina/configs/m2b_trimmed_ad_long.yaml --gpus [0]
"""

import os
import pprint
import random
import sys

import numpy as np
import torch
from torchdrug import core
from torchdrug.utils import comm

sys.path.append("/opt/sina/third_party/biopathnet")
from biopathnet import (  # noqa: F401  (side-effect: class registration)
    dataset,
    layer,
    model,
    task,
    util,
)

# CHANGES/0017 (owner 2026-08-06): brake relaxed — train all epochs.
EARLY_STOP_PATIENCE_DEFAULT = 20  # override via cfg.early_stop_patience


def solver_load(solver, checkpoint, load_optimizer=True):
    if comm.get_rank() == 0:
        logger.warning(f"Load checkpoint from {checkpoint}")
    checkpoint = os.path.expanduser(checkpoint)
    state = torch.load(checkpoint, map_location=solver.device)
    # graphs cause issues on load; strip them (run.py does the same)
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
    if load_optimizer:
        solver.optimizer.load_state_dict(state["optimizer"])
        for state in solver.optimizer.state.values():
            for k, v in state.items():
                if isinstance(v, torch.Tensor):
                    state[k] = v.to(solver.device)
    comm.synchronize()


def train_and_validate(cfg, solver):
    if cfg.train.num_epoch == 0:
        return solver

    patience = int(cfg.get("early_stop_patience", EARLY_STOP_PATIENCE_DEFAULT))
    best_result = float("-inf")
    best_epoch = -1

    for _ in range(cfg.train.num_epoch):
        kwargs = cfg.train.copy()
        kwargs["num_epoch"] = 1
        solver.model.split = "train"
        solver.train(**kwargs)
        solver.save(f"model_epoch_{solver.epoch}.pth")
        solver.model.split = "valid"
        metric = solver.evaluate("valid")
        result = metric[cfg.metric]
        improved = result > best_result
        if improved:
            best_result = result
            best_epoch = solver.epoch
        logger.warning(
            "LONG-RUN epoch %d: valid %s = %.4f (best %.4f @ epoch %d)%s",
            solver.epoch,
            cfg.metric,
            result,
            best_result,
            best_epoch,
            "" if improved else " — no improvement",
        )
        if solver.epoch - best_epoch >= patience:
            logger.warning(
                "LONG-RUN early stop: %d consecutive epochs without validation "
                "improvement (best epoch %d, %.4f)",
                patience,
                best_epoch,
                best_result,
            )
            break

    solver_load(solver, f"model_epoch_{best_epoch}.pth")
    return solver


def test(cfg, solver):
    solver.model.split = "valid"
    solver.evaluate("valid")
    solver.model.split = "test"
    solver.evaluate("test")


if __name__ == "__main__":
    args, vars = util.parse_args()
    cfg = util.load_config(args.config, context=vars)
    working_dir = util.create_working_directory(cfg)
    args.seed = int(args.seed)
    seed_rank = args.seed + int(comm.get_rank())
    random.seed(args.seed)
    np.random.seed(args.seed)
    torch.manual_seed(seed_rank)

    logger = util.get_root_logger()
    logger.warning(f"Working directory: {working_dir}")
    logger.warning(f"Input Seed: {args.seed}")
    logger.warning(f"Set Seed: {seed_rank}")
    if comm.get_rank() == 0:
        logger.warning(f"Config file: {args.config}")
        logger.warning(pprint.pformat(cfg))

    dataset = core.Configurable.load_config_dict(cfg.dataset)
    solver = util.build_solver(cfg, dataset)

    if "checkpoint" in cfg:  # resume: continue training from a saved epoch
        solver_load(solver, cfg.checkpoint)

    train_and_validate(cfg, solver)
    test(cfg, solver)
