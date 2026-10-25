"""eval_exact_check — verify streaming-eval fidelity against exact global metrics.

Context: torchdrug's engine was patched (site-packages) to evaluate per
micro-batch with a triple-weighted running mean ("SINA streaming eval"),
because full-split accumulation OOMs container RAM. Per-triple means are
exact under that weighting ONLY if every batch's metric is a plain mean over
its triples. This script recomputes ranking metrics two ways from a saved
checkpoint on one LODO export:

  1. exact: collect per-triple rankings for the whole split, then global
     mr/mrr/hits@k (authoritative)
  2. streaming: the engine's per-batch-of-4 weighted mean (as patched)

and prints both side by side plus the delta. Rankings themselves are computed
exactly as the task does: rank = sum((pos <= pred) & mask) + 1 (pessimistic,
filtered, type-masked via the task's own target()/predict()).

Usage (biopathnet env, GPU):
  python scripts/pod/eval_exact_check.py --export-dir mondo__0010196_keep_disphe \
      --checkpoint reports/m2b/checkpoints_wave1/mondo__0010196_keep_pass1.pth \
      [--split test] [--eval-batch-size 4]
"""

import argparse
import sys
from pathlib import Path

import torch
from torchdrug import core
from torchdrug import utils as td_utils

REPO_ROOT = Path(__file__).resolve().parent.parent.parent

BIOPATHNET_DIR = Path("/opt/sina/third_party/biopathnet")
if not BIOPATHNET_DIR.is_dir():
    BIOPATHNET_DIR = REPO_ROOT / "third_party" / "biopathnet"
sys.path.append(str(BIOPATHNET_DIR))
from biopathnet import dataset, layer, model, task, util  # noqa: E402, F401

EXPORTS_DIR = Path("/opt/sina/data/processed/m2/lodo_exports")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--export-dir", required=True)
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--split", default="test", choices=["valid", "test"])
    parser.add_argument("--eval-batch-size", type=int, default=4)
    parser.add_argument("-c", "--config", default="configs/m2b_trimmed_ad_long.yaml", type=Path)
    args = parser.parse_args()

    cfg = util.load_config(str(REPO_ROOT / args.config), context={"gpus": [0]})
    cfg.dataset.path = str(EXPORTS_DIR / args.export_dir)
    cfg.engine.gpus = [0]
    cfg.engine.batch_size = args.eval_batch_size

    dataset_obj = core.Configurable.load_config_dict(cfg.dataset)
    solver = util.build_solver(cfg, dataset_obj)
    solver.load(str(args.checkpoint))
    solver.model.split = args.split
    solver.model.eval()

    test_set = getattr(solver, f"{args.split}_set")
    dataloader = torch.utils.data.DataLoader(test_set, args.eval_batch_size)

    all_rankings = []  # per-triple [n, 2]
    streaming_sum: dict[str, float] = {}
    streaming_count = 0
    with torch.no_grad():
        for batch in dataloader:
            batch = td_utils.cuda(batch, device=solver.device)
            pred, target = solver.model.predict_and_target(batch)
            mask, tgt = target
            pos_pred = pred.gather(-1, tgt.unsqueeze(-1).to(pred.device))
            ranking = (
                torch.sum((pos_pred <= pred) & mask.to(pred.device), dim=-1) + 1
            ).float()  # [bs, 2], exact per-triple
            all_rankings.append(ranking.cpu())

            bs = len(pred)
            batch_metric = solver.model.evaluate(pred, target)
            for k, v in batch_metric.items():
                if k in ("mr", "mrr") or k.startswith("hits@"):
                    streaming_sum[k] = streaming_sum.get(k, 0.0) + float(v) * bs
            streaming_count += bs

    ranks = torch.cat(all_rankings)  # [n, 2]
    exact = {
        "mr": float(ranks.mean()),
        "mrr": float((1.0 / ranks).mean()),
    }
    for k in (1, 3, 10, 20, 50, 100):
        exact[f"hits@{k}"] = float((ranks <= k).float().mean())

    streaming = {k: v / max(streaming_count, 1) for k, v in streaming_sum.items()}

    print(f"split={args.split} triples={ranks.shape[0]} batches of {args.eval_batch_size}")
    print(f"{'metric':<10} {'exact':>12} {'streaming':>12} {'delta':>12}")
    for k in ["mr", "mrr", "hits@1", "hits@3", "hits@10", "hits@20", "hits@50", "hits@100"]:
        e = exact.get(k, float("nan"))
        s = streaming.get(k, float("nan"))
        print(f"{k:<10} {e:>12.6f} {s:>12.6f} {e - s:>12.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
