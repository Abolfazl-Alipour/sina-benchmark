#!/usr/bin/env python
# ruff: noqa
"""m2c_arm_e_train — Arm E full-stack runner (C1-C5, v6 §8).

Custom training loop calling NBFNet.forward directly (NOT torchdrug Engine).
Per HANDOFF §4.2:

- Reuse KnowledgeGraphCompletionBiomed.preprocess for fact graph, but CUSTOM loop.
- DIS-GEN: multi-hot BCE over pinned 17,448, per-row edge+inverse removal via edge-ID mask.
- Aux (C4): 1:1 interleave GEN→BPO / GEN→PWY from train1, reoriented, multi-hot, 1.0 weight.
- Val (C3): filtered 1-hot BCE + filtered MRR/Hits@10 over pinned 17,448.
- Early stop: val loss patience 5, min_delta 3% relative, max 20 epochs.

Usage (bpn120 env, GPU):
    P=~/miniconda3/envs/bpn120; OV=~/persist/cuda-12.8-overlay
    PATH="$OV/hostbin:$P/bin:$PATH" CUDA_HOME=$P CXX="$P/bin/x86_64-conda-linux-gnu-g++" \
    TORCH_CUDA_ARCH_LIST="12.0" $P/bin/python -u scripts/pod/m2c_arm_e_train.py \
        --unit mondo__0005180 --batch-size 32 --seed 1337
    --unit mondo__0004975 for AD after PD.

Outputs:
    reports/m2c/checkpoints_completion_ana/<unit>/  (best.pth, last.pth, per-epoch jsonl)
    reports/m2c/arm_e_<unit>.json  (summary)
"""
import argparse
import bisect
import hashlib
import json
import random
import sys
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import torch
import torch.nn.functional as F

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "third_party" / "biopathnet"))
sys.path.insert(0, str(REPO_ROOT / "src"))

from biopathnet import dataset as bpn_dataset, task as bpn_task, model as bpn_model, util  # noqa: E402, F401

from torchdrug import core  # noqa: E402
from torchdrug.layers import functional as td_functional  # noqa: E402

EXPORTS_DIR = REPO_ROOT / "data/processed/m2c/completion_exports_ana"
PINNED_PATH = REPO_ROOT / "data/processed/m2c/pinned_genes_17448.txt"
CKPT_ROOT = REPO_ROOT / "reports/m2c/checkpoints_completion_ana"


def utcnow_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_pinned() -> list[str]:
    txt = PINNED_PATH.read_text(encoding="utf-8").strip().split()
    assert len(txt) == 17448, f"pinned len {len(txt)} != 17448"
    return txt


def build_type_map_from_export(unit_dir: Path):
    """Return dict entity_id -> type_id inferred from entity_types.txt first-appearance.
    Also returns raw mapping for verification vs HANDOFF loader map.
    """
    type_map = {}
    inv_type = {}
    with open(unit_dir / "entity_types.txt") as f:
        for line in f:
            eid, tid = line.strip().split()
            tid = int(tid)
            if eid not in type_map:
                type_map[eid] = tid
            if tid not in inv_type:
                # need name? we just store id->first entity's type token via nodes.parquet ideally
                inv_type[tid] = tid
    return type_map


def discover_relation_ids(rel_vocab: list[str], dataset_obj=None, train_set=None):
    """Heuristic discovery of DIS-GEN, BPO, PWY relation ids.

    Primary: infer DIS-GEN from train2 (all DIS-GEN share one rel id). For BPO/PWY
    use vocab strings when they contain bpo/pwy/pathway, else fall back to
    type-driven inference (done later in main). Returns possibly None for BPO/PWY.
    """
    lower = [r.lower() for r in rel_vocab]
    dis_gen = None
    bpo_rel = None
    pwy_rel = None
    # If dataset provided, infer DIS-GEN directly from train2 relation
    if dataset_obj is not None and train_set is not None:
        try:
            # train_set is Subset of dataset.graph containing train2
            rel_ids = set(int(dataset_obj.graph.edge_list[i, 2].item()) for i in train_set.indices)
            if len(rel_ids) == 1:
                dis_gen = next(iter(rel_ids))
        except Exception:
            pass
    if dis_gen is None:
        for i, r in enumerate(lower):
            if "dis" in r and "gen" in r:
                dis_gen = i
        # Fallback for real vocab where DIS-GEN is ASSOCIATED_WITH
        if dis_gen is None and "associated_with" in lower:
            # Need to disambiguate: many relations use associated_with, but disease-gene is one of them; if we didn't get via train_set, try most frequent disease-gene rel
            pass
    for i, r in enumerate(lower):
        if "bpo" in r or "biological_process" in r or r == "bp":
            if "dis" not in r and ("bpo" in r or "biological" in r):
                bpo_rel = i
        if "pwy" in r or "pathway" in r:
            if "dis" not in r and ("pwy" in r or "pathway" in r):
                pwy_rel = i
    if bpo_rel is None:
        for i, r in enumerate(lower):
            if "go" in r and "gen" in r:
                bpo_rel = i
                break
    if pwy_rel is None:
        for i, r in enumerate(lower):
            if "react" in r:
                pwy_rel = i
                break
    return dis_gen, bpo_rel, pwy_rel


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--unit", default="mondo__0005180", choices=["mondo__0005180", "mondo__0004975"])
    parser.add_argument("--config", type=Path, default=REPO_ROOT / "configs/m2c_arm_e.yaml")
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--seed", type=int, default=1337)
    parser.add_argument("--max-epochs", type=int, default=20)
    parser.add_argument("--patience", type=int, default=5)
    parser.add_argument("--min-delta", type=float, default=0.03, help="relative improvement 3%")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--eval-batch-size", type=int, default=4)
    parser.add_argument("--limit-batches", type=int, default=None, help="for debugging: cap batches per epoch")
    # Probe 2 (owner spec 2026-09-04): ranking-loss variant, opt-in only.
    parser.add_argument("--probe-loss", choices=["bce", "rank"], default="bce",
                        help="bce = frozen default (bit-identical). rank = probe: 1 positive + K masked negatives log-softmax")
    parser.add_argument("--probe-out", type=Path, default=None,
                        help="probe mode output dir (reports/m2c/objective_probe/). When set: probe instrumentation, "
                             "no checkpoint/summary writes, single capped pass only")
    parser.add_argument("--probe-k", default="255",
                        help="rank-loss negatives per query (frozen probe setting 255), or 'all' = full filtered "
                             "softmax over every allowed negative (gathered from the same scores — no extra forward)")
    # CHANGES/0027 (owner directive 2026-09-04): ANA-free KEEP graph + rank-loss epoch.
    parser.add_argument("--exports-dir", type=Path, default=EXPORTS_DIR,
                        help="unit exports root (default: ANA exports; KEEP-only = data/processed/m2c/completion_exports)")
    parser.add_argument("--ckpt-root", type=Path, default=CKPT_ROOT,
                        help="checkpoint root (default: ANA wave dir; override for non-ANA runs)")
    parser.add_argument("--val-every-seconds", type=int, default=0,
                        help="hourly mid-epoch validation (owner 2026-09-04): 0=off (default, frozen behavior); "
                             "3600 = log filtered val loss/MRR/H@10 every hour, logging only")
    parser.add_argument("--width", type=int, default=32,
                        help="model width: input_dim + hidden_dims [w]*6 (frozen 32; width probes only)")
    parser.add_argument("--empty-cache-every", type=int, default=1,
                        help="cuda empty_cache cadence in batches (1 = every batch, frozen behavior; "
                             "higher amortizes the sync when VRAM headroom is large)")
    args = parser.parse_args()
    probe_k_all = str(args.probe_k).lower() == "all"
    probe_k = 255 if probe_k_all else int(args.probe_k)

    device = 0
    if torch.cuda.is_available():
        torch.cuda.init()
        torch.cuda.set_device(device)

    random.seed(args.seed)
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)

    pinned_genes = load_pinned()
    set(pinned_genes)

    unit_dir = args.exports_dir / args.unit
    assert unit_dir.exists(), f"unit dir missing {unit_dir}"

    if args.dry_run:
        print("dry-run config: " + json.dumps(vars(args), sort_keys=True, default=str))
        print(f"dry-run: unit {args.unit} batch {args.batch_size} pinned {len(pinned_genes)}")
        return 0

    # --- dataset & task setup (reuse preprocess for fact graph) ---
    gpus = [device]
    cfg = util.load_config(str(args.config), context={"gpus": gpus})
    cfg.dataset.path = str(unit_dir)
    cfg.engine.gpus = gpus
    cfg.engine.batch_size = args.batch_size
    if args.width != 32:
        # Width probe only (owner 2026-09-04): frozen spec is 32 — never silently widen a wave run
        cfg.task.model.input_dim = args.width
        cfg.task.model.hidden_dims = [args.width] * 6
        print(f"WIDTH OVERRIDE: input_dim/hidden_dims = {args.width} (frozen spec: 32)", flush=True)
    # Ensure preprocessing flags match spec
    dataset_obj = core.Configurable.load_config_dict(cfg.dataset)
    print(f"dataset loaded: {dataset_obj.num_entity} entities, {dataset_obj.num_relation} relations", flush=True)
    # Build solver to get task with fact_graph registered
    # We create task manually to avoid Engine
    task_cfg = cfg.task
    task_cfg.model.num_relation = dataset_obj.num_relation
    task_obj = core.Configurable.load_config_dict(task_cfg)
    # Need train/valid/test splits for preprocess
    train_set, valid_set, test_set = dataset_obj.split()
    print(f"#train2 {len(train_set)} #valid {len(valid_set)} #test {len(test_set)}", flush=True)
    # preprocess registers fact_graph etc
    _train_set_p, _valid_set_p, _test_set_p = task_obj.preprocess(train_set, valid_set, test_set)
    # Move task to device
    task_obj = task_obj.cuda(device) if torch.cuda.is_available() else task_obj
    model = task_obj.model
    model.train()

    # Verify biotype tripwire and pinned universe
    # Use nodes.parquet for biotype but quick check: vocab size
    # Pinned gene indices
    inv_entity = dataset_obj.inv_entity_vocab
    inv_relation = dataset_obj.inv_relation_vocab
    rel_vocab = dataset_obj.relation_vocab if hasattr(dataset_obj, "relation_vocab") else [inv_relation[i] for i in range(dataset_obj.num_relation)]
    # Some datasets store relation_vocab as dict; handle
    if isinstance(rel_vocab, dict):
        rel_vocab = [rel_vocab[i] for i in range(len(rel_vocab))]
    dis_gen_rel, bpo_rel, pwy_rel = discover_relation_ids(rel_vocab, dataset_obj, train_set)
    print(f"relation vocab ({len(rel_vocab)}): {rel_vocab}", flush=True)
    print(f"discovered DIS-GEN={dis_gen_rel} ({rel_vocab[dis_gen_rel] if dis_gen_rel is not None else '?'}) "
          f"BPO={bpo_rel} ({rel_vocab[bpo_rel] if bpo_rel is not None else '?'}) "
          f"PWY={pwy_rel} ({rel_vocab[pwy_rel] if pwy_rel is not None else '?'})", flush=True)
    assert dis_gen_rel is not None, "DIS-GEN relation not found in vocab"
    # Relation disjointness guard (C4)
    assert bpo_rel != dis_gen_rel and pwy_rel != dis_gen_rel, "aux relation overlaps DIS-GEN"

    # Pinned gene indices
    pinned_indices = []
    missing = []
    for g in pinned_genes:
        if g in inv_entity:
            pinned_indices.append(inv_entity[g])
        else:
            missing.append(g)
    assert not missing, f"{len(missing)} pinned genes missing from vocab e.g. {missing[:5]}"
    assert len(pinned_indices) == 17448
    pinned_indices_t = torch.tensor(pinned_indices, dtype=torch.long, device=task_obj.device if hasattr(task_obj, 'device') else torch.device("cuda", device))
    pinned_set_idx = set(pinned_indices)

    # Disease node id
    disease_token = args.unit.replace("__", ":")
    assert disease_token in inv_entity, f"disease {disease_token} not in vocab"
    inv_entity[disease_token]

    # Build per-disease pos sets for DIS-GEN (train2 + valid for val filtering)
    # train2 edges are train_set (Subset)
    def edges_of(subset):
        return [tuple(dataset_obj.graph.edge_list[i].tolist()) for i in subset.indices]

    train2_edges = edges_of(train_set)
    valid_edges = edges_of(valid_set)
    test_edges = edges_of(test_set)

    # Filter to DIS-GEN only
    train2_disgen = [(h, t, r) for h, t, r in train2_edges if r == dis_gen_rel]
    valid_disgen = [(h, t, r) for h, t, r in valid_edges if r == dis_gen_rel]
    test_disgen = [(h, t, r) for h, t, r in test_edges if r == dis_gen_rel]
    print(f"train2 DIS-GEN {len(train2_disgen)} valid DIS-GEN {len(valid_disgen)} test {len(test_edges)}", flush=True)

    # Assert tails ⊆ pinned (C1)
    train2_tails = {t for h, t, r in train2_disgen}
    valid_tails = {t for h, t, r in valid_disgen}
    test_tails = {t for h, t, r in test_disgen}
    assert train2_tails.issubset(pinned_set_idx), "train2 tails not subset of pinned"
    assert valid_tails.issubset(pinned_set_idx), "valid tails not subset of pinned"
    assert test_tails.issubset(pinned_set_idx), "test tails not subset of pinned"

    # Per-disease positives (train2 only)
    pos_by_disease = defaultdict(set)
    for h, t, r in train2_disgen:
        pos_by_disease[h].add(t)
    # Per-disease test positives (probe-mode known-positive mask only; never targets/loss)
    test_pos_by_disease = defaultdict(set)
    for h, t, r in test_disgen:
        test_pos_by_disease[h].add(t)
    # For gene aux, need BPO/PWY pools from train1
    # train1 is not in train_set; fetch from dataset.graph with mask for train1
    # num_samples[0] is train1 count
    num_train1 = dataset_obj.num_samples[0] if hasattr(dataset_obj, "num_samples") else 0
    # Extract train1 edges
    train1_triplets = dataset_obj.graph.edge_list[:num_train1] if num_train1 else []
    # Determine gene type id and BPO/PWY type ids via graph.node_type
    node_type = dataset_obj.graph.node_type  # tensor num_node
    # Infer type ids by looking at pinned genes
    gene_type = int(node_type[pinned_indices[0]].item())
    print(f"gene_type inferred {gene_type}", flush=True)
    # Verify all pinned same type
    assert all(int(node_type[i].item()) == gene_type for i in pinned_indices[:100]), "pinned gene types heterogeneous"

    # BPO/PWY type inference: use aux relation tails if found, else fallback to HANDOFF map
    # Build aux pools
    bpo_tuples = []  # list of (gene_idx, rel, term_idx)
    pwy_tuples = []
    gene_to_bpo = defaultdict(set)
    gene_to_pwy = defaultdict(set)
    bpo_term_set = set()
    pwy_term_set = set()
    fallback_bpo_type = 1
    fallback_pwy_type = 7
    for idx in range(num_train1):
        h, t, r = train1_triplets[idx].tolist() if isinstance(train1_triplets, torch.Tensor) else tuple(train1_triplets[idx])
        ht = int(node_type[h].item())
        tt = int(node_type[t].item())
        # Check if this edge involves gene + BPO/PWY
        is_bpo_edge = (r == bpo_rel) if bpo_rel is not None else False
        is_pwy_edge = (r == pwy_rel) if pwy_rel is not None else False
        # If relation ids not discovered, fallback to type-based detection
        if bpo_rel is None or pwy_rel is None:
            # use type ids
            if (ht == gene_type and tt == fallback_bpo_type) or (tt == gene_type and ht == fallback_bpo_type):
                is_bpo_edge = True
            if (ht == gene_type and tt == fallback_pwy_type) or (tt == gene_type and ht == fallback_pwy_type):
                is_pwy_edge = True
        if is_bpo_edge:
            # original orientation in train1 is typically BPO/PWY -> GEN (e.g. go:... INTERACTS_WITH ensembl:...)
            # For query we reorient as GEN head -> term tail, but for masking we need original (h,t,r) as stored
            if ht == gene_type and tt != gene_type:
                g, term = h, t
                orig = (h, t, r)
            elif tt == gene_type and ht != gene_type:
                g, term = t, h
                orig = (h, t, r)
            else:
                continue
            # store (orig_h, orig_t, r, gene, term) — orig for masking, gene/term for query
            bpo_tuples.append((orig[0], orig[1], r, g, term))
            gene_to_bpo[g].add(term)
            bpo_term_set.add(term)
        if is_pwy_edge:
            if ht == gene_type and tt != gene_type:
                g, term = h, t
                orig = (h, t, r)
            elif tt == gene_type and ht != gene_type:
                g, term = t, h
                orig = (h, t, r)
            else:
                continue
            pwy_tuples.append((orig[0], orig[1], r, g, term))
            gene_to_pwy[g].add(term)
            pwy_term_set.add(term)

    print(f"aux pools: BPO {len(bpo_tuples)} (genes {len(gene_to_bpo)} terms {len(bpo_term_set)}) "
          f"PWY {len(pwy_tuples)} (genes {len(gene_to_pwy)} terms {len(pwy_term_set)})", flush=True)
    # If pools empty, try generic fallback: any train1 edge where one endpoint is gene and other is not disease/drug/anatomy etc and relation not DIS-GEN
    if not bpo_tuples and not pwy_tuples:
        print("WARNING: aux pools empty — scanning for any gene-* train1 edges", flush=True)
        # collect all distinct type pairs for train1
        type_pairs = set()
        for idx in range(min(num_train1, 5000)):
            h, t, r = train1_triplets[idx].tolist()
            type_pairs.add((int(node_type[h].item()), int(node_type[t].item()), r))
        print(f"sample type pairs {list(type_pairs)[:10]}", flush=True)

    # Build vocab lists for aux scoring: BPO vocab = sorted term indices, PWY vocab = sorted
    bpo_vocab = sorted(bpo_term_set)
    pwy_vocab = sorted(pwy_term_set)
    # If still empty, fallback to HANDOFF type ids: collect all nodes of those types
    if not bpo_vocab:
        bpo_vocab = [i for i in range(dataset_obj.num_entity) if int(node_type[i].item()) == fallback_bpo_type]
        print(f"fallback BPO vocab from node_type {len(bpo_vocab)}", flush=True)
    if not pwy_vocab:
        pwy_vocab = [i for i in range(dataset_obj.num_entity) if int(node_type[i].item()) == fallback_pwy_type]
        print(f"fallback PWY vocab from node_type {len(pwy_vocab)}", flush=True)

    bpo_vocab_t = torch.tensor(bpo_vocab, dtype=torch.long, device=task_obj.device if hasattr(task_obj, 'device') else torch.device("cuda", device))
    pwy_vocab_t = torch.tensor(pwy_vocab, dtype=torch.long, device=task_obj.device if hasattr(task_obj, 'device') else torch.device("cuda", device))

    # Mapping term -> position in vocab for target construction
    bpo_pos_map = {term: idx for idx, term in enumerate(bpo_vocab)}
    pwy_pos_map = {term: idx for idx, term in enumerate(pwy_vocab)}
    # For pinned genes, map gene idx -> position in pinned vocab
    pinned_pos_map = {g: i for i, g in enumerate(pinned_indices)}

    # Optimizer
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
    # Fact graph (message passing graph) - use task's fact_graph
    fact_graph = task_obj.fact_graph
    # Ensure on device
    if torch.cuda.is_available():
        fact_graph = fact_graph.cuda(device)

    # Helper: mask edges per batch (remove queried edge IDs)
    def masked_fact_graph(batch_h, batch_t, batch_r, strict=True):
        B = batch_h.shape[0]
        pattern = torch.stack([batch_h, batch_t, batch_r], dim=-1)
        edge_index, _ = fact_graph.match(pattern)
        if strict:
            assert edge_index.numel() == B, f"mask: expected {B} edge hits, got {edge_index.numel()}"
        else:
            if edge_index.numel() != B:
                # For aux where orientation or relation may vary, fallback to no masking and warn
                print(f"WARNING aux mask: expected {B} hits got {edge_index.numel()}, skipping mask", flush=True)
                return fact_graph
        edge_mask = ~td_functional.as_mask(edge_index, fact_graph.num_edge)
        masked = fact_graph.edge_mask(edge_mask)
        if strict:
            assert masked.num_edge == fact_graph.num_edge - B
        return masked

    def run_validation():
        """Filtered 1-hot BCE + filtered MRR/Hits@10 over pinned 17,448 (C3 evaluator).
        Scoring only — val edges are not in the message graph, per-row removal is training-only.
        Restores model.train() before returning."""
        model.eval()
        val_losses = []
        val_ranks = []
        with torch.no_grad():
            for (h, t, r) in valid_disgen:
                bh = torch.tensor([h], dtype=torch.long, device=fact_graph.device)
                bt = torch.tensor([t], dtype=torch.long, device=fact_graph.device)
                br = torch.tensor([r], dtype=torch.long, device=fact_graph.device)
                h_idx = bh.unsqueeze(1).expand(-1, len(pinned_indices))
                r_idx = br.unsqueeze(1).expand(-1, len(pinned_indices))
                t_idx = pinned_indices_t.unsqueeze(0)
                scores = model(fact_graph, h_idx, t_idx, r_idx, all_loss=None, metric=None).squeeze(0)  # [17448]
                # Filtered 1-hot BCE: exclude train2 ∪ valid positives of this disease (except query)
                disease = h
                train_pos = pos_by_disease.get(disease, set())
                valid_pos = valid_pos_by_disease.get(disease, set())
                excluded = (train_pos | valid_pos) - {t}
                t_pos = pinned_pos_map[t]
                filt_scores = []
                filt_targets = []
                for idx, gene_idx in enumerate(pinned_indices):
                    if gene_idx in excluded:
                        continue
                    filt_scores.append(scores[idx].item())
                    filt_targets.append(1.0 if gene_idx == t else 0.0)
                filt_scores_t = torch.tensor(filt_scores, dtype=torch.float32)
                filt_targets_t = torch.tensor(filt_targets, dtype=torch.float32)
                loss_val = F.binary_cross_entropy_with_logits(filt_scores_t, filt_targets_t)
                val_losses.append(float(loss_val.item()))
                # Filtered ranking, optimistic ties (bisect_left on -scores)
                sorted_filt = sorted(filt_scores, reverse=True)
                true_score = scores[t_pos].item()
                neg_sorted = [-s for s in sorted_filt]
                rank = bisect.bisect_left(neg_sorted, -true_score) + 1
                val_ranks.append(rank)
        val_loss = float(np.mean(val_losses)) if val_losses else 0.0
        val_mrr = float(np.mean([1.0 / r for r in val_ranks])) if val_ranks else 0.0
        val_h10 = float(np.mean([1.0 if r <= 10 else 0.0 for r in val_ranks])) if val_ranks else 0.0
        model.train()
        return val_loss, val_mrr, val_h10

    # Training data: DIS-GEN positives list
    # Shuffle each epoch
    disgen_list = train2_disgen  # list of (h,t,r)
    # For validation, prepare structures
    # Valid edges for loss: list of (h,t,r)
    # Need per-disease pos sets for filtering (train2 + valid)
    valid_pos_by_disease = defaultdict(set)
    for h, t, r in valid_disgen:
        valid_pos_by_disease[h].add(t)
    {d: pos_by_disease[d] | valid_pos_by_disease.get(d, set()) for d in pos_by_disease}
    # For filtered ranking, also include train2+valid positives per disease

    # ---- Probe 2 instrumentation (owner spec 2026-09-04) — active ONLY with --probe-out ----
    probe_mode = args.probe_out is not None
    PROBE_K = 255
    PROBE_TEMPERATURE = 1.0
    probe_gen = torch.Generator().manual_seed(args.seed)
    probe_ranks: list = []
    probe_top10_degs: list = []
    probe_allowed_by_disease: dict = {}
    probe_audit: dict = {}
    probe_jsonl = None
    pinned_deg = np.zeros(len(pinned_indices), dtype=np.float64)
    if probe_mode or args.probe_loss == "rank":
        if probe_mode:
            args.probe_out.mkdir(parents=True, exist_ok=True)
            probe_jsonl = args.probe_out / f"probe_{args.probe_loss}_{args.unit}_b{args.batch_size}_s{args.seed}.jsonl"
        # Known-positive mask per disease: train2 ∪ valid ∪ test tails (C1 pinned positions).
        pinned_pos_set = set(range(len(pinned_indices)))
        for d in pos_by_disease:
            known_genes = pos_by_disease[d] | valid_pos_by_disease.get(d, set()) | test_pos_by_disease.get(d, set())
            known_pos = {pinned_pos_map[g] for g in known_genes if g in pinned_pos_map}
            assert len(known_pos) == len(known_genes), "known positive outside pinned list (C1 violated)"
            allowed = sorted(pinned_pos_set - known_pos)
            assert len(allowed) == len(pinned_indices) - len(known_pos), \
                f"mask count {len(allowed)} != {len(pinned_indices)} - {len(known_pos)}"
            probe_allowed_by_disease[d] = torch.tensor(allowed, dtype=torch.long)
            # Contamination audit: what fraction of UNMASKED samples would hit known positives
            audit_gen = torch.Generator().manual_seed(args.seed + 1)
            hits, draws = 0, 0
            for _ in range(200):
                perm = torch.randperm(len(pinned_indices), generator=audit_gen)[:PROBE_K]
                hits += sum(1 for p in perm.tolist() if p in known_pos)
                draws += PROBE_K
            probe_audit[str(d)] = {
                "known_positives_masked": len(known_pos),
                "allowed_negatives": len(allowed),
                "expected_unmasked_contamination": round(len(known_pos) / len(pinned_indices), 4),
                "empirical_unmasked_contamination": round(hits / draws, 4),
                "audit_draws": draws,
            }
            print(f"probe mask disease {d}: known {len(known_pos)} masked, allowed {len(allowed)}; "
                  f"unmasked contamination expected {len(known_pos) / len(pinned_indices):.4f} "
                  f"empirical {hits / draws:.4f}", flush=True)
        # Node degrees (popularity-bias audit) from the full export graph
        el = dataset_obj.graph.edge_list
        el_cpu = el.cpu() if torch.is_tensor(el) else torch.tensor(el)
        deg = (np.bincount(el_cpu[:, 0].numpy(), minlength=dataset_obj.num_entity)
               + np.bincount(el_cpu[:, 1].numpy(), minlength=dataset_obj.num_entity))
        pinned_deg = deg[np.asarray(pinned_indices)]

    # Checkpoint paths
    ckpt_dir = args.ckpt_root / args.unit
    ckpt_dir.mkdir(parents=True, exist_ok=True)
    best_path = ckpt_dir / "best.pth"
    last_path = ckpt_dir / "last.pth"
    jsonl_path = ckpt_dir / f"train_log_{args.unit}.jsonl"
    summary_path = REPO_ROOT / f"reports/m2c/arm_e_{args.unit}.json"

    start_epoch = 1
    best_val_loss = float("inf")
    best_epoch = 0
    patience_counter = 0

    # Resume
    if args.resume and last_path.exists():
        try:
            ckpt = torch.load(last_path, map_location="cpu")
            model.load_state_dict(ckpt["model"])
            optimizer.load_state_dict(ckpt["optimizer"])
            start_epoch = ckpt.get("epoch", 0) + 1
            best_val_loss = ckpt.get("best_val_loss", float("inf"))
            best_epoch = ckpt.get("best_epoch", 0)
            patience_counter = ckpt.get("patience_counter", 0)
            print(f"resumed from epoch {start_epoch} best {best_val_loss:.4f} best_epoch {best_epoch}", flush=True)
        except Exception as e:
            print(f"resume failed: {e}", flush=True)

    # RSS guard
    try:
        import psutil
        has_psutil = True
    except ImportError:
        has_psutil = False

    def check_rss():
        if has_psutil:
            rss_gb = psutil.Process().memory_info().rss / 1e9
            for child in psutil.Process().children(recursive=True):
                with contextlib.suppress(Exception):
                    rss_gb += child.memory_info().rss / 1e9
            if rss_gb > 40:
                raise RuntimeError(f"RSS {rss_gb:.1f}GB exceeds 40GB guard")

    # Per-epoch loop
    # Precompute pinned gene idx tensor for scoring: we need h_index [B, 17448], t_index [B,17448], r_index [B,17448]
    # For DIS-GEN, t_index is pinned_vocab per row; h_index is disease repeated; r_index is dis_gen_rel repeated.
    # For aux, similarly with bpo/pwy vocab.

    num_batches_per_epoch = (len(disgen_list) + args.batch_size - 1) // args.batch_size
    print(f"batches per epoch {num_batches_per_epoch} (disgen {len(disgen_list)} bs {args.batch_size})", flush=True)

    for epoch in range(start_epoch, args.max_epochs + 1):
        model.train()
        epoch_start = time.perf_counter()
        last_val_ts = time.perf_counter()  # for --val-every-seconds hourly mid-epoch validation
        if torch.cuda.is_available():
            torch.cuda.reset_peak_memory_stats(device)
        train_losses = []
        batch_times = []

        # Shuffle DIS-GEN
        random.shuffle(disgen_list)
        # Shuffle aux pools
        random.shuffle(bpo_tuples)
        random.shuffle(pwy_tuples)
        bpo_ptr = 0
        pwy_ptr = 0

        for batch_idx in range(num_batches_per_epoch):
            if args.limit_batches and batch_idx >= args.limit_batches:
                break
            check_rss()
            b_start = batch_idx * args.batch_size
            b_end = min(b_start + args.batch_size, len(disgen_list))
            batch = disgen_list[b_start:b_end]
            B = len(batch)
            if B == 0:
                continue
            batch_h = torch.tensor([h for h, t, r in batch], dtype=torch.long, device=fact_graph.device)
            batch_t = torch.tensor([t for h, t, r in batch], dtype=torch.long, device=fact_graph.device)
            batch_r = torch.tensor([r for h, t, r in batch], dtype=torch.long, device=fact_graph.device)

            # Build masked graph for DIS-GEN
            masked_graph = masked_fact_graph(batch_h, batch_t, batch_r)

            # Build scoring tensors for DIS-GEN
            # h_index [B,17448], t_index [B,17448], r_index [B,17448]
            h_index = batch_h.unsqueeze(1).expand(-1, len(pinned_indices))
            r_index = batch_r.unsqueeze(1).expand(-1, len(pinned_indices))
            t_index = pinned_indices_t.unsqueeze(0).expand(B, -1)

            # Forward
            t0 = time.perf_counter()
            scores = model(masked_graph, h_index, t_index, r_index, all_loss=None, metric=None)  # [B,17448]
            if probe_mode or args.probe_loss == "rank":
                row_pos = torch.tensor([pinned_pos_map[t] for h, t, r in batch], dtype=torch.long)  # CPU [B]
            if args.probe_loss == "bce":
                # Build multi-hot target
                target = torch.zeros_like(scores)
                for i, (h, t, r) in enumerate(batch):
                    positives = pos_by_disease[h]  # set of gene idx
                    # map to positions in pinned vocab
                    for gene_idx in positives:
                        pos = pinned_pos_map.get(gene_idx)
                        if pos is not None:
                            target[i, pos] = 1.0
                # BCE
                loss_dis = F.binary_cross_entropy_with_logits(scores, target)
            else:
                # Probe 2 rank variant: 1 held-out positive (row's own target tail) + PROBE_K
                # seeded negatives sampled WITHOUT replacement from pinned list EXCLUDING all
                # known positives (train2+valid+test) of the row's disease. Gather down the
                # already-computed scores (NO rescoring); log-softmax temperature 1.0; NLL of
                # the positive placed at index 0.
                rows_by_disease = defaultdict(list)
                for i, (h, t, r) in enumerate(batch):
                    rows_by_disease[h].append(i)
                if probe_k_all:
                    # K=all: full filtered softmax over every allowed negative per disease
                    # (gathered from the same scores — no extra forward; CHANGES/0027).
                    loss_terms = []
                    for h, rows in rows_by_disease.items():
                        allowed = probe_allowed_by_disease[h]
                        sel_g = torch.cat([row_pos[rows].unsqueeze(1),
                                           allowed.unsqueeze(0).expand(len(rows), -1)], dim=1)
                        gathered_g = torch.gather(scores[rows], 1, sel_g.to(scores.device))
                        logp_g = F.log_softmax(gathered_g / PROBE_TEMPERATURE, dim=1)
                        loss_g = F.nll_loss(logp_g, torch.zeros(len(rows), dtype=torch.long,
                                                                device=scores.device))
                        loss_terms.append(loss_g * len(rows))
                    loss_dis = torch.stack(loss_terms).sum() / B
                    sel = gathered = logp = None
                else:
                    sel = torch.empty((B, 1 + probe_k), dtype=torch.long)
                    sel[:, 0] = row_pos
                    for h, rows in rows_by_disease.items():
                        allowed = probe_allowed_by_disease[h]
                        # vectorized per-row without-replacement: argsort of seeded random keys
                        keys = torch.rand(len(rows), allowed.numel(), generator=probe_gen)
                        sel[rows, 1:] = allowed[keys.argsort(dim=1)[:, :probe_k]]
                    assert not (sel[:, 1:] == sel[:, 0:1]).any(), "positive leaked into negatives"
                    sel = sel.to(scores.device)  # [B, 1+K] pinned positions
                    gathered = torch.gather(scores, 1, sel)
                    logp = F.log_softmax(gathered / PROBE_TEMPERATURE, dim=1)
                    loss_dis = F.nll_loss(logp, torch.zeros(B, dtype=torch.long, device=scores.device))

            # Probe instrumentation (both variants): filtered running rank of the queried
            # positive from the SAME forward scores (no extra forward). Filter = train2 ∪
            # valid positives of the row's disease minus the query itself (C3, matches val).
            if probe_mode:
                with torch.no_grad():
                    excl_mask = torch.zeros_like(scores, dtype=torch.bool)
                    for i, (h, t, r) in enumerate(batch):
                        excl = (pos_by_disease[h] | valid_pos_by_disease.get(h, set())) - {t}
                        pos_list = [pinned_pos_map[g] for g in excl if g in pinned_pos_map]
                        if pos_list:
                            excl_mask[i, pos_list] = True
                    filt_sc = scores.masked_fill(excl_mask, float("-inf"))
                    true_sc = scores.gather(1, row_pos.to(scores.device).unsqueeze(1))
                    ranks = (filt_sc > true_sc).sum(dim=1) + 1  # optimistic rank
                    probe_ranks.extend(ranks.cpu().tolist())
                    top10_pos = torch.topk(filt_sc, k=10, dim=1).indices.cpu().numpy()
                    probe_top10_degs.extend(pinned_deg[top10_pos].ravel().tolist())
                    del excl_mask, filt_sc, true_sc, ranks, top10_pos

            # Free DIS-GEN forward graph via sequential backward: backward DIS-GEN before AUX forward
            optimizer.zero_grad()
            loss_dis.backward()
            # keep loss_dis value for logging, free graph tensors
            loss_dis_val = float(loss_dis.item())
            if args.probe_loss == "bce":
                del scores, target
            else:
                del scores, sel, gathered, logp
            if torch.cuda.is_available() and (batch_idx + 1) % args.empty_cache_every == 0:
                torch.cuda.empty_cache()

            # Aux 1:1 interleave — forward AUX after DIS-GEN backward to avoid 2× graph + inplace version 2
            aux_batch = None
            aux_vocab = None
            aux_vocab_t = None
            aux_pos_map = None
            gene_to_aux = None
            loss = None
            loss_aux_val = None
            if bpo_tuples or pwy_tuples:
                use_bpo = (batch_idx % 2 == 0 and bpo_ptr < len(bpo_tuples)) or pwy_ptr >= len(pwy_tuples)
                if use_bpo and bpo_ptr < len(bpo_tuples):
                    aux_slice = bpo_tuples[bpo_ptr:bpo_ptr + B]
                    bpo_ptr += len(aux_slice)
                    if len(aux_slice) < B and bpo_ptr >= len(bpo_tuples):
                        bpo_ptr = 0
                    aux_batch = aux_slice
                    aux_vocab = bpo_vocab
                    aux_vocab_t = bpo_vocab_t
                    aux_pos_map = bpo_pos_map
                    gene_to_aux = gene_to_bpo
                else:
                    aux_slice = pwy_tuples[pwy_ptr:pwy_ptr + B]
                    pwy_ptr += len(aux_slice)
                    if len(aux_slice) < B and pwy_ptr >= len(pwy_tuples):
                        pwy_ptr = 0
                    aux_batch = aux_slice
                    aux_vocab = pwy_vocab
                    aux_vocab_t = pwy_vocab_t
                    aux_pos_map = pwy_pos_map
                    gene_to_aux = gene_to_pwy
                if aux_batch:
                    Ab = len(aux_batch)
                    ah_orig = torch.tensor([h for h, t, r, _, _ in aux_batch], dtype=torch.long, device=fact_graph.device)
                    at_orig = torch.tensor([t for h, t, r, _, _ in aux_batch], dtype=torch.long, device=fact_graph.device)
                    ar = torch.tensor([r for h, t, r, _, _ in aux_batch], dtype=torch.long, device=fact_graph.device)
                    ah = torch.tensor([g for h, t, r, g, _ in aux_batch], dtype=torch.long, device=fact_graph.device)
                    masked_aux = masked_fact_graph(ah_orig, at_orig, ar, strict=False)
                    ah_index = ah.unsqueeze(1).expand(-1, len(aux_vocab))
                    ar_index = ar.unsqueeze(1).expand(-1, len(aux_vocab))
                    at_index = aux_vocab_t.unsqueeze(0).expand(Ab, -1)
                    aux_scores = model(masked_aux, ah_index, at_index, ar_index, all_loss=None, metric=None)
                    aux_target = torch.zeros_like(aux_scores)
                    for i, (h, t, r, g, term) in enumerate(aux_batch):
                        positives = gene_to_aux.get(g, set())
                        for term_pos in positives:
                            pos = aux_pos_map.get(term_pos)
                            if pos is not None:
                                aux_target[i, pos] = 1.0
                    loss_aux = F.binary_cross_entropy_with_logits(aux_scores, aux_target)
                    loss_aux.backward()
                    optimizer.step()
                    loss_aux_val = float(loss_aux.item())
                    loss = loss_dis_val + loss_aux_val
                    del loss_dis, loss_aux, aux_scores, aux_target
                else:
                    optimizer.step()
                    loss = loss_dis_val
                    del loss_dis
            else:
                optimizer.step()
                loss = loss_dis_val
                del loss_dis
            # need loss as float for train_losses
            if isinstance(loss, float):
                loss_val = loss
            else:
                loss_val = float(loss.item()) if hasattr(loss, "item") else float(loss)
                del loss

            train_losses.append(float(loss_val))
            batch_times.append((time.perf_counter() - t0) * 1000)
            if batch_idx % 250 == 0:
                run_avg = float(np.mean(train_losses[-250:]))
                pct = 100.0 * (batch_idx + 1) / num_batches_per_epoch
                el = time.perf_counter() - epoch_start
                eta_h = (el / (batch_idx + 1)) * (num_batches_per_epoch - batch_idx - 1) / 3600
                print(f"  epoch {epoch} batch {batch_idx + 1}/{num_batches_per_epoch} ({pct:.1f}%) "
                      f"loss250 {run_avg:.4f} eta {eta_h:.1f}h", flush=True)

            if probe_mode and (batch_idx + 1) % 50 == 0:
                win_ranks = probe_ranks[-(50 * args.batch_size):]
                win_degs = probe_top10_degs[-(50 * args.batch_size * 10):]
                run_mrr = float(np.mean([1.0 / r for r in probe_ranks])) if probe_ranks else 0.0
                entry = {
                    "variant": args.probe_loss,
                    "batch": batch_idx + 1,
                    "loss50": float(np.mean(train_losses[-50:])),
                    "window_train_filtered_mrr": float(np.mean([1.0 / r for r in win_ranks])) if win_ranks else 0.0,
                    "running_train_filtered_mrr": run_mrr,
                    "top10_degree_mean": float(np.mean(win_degs)) if win_degs else 0.0,
                    "top10_degree_median": float(np.median(win_degs)) if win_degs else 0.0,
                    "ms_per_batch50": float(np.mean(batch_times[-50:])),
                    "peak_reserved_vram_gb": round(torch.cuda.max_memory_reserved(device) / 1e9, 2)
                    if torch.cuda.is_available() else 0.0,
                    "elapsed_s": round(time.perf_counter() - epoch_start, 1),
                    "wall_iso": utcnow_iso(),
                }
                with open(probe_jsonl, "a", encoding="utf-8") as f:
                    f.write(json.dumps(entry) + "\n")
                print(f"  probe[{args.probe_loss}] batch {batch_idx + 1} loss50 {entry['loss50']:.5f} "
                      f"winMRR {entry['window_train_filtered_mrr']:.5f} runMRR {run_mrr:.5f} "
                      f"top10deg {entry['top10_degree_mean']:.0f} "
                      f"res {entry['peak_reserved_vram_gb']}GB", flush=True)

            # Hourly mid-epoch validation (owner 2026-09-04): early-signal curve, LOGGING ONLY —
            # checkpointing and early-stop remain epoch-boundary decisions (frozen discipline).
            if args.val_every_seconds > 0 and (time.perf_counter() - last_val_ts) >= args.val_every_seconds:
                hv_loss, hv_mrr, hv_h10 = run_validation()
                last_val_ts = time.perf_counter()
                hv_entry = {
                    "kind": "hourly_val",
                    "epoch": epoch,
                    "batch": batch_idx + 1,
                    "val_loss": hv_loss,
                    "val_mrr": hv_mrr,
                    "val_hits10": hv_h10,
                    "elapsed_s": round(time.perf_counter() - epoch_start, 1),
                    "wall_iso": utcnow_iso(),
                }
                hv_target = probe_jsonl if probe_mode else jsonl_path
                with open(hv_target, "a", encoding="utf-8") as f:
                    f.write(json.dumps(hv_entry) + "\n")
                print(f"  [hourly-val] epoch {epoch} batch {batch_idx + 1} val {hv_loss:.6f} "
                      f"mrr {hv_mrr:.5f} h10 {hv_h10:.3f}", flush=True)

        # End of epoch training
        avg_train_loss = float(np.mean(train_losses)) if train_losses else 0.0
        ms_per_batch = float(np.mean(batch_times)) if batch_times else 0.0
        epoch_train_sec = time.perf_counter() - epoch_start
        peak_vram_gb = round(torch.cuda.max_memory_allocated(device) / 1e9, 2) if torch.cuda.is_available() else 0.0

        # Validation: filtered 1-hot BCE + filtered MRR/Hits@10 (C3 evaluator)
        val_loss, val_mrr, val_h10 = run_validation()

        # Log
        log_entry = {
            "epoch": epoch,
            "train_loss": avg_train_loss,
            "val_loss": val_loss,
            "val_mrr": val_mrr,
            "val_hits10": val_h10,
            "peak_vram_gb": peak_vram_gb,
            "ms_per_batch": ms_per_batch,
            "epoch_seconds": round(epoch_train_sec, 1),
            "num_train_batches": len(train_losses),
            "wall_iso": utcnow_iso(),
        }
        print(f"epoch {epoch}: train {avg_train_loss:.4f} val {val_loss:.4f} mrr {val_mrr:.4f} h10 {val_h10:.3f} "
              f"vram {peak_vram_gb}GB {ms_per_batch:.1f}ms/batch ({epoch_train_sec:.1f}s)", flush=True)
        if probe_mode:
            # Probe 2: single capped pass, no checkpoint/summary writes outside probe-out dir.
            probe_summary = {
                "variant": args.probe_loss,
                "unit": args.unit,
                "seed": args.seed,
                "batch_size": args.batch_size,
                "limit_batches": args.limit_batches,
                "probe_k": args.probe_k,
                "temperature": PROBE_TEMPERATURE,
                "contamination_audit": probe_audit,
                "train_loss_avg": avg_train_loss,
                "running_train_filtered_mrr": float(np.mean([1.0 / r for r in probe_ranks])) if probe_ranks else 0.0,
                "val_loss": val_loss,
                "val_mrr": val_mrr,
                "val_hits10": val_h10,
                "peak_vram_gb": peak_vram_gb,
                "peak_reserved_vram_gb": round(torch.cuda.max_memory_reserved(device) / 1e9, 2)
                if torch.cuda.is_available() else 0.0,
                "ms_per_batch": ms_per_batch,
                "epoch_seconds": round(epoch_train_sec, 1),
                "jsonl": str(probe_jsonl),
                "wall_iso": utcnow_iso(),
            }
            probe_summary_path = args.probe_out / f"probe_summary_{args.probe_loss}_{args.unit}_b{args.batch_size}_s{args.seed}.json"
            probe_summary_path.write_text(json.dumps(probe_summary, indent=2) + "\n", encoding="utf-8")
            print(f"probe[{args.probe_loss}] done -> {probe_summary_path} (no checkpoints written)", flush=True)
            break

        # Append jsonl
        with open(jsonl_path, "a", encoding="utf-8") as f:
            f.write(json.dumps(log_entry) + "\n")

        # Checkpoint
        ckpt = {
            "epoch": epoch,
            "model": model.state_dict(),
            "optimizer": optimizer.state_dict(),
            "best_val_loss": best_val_loss,
            "best_epoch": best_epoch,
            "patience_counter": patience_counter,
            "config": {"unit": args.unit, "batch_size": args.batch_size, "seed": args.seed},
        }
        torch.save(ckpt, last_path)
        # Owner request 2026-09-02: keep EVERY epoch's checkpoint (epoch_001.pth, ...)
        torch.save(ckpt, ckpt_dir / f"epoch_{epoch:03d}.pth")
        # Best tracking
        if val_loss < best_val_loss * (1 - args.min_delta):
            best_val_loss = val_loss
            best_epoch = epoch
            patience_counter = 0
            torch.save(ckpt, best_path)
            print(f"  -> new best val {best_val_loss:.4f} at epoch {epoch}", flush=True)
        else:
            patience_counter += 1
            print(f"  -> no improvement patience {patience_counter}/{args.patience}", flush=True)

        # Early stop check
        if patience_counter >= args.patience:
            print(f"early stop at epoch {epoch} (patience {args.patience})", flush=True)
            break
        if epoch >= args.max_epochs:
            print(f"max epochs {args.max_epochs} reached", flush=True)
            break

    # Final summary (probe mode writes nothing outside --probe-out)
    if not probe_mode:
        summary = {
            "unit": args.unit,
            "seed": args.seed,
            "batch_size": args.batch_size,
            "best_epoch": best_epoch,
            "best_val_loss": best_val_loss,
            "patience": args.patience,
            "min_delta": args.min_delta,
            "max_epochs": args.max_epochs,
            "finished_at": utcnow_iso(),
            "checkpoints": {"best": str(best_path), "last": str(last_path), "log": str(jsonl_path)},
        }
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        print(f"done unit {args.unit} best epoch {best_epoch} val {best_val_loss:.4f} -> {summary_path}", flush=True)
    else:
        print(f"done probe[{args.probe_loss}] unit {args.unit} (probe mode — no wave summary written)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
