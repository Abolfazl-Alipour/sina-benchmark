"""PyKEEN KGE baseline wrapper (spec §6.3).

Trains deterministic knowledge-graph embeddings on a training edge set and
scores candidate tails for (query, relation) pairs in the same session.
"""

from pathlib import Path

import numpy as np
import polars as pl
import torch
from pykeen.pipeline import pipeline
from pykeen.predict import predict_target
from pykeen.triples import TriplesFactory


def _to_triples(edges: pl.DataFrame) -> np.ndarray:
    return np.array(
        [
            (r["src_id"], r["relation"], r["dst_id"])
            for r in edges.select("src_id", "relation", "dst_id").to_dicts()
        ],
        dtype=str,
    )


def train_and_score_kge(
    train_edges: pl.DataFrame,
    query_nodes: list[str],
    candidate_nodes: list[str],
    relation: str,
    model: str = "TransE",
    embedding_dim: int = 256,
    epochs: int = 50,
    seed: int = 1337,
    output_dir: Path | None = None,
    device: str = "cpu",
) -> pl.DataFrame:
    """Train a PyKEEN model and score candidate tails for each query node."""
    triples = _to_triples(train_edges)
    factory = TriplesFactory.from_labeled_triples(triples)
    result = pipeline(
        training=factory,
        testing=factory,  # not used for training; we score manually
        model=model,
        model_kwargs={"embedding_dim": embedding_dim},
        training_loop="sLCWA",
        optimizer="Adam",
        optimizer_kwargs={"lr": 0.001},
        loss="softplus",
        training_kwargs={"num_epochs": epochs, "use_tqdm": False},
        random_seed=seed,
        device=device,
    )
    if output_dir is not None:
        result.save_to_directory(output_dir)

    model_obj = result.model
    model_obj.eval()
    rows: list[tuple[str, str, float]] = []
    with torch.no_grad():
        for query in query_nodes:
            try:
                pred = predict_target(
                    model=model_obj,
                    triples_factory=factory,
                    head=query,
                    relation=relation,
                    targets=candidate_nodes,
                )
            except KeyError:
                # query or relation unseen in training
                continue
            for _, r in pred.df.iterrows():
                rows.append((query, r["tail_label"], float(r["score"])))
    return pl.DataFrame(
        {
            "query_id": [r[0] for r in rows],
            "candidate_id": [r[1] for r in rows],
            "score": [r[2] for r in rows],
        },
        schema={"query_id": pl.String, "candidate_id": pl.String, "score": pl.Float64},
    )
