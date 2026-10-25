#!/usr/bin/env python
"""m2_monogenic_rank — monogenic-gene rank extraction for monogenic panel (A8).

For each monogenic unit (primary + secondary panel, n≈10), report exact rank of
the frozen monogenic gene. Panel membership = monogenic units only (trait nodes
like skin aging never panel members). Gene list frozen with ClinGen/OMIM.

Usage: python scripts/m2_monogenic_rank.py --config configs/m2_monogenic_rank.yaml [--dry-run]
"""

import argparse
from pathlib import Path

import polars as pl
from pydantic import BaseModel

from sina.config import load_config
from sina.logging import configure_logging, get_logger

log = get_logger(__name__)


class MonogenicConfig(BaseModel):
    lodo_exports_dir: Path
    monogenic_genes: dict[str, str]  # curie -> ensembl:ENSG
    output_path: Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    configure_logging()
    cfg = load_config(args.config, MonogenicConfig)
    if args.dry_run:
        print(f"would rank monogenic genes for {len(cfg.monogenic_genes)} panel units")
        return 0
    rows = []
    for curie, gene in cfg.monogenic_genes.items():
        export_dir = cfg.lodo_exports_dir / curie.replace(":", "__")
        if not export_dir.exists():
            rows.append({"unit": curie, "monogenic_gene": gene, "rank": None, "note": "no export"})
            continue
        # Placeholder: rank would be computed from model scores; here we record gene
        rows.append(
            {
                "unit": curie,
                "monogenic_gene": gene,
                "rank": None,
                "note": "model not run — rank pending",
            }
        )
    pl.DataFrame(rows).write_parquet(cfg.output_path)
    log.info("monogenic_rank_placeholder", units=len(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
