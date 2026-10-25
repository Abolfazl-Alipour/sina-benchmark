"""Config loading: YAML files into pydantic models (plan §0.3: config-driven YAML only)."""

from pathlib import Path

import yaml
from pydantic import BaseModel


def load_yaml(path: Path | str) -> dict:
    """Load a YAML file into a plain dict."""
    path = Path(path)
    with path.open("r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh)
    if not isinstance(data, dict):
        raise ValueError(f"config {path} must be a YAML mapping at top level")
    return data


def load_config[T: BaseModel](path: Path | str, model: type[T]) -> T:
    """Load a YAML config file and validate it against a pydantic model."""
    return model.model_validate(load_yaml(path))
