"""Caching utilities."""

import functools
import hashlib
import json
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd

from .logging import get_logger

logger = get_logger()


def compute_cache_key(*args: Any, **kwargs: Any) -> str:
    """Compute deterministic cache key from inputs.

    Args:
        *args: Positional arguments
        **kwargs: Keyword arguments

    Returns:
        12-character hash string
    """
    key_data = {
        "args": [_serialize_arg(a) for a in args],
        "kwargs": {k: _serialize_arg(v) for k, v in sorted(kwargs.items())},
    }
    key_str = json.dumps(key_data, sort_keys=True)
    return hashlib.md5(key_str.encode()).hexdigest()[:12]


def _serialize_arg(arg: Any) -> str:
    """Serialize argument for hashing."""
    if hasattr(arg, "shape"):  # numpy array or pandas
        return f"{type(arg).__name__}:{arg.shape}"
    if isinstance(arg, Path):
        return str(arg)
    return str(arg)


class CacheEncoder(json.JSONEncoder):
    """JSON encoder that handles numpy arrays and pandas objects."""

    def default(self, obj: Any) -> Any:
        if isinstance(obj, np.ndarray):
            return {"__numpy__": True, "data": obj.tolist(), "dtype": str(obj.dtype)}
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, pd.DataFrame):
            return {
                "__dataframe__": True,
                "data": obj.to_dict(orient="split"),
            }
        if isinstance(obj, pd.Series):
            return {
                "__series__": True,
                "data": obj.to_dict(),
                "name": obj.name,
            }
        return super().default(obj)


def cache_decode(obj: dict) -> Any:
    """Decode cached JSON back to numpy/pandas objects."""
    if isinstance(obj, dict):
        if obj.get("__numpy__"):
            return np.array(obj["data"], dtype=obj["dtype"])
        if obj.get("__dataframe__"):
            return pd.DataFrame(**obj["data"])
        if obj.get("__series__"):
            return pd.Series(obj["data"], name=obj["name"])
    return obj


def cached(cache_dir: str = "cache") -> Callable:
    """Decorator for caching function results using JSON.

    Args:
        cache_dir: Directory for cache files

    Returns:
        Decorated function
    """

    def decorator(func: Callable) -> Callable:
        @functools.wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            key = compute_cache_key(func.__name__, *args, **kwargs)
            cache_path = Path(cache_dir) / f"{func.__name__}_{key}.json"

            if cache_path.exists():
                logger.info(f"Loading cached: {cache_path.name}")
                with open(cache_path) as f:
                    data = json.load(f, object_hook=cache_decode)
                return data

            result = func(*args, **kwargs)

            cache_path.parent.mkdir(parents=True, exist_ok=True)
            with open(cache_path, "w") as f:
                json.dump(result, f, cls=CacheEncoder)
            logger.info(f"Cached: {cache_path.name}")

            return result

        return wrapper

    return decorator


def clear_cache(cache_dir: str = "cache", pattern: str = "*") -> int:
    """Clear cache files.

    Args:
        cache_dir: Directory containing cache files
        pattern: Glob pattern for files to remove

    Returns:
        Number of files removed
    """
    cache_path = Path(cache_dir)
    if not cache_path.exists():
        return 0

    removed = 0
    for f in cache_path.glob(f"{pattern}.json"):
        f.unlink()
        removed += 1
        logger.info(f"Removed cache: {f.name}")

    return removed
