"""Reproducibility utilities."""

import random
import subprocess

import numpy as np


def set_seeds(seed: int = 42) -> None:
    """Set random seeds for reproducibility.

    Args:
        seed: Random seed value
    """
    random.seed(seed)
    np.random.seed(seed)

    try:
        import torch

        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
            torch.backends.cudnn.deterministic = True
            torch.backends.cudnn.benchmark = False
    except ImportError:
        pass


def get_git_hash() -> str | None:
    """Get current git commit hash.

    Returns:
        Short git hash or None if not in a git repo
    """
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        )
        return result.stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


def get_run_metadata(seed: int, config_path: str) -> dict:
    """Get metadata for a pipeline run.

    Args:
        seed: Random seed used
        config_path: Path to config file

    Returns:
        Dictionary with run metadata
    """
    from datetime import datetime

    return {
        "timestamp": datetime.now().isoformat(),
        "git_hash": get_git_hash(),
        "random_seed": seed,
        "config_path": str(config_path),
    }
