"""GWAS analysis module."""

from .runner import GWASRunner, GWASResult, GWASError
from .stats import (
    calculate_lambda_gc,
    calculate_prps,
    bonferroni_correction,
    fdr_correction,
)

__all__ = [
    "GWASRunner",
    "GWASResult",
    "GWASError",
    "calculate_lambda_gc",
    "calculate_prps",
    "bonferroni_correction",
    "fdr_correction",
]
