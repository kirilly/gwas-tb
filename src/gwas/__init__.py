"""GWAS analysis module."""

from .runner import GWASError, GWASResult, GWASRunner
from .stats import (
    bonferroni_correction,
    calculate_lambda_gc,
    calculate_prps,
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
