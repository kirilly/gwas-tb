"""Tests for GWAS module."""

import numpy as np
import pytest

from src.gwas.stats import (
    calculate_lambda_gc,
    calculate_prps,
    bonferroni_correction,
    fdr_correction,
    odds_ratio_ci,
)


def test_lambda_gc_null():
    """Test lambda GC with null distribution."""
    np.random.seed(42)
    p_values = np.random.uniform(0, 1, 10000)

    lambda_gc = calculate_lambda_gc(p_values)

    # Should be close to 1.0 for uniform distribution
    assert 0.9 < lambda_gc < 1.1


def test_lambda_gc_inflated():
    """Test lambda GC with inflated p-values."""
    np.random.seed(42)
    # Biased toward small p-values
    p_values = np.random.uniform(0, 0.2, 10000)

    lambda_gc = calculate_lambda_gc(p_values)

    # Should be > 1.0 for inflated distribution
    assert lambda_gc > 1.5


def test_lambda_gc_handles_nan():
    """Test lambda GC handles NaN values."""
    p_values = np.array([0.1, 0.5, np.nan, 0.01, 0.9])

    lambda_gc = calculate_lambda_gc(p_values)

    assert not np.isnan(lambda_gc)


def test_calculate_prps(mock_kinship):
    """Test PRPS calculation."""
    np.random.seed(42)
    n = mock_kinship.shape[0]
    snp = np.random.randint(0, 2, n)

    prps = calculate_prps(snp, mock_kinship)

    assert 0 <= prps <= 1


def test_calculate_prps_monomorphic(mock_kinship):
    """Test PRPS with monomorphic variant."""
    n = mock_kinship.shape[0]
    snp = np.zeros(n)

    prps = calculate_prps(snp, mock_kinship)

    assert prps == 0.0


def test_bonferroni_correction():
    """Test Bonferroni correction."""
    p_values = np.array([0.001, 0.01, 0.05, 0.1])

    adjusted, threshold = bonferroni_correction(p_values, alpha=0.05)

    assert threshold == 0.05 / 4
    assert all(adjusted >= p_values)
    assert all(adjusted <= 1.0)


def test_fdr_correction():
    """Test FDR correction."""
    p_values = np.array([0.001, 0.01, 0.02, 0.5])

    adjusted, rejected = fdr_correction(p_values, alpha=0.05)

    assert len(adjusted) == len(p_values)
    assert len(rejected) == len(p_values)
    assert isinstance(rejected[0], (bool, np.bool_))


def test_odds_ratio_ci():
    """Test odds ratio calculation."""
    # 2x2 table
    n11, n10, n01, n00 = 50, 30, 20, 100

    odds_ratio, ci_lower, ci_upper = odds_ratio_ci(n11, n10, n01, n00)

    expected_or = (50 * 100) / (30 * 20)
    assert abs(odds_ratio - expected_or) < 0.01
    assert ci_lower < odds_ratio < ci_upper


def test_odds_ratio_ci_zero_cell():
    """Test odds ratio with zero cell (Haldane correction)."""
    n11, n10, n01, n00 = 0, 30, 20, 100

    odds_ratio, ci_lower, ci_upper = odds_ratio_ci(n11, n10, n01, n00)

    # Should not fail with zero cell
    assert not np.isnan(odds_ratio)
    assert ci_lower < odds_ratio < ci_upper
