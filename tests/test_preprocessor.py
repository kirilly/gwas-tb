"""Tests for preprocessor module."""

import numpy as np
import pandas as pd
import pytest

from src.config import QCConfig
from src.data import Preprocessor


@pytest.fixture
def preprocessor():
    """Create preprocessor with default config."""
    config = QCConfig(min_maf=0.05, max_missing_rate=0.1, min_samples_per_class=5)
    return Preprocessor(config)


def test_calculate_maf(preprocessor, small_snp_matrix):
    """Test MAF calculation."""
    maf = preprocessor.calculate_maf(small_snp_matrix)

    assert len(maf) == small_snp_matrix.shape[1]
    assert all(0 <= m <= 0.5 for m in maf)


def test_filter_by_maf(preprocessor, small_snp_matrix):
    """Test MAF filtering."""
    filtered, removed = preprocessor.filter_by_maf(small_snp_matrix, min_maf=0.1)

    assert filtered.shape[1] <= small_snp_matrix.shape[1]
    assert len(removed) == small_snp_matrix.shape[1] - filtered.shape[1]


def test_filter_by_missingness_samples(preprocessor, small_snp_matrix):
    """Test sample missingness filtering."""
    filtered, removed = preprocessor.filter_by_missingness(
        small_snp_matrix, max_missing=0.05, axis=0
    )

    assert filtered.shape[0] <= small_snp_matrix.shape[0]


def test_filter_by_missingness_variants(preprocessor, small_snp_matrix):
    """Test variant missingness filtering."""
    filtered, removed = preprocessor.filter_by_missingness(
        small_snp_matrix, max_missing=0.05, axis=1
    )

    assert filtered.shape[1] <= small_snp_matrix.shape[1]


def test_remove_monomorphic(preprocessor):
    """Test removal of monomorphic variants."""
    # Create matrix with some monomorphic columns
    data = np.array(
        [
            [0, 0, 1, 2],
            [0, 1, 1, 2],
            [0, 2, 1, 0],
        ]
    )
    df = pd.DataFrame(data, columns=["mono1", "poly1", "mono2", "poly2"])

    filtered, removed = preprocessor.remove_monomorphic(df)

    assert "mono1" in removed  # All zeros
    assert "mono2" in removed  # All ones
    assert "poly1" in filtered.columns
    assert "poly2" in filtered.columns


def test_run_qc(preprocessor, small_snp_matrix, mock_phenotypes):
    """Test full QC pipeline."""
    snps, pheno, report = preprocessor.run_qc(small_snp_matrix, mock_phenotypes)

    assert report.n_samples_before == 100
    assert report.n_variants_before == 500
    assert report.n_samples_after <= report.n_samples_before
    assert report.n_variants_after <= report.n_variants_before


def test_encode_missing(preprocessor, small_snp_matrix):
    """Test missing value encoding."""
    encoded = preprocessor.encode_missing(small_snp_matrix, strategy="mean")

    # Should have no -1 or NaN values
    assert (encoded == -1).sum().sum() == 0
    assert not encoded.isna().any().any()
