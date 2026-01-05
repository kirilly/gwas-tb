"""Pytest fixtures for TB GWAS tests."""

import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def small_snp_matrix() -> pd.DataFrame:
    """Generate small SNP matrix for testing."""
    np.random.seed(42)
    n_samples = 100
    n_variants = 500

    # Generate random genotypes (0, 1, 2)
    data = np.random.randint(0, 3, (n_samples, n_variants))

    # Add some missing values
    missing_mask = np.random.random((n_samples, n_variants)) < 0.02
    data = data.astype(float)
    data[missing_mask] = -1

    return pd.DataFrame(
        data,
        index=[f"sample_{i}" for i in range(n_samples)],
        columns=[f"var_{i}" for i in range(n_variants)],
    )


@pytest.fixture
def mock_phenotype() -> pd.Series:
    """Generate mock binary phenotype."""
    np.random.seed(42)
    n_samples = 100

    return pd.Series(
        np.random.randint(0, 2, n_samples),
        index=[f"sample_{i}" for i in range(n_samples)],
        name="RIF",
    )


@pytest.fixture
def mock_phenotypes() -> pd.DataFrame:
    """Generate mock phenotypes for multiple drugs."""
    np.random.seed(42)
    n_samples = 100
    drugs = ["RIF", "INH", "FQ"]

    data = np.random.randint(0, 2, (n_samples, len(drugs)))

    return pd.DataFrame(
        data,
        index=[f"sample_{i}" for i in range(n_samples)],
        columns=drugs,
    )


@pytest.fixture
def mock_kinship() -> np.ndarray:
    """Generate mock kinship matrix."""
    np.random.seed(42)
    n_samples = 100

    # Generate random symmetric PSD matrix
    A = np.random.randn(n_samples, n_samples)
    kinship = A @ A.T / n_samples

    # Ensure diagonal is 1
    d = np.sqrt(np.diag(kinship))
    kinship = kinship / np.outer(d, d)

    return kinship


@pytest.fixture
def temp_config_file(tmp_path) -> str:
    """Create temporary config file."""
    config_content = """
project:
  name: "test-gwas"
  random_seed: 42
  n_jobs: 2
  log_level: "DEBUG"

data:
  snp_matrix: "test_snps.csv"
  phenotypes: "test_pheno.csv"
  who_catalogue: "who.csv"

qc:
  min_maf: 0.01
  max_missing_rate: 0.1
  min_samples_per_class: 5

phylogeny:
  method: "fasttree"
  kinship_method: "snp"
  n_threads: 2

gwas:
  method: "lmm"
  covariates: []
  p_threshold: 0.001
  fdr_threshold: 0.1

drugs:
  order: ["RIF", "INH"]
  conditional: false

output:
  results_dir: "results/"
  cache_dir: "cache/"
  figures_dir: "figures/"
  logs_dir: "logs/"
"""
    config_file = tmp_path / "test_config.yaml"
    config_file.write_text(config_content)
    return str(config_file)
