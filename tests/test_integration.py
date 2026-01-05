"""Integration tests for GWAS pipeline.

These tests run the full pipeline on small synthetic data to verify
end-to-end functionality without requiring real data.
"""

import numpy as np
import pandas as pd
import pytest

from src.config import GWASConfig, QCConfig
from src.data.preprocessor import Preprocessor
from src.gwas.runner import GWASRunner
from src.gwas.stats import calculate_lambda_gc, calculate_prps_batch


@pytest.fixture
def integration_data() -> tuple[pd.DataFrame, pd.Series, np.ndarray]:
    """Generate synthetic data for integration testing.

    Creates 100 samples × 100 variants with:
    - 5 causal variants with strong effect
    - Realistic kinship structure
    - Binary phenotype correlated with causal variants
    """
    np.random.seed(42)
    n_samples = 100
    n_variants = 100
    n_causal = 5

    # Generate SNP matrix (0, 1, 2 encoding)
    snps = np.random.binomial(2, 0.3, (n_samples, n_variants)).astype(float)

    # Generate kinship (population structure)
    # Simulate 3 subpopulations
    pop_labels = np.repeat([0, 1, 2], [40, 30, 30])
    kinship = np.zeros((n_samples, n_samples))
    for i in range(n_samples):
        for j in range(n_samples):
            if pop_labels[i] == pop_labels[j]:
                kinship[i, j] = 0.1  # Within-population relatedness
            kinship[i, j] += 0.05  # Baseline relatedness
    np.fill_diagonal(kinship, 1.0)

    # Generate phenotype correlated with causal variants
    causal_idx = np.random.choice(n_variants, n_causal, replace=False)
    genetic_effect = snps[:, causal_idx].sum(axis=1)
    # Add noise and convert to binary
    liability = genetic_effect + np.random.normal(0, 2, n_samples)
    phenotype = (liability > np.median(liability)).astype(int)

    # Create DataFrames
    sample_ids = [f"sample_{i}" for i in range(n_samples)]
    variant_ids = [f"{1000000 + i * 100}_A_G" for i in range(n_variants)]

    snps_df = pd.DataFrame(snps, index=sample_ids, columns=variant_ids)
    pheno_series = pd.Series(phenotype, index=sample_ids, name="TestDrug")

    return snps_df, pheno_series, kinship


class TestGWASIntegration:
    """Integration tests for the GWAS pipeline."""

    def test_gwas_pipeline_lmm(self, integration_data: tuple) -> None:
        """Test full LMM GWAS pipeline runs without errors."""
        snps, phenotype, kinship = integration_data

        # Configure and run GWAS
        config = GWASConfig(method="lmm", min_maf=0.01, p_threshold=0.05)
        runner = GWASRunner(config, n_jobs=1)

        result = runner.run(
            snps=snps,
            phenotype=phenotype,
            kinship=kinship,
            drug_name="TestDrug",
        )

        # Verify result structure
        assert result is not None
        assert hasattr(result, "variants")
        assert hasattr(result, "lambda_gc")
        assert hasattr(result, "n_significant")

        # Verify variants DataFrame
        assert len(result.variants) == snps.shape[1]
        assert "variant" in result.variants.columns
        assert "p_value" in result.variants.columns
        assert "prps" in result.variants.columns

        # Verify p-values are valid
        p_values = result.variants["p_value"].values
        assert np.all(p_values >= 0)
        assert np.all(p_values <= 1)
        assert not np.any(np.isnan(p_values))

        # Lambda GC should be reasonable (not extremely inflated)
        assert 0.0 < result.lambda_gc < 3.0

    def test_gwas_pipeline_simple(self, integration_data: tuple) -> None:
        """Test simple (non-LMM) GWAS pipeline."""
        snps, phenotype, kinship = integration_data

        # Simple method still needs kinship for sample alignment
        config = GWASConfig(method="simple", min_maf=0.01)
        runner = GWASRunner(config, n_jobs=1)

        result = runner.run(
            snps=snps,
            phenotype=phenotype,
            kinship=kinship,
            drug_name="TestDrug",
        )

        assert result is not None
        assert len(result.variants) == snps.shape[1]

    def test_preprocessor_qc(self, integration_data: tuple) -> None:
        """Test QC preprocessing pipeline."""
        snps, phenotype, _ = integration_data

        # Add some variants that should be filtered
        snps_with_issues = snps.copy()
        # Add monomorphic variant
        snps_with_issues["mono_var"] = 0
        # Add low MAF variant
        snps_with_issues["low_maf_var"] = 0
        snps_with_issues.iloc[0, -1] = 1  # Only one sample has the variant

        # Create phenotypes DataFrame (preprocessor expects DataFrame, not Series)
        phenotypes_df = pd.DataFrame({"TestDrug": phenotype})

        qc_config = QCConfig(min_maf=0.05, max_missing_rate=0.1)
        preprocessor = Preprocessor(qc_config)
        snps_qc, pheno_qc, report = preprocessor.run_qc(snps_with_issues, phenotypes_df)

        # Verify QC removed problematic variants
        assert "mono_var" not in snps_qc.columns
        assert "low_maf_var" not in snps_qc.columns
        assert report.n_variants_before > report.n_variants_after

    def test_prps_batch_calculation(self, integration_data: tuple) -> None:
        """Test batch PRPS calculation on integration data."""
        snps, _, kinship = integration_data

        # Compute eigendecomposition
        eigenvalues, eigenvectors = np.linalg.eigh(kinship)

        # Calculate PRPS for all variants
        prps_scores = calculate_prps_batch(
            snps.values,
            eigenvectors,
            eigenvalues,
        )

        # Verify output
        assert len(prps_scores) == snps.shape[1]
        assert np.all(prps_scores >= 0)
        assert np.all(prps_scores <= 1)
        assert not np.any(np.isnan(prps_scores))

    def test_lambda_gc_on_synthetic_data(self, integration_data: tuple) -> None:
        """Test lambda GC calculation on synthetic data."""
        # Generate null p-values (uniform under null)
        np.random.seed(123)
        null_pvalues = np.random.uniform(0, 1, 1000)

        lambda_gc = calculate_lambda_gc(null_pvalues)

        # Under null, lambda should be close to 1.0
        assert 0.8 < lambda_gc < 1.2


class TestEndToEnd:
    """End-to-end tests simulating real usage."""

    def test_minimal_gwas_run(self) -> None:
        """Minimal test: 50 samples, 50 variants, runs in <5 seconds."""
        np.random.seed(999)

        # Minimal data
        n_samples, n_variants = 50, 50
        snps = pd.DataFrame(
            np.random.binomial(2, 0.3, (n_samples, n_variants)),
            index=[f"s{i}" for i in range(n_samples)],
            columns=[f"v{i}" for i in range(n_variants)],
        ).astype(float)

        phenotype = pd.Series(
            np.random.randint(0, 2, n_samples),
            index=snps.index,
            name="Drug",
        )

        # Simple kinship
        kinship = np.eye(n_samples) * 0.9 + 0.1

        # Run GWAS
        config = GWASConfig(method="lmm", min_maf=0.01)
        runner = GWASRunner(config, n_jobs=1)
        result = runner.run(snps, phenotype, kinship, drug_name="Drug")

        # Basic assertions
        assert result is not None
        assert len(result.variants) == n_variants
        assert result.drug == "Drug"
