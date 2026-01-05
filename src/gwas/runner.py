"""GWAS analysis runner."""

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd
from tqdm import tqdm

from src.config import GWASConfig
from src.utils import get_logger
from .stats import calculate_lambda_gc, calculate_prps, fdr_correction, odds_ratio_ci

logger = get_logger()


class GWASError(Exception):
    """Exception for GWAS-related errors."""

    pass


@dataclass
class GWASResult:
    """GWAS analysis result."""

    drug: str
    variants: pd.DataFrame  # variant, gene, p_value, OR, CI_low, CI_high, prps
    lambda_gc: float
    n_significant: int
    n_samples: int
    n_variants: int
    covariates_used: list[str] = field(default_factory=list)

    def get_significant(self, threshold: float = 5e-8) -> pd.DataFrame:
        """Get significant variants below threshold."""
        return self.variants[self.variants["p_value"] < threshold].copy()

    def get_top_hits(self, n: int = 20) -> pd.DataFrame:
        """Get top N hits by p-value."""
        return self.variants.nsmallest(n, "p_value").copy()


class GWASRunner:
    """Run genome-wide association analysis."""

    def __init__(self, config: GWASConfig) -> None:
        """Initialize GWASRunner.

        Args:
            config: GWAS configuration
        """
        self.config = config

    def run(
        self,
        snps: pd.DataFrame,
        phenotype: pd.Series,
        kinship: np.ndarray,
        covariates: Optional[pd.DataFrame] = None,
        drug_name: str = "unknown",
    ) -> GWASResult:
        """Run GWAS analysis.

        Args:
            snps: SNP matrix [n_samples, n_variants]
            phenotype: Binary phenotype [n_samples]
            kinship: Kinship matrix [n_samples, n_samples]
            covariates: Optional covariates [n_samples, n_covariates]
            drug_name: Name of drug being analyzed

        Returns:
            GWASResult with per-variant statistics
        """
        logger.info(f"Running GWAS for {drug_name}")
        logger.info(f"Samples: {snps.shape[0]}, Variants: {snps.shape[1]}")

        # Align data
        common = snps.index.intersection(phenotype.index)
        if len(common) < 50:
            raise GWASError(f"Insufficient samples: {len(common)}")

        snps_aligned = snps.loc[common]
        pheno_aligned = phenotype.loc[common]
        kinship_aligned = self._align_kinship(kinship, snps.index, common)

        if covariates is not None:
            covariates = covariates.loc[common]

        # Run analysis
        if self.config.method == "lmm":
            results = self._run_lmm(snps_aligned, pheno_aligned, kinship_aligned, covariates)
        else:
            results = self._run_simple(snps_aligned, pheno_aligned, covariates)

        # Calculate PRPS for each variant
        logger.info("Calculating PRPS scores")
        prps_scores = []
        for var in tqdm(snps_aligned.columns, desc="PRPS", disable=len(snps_aligned.columns) < 100):
            prps = calculate_prps(snps_aligned[var].values, kinship_aligned)
            prps_scores.append(prps)
        results["prps"] = prps_scores

        # Calculate lambda GC
        lambda_gc = calculate_lambda_gc(results["p_value"].values)
        logger.info(f"Lambda GC: {lambda_gc:.3f}")

        if lambda_gc > 1.1:
            logger.warning(f"Inflated Lambda GC ({lambda_gc:.3f}) - check population structure")

        # Apply FDR correction
        results["p_adjusted"], results["significant"] = fdr_correction(
            results["p_value"].values, alpha=self.config.fdr_threshold
        )

        n_significant = (results["p_value"] < self.config.p_threshold).sum()
        logger.info(f"Significant variants (p < {self.config.p_threshold}): {n_significant}")

        return GWASResult(
            drug=drug_name,
            variants=results,
            lambda_gc=lambda_gc,
            n_significant=n_significant,
            n_samples=len(common),
            n_variants=snps_aligned.shape[1],
            covariates_used=list(covariates.columns) if covariates is not None else [],
        )

    def _run_lmm(
        self,
        snps: pd.DataFrame,
        phenotype: pd.Series,
        kinship: np.ndarray,
        covariates: Optional[pd.DataFrame],
    ) -> pd.DataFrame:
        """Run Linear Mixed Model GWAS using pyseer-style approach."""
        from scipy import stats
        from scipy.linalg import cho_factor, cho_solve

        n_samples = len(phenotype)
        n_variants = snps.shape[1]

        # Prepare response
        y = phenotype.values.astype(float)

        # Prepare fixed effects design matrix
        if covariates is not None and len(covariates.columns) > 0:
            X_cov = covariates.values
        else:
            X_cov = np.ones((n_samples, 1))

        # Eigendecomposition of kinship for efficient LMM
        eigenvalues, eigenvectors = np.linalg.eigh(kinship)
        eigenvalues = np.maximum(eigenvalues, 1e-10)  # Ensure positive

        # Transform data
        y_rotated = eigenvectors.T @ y
        X_cov_rotated = eigenvectors.T @ X_cov

        # Results storage
        results = []

        logger.info("Running LMM GWAS")
        for i, var_id in enumerate(tqdm(snps.columns, desc="GWAS", disable=n_variants < 100)):
            snp_vec = snps.iloc[:, i].values.astype(float)

            # Skip if monomorphic
            if np.var(snp_vec) == 0:
                results.append({
                    "variant": var_id,
                    "p_value": np.nan,
                    "odds_ratio": np.nan,
                    "ci_lower": np.nan,
                    "ci_upper": np.nan,
                    "beta": np.nan,
                    "se": np.nan,
                })
                continue

            # Rotate SNP
            snp_rotated = eigenvectors.T @ snp_vec

            # Build design matrix
            X = np.column_stack([X_cov_rotated, snp_rotated])

            # Weighted least squares (weights = 1 / eigenvalues)
            try:
                weights = 1.0 / eigenvalues
                W = np.diag(weights)

                # Solve weighted least squares
                XtWX = X.T @ W @ X
                XtWy = X.T @ W @ y_rotated

                # Use robust solve
                beta = np.linalg.lstsq(XtWX, XtWy, rcond=None)[0]

                # Residuals and variance
                residuals = y_rotated - X @ beta
                sigma2 = np.sum(weights * residuals**2) / (n_samples - X.shape[1])

                # Standard errors
                try:
                    cov_beta = sigma2 * np.linalg.inv(XtWX)
                    se = np.sqrt(np.diag(cov_beta))
                except np.linalg.LinAlgError:
                    se = np.full(len(beta), np.nan)

                # SNP effect is last coefficient
                snp_beta = beta[-1]
                snp_se = se[-1]

                # Wald test
                if snp_se > 0:
                    z = snp_beta / snp_se
                    p_value = 2 * stats.norm.sf(abs(z))
                else:
                    p_value = np.nan

                # Odds ratio (approximate for binary trait)
                odds_ratio = np.exp(snp_beta)
                ci_lower = np.exp(snp_beta - 1.96 * snp_se)
                ci_upper = np.exp(snp_beta + 1.96 * snp_se)

            except Exception:
                p_value = np.nan
                snp_beta = np.nan
                snp_se = np.nan
                odds_ratio = np.nan
                ci_lower = np.nan
                ci_upper = np.nan

            results.append({
                "variant": var_id,
                "p_value": p_value,
                "odds_ratio": odds_ratio,
                "ci_lower": ci_lower,
                "ci_upper": ci_upper,
                "beta": snp_beta,
                "se": snp_se,
            })

        return pd.DataFrame(results)

    def _run_simple(
        self,
        snps: pd.DataFrame,
        phenotype: pd.Series,
        covariates: Optional[pd.DataFrame],
    ) -> pd.DataFrame:
        """Run simple logistic regression GWAS (no kinship correction)."""
        from scipy import stats
        import statsmodels.api as sm

        results = []

        logger.info("Running simple GWAS (no kinship)")
        for var_id in tqdm(snps.columns, desc="GWAS"):
            snp_vec = snps[var_id].values

            if np.var(snp_vec) == 0:
                results.append({
                    "variant": var_id,
                    "p_value": np.nan,
                    "odds_ratio": np.nan,
                    "ci_lower": np.nan,
                    "ci_upper": np.nan,
                    "beta": np.nan,
                    "se": np.nan,
                })
                continue

            try:
                # Build design matrix
                if covariates is not None:
                    X = np.column_stack([snp_vec, covariates.values])
                else:
                    X = snp_vec.reshape(-1, 1)
                X = sm.add_constant(X)

                # Logistic regression
                model = sm.Logit(phenotype.values, X)
                result = model.fit(disp=0, method="bfgs", maxiter=100)

                # SNP coefficient is at index 1
                snp_beta = result.params[1]
                snp_se = result.bse[1]
                p_value = result.pvalues[1]

                odds_ratio = np.exp(snp_beta)
                ci_lower = np.exp(snp_beta - 1.96 * snp_se)
                ci_upper = np.exp(snp_beta + 1.96 * snp_se)

            except Exception:
                p_value = np.nan
                snp_beta = np.nan
                snp_se = np.nan
                odds_ratio = np.nan
                ci_lower = np.nan
                ci_upper = np.nan

            results.append({
                "variant": var_id,
                "p_value": p_value,
                "odds_ratio": odds_ratio,
                "ci_lower": ci_lower,
                "ci_upper": ci_upper,
                "beta": snp_beta,
                "se": snp_se,
            })

        return pd.DataFrame(results)

    def _align_kinship(
        self,
        kinship: np.ndarray,
        original_index: pd.Index,
        target_index: pd.Index,
    ) -> np.ndarray:
        """Align kinship matrix to target sample order."""
        # Get positions of target samples in original
        positions = [original_index.get_loc(s) for s in target_index]
        return kinship[np.ix_(positions, positions)]

    def run_conditional(
        self,
        snps: pd.DataFrame,
        phenotypes: pd.DataFrame,
        kinship: np.ndarray,
        drug_order: list[str],
    ) -> dict[str, GWASResult]:
        """Run conditional multi-drug GWAS.

        For each drug (in order):
        1. Run GWAS with previous drugs' top hits as covariates
        2. Store drug-specific results

        Args:
            snps: SNP matrix
            phenotypes: Phenotype DataFrame with drug columns
            kinship: Kinship matrix
            drug_order: Order of drugs to analyze

        Returns:
            Dict mapping drug name to GWASResult
        """
        results = {}
        accumulated_covariates: list[str] = []

        for drug in drug_order:
            if drug not in phenotypes.columns:
                logger.warning(f"Drug {drug} not in phenotypes, skipping")
                continue

            logger.info(f"Conditional GWAS for {drug}")

            # Build covariate matrix from previous drugs' top hits
            if accumulated_covariates:
                covariates = snps[accumulated_covariates].copy()
            else:
                covariates = None

            # Run GWAS
            result = self.run(
                snps=snps,
                phenotype=phenotypes[drug],
                kinship=kinship,
                covariates=covariates,
                drug_name=drug,
            )
            results[drug] = result

            # Add top hits to covariates for next drug
            top_hits = result.get_significant(self.config.p_threshold)
            if len(top_hits) > 0:
                new_covars = top_hits["variant"].tolist()[:10]  # Max 10 per drug
                accumulated_covariates.extend(new_covars)
                logger.info(f"Added {len(new_covars)} covariates from {drug}")

        return results
