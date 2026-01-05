"""GWAS analysis runner."""

import os
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from tqdm import tqdm

from src.config import GWASConfig
from src.utils import get_logger
from .stats import calculate_lambda_gc, calculate_prps, calculate_prps_batch, fdr_correction, odds_ratio_ci

logger = get_logger()

PYSEER_AVAILABLE = False
try:
    result = subprocess.run(["pyseer", "--version"], capture_output=True, text=True)
    if result.returncode == 0:
        PYSEER_AVAILABLE = True
        logger.info(f"pyseer available: {result.stdout.strip()}")
except FileNotFoundError:
    pass


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

    def __init__(self, config: GWASConfig, n_jobs: int = 1) -> None:
        """Initialize GWASRunner.

        Args:
            config: GWAS configuration
            n_jobs: Number of parallel jobs for pyseer (default 1)
        """
        self.config = config
        self.n_jobs = n_jobs

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

        # Run analysis - returns results and optionally eigendecomposition
        eigenvalues = None
        eigenvectors = None

        if self.config.method == "lmm":
            if PYSEER_AVAILABLE:
                logger.info(f"Using pyseer Python API (no I/O overhead)")
                results, eigenvalues, eigenvectors = self._run_pyseer_api(
                    snps_aligned, pheno_aligned, kinship_aligned, covariates
                )
            else:
                logger.info("Using native LMM (serial)")
                results, eigenvalues, eigenvectors = self._run_lmm(
                    snps_aligned, pheno_aligned, kinship_aligned, covariates
                )
        else:
            results = self._run_simple(snps_aligned, pheno_aligned, covariates)

        # Calculate PRPS - use optimized eigenspace version if we have eigendecomposition
        if eigenvalues is not None and eigenvectors is not None:
            logger.info("Calculating PRPS scores (optimized eigenspace batch)")
            prps_scores = calculate_prps_batch(
                snps_aligned.values, eigenvectors, eigenvalues
            )
        else:
            # Fallback to naive O(n²) version with parallelization
            logger.info(f"Calculating PRPS scores with {self.n_jobs} jobs (naive method)")
            if self.n_jobs > 1:
                from joblib import Parallel, delayed
                prps_scores = Parallel(n_jobs=self.n_jobs)(
                    delayed(calculate_prps)(snps_aligned[var].values, kinship_aligned)
                    for var in tqdm(snps_aligned.columns, desc="PRPS", disable=len(snps_aligned.columns) < 100)
                )
            else:
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
    ) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
        """Run Linear Mixed Model GWAS using pyseer-style approach.

        Returns:
            Tuple of (results_df, eigenvalues, eigenvectors)
        """
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
        logger.info("Computing eigendecomposition of kinship matrix...")
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

        return pd.DataFrame(results), eigenvalues, eigenvectors

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

    def _run_pyseer_api(
        self,
        snps: pd.DataFrame,
        phenotype: pd.Series,
        kinship: np.ndarray,
        covariates: Optional[pd.DataFrame],
    ) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
        """Run GWAS using pyseer Python API (no file I/O overhead).

        Returns:
            Tuple of (results_df, eigenvalues, eigenvectors)
        """
        from pyseer.lmm import lmm_cov, fit_lmm_block
        from joblib import Parallel, delayed

        n_samples = len(phenotype)
        n_variants = snps.shape[1]

        # Prepare phenotype as column vector
        y = np.reshape(phenotype.values.astype(float), (-1, 1))

        # Prepare covariates (add intercept)
        if covariates is not None and len(covariates.columns) > 0:
            X = np.c_[covariates.values, np.ones((n_samples, 1))]
        else:
            X = np.ones((n_samples, 1))

        # Initialize LMM with kinship matrix directly (no file I/O!)
        logger.info("Initializing LMM (eigendecomposition)...")
        lmm = lmm_cov(X=X, Y=y, K=kinship, regressX=True, inplace=False)

        # Find optimal heritability (this triggers eigendecomposition)
        logger.info("Finding optimal h2...")
        result = lmm.findH2()
        h2 = result['h2']
        logger.info(f"Heritability h2 = {h2:.3f}")

        # Extract eigendecomposition from pyseer LMM object for PRPS reuse
        # pyseer populates lmm.U (eigenvectors) and lmm.S (eigenvalues) during findH2()
        eigenvectors = lmm.U
        eigenvalues = lmm.S

        # Process variants in blocks for parallelization
        block_size = max(100, n_variants // (self.n_jobs * 2))
        n_blocks = (n_variants + block_size - 1) // block_size

        logger.info(f"Processing {n_variants} variants in {n_blocks} blocks (size={block_size}) with {self.n_jobs} jobs")

        def process_block(start_idx: int, end_idx: int) -> dict:
            """Process a block of variants."""
            variant_block = snps.iloc[:, start_idx:end_idx].values.astype(float)
            return fit_lmm_block(lmm, h2, variant_block)

        # Parallel block processing
        block_ranges = [(i * block_size, min((i + 1) * block_size, n_variants))
                        for i in range(n_blocks)]

        if self.n_jobs > 1:
            block_results = Parallel(n_jobs=self.n_jobs, prefer="threads")(
                delayed(process_block)(start, end) for start, end in block_ranges
            )
        else:
            block_results = [process_block(start, end) for start, end in block_ranges]

        # Combine results
        all_p_values = np.concatenate([r['p_values'] for r in block_results])
        all_betas = np.concatenate([r['beta'] for r in block_results])
        all_ses = np.concatenate([r['bse'] for r in block_results])

        # Build results DataFrame
        results = pd.DataFrame({
            "variant": snps.columns,
            "p_value": all_p_values,
            "beta": all_betas,
            "se": all_ses,
            "odds_ratio": np.exp(all_betas),
            "ci_lower": np.exp(all_betas - 1.96 * all_ses),
            "ci_upper": np.exp(all_betas + 1.96 * all_ses),
        })

        logger.info(f"pyseer API returned {len(results)} variants")
        return results, eigenvalues, eigenvectors

    def _run_pyseer_cli(
        self,
        snps: pd.DataFrame,
        phenotype: pd.Series,
        kinship: np.ndarray,
        covariates: Optional[pd.DataFrame],
    ) -> pd.DataFrame:
        """Run GWAS using pyseer CLI (kept as fallback)."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # 1. Write phenotypes (tab-separated: sample, phenotype)
            pheno_file = tmpdir / "phenotypes.tsv"
            pheno_df = pd.DataFrame({"samples": phenotype.index, "phenotype": phenotype.values})
            pheno_df.to_csv(pheno_file, sep="\t", index=False)

            # 2. Write variants in Rtab format (variant as row, samples as columns)
            pres_file = tmpdir / "variants.Rtab"
            snps_T = snps.T
            snps_T.index.name = "Gene"
            snps_T.to_csv(pres_file, sep="\t")

            # 3. Write similarity matrix (kinship)
            sim_file = tmpdir / "similarity.tsv"
            sim_df = pd.DataFrame(kinship, index=phenotype.index, columns=phenotype.index)
            sim_df.to_csv(sim_file, sep="\t")

            # 4. Run pyseer
            output_file = tmpdir / "results.tsv"
            # Use smaller block_size for better parallelization
            # Default 3000 gives only ~2 blocks for 6K variants
            # 500 gives ~14 blocks → real parallel speedup
            block_size = min(500, max(100, snps.shape[1] // (self.n_jobs * 2)))

            cmd = [
                "pyseer",
                "--phenotypes", str(pheno_file),
                "--pres", str(pres_file),
                "--similarity", str(sim_file),
                "--lmm",
                "--cpu", str(self.n_jobs),
                "--block_size", str(block_size),
                "--min-af", str(self.config.min_maf),
                "--max-af", str(1 - self.config.min_maf),
            ]

            logger.info(f"Running pyseer: {' '.join(cmd)}")
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                cwd=str(tmpdir),
            )

            if result.returncode != 0:
                logger.error(f"pyseer failed: {result.stderr}")
                raise GWASError(f"pyseer failed: {result.stderr}")

            # 5. Parse output
            # pyseer outputs to stdout in format:
            # variant af  filter-pvalue   lrt-pvalue  beta    beta-std-err    intercept   ...
            lines = result.stdout.strip().split("\n")
            if not lines:
                raise GWASError("pyseer returned empty output")

            header = lines[0].split("\t")
            data = []
            for line in lines[1:]:
                fields = line.split("\t")
                if len(fields) >= 6:
                    try:
                        data.append({
                            "variant": fields[0],
                            "af": float(fields[1]) if fields[1] != "NA" else np.nan,
                            "p_value": float(fields[3]) if fields[3] != "NA" else np.nan,  # lrt-pvalue
                            "beta": float(fields[4]) if fields[4] != "NA" else np.nan,
                            "se": float(fields[5]) if fields[5] != "NA" else np.nan,
                        })
                    except (ValueError, IndexError) as e:
                        logger.warning(f"Failed to parse line: {line[:50]}... ({e})")

            results = pd.DataFrame(data)

            # Calculate odds ratios from beta
            results["odds_ratio"] = np.exp(results["beta"])
            results["ci_lower"] = np.exp(results["beta"] - 1.96 * results["se"])
            results["ci_upper"] = np.exp(results["beta"] + 1.96 * results["se"])

            logger.info(f"pyseer returned {len(results)} variants")
            return results

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
