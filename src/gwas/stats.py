"""Statistical functions for GWAS analysis."""

import numpy as np
from scipy import stats
from typing import Optional


def calculate_lambda_gc(p_values: np.ndarray) -> float:
    """Calculate genomic inflation factor (Lambda GC).

    Lambda = median(chi2_observed) / median(chi2_expected)

    Interpretation:
        < 1.0: Deflated (possibly underpowered)
        1.0: Perfect
        1.0-1.1: Acceptable
        > 1.1: Inflated (population structure not corrected)
        > 1.5: Severely inflated (results unreliable)

    Args:
        p_values: Array of p-values from GWAS

    Returns:
        Lambda GC value
    """
    # Remove NaN and zero p-values
    p_values = np.asarray(p_values)
    p_values = p_values[~np.isnan(p_values)]
    p_values = p_values[p_values > 0]

    if len(p_values) == 0:
        return np.nan

    # Convert p-values to chi-squared statistics (1 df)
    chi2_observed = stats.chi2.ppf(1 - p_values, df=1)

    # Expected median chi2 for 1 df
    chi2_expected = stats.chi2.ppf(0.5, df=1)

    lambda_gc = np.median(chi2_observed) / chi2_expected
    return float(lambda_gc)


def calculate_prps(
    snp: np.ndarray,
    kinship: np.ndarray,
) -> float:
    """Calculate Phylogeny-Related Parallelism Score.

    PRPS measures how well a variant follows the phylogeny.
    High PRPS (>0.5): Variant follows phylogeny (possibly hitchhiking)
    Low PRPS (<0.3): Convergent evolution (likely under selection)

    Args:
        snp: Genotype vector (0/1/2) for one variant
        kinship: Kinship matrix [n_samples, n_samples]

    Returns:
        PRPS score in [0, 1]
    """
    snp = np.asarray(snp).flatten()
    n = len(snp)

    if kinship.shape[0] != n:
        raise ValueError(f"Kinship size {kinship.shape[0]} != SNP length {n}")

    # Calculate expected covariance under phylogeny
    snp_centered = snp - snp.mean()
    var_snp = np.var(snp)

    if var_snp == 0:
        return 0.0

    # PRPS = correlation between genotype similarity and kinship
    # Create genotype similarity matrix
    geno_similarity = 1 - np.abs(np.subtract.outer(snp, snp)) / 2

    # Flatten upper triangle for correlation
    upper_idx = np.triu_indices(n, k=1)
    geno_flat = geno_similarity[upper_idx]
    kinship_flat = kinship[upper_idx]

    # Calculate correlation
    if np.std(geno_flat) == 0 or np.std(kinship_flat) == 0:
        return 0.0

    prps = np.corrcoef(geno_flat, kinship_flat)[0, 1]

    # Ensure in [0, 1]
    prps = max(0, min(1, (prps + 1) / 2))

    return float(prps)


def bonferroni_correction(
    p_values: np.ndarray,
    alpha: float = 0.05,
) -> tuple[np.ndarray, float]:
    """Apply Bonferroni multiple testing correction.

    Args:
        p_values: Array of p-values
        alpha: Family-wise error rate

    Returns:
        Tuple of (adjusted_p_values, significance_threshold)
    """
    p_values = np.asarray(p_values)
    n_tests = len(p_values[~np.isnan(p_values)])

    threshold = alpha / n_tests
    adjusted = np.minimum(p_values * n_tests, 1.0)

    return adjusted, threshold


def fdr_correction(
    p_values: np.ndarray,
    alpha: float = 0.05,
    method: str = "bh",
) -> tuple[np.ndarray, np.ndarray]:
    """Apply FDR multiple testing correction.

    Args:
        p_values: Array of p-values
        alpha: False discovery rate
        method: "bh" for Benjamini-Hochberg

    Returns:
        Tuple of (adjusted_p_values, rejection_mask)
    """
    from statsmodels.stats.multitest import multipletests

    p_values = np.asarray(p_values)

    # Handle NaN
    valid_mask = ~np.isnan(p_values)
    valid_p = p_values[valid_mask]

    if len(valid_p) == 0:
        return p_values, np.zeros(len(p_values), dtype=bool)

    reject, adjusted_valid, _, _ = multipletests(
        valid_p, alpha=alpha, method="fdr_bh"
    )

    # Reconstruct full arrays
    adjusted = np.full_like(p_values, np.nan)
    adjusted[valid_mask] = adjusted_valid

    rejection = np.zeros(len(p_values), dtype=bool)
    rejection[valid_mask] = reject

    return adjusted, rejection


def odds_ratio_ci(
    n11: int, n10: int, n01: int, n00: int,
    alpha: float = 0.05,
) -> tuple[float, float, float]:
    """Calculate odds ratio and confidence interval.

    Args:
        n11: Count of cases with variant
        n10: Count of cases without variant
        n01: Count of controls with variant
        n00: Count of controls without variant
        alpha: Significance level for CI

    Returns:
        Tuple of (odds_ratio, ci_lower, ci_upper)
    """
    # Add 0.5 to avoid division by zero (Haldane correction)
    if n11 == 0 or n10 == 0 or n01 == 0 or n00 == 0:
        n11 += 0.5
        n10 += 0.5
        n01 += 0.5
        n00 += 0.5

    odds_ratio = (n11 * n00) / (n10 * n01)

    # Log odds ratio and SE
    log_or = np.log(odds_ratio)
    se = np.sqrt(1/n11 + 1/n10 + 1/n01 + 1/n00)

    # CI
    z = stats.norm.ppf(1 - alpha/2)
    ci_lower = np.exp(log_or - z * se)
    ci_upper = np.exp(log_or + z * se)

    return float(odds_ratio), float(ci_lower), float(ci_upper)


def wald_test(
    coef: float,
    se: float,
) -> float:
    """Wald test for coefficient significance.

    Args:
        coef: Coefficient estimate
        se: Standard error

    Returns:
        Two-tailed p-value
    """
    if se == 0:
        return 0.0 if coef != 0 else 1.0

    z = coef / se
    p_value = 2 * stats.norm.sf(abs(z))
    return float(p_value)
