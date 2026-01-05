"""Data preprocessing and QC module."""

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from src.config import QCConfig
from src.utils import get_logger

logger = get_logger()


@dataclass
class QCReport:
    """Quality control report."""

    n_samples_before: int
    n_samples_after: int
    n_variants_before: int
    n_variants_after: int
    removed_samples: list[str] = field(default_factory=list)
    removed_variants: list[str] = field(default_factory=list)
    maf_distribution: np.ndarray = field(default_factory=lambda: np.array([]))

    def summary(self) -> str:
        """Generate summary string."""
        return (
            f"QC Report:\n"
            f"  Samples: {self.n_samples_before} -> {self.n_samples_after} "
            f"(removed {self.n_samples_before - self.n_samples_after})\n"
            f"  Variants: {self.n_variants_before} -> {self.n_variants_after} "
            f"(removed {self.n_variants_before - self.n_variants_after})"
        )


class Preprocessor:
    """Data preprocessing and quality control."""

    def __init__(self, config: QCConfig) -> None:
        """Initialize Preprocessor.

        Args:
            config: QC configuration
        """
        self.config = config

    def run_qc(
        self,
        snps: pd.DataFrame,
        phenotypes: pd.DataFrame,
    ) -> tuple[pd.DataFrame, pd.DataFrame, QCReport]:
        """Run full QC pipeline.

        Steps:
        1. Filter samples by missingness
        2. Filter variants by missingness
        3. Filter variants by MAF
        4. Remove monomorphic variants
        5. Align samples between snps and phenotypes

        Args:
            snps: SNP matrix
            phenotypes: Phenotype data

        Returns:
            Tuple of (filtered_snps, filtered_phenotypes, qc_report)
        """
        n_samples_before = snps.shape[0]
        n_variants_before = snps.shape[1]
        removed_samples: list[str] = []
        removed_variants: list[str] = []

        logger.info("Starting QC pipeline")

        # Step 1: Filter samples by missingness
        snps, removed = self.filter_by_missingness(
            snps, max_missing=self.config.max_missing_rate, axis=0
        )
        removed_samples.extend(removed)

        # Step 2: Filter variants by missingness
        snps, removed = self.filter_by_missingness(
            snps, max_missing=self.config.max_missing_rate, axis=1
        )
        removed_variants.extend(removed)

        # Step 3: Filter variants by MAF
        snps, removed = self.filter_by_maf(snps, min_maf=self.config.min_maf)
        removed_variants.extend(removed)

        # Step 4: Remove monomorphic variants
        snps, removed = self.remove_monomorphic(snps)
        removed_variants.extend(removed)

        # Step 5: Align samples
        common_samples = snps.index.intersection(phenotypes.index)
        snps = snps.loc[common_samples]
        phenotypes = phenotypes.loc[common_samples]

        # Calculate MAF distribution for report
        maf = self.calculate_maf(snps)

        report = QCReport(
            n_samples_before=n_samples_before,
            n_samples_after=snps.shape[0],
            n_variants_before=n_variants_before,
            n_variants_after=snps.shape[1],
            removed_samples=removed_samples,
            removed_variants=removed_variants,
            maf_distribution=maf.values,
        )

        logger.info(report.summary())
        return snps, phenotypes, report

    def filter_by_maf(
        self,
        snps: pd.DataFrame,
        min_maf: float = 0.01,
    ) -> tuple[pd.DataFrame, list[str]]:
        """Filter variants by Minor Allele Frequency.

        Args:
            snps: SNP matrix
            min_maf: Minimum MAF threshold

        Returns:
            Tuple of (filtered_snps, removed_variant_ids)
        """
        maf = self.calculate_maf(snps)
        keep = (maf >= min_maf) & (maf <= 1 - min_maf)
        removed = snps.columns[~keep].tolist()

        logger.info(f"MAF filter: removed {len(removed)} variants (MAF < {min_maf})")
        return snps.loc[:, keep], removed

    def filter_by_missingness(
        self,
        snps: pd.DataFrame,
        max_missing: float = 0.05,
        axis: int = 0,
    ) -> tuple[pd.DataFrame, list[str]]:
        """Filter by missingness rate.

        Args:
            snps: SNP matrix
            max_missing: Maximum allowed missing rate
            axis: 0 for samples, 1 for variants

        Returns:
            Tuple of (filtered_snps, removed_ids)
        """
        # Calculate missing rate (count -1 and NaN as missing)
        missing_mask = (snps == -1) | snps.isna()
        missing_rate = missing_mask.mean(axis=1 - axis)

        keep = missing_rate <= max_missing
        if axis == 0:
            removed = snps.index[~keep].tolist()
            filtered = snps.loc[keep]
            logger.info(f"Sample missingness filter: removed {len(removed)} samples")
        else:
            removed = snps.columns[~keep].tolist()
            filtered = snps.loc[:, keep]
            logger.info(f"Variant missingness filter: removed {len(removed)} variants")

        return filtered, removed

    def remove_monomorphic(
        self,
        snps: pd.DataFrame,
    ) -> tuple[pd.DataFrame, list[str]]:
        """Remove monomorphic variants (no variation).

        Args:
            snps: SNP matrix

        Returns:
            Tuple of (filtered_snps, removed_variant_ids)
        """
        # A variant is monomorphic if all non-missing values are the same
        n_unique = snps.apply(lambda x: x[x >= 0].nunique(), axis=0)
        keep = n_unique > 1
        removed = snps.columns[~keep].tolist()

        logger.info(f"Removed {len(removed)} monomorphic variants")
        return snps.loc[:, keep], removed

    def calculate_maf(self, snps: pd.DataFrame) -> pd.Series:
        """Calculate Minor Allele Frequency for each variant.

        Args:
            snps: SNP matrix with values 0, 1, 2 (-1 for missing)

        Returns:
            Series of MAF values indexed by variant ID
        """
        # Replace -1 with NaN for calculation
        snps_clean = snps.replace(-1, np.nan)

        # Calculate allele frequency
        # For diploid (0/1/2): AF = sum / (2 * n) = mean / 2
        # For haploid (0/1): AF = mean
        max_val = snps_clean.max().max()
        if max_val > 1:
            # Diploid encoding (0/1/2)
            af = snps_clean.mean(axis=0) / 2
        else:
            # Haploid encoding (0/1)
            af = snps_clean.mean(axis=0)

        # MAF is the smaller of af and 1-af
        maf = pd.concat([af, 1 - af], axis=1).min(axis=1)
        return maf

    def encode_missing(
        self,
        snps: pd.DataFrame,
        strategy: str = "mean",
    ) -> pd.DataFrame:
        """Encode missing values.

        Args:
            snps: SNP matrix with -1 or NaN for missing
            strategy: Encoding strategy ("mean", "mode", "zero")

        Returns:
            SNP matrix with missing values encoded
        """
        snps = snps.replace(-1, np.nan)

        if strategy == "mean":
            snps = snps.fillna(snps.mean())
        elif strategy == "mode":
            snps = snps.fillna(snps.mode().iloc[0])
        elif strategy == "zero":
            snps = snps.fillna(0)
        else:
            raise ValueError(f"Unknown encoding strategy: {strategy}")

        return snps
