"""Data loading module."""

from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from src.config import DataConfig
from src.utils import get_logger

logger = get_logger()


class DataValidationError(Exception):
    """Exception raised for data validation errors."""

    pass


class DataLoader:
    """Load and validate input data for GWAS analysis."""

    def __init__(self, config: DataConfig) -> None:
        """Initialize DataLoader.

        Args:
            config: Data configuration
        """
        self.config = config

    def load_snp_matrix(
        self,
        path: Path | None = None,
        n_samples: int | None = None,
        n_variants: int | None = None,
    ) -> pd.DataFrame:
        """Load SNP matrix from file.

        Args:
            path: Path to SNP matrix (HDF5 or CSV). Uses config if None.
            n_samples: Limit number of samples (for debugging)
            n_variants: Limit number of variants (for debugging)

        Returns:
            DataFrame with shape [n_samples, n_variants]
            Index: sample IDs
            Columns: variant IDs
            Values: 0, 1, 2 (genotype) or -1 (missing)

        Raises:
            DataValidationError: If data fails validation
        """
        path = Path(path) if path else Path(self.config.snp_matrix)
        logger.info(f"Loading SNP matrix from {path}")

        if path.suffix in (".h5", ".hdf5"):
            snps = self._load_hdf5(path, n_samples, n_variants)
        elif path.suffix == ".csv":
            snps = self._load_csv(path, n_samples, n_variants)
        else:
            raise DataValidationError(f"Unsupported file format: {path.suffix}")

        errors = self._validate_snp_matrix(snps)
        if errors:
            raise DataValidationError(f"SNP matrix validation failed: {errors}")

        logger.info(f"Loaded SNP matrix: {snps.shape[0]} samples x {snps.shape[1]} variants")
        return snps

    def _load_hdf5(
        self,
        path: Path,
        n_samples: int | None,
        n_variants: int | None,
    ) -> pd.DataFrame:
        """Load SNP matrix from HDF5 file."""
        with h5py.File(path, "r") as f:
            # Get dataset names
            if "snps" in f:
                data_key = "snps"
            elif "genotypes" in f:
                data_key = "genotypes"
            else:
                data_key = list(f.keys())[0]

            dataset = f[data_key]

            # Get dimensions
            n_s = n_samples if n_samples else dataset.shape[0]
            n_v = n_variants if n_variants else dataset.shape[1]

            # Load data
            data = dataset[:n_s, :n_v]

            # Get sample and variant IDs if available
            if "sample_ids" in f:
                sample_ids = [
                    s.decode() if isinstance(s, bytes) else s for s in f["sample_ids"][:n_s]
                ]
            else:
                sample_ids = [f"sample_{i}" for i in range(n_s)]

            if "variant_ids" in f:
                variant_ids = [
                    v.decode() if isinstance(v, bytes) else v for v in f["variant_ids"][:n_v]
                ]
            else:
                variant_ids = [f"var_{i}" for i in range(n_v)]

        return pd.DataFrame(data, index=sample_ids, columns=variant_ids)

    def _load_csv(
        self,
        path: Path,
        n_samples: int | None,
        n_variants: int | None,
    ) -> pd.DataFrame:
        """Load SNP matrix from CSV file."""
        df = pd.read_csv(path, index_col=0, nrows=n_samples)
        if n_variants:
            df = df.iloc[:, :n_variants]
        return df

    def _validate_snp_matrix(self, snps: pd.DataFrame) -> list[str]:
        """Validate SNP matrix.

        Checks:
        - No duplicate sample IDs
        - No duplicate variant IDs
        - Values in {-1, 0, 1, 2}
        - At least 50 samples
        - At least 100 variants
        """
        errors = []

        if snps.index.duplicated().any():
            errors.append("Duplicate sample IDs found")

        if snps.columns.duplicated().any():
            errors.append("Duplicate variant IDs found")

        unique_values = set(snps.values.flatten())
        invalid_values = unique_values - {-1, 0, 1, 2, np.nan}
        # Allow NaN and convert to -1
        valid_set = {-1, 0, 1, 2}
        for v in unique_values:
            if pd.isna(v):
                continue
            if v not in valid_set:
                errors.append(f"Invalid genotype values found: {invalid_values}")
                break

        if snps.shape[0] < 50:
            errors.append(f"Insufficient samples: {snps.shape[0]} < 50")

        if snps.shape[1] < 100:
            errors.append(f"Insufficient variants: {snps.shape[1]} < 100")

        return errors

    def load_phenotypes(
        self,
        path: Path | None = None,
        drugs: list[str] | None = None,
    ) -> pd.DataFrame:
        """Load phenotype data.

        Args:
            path: Path to phenotype file. Uses config if None.
            drugs: List of drugs to load. Loads all if None.

        Returns:
            DataFrame with shape [n_samples, n_drugs]
            Index: sample IDs
            Columns: drug names
            Values: 0 (susceptible), 1 (resistant), NaN (missing)

        Raises:
            DataValidationError: If data fails validation
        """
        path = Path(path) if path else Path(self.config.phenotypes)
        logger.info(f"Loading phenotypes from {path}")

        df = pd.read_csv(path, index_col=0)

        if drugs:
            missing = set(drugs) - set(df.columns)
            if missing:
                raise DataValidationError(f"Missing drug columns: {missing}")
            df = df[drugs]

        errors = self._validate_phenotypes(df)
        if errors:
            raise DataValidationError(f"Phenotype validation failed: {errors}")

        logger.info(f"Loaded phenotypes: {df.shape[0]} samples x {df.shape[1]} drugs")
        return df

    def _validate_phenotypes(self, phenotypes: pd.DataFrame) -> list[str]:
        """Validate phenotype data."""
        errors = []

        if phenotypes.index.duplicated().any():
            errors.append("Duplicate sample IDs found")

        # Check values are 0, 1, or NaN
        for col in phenotypes.columns:
            unique = set(phenotypes[col].dropna().unique())
            if not unique.issubset({0, 1, 0.0, 1.0}):
                errors.append(f"Invalid phenotype values in {col}: {unique}")

        return errors

    def load_who_catalogue(
        self,
        drug: str | None = None,
    ) -> pd.DataFrame:
        """Load WHO mutation catalogue.

        Args:
            drug: Filter to specific drug. Returns all if None.

        Returns:
            DataFrame with WHO catalogue data
        """
        path = Path(self.config.who_catalogue)
        logger.info(f"Loading WHO catalogue from {path}")

        if path.suffix in (".xlsx", ".xls"):
            df = pd.read_excel(path)
        else:
            df = pd.read_csv(path)

        if drug and "drug" in df.columns:
            df = df[df["drug"].str.upper() == drug.upper()]

        logger.info(f"Loaded WHO catalogue: {len(df)} entries")
        return df

    def align_samples(
        self,
        snps: pd.DataFrame,
        phenotypes: pd.DataFrame,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Align samples between SNP matrix and phenotypes.

        Args:
            snps: SNP matrix
            phenotypes: Phenotype data

        Returns:
            Tuple of aligned (snps, phenotypes) with common samples
        """
        common = snps.index.intersection(phenotypes.index)
        if len(common) == 0:
            raise DataValidationError("No common samples between SNPs and phenotypes")

        n_snp_only = len(snps.index) - len(common)
        n_pheno_only = len(phenotypes.index) - len(common)

        if n_snp_only > 0 or n_pheno_only > 0:
            logger.warning(f"Sample mismatch: {n_snp_only} SNP-only, {n_pheno_only} phenotype-only")

        logger.info(f"Aligned {len(common)} common samples")
        return snps.loc[common], phenotypes.loc[common]
