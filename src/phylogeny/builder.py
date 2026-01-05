"""Phylogeny building and kinship computation."""

import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo

from src.config import PhyloConfig
from src.utils import get_logger

logger = get_logger()


class PhylogenyError(Exception):
    """Exception for phylogeny-related errors."""

    pass


class PhylogenyBuilder:
    """Build phylogenetic trees and compute kinship matrices."""

    def __init__(self, config: PhyloConfig) -> None:
        """Initialize PhylogenyBuilder.

        Args:
            config: Phylogeny configuration
        """
        self.config = config

    def build_tree(
        self,
        alignment_path: Path,
        output_path: Path | None = None,
    ) -> Path:
        """Build phylogenetic tree from alignment.

        Args:
            alignment_path: Path to FASTA alignment
            output_path: Where to save tree. Uses temp file if None.

        Returns:
            Path to Newick tree file

        Raises:
            PhylogenyError: If tree building fails
        """
        alignment_path = Path(alignment_path)
        if not alignment_path.exists():
            raise PhylogenyError(f"Alignment not found: {alignment_path}")

        if output_path is None:
            output_path = Path(tempfile.mktemp(suffix=".nwk"))

        logger.info(f"Building tree from {alignment_path} using {self.config.method}")

        if self.config.method == "fasttree":
            self._run_fasttree(alignment_path, output_path)
        elif self.config.method == "iqtree":
            self._run_iqtree(alignment_path, output_path)
        else:
            raise PhylogenyError(f"Unknown method: {self.config.method}")

        logger.info(f"Tree saved to {output_path}")
        return output_path

    def _run_fasttree(self, alignment: Path, output: Path) -> None:
        """Run FastTree."""
        cmd = ["fasttree", "-nt", "-gtr", "-fastest", str(alignment)]
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True,
            )
            output.write_text(result.stdout)
        except subprocess.CalledProcessError as e:
            raise PhylogenyError(f"FastTree failed: {e.stderr}")
        except FileNotFoundError:
            raise PhylogenyError("FastTree not found. Install via: pixi add fasttree")

    def _run_iqtree(self, alignment: Path, output: Path) -> None:
        """Run IQ-TREE."""
        cmd = [
            "iqtree",
            "-s", str(alignment),
            "-m", "GTR+G",
            "-nt", str(self.config.n_threads),
            "-pre", str(output.with_suffix("")),
            "-fast",
        ]
        try:
            subprocess.run(cmd, capture_output=True, text=True, check=True)
            # IQ-TREE outputs to .treefile
            treefile = output.with_suffix(".treefile")
            if treefile.exists():
                treefile.rename(output)
        except subprocess.CalledProcessError as e:
            raise PhylogenyError(f"IQ-TREE failed: {e.stderr}")
        except FileNotFoundError:
            raise PhylogenyError("IQ-TREE not found. Install via: pixi add iqtree")

    def compute_kinship(
        self,
        tree_path: Path | None = None,
        snps: pd.DataFrame | None = None,
        method: str | None = None,
    ) -> np.ndarray:
        """Compute kinship matrix.

        Args:
            tree_path: Path to Newick tree (for tree-based kinship)
            snps: SNP matrix (for SNP-based kinship)
            method: Override config method ("tree" or "snp")

        Returns:
            Kinship matrix [n_samples, n_samples]
            Guaranteed positive semi-definite

        Raises:
            PhylogenyError: If computation fails
        """
        method = method or self.config.kinship_method

        if method == "tree":
            if tree_path is None:
                raise PhylogenyError("Tree path required for tree-based kinship")
            kinship = self._kinship_from_tree(tree_path)
        elif method == "snp":
            if snps is None:
                raise PhylogenyError("SNP matrix required for SNP-based kinship")
            kinship = self._kinship_from_snps(snps)
        else:
            raise PhylogenyError(f"Unknown kinship method: {method}")

        # Ensure PSD
        kinship = self._ensure_psd(kinship)

        logger.info(f"Computed kinship matrix: {kinship.shape}")
        return kinship

    def _kinship_from_tree(self, tree_path: Path) -> np.ndarray:
        """Compute kinship from phylogenetic tree using patristic distances."""
        tree_path = Path(tree_path)
        if not tree_path.exists():
            raise PhylogenyError(f"Tree not found: {tree_path}")

        tree = Phylo.read(tree_path, "newick")

        # Get all terminal nodes (samples)
        terminals = tree.get_terminals()
        n = len(terminals)
        names = [t.name for t in terminals]

        logger.info(f"Computing patristic distances for {n} samples")

        # Compute pairwise distances
        distances = np.zeros((n, n))
        for i, t1 in enumerate(terminals):
            for j, t2 in enumerate(terminals):
                if i <= j:
                    d = tree.distance(t1, t2)
                    distances[i, j] = d
                    distances[j, i] = d

        # Convert distance to similarity (kinship)
        # Using exponential decay: K = exp(-d / median_d)
        median_d = np.median(distances[distances > 0])
        kinship = np.exp(-distances / median_d)

        return kinship

    def _kinship_from_snps(self, snps: pd.DataFrame) -> np.ndarray:
        """Compute kinship from SNP matrix (IBS-based).

        Uses the realized relationship matrix (GRM) approach.
        """
        # Replace missing values with mean
        X = snps.replace(-1, np.nan).values
        col_means = np.nanmean(X, axis=0)
        X = np.where(np.isnan(X), col_means, X)

        # Center and scale
        X_centered = X - X.mean(axis=0)
        std = X.std(axis=0)
        std[std == 0] = 1  # Avoid division by zero
        X_scaled = X_centered / std

        # Compute GRM: K = X @ X.T / p
        n, p = X_scaled.shape
        kinship = X_scaled @ X_scaled.T / p

        return kinship

    def _ensure_psd(self, matrix: np.ndarray, eps: float = 1e-6) -> np.ndarray:
        """Ensure matrix is positive semi-definite.

        Args:
            matrix: Square symmetric matrix
            eps: Small value to add to negative eigenvalues

        Returns:
            PSD matrix
        """
        # Eigendecomposition
        eigenvalues, eigenvectors = np.linalg.eigh(matrix)

        # Fix negative eigenvalues
        eigenvalues = np.maximum(eigenvalues, eps)

        # Reconstruct
        psd = eigenvectors @ np.diag(eigenvalues) @ eigenvectors.T

        # Ensure symmetry
        psd = (psd + psd.T) / 2

        return psd

    def validate_kinship(self, kinship: np.ndarray) -> bool:
        """Check if kinship matrix is valid.

        Args:
            kinship: Kinship matrix to validate

        Returns:
            True if valid, False otherwise
        """
        # Check square
        if kinship.ndim != 2 or kinship.shape[0] != kinship.shape[1]:
            logger.warning("Kinship matrix is not square")
            return False

        # Check symmetric
        if not np.allclose(kinship, kinship.T):
            logger.warning("Kinship matrix is not symmetric")
            return False

        # Check PSD (all eigenvalues >= 0)
        eigenvalues = np.linalg.eigvalsh(kinship)
        if np.any(eigenvalues < -1e-6):
            logger.warning(f"Kinship matrix is not PSD: min eigenvalue = {eigenvalues.min()}")
            return False

        return True

    def get_sample_order(self, tree_path: Path) -> list[str]:
        """Get sample order from tree.

        Args:
            tree_path: Path to Newick tree

        Returns:
            List of sample names in tree order
        """
        tree = Phylo.read(tree_path, "newick")
        return [t.name for t in tree.get_terminals()]
