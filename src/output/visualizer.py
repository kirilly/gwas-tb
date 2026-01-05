"""Visualization module for GWAS results."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

from src.config import OutputConfig
from src.gwas import GWASResult
from src.utils import get_logger

logger = get_logger()

# Set style
plt.style.use("seaborn-v0_8-whitegrid")


class Visualizer:
    """Generate GWAS visualizations."""

    def __init__(self, config: OutputConfig) -> None:
        """Initialize Visualizer.

        Args:
            config: Output configuration
        """
        self.config = config
        self.figures_dir = Path(config.figures_dir)
        self.figures_dir.mkdir(parents=True, exist_ok=True)

    def manhattan_plot(
        self,
        results: GWASResult,
        output_path: Path | None = None,
        highlight_genes: list[str] | None = None,
        threshold: float = 5e-8,
    ) -> Path:
        """Generate Manhattan plot.

        Args:
            results: GWAS results
            output_path: Output path for figure
            highlight_genes: Genes to highlight
            threshold: Significance threshold line

        Returns:
            Path to saved figure
        """
        if output_path is None:
            output_path = self.figures_dir / f"{results.drug}_manhattan.png"

        df = results.variants.copy()
        df["-log10(p)"] = -np.log10(df["p_value"].clip(lower=1e-300))

        # Extract position from variant ID if available
        # Assuming format like "pos_123456" or just numeric
        df["position"] = range(len(df))

        fig, ax = plt.subplots(figsize=(14, 6))

        # Plot points
        colors = ["#1f77b4", "#aec7e8"]
        ax.scatter(
            df["position"],
            df["-log10(p)"],
            c=[colors[i % 2] for i in range(len(df))],
            s=10,
            alpha=0.7,
        )

        # Significance threshold
        threshold_y = -np.log10(threshold)
        ax.axhline(y=threshold_y, color="red", linestyle="--", label=f"p = {threshold:.0e}")

        # Suggestive threshold
        suggestive_y = -np.log10(1e-5)
        ax.axhline(y=suggestive_y, color="orange", linestyle=":", alpha=0.7)

        # Highlight significant variants
        significant = df[df["p_value"] < threshold]
        if len(significant) > 0:
            ax.scatter(
                significant["position"],
                significant["-log10(p)"],
                c="red",
                s=30,
                zorder=5,
                label=f"Significant (n={len(significant)})",
            )

        ax.set_xlabel("Variant Position")
        ax.set_ylabel("-log10(p-value)")
        ax.set_title(f"Manhattan Plot - {results.drug}")
        ax.legend(loc="upper right")

        # Add lambda GC annotation
        ax.text(
            0.02, 0.98,
            f"Lambda GC = {results.lambda_gc:.3f}",
            transform=ax.transAxes,
            fontsize=10,
            verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
        )

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()

        logger.info(f"Manhattan plot saved to {output_path}")
        return output_path

    def qq_plot(
        self,
        results: GWASResult,
        output_path: Path | None = None,
    ) -> Path:
        """Generate QQ plot of p-values.

        Args:
            results: GWAS results
            output_path: Output path for figure

        Returns:
            Path to saved figure
        """
        if output_path is None:
            output_path = self.figures_dir / f"{results.drug}_qq.png"

        p_values = results.variants["p_value"].dropna().values
        p_values = p_values[p_values > 0]

        n = len(p_values)
        expected = -np.log10(np.arange(1, n + 1) / (n + 1))
        observed = -np.log10(np.sort(p_values))

        fig, ax = plt.subplots(figsize=(8, 8))

        # Plot points
        ax.scatter(expected, observed, s=10, alpha=0.7, c="#1f77b4")

        # Diagonal line
        max_val = max(expected.max(), observed.max())
        ax.plot([0, max_val], [0, max_val], "r--", label="Expected")

        # Confidence bands (95%)
        ci_lower = stats.beta.ppf(0.025, np.arange(1, n + 1), n - np.arange(1, n + 1) + 1)
        ci_upper = stats.beta.ppf(0.975, np.arange(1, n + 1), n - np.arange(1, n + 1) + 1)
        ax.fill_between(
            expected,
            -np.log10(ci_upper),
            -np.log10(ci_lower),
            alpha=0.2,
            color="gray",
            label="95% CI",
        )

        ax.set_xlabel("Expected -log10(p)")
        ax.set_ylabel("Observed -log10(p)")
        ax.set_title(f"QQ Plot - {results.drug}")

        # Lambda annotation
        ax.text(
            0.05, 0.95,
            f"Lambda GC = {results.lambda_gc:.3f}",
            transform=ax.transAxes,
            fontsize=12,
            verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
        )

        ax.legend(loc="lower right")
        ax.set_aspect("equal", adjustable="box")

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()

        logger.info(f"QQ plot saved to {output_path}")
        return output_path

    def pca_plot(
        self,
        snps: pd.DataFrame,
        labels: pd.Series | None = None,
        output_path: Path | None = None,
        n_components: int = 2,
    ) -> Path:
        """Generate PCA plot of samples.

        Args:
            snps: SNP matrix
            labels: Optional sample labels for coloring
            output_path: Output path for figure
            n_components: Number of PC components to show

        Returns:
            Path to saved figure
        """
        from sklearn.decomposition import PCA
        from sklearn.preprocessing import StandardScaler

        if output_path is None:
            output_path = self.figures_dir / "pca.png"

        # Prepare data
        X = snps.replace(-1, np.nan).fillna(snps.replace(-1, np.nan).mean()).values
        X = StandardScaler().fit_transform(X)

        # PCA
        pca = PCA(n_components=n_components)
        pcs = pca.fit_transform(X)

        fig, ax = plt.subplots(figsize=(10, 8))

        if labels is not None:
            labels_aligned = labels.loc[snps.index]
            scatter = ax.scatter(
                pcs[:, 0], pcs[:, 1],
                c=labels_aligned.values,
                cmap="coolwarm",
                s=30,
                alpha=0.7,
            )
            plt.colorbar(scatter, ax=ax, label="Phenotype")
        else:
            ax.scatter(pcs[:, 0], pcs[:, 1], s=30, alpha=0.7)

        ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
        ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
        ax.set_title("PCA of Samples")

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()

        logger.info(f"PCA plot saved to {output_path}")
        return output_path

    def cross_drug_heatmap(
        self,
        results: dict[str, GWASResult],
        output_path: Path | None = None,
        top_n: int = 50,
    ) -> Path:
        """Generate heatmap of cross-drug associations.

        Args:
            results: Dict of drug -> GWASResult
            output_path: Output path for figure
            top_n: Number of top variants per drug

        Returns:
            Path to saved figure
        """
        if output_path is None:
            output_path = self.figures_dir / "cross_drug_heatmap.png"

        drugs = list(results.keys())

        # Get union of top variants
        all_variants = set()
        for drug, result in results.items():
            top = result.get_top_hits(top_n)
            all_variants.update(top["variant"].tolist())

        all_variants = sorted(all_variants)

        # Build matrix of -log10(p)
        matrix = np.zeros((len(all_variants), len(drugs)))
        for j, drug in enumerate(drugs):
            var_df = results[drug].variants.set_index("variant")
            for i, var in enumerate(all_variants):
                if var in var_df.index:
                    p = var_df.loc[var, "p_value"]
                    matrix[i, j] = -np.log10(max(p, 1e-300))
                else:
                    matrix[i, j] = 0

        fig, ax = plt.subplots(figsize=(10, max(8, len(all_variants) * 0.2)))

        sns.heatmap(
            matrix,
            xticklabels=drugs,
            yticklabels=all_variants if len(all_variants) <= 30 else False,
            cmap="YlOrRd",
            ax=ax,
            cbar_kws={"label": "-log10(p)"},
        )

        ax.set_title("Cross-Drug Association Heatmap")
        ax.set_xlabel("Drug")
        ax.set_ylabel("Variant")

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()

        logger.info(f"Cross-drug heatmap saved to {output_path}")
        return output_path

    def maf_distribution(
        self,
        maf: np.ndarray,
        output_path: Path | None = None,
    ) -> Path:
        """Plot MAF distribution.

        Args:
            maf: Array of MAF values
            output_path: Output path for figure

        Returns:
            Path to saved figure
        """
        if output_path is None:
            output_path = self.figures_dir / "maf_distribution.png"

        fig, ax = plt.subplots(figsize=(10, 6))

        ax.hist(maf, bins=50, edgecolor="black", alpha=0.7)
        ax.axvline(x=0.01, color="red", linestyle="--", label="MAF = 0.01")
        ax.axvline(x=0.05, color="orange", linestyle="--", label="MAF = 0.05")

        ax.set_xlabel("Minor Allele Frequency")
        ax.set_ylabel("Count")
        ax.set_title("MAF Distribution")
        ax.legend()

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()

        logger.info(f"MAF distribution plot saved to {output_path}")
        return output_path
