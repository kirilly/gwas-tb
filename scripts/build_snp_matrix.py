#!/usr/bin/env python
"""Build SNP matrix from per-sample variant files."""

import argparse
import sys
from pathlib import Path
from collections import defaultdict
from typing import Optional

import numpy as np
import pandas as pd
from tqdm import tqdm

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build SNP matrix from variant files")
    parser.add_argument(
        "--variants-dir",
        type=Path,
        default=Path("data/raw/research-code/db/nucl_data"),
        help="Directory containing .variants files",
    )
    parser.add_argument(
        "--pheno-dir",
        type=Path,
        default=Path("data/raw/research-code/db/pheno"),
        help="Directory containing .pheno files",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/processed"),
        help="Output directory",
    )
    parser.add_argument(
        "--drug",
        type=str,
        default="Rifampicin",
        help="Drug to filter samples by (default: Rifampicin)",
    )
    parser.add_argument(
        "--max-samples",
        type=int,
        default=None,
        help="Maximum samples to process (for testing)",
    )
    parser.add_argument(
        "--min-maf",
        type=float,
        default=0.01,
        help="Minimum minor allele frequency",
    )
    parser.add_argument(
        "--snps-only",
        action="store_true",
        help="Only include SNPs (exclude indels)",
    )
    return parser.parse_args()


def load_phenotypes(pheno_dir: Path, drug: str) -> dict[str, int]:
    """Load phenotype data for a drug.

    Returns:
        Dict mapping sample_id -> phenotype (0 or 1)
    """
    pheno_file = pheno_dir / f"{drug}.pheno"
    if not pheno_file.exists():
        raise FileNotFoundError(f"Phenotype file not found: {pheno_file}")

    phenotypes = {}
    with open(pheno_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                sample_id = parts[0]
                pheno = int(parts[1])
                phenotypes[sample_id] = pheno

    print(f"Loaded {len(phenotypes)} phenotypes for {drug}")
    return phenotypes


def get_sample_files(variants_dir: Path, phenotypes: dict[str, int]) -> list[tuple[str, Path]]:
    """Get variant files for samples with phenotypes.

    Returns:
        List of (sample_id, file_path) tuples
    """
    sample_files = []
    variant_files_by_stem = {f.stem: f for f in variants_dir.glob("*.variants")}
    print(f"Found {len(variant_files_by_stem)} variant files")

    for sample_id in phenotypes:
        if sample_id in variant_files_by_stem:
            sample_files.append((sample_id, variant_files_by_stem[sample_id]))

    print(f"Found {len(sample_files)} samples with both variants and phenotypes")
    return sample_files


def collect_variants(
    sample_files: list[tuple[str, Path]],
    snps_only: bool = False,
    max_samples: Optional[int] = None,
) -> tuple[dict[str, set[str]], set[str]]:
    """First pass: collect all variants and which samples have them.

    Returns:
        Tuple of (variant_samples dict, all_samples set)
    """
    variant_samples: dict[str, set[str]] = defaultdict(set)
    all_samples = set()

    files_to_process = sample_files[:max_samples] if max_samples else sample_files

    print("Pass 1: Collecting variants...")
    for sample_id, variant_file in tqdm(files_to_process):
        all_samples.add(sample_id)

        with open(variant_file) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    pos, ref, alt = parts[0], parts[1], parts[2]

                    # Filter for SNPs only if requested
                    if snps_only and (len(ref) != 1 or len(alt) != 1):
                        continue

                    # Create variant ID
                    variant_id = f"{pos}_{ref}_{alt}"
                    variant_samples[variant_id].add(sample_id)

    print(f"Found {len(variant_samples)} unique variants in {len(all_samples)} samples")
    return variant_samples, all_samples


def filter_by_maf(
    variant_samples: dict[str, set[str]],
    n_samples: int,
    min_maf: float,
) -> list[str]:
    """Filter variants by minor allele frequency.

    Returns:
        List of variant IDs passing MAF filter
    """
    passing_variants = []

    for variant_id, samples in variant_samples.items():
        af = len(samples) / n_samples
        maf = min(af, 1 - af)
        if maf >= min_maf:
            passing_variants.append(variant_id)

    print(f"Variants passing MAF >= {min_maf}: {len(passing_variants)}")
    return passing_variants


def build_matrix(
    sample_files: list[tuple[str, Path]],
    variants: list[str],
    max_samples: Optional[int] = None,
) -> pd.DataFrame:
    """Build SNP matrix.

    Returns:
        DataFrame with samples as rows, variants as columns
    """
    variant_set = set(variants)
    variant_to_idx = {v: i for i, v in enumerate(variants)}

    files_to_process = sample_files[:max_samples] if max_samples else sample_files
    n_samples = len(files_to_process)
    n_variants = len(variants)

    print(f"Pass 2: Building {n_samples} x {n_variants} matrix...")

    # Initialize matrix with zeros
    matrix = np.zeros((n_samples, n_variants), dtype=np.int8)
    sample_ids = []

    for i, (sample_id, variant_file) in enumerate(tqdm(files_to_process)):
        sample_ids.append(sample_id)

        with open(variant_file) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    pos, ref, alt = parts[0], parts[1], parts[2]
                    variant_id = f"{pos}_{ref}_{alt}"

                    if variant_id in variant_to_idx:
                        matrix[i, variant_to_idx[variant_id]] = 1

    df = pd.DataFrame(matrix, index=sample_ids, columns=variants)
    return df


def main() -> int:
    args = parse_args()

    # Ensure output directory exists
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load phenotypes
    phenotypes = load_phenotypes(args.pheno_dir, args.drug)

    # Get sample files
    sample_files = get_sample_files(args.variants_dir, phenotypes)

    if len(sample_files) == 0:
        print("ERROR: No samples found with both variants and phenotypes")
        return 1

    # Collect all variants (first pass)
    variant_samples, all_samples = collect_variants(
        sample_files,
        snps_only=args.snps_only,
        max_samples=args.max_samples,
    )

    # Filter by MAF
    filtered_variants = filter_by_maf(
        variant_samples,
        len(all_samples),
        args.min_maf,
    )

    # Sort variants by position
    filtered_variants.sort(key=lambda x: int(x.split("_")[0]))

    # Build matrix (second pass)
    snp_matrix = build_matrix(
        sample_files,
        filtered_variants,
        max_samples=args.max_samples,
    )

    # Create phenotype DataFrame
    pheno_df = pd.DataFrame(
        {args.drug: [phenotypes[s] for s in snp_matrix.index]},
        index=snp_matrix.index,
    )

    # Save outputs
    snp_file = args.output_dir / "snps.csv"
    pheno_file = args.output_dir / "phenotypes.csv"

    print(f"Saving SNP matrix to {snp_file}...")
    snp_matrix.to_csv(snp_file)

    print(f"Saving phenotypes to {pheno_file}...")
    pheno_df.to_csv(pheno_file)

    # Print summary
    print("\n=== Summary ===")
    print(f"Samples: {snp_matrix.shape[0]}")
    print(f"Variants: {snp_matrix.shape[1]}")
    print(f"Phenotype distribution:")
    print(f"  Susceptible (0): {(pheno_df[args.drug] == 0).sum()}")
    print(f"  Resistant (1): {(pheno_df[args.drug] == 1).sum()}")
    print(f"Output files:")
    print(f"  {snp_file}")
    print(f"  {pheno_file}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
