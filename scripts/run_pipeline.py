#!/usr/bin/env python
"""Main entry point for TB GWAS pipeline."""

import argparse
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.config import Config
from src.utils import setup_logging, set_seeds, get_run_metadata


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="TB GWAS Pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--config", "-c",
        type=Path,
        default=Path("config/base.yaml"),
        help="Path to configuration file",
    )
    parser.add_argument(
        "--phase",
        choices=["data", "phylo", "gwas", "abess", "epistasis", "all"],
        default="all",
        help="Pipeline phase to run",
    )
    parser.add_argument(
        "--drug",
        type=str,
        default=None,
        help="Single drug to analyze (default: all in config)",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable debug mode",
    )
    return parser.parse_args()


def main() -> int:
    """Main entry point."""
    args = parse_args()

    # Load config
    config = Config.from_yaml(args.config)

    # Override for debug
    if args.debug:
        config.project.log_level = "DEBUG"

    # Setup logging
    logger = setup_logging(
        level=config.project.log_level,
        log_dir=config.output.logs_dir,
        name=config.project.name,
    )

    # Set seeds
    set_seeds(config.project.random_seed)
    logger.info(f"Random seed: {config.project.random_seed}")

    # Validate config
    errors = config.validate()
    if errors:
        logger.error(f"Configuration errors: {errors}")
        return 1

    # Log run metadata
    metadata = get_run_metadata(config.project.random_seed, str(args.config))
    logger.info(f"Run metadata: {metadata}")

    # Run pipeline
    logger.info(f"Starting pipeline phase: {args.phase}")

    try:
        if args.phase in ("data", "all"):
            run_data_phase(config, logger)

        if args.phase in ("phylo", "all"):
            run_phylo_phase(config, logger)

        if args.phase in ("gwas", "all"):
            run_gwas_phase(config, logger, drug=args.drug)

        logger.info("Pipeline completed successfully")
        return 0

    except Exception as e:
        logger.exception(f"Pipeline failed: {e}")
        return 1


def run_data_phase(config: Config, logger) -> None:
    """Run data loading and QC phase."""
    from src.data import DataLoader, Preprocessor

    logger.info("=== Data Phase ===")

    loader = DataLoader(config.data)
    preprocessor = Preprocessor(config.qc)

    # Load data
    snps = loader.load_snp_matrix()
    phenotypes = loader.load_phenotypes(drugs=config.drugs.order)

    # Run QC
    snps_filtered, phenotypes_filtered, qc_report = preprocessor.run_qc(snps, phenotypes)

    logger.info(qc_report.summary())

    # Save filtered data
    results_dir = Path(config.output.results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    snps_filtered.to_csv(results_dir / "snps_filtered.csv")
    phenotypes_filtered.to_csv(results_dir / "phenotypes_filtered.csv")

    logger.info(f"Filtered data saved to {results_dir}")


def run_phylo_phase(config: Config, logger) -> None:
    """Run phylogeny building phase."""
    from src.phylogeny import PhylogenyBuilder

    logger.info("=== Phylogeny Phase ===")

    builder = PhylogenyBuilder(config.phylogeny)

    # Build or load tree
    if config.data.tree:
        tree_path = Path(config.data.tree)
        logger.info(f"Using pre-built tree: {tree_path}")
    else:
        tree_path = builder.build_tree(
            Path(config.data.alignment),
            output_path=Path(config.output.cache_dir) / "tree.nwk",
        )

    # Compute kinship
    kinship = builder.compute_kinship(tree_path=tree_path)

    if not builder.validate_kinship(kinship):
        logger.warning("Kinship matrix validation failed")

    # Save kinship
    import numpy as np
    cache_dir = Path(config.output.cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    np.save(cache_dir / "kinship.npy", kinship)

    logger.info(f"Kinship matrix saved to {cache_dir / 'kinship.npy'}")


def run_gwas_phase(config: Config, logger, drug: str = None) -> None:
    """Run GWAS analysis phase."""
    import numpy as np
    import pandas as pd
    from src.gwas import GWASRunner
    from src.output import Visualizer

    logger.info("=== GWAS Phase ===")

    # Load filtered data
    results_dir = Path(config.output.results_dir)
    snps = pd.read_csv(results_dir / "snps_filtered.csv", index_col=0)
    phenotypes = pd.read_csv(results_dir / "phenotypes_filtered.csv", index_col=0)

    # Load kinship
    cache_dir = Path(config.output.cache_dir)
    kinship = np.load(cache_dir / "kinship.npy")

    runner = GWASRunner(config.gwas)
    visualizer = Visualizer(config.output)

    # Determine drugs to analyze
    drugs = [drug] if drug else config.drugs.order

    for drug_name in drugs:
        if drug_name not in phenotypes.columns:
            logger.warning(f"Drug {drug_name} not in phenotypes, skipping")
            continue

        logger.info(f"Running GWAS for {drug_name}")

        result = runner.run(
            snps=snps,
            phenotype=phenotypes[drug_name],
            kinship=kinship,
            drug_name=drug_name,
        )

        # Generate plots
        visualizer.manhattan_plot(result)
        visualizer.qq_plot(result)

        # Save results
        result.variants.to_csv(results_dir / f"{drug_name}_gwas.csv", index=False)
        logger.info(f"Results saved for {drug_name}")

        # Log summary
        logger.info(f"{drug_name} summary:")
        logger.info(f"  Lambda GC: {result.lambda_gc:.3f}")
        logger.info(f"  Significant: {result.n_significant}")
        logger.info(f"  Top hits: {result.get_top_hits(5)['variant'].tolist()}")


if __name__ == "__main__":
    sys.exit(main())
