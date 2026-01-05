"""Configuration dataclasses and YAML loading."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml


@dataclass
class ProjectConfig:
    """Project-level configuration."""

    name: str = "tb-gwas"
    random_seed: int = 42
    n_jobs: int = 4
    log_level: str = "INFO"


@dataclass
class DataConfig:
    """Data paths configuration."""

    snp_matrix: str = "data/processed/snps.h5"
    phenotypes: str = "data/processed/phenotypes.csv"
    who_catalogue: str = "data/external/who_catalogue_2024.xlsx"
    alignment: str = "data/processed/alignment.fasta"
    tree: str = ""  # Optional pre-built tree


@dataclass
class QCConfig:
    """Quality control parameters."""

    min_maf: float = 0.01
    max_missing_rate: float = 0.05
    min_samples_per_class: int = 20


@dataclass
class PhyloConfig:
    """Phylogeny building configuration."""

    method: str = "fasttree"  # or "iqtree"
    kinship_method: str = "tree"  # or "snp"
    n_threads: int = 4


@dataclass
class GWASConfig:
    """GWAS analysis configuration."""

    method: str = "lmm"
    covariates: list[str] = field(default_factory=list)
    p_threshold: float = 5e-8
    fdr_threshold: float = 0.05


@dataclass
class DrugConfig:
    """Drug analysis configuration."""

    order: list[str] = field(
        default_factory=lambda: ["RIF", "INH", "FQ", "STR", "EMB", "PZA"]
    )
    conditional: bool = True


@dataclass
class ABESSConfig:
    """ABESS feature selection configuration."""

    max_features: int = 50
    n_iterations: int = 2
    cv_folds: int = 5


@dataclass
class EpistasisConfig:
    """Epistasis analysis configuration."""

    known_pairs: list[list[str]] = field(
        default_factory=lambda: [
            ["rpoB", "rpoC"],
            ["rpoB", "rpoA"],
            ["katG", "ahpC"],
            ["gyrA", "rpoC"],
        ]
    )
    p_threshold: float = 0.001


@dataclass
class OutputConfig:
    """Output paths configuration."""

    results_dir: str = "results/"
    cache_dir: str = "cache/"
    figures_dir: str = "results/figures/"
    logs_dir: str = "logs/"


@dataclass
class Config:
    """Main configuration container."""

    project: ProjectConfig = field(default_factory=ProjectConfig)
    data: DataConfig = field(default_factory=DataConfig)
    qc: QCConfig = field(default_factory=QCConfig)
    phylogeny: PhyloConfig = field(default_factory=PhyloConfig)
    gwas: GWASConfig = field(default_factory=GWASConfig)
    drugs: DrugConfig = field(default_factory=DrugConfig)
    abess: ABESSConfig = field(default_factory=ABESSConfig)
    epistasis: EpistasisConfig = field(default_factory=EpistasisConfig)
    output: OutputConfig = field(default_factory=OutputConfig)

    @classmethod
    def from_yaml(cls, path: Path | str) -> "Config":
        """Load configuration from YAML file."""
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"Config file not found: {path}")

        with open(path) as f:
            data = yaml.safe_load(f)

        return cls._from_dict(data)

    @classmethod
    def _from_dict(cls, data: dict[str, Any]) -> "Config":
        """Create Config from dictionary."""
        return cls(
            project=ProjectConfig(**data.get("project", {})),
            data=DataConfig(**data.get("data", {})),
            qc=QCConfig(**data.get("qc", {})),
            phylogeny=PhyloConfig(**data.get("phylogeny", {})),
            gwas=GWASConfig(**data.get("gwas", {})),
            drugs=DrugConfig(**data.get("drugs", {})),
            abess=ABESSConfig(**data.get("abess", {})),
            epistasis=EpistasisConfig(**data.get("epistasis", {})),
            output=OutputConfig(**data.get("output", {})),
        )

    def validate(self) -> list[str]:
        """Validate configuration, return list of errors."""
        errors = []

        # Check MAF threshold
        if not 0 < self.qc.min_maf < 0.5:
            errors.append(f"Invalid MAF threshold: {self.qc.min_maf}")

        # Check missing rate threshold
        if not 0 < self.qc.max_missing_rate < 1:
            errors.append(f"Invalid missing rate threshold: {self.qc.max_missing_rate}")

        # Check p-value thresholds
        if not 0 < self.gwas.p_threshold < 1:
            errors.append(f"Invalid p-value threshold: {self.gwas.p_threshold}")

        # Check phylogeny method
        if self.phylogeny.method not in ("fasttree", "iqtree"):
            errors.append(f"Invalid phylogeny method: {self.phylogeny.method}")

        # Check kinship method
        if self.phylogeny.kinship_method not in ("tree", "snp"):
            errors.append(f"Invalid kinship method: {self.phylogeny.kinship_method}")

        return errors

    def to_dict(self) -> dict[str, Any]:
        """Convert config to dictionary."""
        import dataclasses

        def to_dict_recursive(obj: Any) -> Any:
            if dataclasses.is_dataclass(obj):
                return {k: to_dict_recursive(v) for k, v in dataclasses.asdict(obj).items()}
            return obj

        return to_dict_recursive(self)
