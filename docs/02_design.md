# TB GWAS: Design Specification

## Architecture Overview

### Design Principles

1. **Modularity**: Each component is independent and testable
2. **Configuration-driven**: Behavior controlled via YAML, not code changes
3. **Fail-fast**: Validate inputs early, provide clear error messages
4. **Cache-friendly**: Expensive computations cached automatically
5. **Reproducible**: Fixed seeds, versioned environment, logged parameters

### High-Level Architecture

```
┌─────────────────────────────────────────────────────────────────────────┐
│                              CLI / Runner                                │
│                         (scripts/run_pipeline.py)                        │
└─────────────────────────────────────────────────────────────────────────┘
                                     │
                                     ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                           Pipeline Orchestrator                          │
│                          (src/pipeline/runner.py)                        │
│  - Loads config                                                          │
│  - Manages execution order                                               │
│  - Handles caching                                                       │
│  - Coordinates checkpoints                                               │
└─────────────────────────────────────────────────────────────────────────┘
                                     │
         ┌───────────────────────────┼───────────────────────────┐
         │                           │                           │
         ▼                           ▼                           ▼
┌─────────────────┐       ┌─────────────────┐       ┌─────────────────┐
│   Data Layer    │       │  Analysis Layer │       │  Output Layer   │
│                 │       │                 │       │                 │
│ • DataLoader    │──────▶│ • PhyloBuilder  │──────▶│ • Visualizer    │
│ • Preprocessor  │       │ • GWASRunner    │       │ • Reporter      │
│ • Validator     │       │ • FeatureSelect │       │ • Exporter      │
│                 │       │ • Epistasis     │       │                 │
└─────────────────┘       │ • Model         │       └─────────────────┘
                          └─────────────────┘
```

---

## Component Design

### 1. Configuration System

**Location**: `src/config/`

**Purpose**: Centralized configuration management

```yaml
# config/base.yaml
project:
  name: "tb-gwas"
  random_seed: 42
  n_jobs: 4
  log_level: "INFO"

data:
  snp_matrix: "data/processed/snps.h5"
  phenotypes: "data/processed/phenotypes.csv"
  who_catalogue: "data/external/who_catalogue_2024.xlsx"
  
qc:
  min_maf: 0.01
  max_missing_rate: 0.05
  min_samples_per_class: 20

phylogeny:
  method: "fasttree"  # or "iqtree"
  kinship_method: "tree"  # or "snp"
  n_threads: 4

gwas:
  method: "lmm"
  covariates: []
  p_threshold: 5e-8
  fdr_threshold: 0.05

drugs:
  order: ["RIF", "INH", "FQ", "STR", "EMB", "PZA"]
  conditional: true  # Include previous drugs' hits as covariates

abess:
  max_features: 50
  n_iterations: 2
  cv_folds: 5

epistasis:
  known_pairs:
    - ["rpoB", "rpoC"]
    - ["rpoB", "rpoA"]
    - ["katG", "ahpC"]
    - ["gyrA", "rpoC"]
  p_threshold: 0.001

output:
  results_dir: "results/"
  cache_dir: "cache/"
  figures_dir: "results/figures/"
```

**Interface**:

```python
@dataclass
class Config:
    project: ProjectConfig
    data: DataConfig
    qc: QCConfig
    phylogeny: PhyloConfig
    gwas: GWASConfig
    drugs: DrugConfig
    abess: ABESSConfig
    epistasis: EpistasisConfig
    output: OutputConfig
    
    @classmethod
    def from_yaml(cls, path: Path) -> "Config": ...
    
    def validate(self) -> List[str]:
        """Return list of validation errors, empty if valid."""
```

---

### 2. Data Layer

**Location**: `src/data/`

#### 2.1 DataLoader

**Purpose**: Load and validate input data

**Interface**:

```python
class DataLoader:
    def __init__(self, config: DataConfig): ...
    
    def load_snp_matrix(
        self, 
        path: Path = None,
        n_samples: int = None,  # For debugging
        n_variants: int = None
    ) -> pd.DataFrame:
        """
        Load SNP matrix.
        
        Returns:
            DataFrame with shape [n_samples, n_variants]
            Index: sample IDs
            Columns: variant IDs
            Values: 0, 1, 2 (genotype) or -1 (missing)
        
        Raises:
            DataValidationError: If data fails validation
        """
    
    def load_phenotypes(
        self,
        path: Path = None,
        drugs: List[str] = None
    ) -> pd.DataFrame:
        """
        Load phenotype data.
        
        Returns:
            DataFrame with shape [n_samples, n_drugs]
            Index: sample IDs
            Columns: drug names
            Values: 0 (susceptible), 1 (resistant), NaN (missing)
        """
    
    def load_who_catalogue(
        self,
        drug: str = None
    ) -> pd.DataFrame:
        """Load WHO mutation catalogue."""
```

#### 2.2 Preprocessor

**Purpose**: Quality control and filtering

**Interface**:

```python
class Preprocessor:
    def __init__(self, config: QCConfig): ...
    
    def run_qc(
        self,
        snps: pd.DataFrame,
        phenotypes: pd.DataFrame
    ) -> Tuple[pd.DataFrame, pd.DataFrame, QCReport]:
        """
        Run full QC pipeline.
        
        Steps:
        1. Filter samples by missingness
        2. Filter variants by missingness
        3. Filter variants by MAF
        4. Remove monomorphic variants
        5. Align samples between snps and phenotypes
        
        Returns:
            filtered_snps, filtered_phenotypes, qc_report
        """
    
    def filter_by_maf(
        self,
        snps: pd.DataFrame,
        min_maf: float = 0.01
    ) -> pd.DataFrame: ...
    
    def filter_by_missingness(
        self,
        snps: pd.DataFrame,
        max_missing: float = 0.05,
        axis: int = 0  # 0=samples, 1=variants
    ) -> pd.DataFrame: ...

@dataclass
class QCReport:
    n_samples_before: int
    n_samples_after: int
    n_variants_before: int
    n_variants_after: int
    removed_samples: List[str]
    removed_variants: List[str]
    maf_distribution: np.ndarray
```

---

### 3. Phylogeny Module

**Location**: `src/phylogeny/`

**Purpose**: Build phylogenetic tree and compute kinship matrix

**Interface**:

```python
class PhylogenyBuilder:
    def __init__(self, config: PhyloConfig): ...
    
    def build_tree(
        self,
        alignment_path: Path,
        output_path: Path = None
    ) -> Path:
        """
        Build phylogenetic tree.
        
        Args:
            alignment_path: Path to FASTA alignment
            output_path: Where to save tree (default: cache)
        
        Returns:
            Path to Newick tree file
        
        Methods:
            - fasttree: Fast, ~30 min for 500 samples
            - iqtree: Accurate, ~4 hours for 500 samples
        """
    
    def compute_kinship(
        self,
        tree_path: Path = None,
        snps: pd.DataFrame = None,
        method: str = None  # "tree" or "snp"
    ) -> np.ndarray:
        """
        Compute kinship matrix.
        
        Args:
            tree_path: Path to Newick tree (for tree-based kinship)
            snps: SNP matrix (for SNP-based kinship)
            method: Override config method
        
        Returns:
            Kinship matrix [n_samples, n_samples]
            Guaranteed positive semi-definite
        """
    
    def validate_kinship(
        self,
        kinship: np.ndarray
    ) -> bool:
        """Check kinship is valid (PSD, correct dimensions)."""
```

**Data Flow**:

```
alignment.fasta
       │
       ▼ (FastTree/IQ-TREE)
  phylogeny.nwk
       │
       ▼ (pyseer phylogeny_distance)
  kinship_matrix.npy
```

---

### 4. GWAS Module

**Location**: `src/gwas/`

#### 4.1 GWASRunner

**Purpose**: Run genome-wide association analysis

**Interface**:

```python
class GWASRunner:
    def __init__(self, config: GWASConfig): ...
    
    def run(
        self,
        snps: pd.DataFrame,
        phenotype: pd.Series,
        kinship: np.ndarray,
        covariates: pd.DataFrame = None
    ) -> GWASResult:
        """
        Run GWAS analysis.
        
        Args:
            snps: [n_samples, n_variants]
            phenotype: [n_samples] binary
            kinship: [n_samples, n_samples]
            covariates: [n_samples, n_covariates] optional
        
        Returns:
            GWASResult with per-variant statistics
        """
    
    def run_conditional(
        self,
        snps: pd.DataFrame,
        phenotypes: pd.DataFrame,
        kinship: np.ndarray,
        drug_order: List[str]
    ) -> Dict[str, GWASResult]:
        """
        Run conditional multi-drug GWAS.
        
        For each drug (in order):
        1. Run GWAS with previous drugs' top hits as covariates
        2. Filter by PRPS
        3. Store drug-specific results
        
        Returns:
            Dict mapping drug name to GWASResult
        """

@dataclass
class GWASResult:
    drug: str
    variants: pd.DataFrame  # variant, gene, p_value, OR, CI_low, CI_high, prps
    lambda_gc: float
    n_significant: int
    n_samples: int
    n_variants: int
    covariates_used: List[str]
```

#### 4.2 StatisticalTests

**Purpose**: Core statistical functions

```python
def calculate_lambda_gc(p_values: np.ndarray) -> float:
    """
    Calculate genomic inflation factor.
    
    Lambda = median(chi2_observed) / median(chi2_expected)
    
    Interpretation:
        < 1.0: Deflated (possibly underpowered)
        1.0: Perfect
        1.0-1.1: Acceptable
        > 1.1: Inflated (population structure not corrected)
        > 1.5: Severely inflated (results unreliable)
    """

def calculate_prps(
    snp: np.ndarray,
    kinship: np.ndarray
) -> float:
    """
    Calculate Phylogeny-Related Parallelism Score.
    
    High PRPS (>0.5): Variant follows phylogeny (possibly hitchhiking)
    Low PRPS (<0.3): Convergent evolution (likely under selection)
    """

def bonferroni_correction(
    p_values: np.ndarray,
    alpha: float = 0.05
) -> Tuple[np.ndarray, float]:
    """Return adjusted p-values and threshold."""

def fdr_correction(
    p_values: np.ndarray,
    alpha: float = 0.05,
    method: str = "bh"  # Benjamini-Hochberg
) -> Tuple[np.ndarray, np.ndarray]:
    """Return adjusted p-values and rejection mask."""
```

---

### 5. Feature Selection Module

**Location**: `src/feature_selection/`

**Purpose**: ABESS-based feature selection

**Interface**:

```python
class ABESSSelector:
    def __init__(self, config: ABESSConfig): ...
    
    def select(
        self,
        snps: pd.DataFrame,
        phenotype: pd.Series,
        n_iterations: int = 2
    ) -> ABESSResult:
        """
        Run iterative ABESS selection.
        
        Iteration 1: Find strong associations
        Iteration 2+: Remove explained samples, find weak associations
        
        Returns:
            ABESSResult with selected features and CV metrics
        """
    
    def cross_validate(
        self,
        snps: pd.DataFrame,
        phenotype: pd.Series,
        selected_features: List[str]
    ) -> CVMetrics:
        """Evaluate selected features via cross-validation."""

@dataclass
class ABESSResult:
    drug: str
    selected_features: List[str]
    coefficients: Dict[str, float]
    cv_auc: float
    cv_auc_std: float
    iterations: List[IterationResult]

@dataclass
class IterationResult:
    iteration: int
    n_samples: int
    n_resistant: int
    selected_features: List[str]
    explained_samples: int
```

---

### 6. Epistasis Module

**Location**: `src/epistasis/`

**Purpose**: Detect genetic interactions

**Interface**:

```python
class EpistasisAnalyzer:
    def __init__(self, config: EpistasisConfig): ...
    
    def test_known_pairs(
        self,
        snps: pd.DataFrame,
        phenotype: pd.Series,
        kinship: np.ndarray
    ) -> List[EpistasisResult]:
        """
        Test known epistatic pairs from config.
        
        For each pair:
        1. Get all SNPs in gene1 and gene2
        2. Test pairwise interactions
        3. Report significant interactions
        """
    
    def mediation_analysis(
        self,
        treatment: pd.Series,
        mediator: pd.Series,
        outcome: pd.Series
    ) -> MediationResult:
        """
        Test if mediator explains treatment effect.
        
        Example: Does rpoC mediate rpoB effect?
        """
    
    def detect_sign_epistasis(
        self,
        snp1: pd.Series,
        snp2: pd.Series,
        phenotype: pd.Series
    ) -> SignEpistasisResult:
        """
        Detect when two beneficial mutations are deleterious together.
        """

@dataclass
class EpistasisResult:
    gene1: str
    gene2: str
    snp1: str
    snp2: str
    p_interaction: float
    coef_interaction: float
    odds_ratio: float
    ci_lower: float
    ci_upper: float
    type: str  # "synergistic" or "antagonistic"

@dataclass
class MediationResult:
    treatment: str
    mediator: str
    total_effect: float
    direct_effect: float
    indirect_effect: float
    proportion_mediated: float
    p_sobel: float
    significant: bool

@dataclass
class SignEpistasisResult:
    snp1: str
    snp2: str
    effect1: float  # Effect of snp1 alone
    effect2: float  # Effect of snp2 alone
    effect_both: float  # Effect when both present
    expected_both: float  # Expected additive effect
    epistasis: float  # Deviation from additivity
    is_sign_epistasis: bool
```

---

### 7. Validation Module

**Location**: `src/validation/`

**Purpose**: Validate results against WHO catalogue

**Interface**:

```python
class WHOValidator:
    def __init__(self, catalogue_path: Path): ...
    
    def validate(
        self,
        results: GWASResult,
        drug: str,
        p_threshold: float = 5e-8
    ) -> ValidationReport:
        """
        Compare GWAS results with WHO catalogue.
        
        Returns:
            ValidationReport with TP, FP, FN, sensitivity, PPV
        """
    
    def cross_drug_contamination(
        self,
        results: Dict[str, GWASResult],
        p_threshold: float = 5e-8
    ) -> ContaminationReport:
        """
        Check for cross-drug associations.
        
        Reports cases where:
        - RIF-specific genes appear in INH results
        - INH-specific genes appear in RIF results
        - etc.
        """

@dataclass
class ValidationReport:
    drug: str
    n_significant: int
    who_genes: List[str]
    found_genes: List[str]
    true_positives: List[str]
    false_positives: List[str]
    false_negatives: List[str]
    sensitivity: float
    ppv: float
    f1_score: float
```

---

### 8. Output Module

**Location**: `src/output/`

#### 8.1 Visualizer

```python
class Visualizer:
    def __init__(self, config: OutputConfig): ...
    
    def manhattan_plot(
        self,
        results: GWASResult,
        output_path: Path,
        highlight_genes: List[str] = None
    ) -> Path: ...
    
    def qq_plot(
        self,
        results: GWASResult,
        output_path: Path
    ) -> Path: ...
    
    def pca_plot(
        self,
        snps: pd.DataFrame,
        labels: pd.Series,
        output_path: Path
    ) -> Path: ...
    
    def cross_drug_heatmap(
        self,
        results: Dict[str, GWASResult],
        output_path: Path
    ) -> Path: ...
```

#### 8.2 Reporter

```python
class Reporter:
    def __init__(self, config: OutputConfig): ...
    
    def generate_report(
        self,
        gwas_results: Dict[str, GWASResult],
        validation_reports: Dict[str, ValidationReport],
        epistasis_results: List[EpistasisResult],
        abess_results: Dict[str, ABESSResult]
    ) -> Report:
        """Generate comprehensive analysis report."""
    
    def export_results(
        self,
        report: Report,
        format: str = "all"  # "csv", "json", "html", "all"
    ) -> List[Path]: ...
```

---

### 9. Pipeline Orchestrator

**Location**: `src/pipeline/`

**Purpose**: Coordinate all components

```python
class PipelineRunner:
    def __init__(self, config_path: Path): ...
    
    def run(
        self,
        phases: List[str] = None  # ["data", "phylo", "gwas", "abess", "epistasis"]
    ) -> PipelineResult:
        """
        Run full or partial pipeline.
        
        Features:
        - Automatic caching of intermediate results
        - Resume from checkpoint
        - Parallel execution where possible
        - Progress logging
        """
    
    def run_from_checkpoint(
        self,
        checkpoint: str
    ) -> PipelineResult:
        """Resume from saved checkpoint."""
    
    def validate_checkpoint(
        self,
        checkpoint: str,
        bioinformatician_approved: bool = False
    ) -> bool:
        """
        Validate checkpoint results.
        
        If bioinformatician_approved=False, only automated checks.
        If True, record approval and proceed.
        """
```

---

## Caching Strategy

### Cache Structure

```
cache/
├── data/
│   ├── snps_filtered_{hash}.h5
│   └── phenotypes_filtered_{hash}.csv
├── phylogeny/
│   ├── tree_{method}_{hash}.nwk
│   └── kinship_{method}_{hash}.npy
├── gwas/
│   ├── {drug}_gwas_{hash}.parquet
│   └── {drug}_conditional_{hash}.parquet
├── abess/
│   └── {drug}_features_{hash}.json
└── checkpoints/
    ├── checkpoint_phase1.pkl
    └── checkpoint_phase2.pkl
```

### Cache Key Generation

```python
def compute_cache_key(*args, **kwargs) -> str:
    """
    Compute deterministic cache key from inputs.
    
    Uses:
    - Config hash
    - Input data hash (fast, samples subset)
    - Function name
    - Parameters
    """
    import hashlib
    import json
    
    key_data = {
        "args": [str(a) for a in args],
        "kwargs": {k: str(v) for k, v in sorted(kwargs.items())}
    }
    key_str = json.dumps(key_data, sort_keys=True)
    return hashlib.md5(key_str.encode()).hexdigest()[:12]
```

### Cache Decorator

```python
def cached(cache_dir: str = "cache"):
    """Decorator for caching function results."""
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            key = compute_cache_key(func.__name__, *args, **kwargs)
            cache_path = Path(cache_dir) / f"{func.__name__}_{key}.pkl"
            
            if cache_path.exists():
                logger.info(f"Loading cached: {cache_path}")
                return pickle.load(open(cache_path, "rb"))
            
            result = func(*args, **kwargs)
            
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            pickle.dump(result, open(cache_path, "wb"))
            logger.info(f"Cached: {cache_path}")
            
            return result
        return wrapper
    return decorator
```

---

## Error Handling

### Exception Hierarchy

```python
class TBGWASError(Exception):
    """Base exception for all TB GWAS errors."""

class DataValidationError(TBGWASError):
    """Invalid input data."""

class ConfigurationError(TBGWASError):
    """Invalid configuration."""

class PhylogenyError(TBGWASError):
    """Error building phylogeny or kinship."""

class GWASError(TBGWASError):
    """Error in GWAS analysis."""

class InsufficientSamplesError(GWASError):
    """Not enough samples for analysis."""

class ConvergenceError(GWASError):
    """Model failed to converge."""
```

### Validation Points

```python
def validate_snp_matrix(snps: pd.DataFrame) -> List[str]:
    """
    Validate SNP matrix.
    
    Checks:
    - No duplicate sample IDs
    - No duplicate variant IDs
    - Values in {0, 1, 2, -1}
    - At least 50 samples
    - At least 100 variants
    - Missing rate < 50%
    
    Returns:
        List of error messages (empty if valid)
    """
```

---

## Logging Strategy

### Log Levels

| Level | Usage |
|-------|-------|
| DEBUG | Detailed execution trace |
| INFO | Progress, key metrics, cache hits |
| WARNING | Non-fatal issues (high lambda, low samples) |
| ERROR | Failures that stop execution |

### Log Format

```python
# Console: concise
"2026-01-05 10:15:32 | INFO | Loading SNP matrix: 500 samples × 50000 variants"

# File: detailed
"2026-01-05 10:15:32 | INFO | data.loader | Loading SNP matrix: 500 samples × 50000 variants | memory=2.1GB"
```

### Structured Logging for Metrics

```python
logger.info("gwas_complete", extra={
    "drug": "RIF",
    "n_significant": 25,
    "lambda_gc": 1.05,
    "runtime_seconds": 342,
    "memory_peak_gb": 4.2
})
```

---

## Testing Strategy

### Unit Tests

```python
# tests/test_gwas.py
def test_lambda_gc_calculation():
    """Test genomic inflation calculation."""
    # Null distribution should give lambda ≈ 1.0
    p_values = np.random.uniform(0, 1, 10000)
    lambda_gc = calculate_lambda_gc(p_values)
    assert 0.95 < lambda_gc < 1.05

def test_lambda_gc_inflated():
    """Test detection of inflated p-values."""
    # Artificially inflated
    p_values = np.random.uniform(0, 0.1, 10000)  # Too many small p-values
    lambda_gc = calculate_lambda_gc(p_values)
    assert lambda_gc > 1.5
```

### Integration Tests

```python
# tests/test_pipeline.py
def test_full_pipeline_small():
    """Test full pipeline on small dataset."""
    result = PipelineRunner("config/test.yaml").run()
    
    assert result.gwas_results["RIF"].lambda_gc < 1.2
    assert "rpoB" in result.validation_reports["RIF"].found_genes
```

### Fixtures

```python
# tests/conftest.py
@pytest.fixture
def small_snp_matrix():
    """Generate small SNP matrix for testing."""
    np.random.seed(42)
    return pd.DataFrame(
        np.random.randint(0, 3, (100, 500)),
        index=[f"sample_{i}" for i in range(100)],
        columns=[f"var_{i}" for i in range(500)]
    )

@pytest.fixture
def mock_phenotype():
    """Generate mock phenotype."""
    np.random.seed(42)
    return pd.Series(
        np.random.randint(0, 2, 100),
        index=[f"sample_{i}" for i in range(100)],
        name="RIF"
    )
```

---

## File Structure

```
tb-gwas/
├── pixi.toml                      # Environment
├── pyproject.toml                 # Package metadata
├── README.md
│
├── config/
│   ├── base.yaml                  # Default config
│   ├── debug.yaml                 # Debug mode (small data)
│   └── production.yaml            # Full analysis
│
├── src/
│   ├── __init__.py
│   ├── config/
│   │   ├── __init__.py
│   │   └── config.py              # Config dataclasses
│   ├── data/
│   │   ├── __init__.py
│   │   ├── loader.py              # DataLoader
│   │   └── preprocessor.py        # Preprocessor
│   ├── phylogeny/
│   │   ├── __init__.py
│   │   └── builder.py             # PhylogenyBuilder
│   ├── gwas/
│   │   ├── __init__.py
│   │   ├── runner.py              # GWASRunner
│   │   └── stats.py               # Statistical functions
│   ├── feature_selection/
│   │   ├── __init__.py
│   │   └── abess.py               # ABESSSelector
│   ├── epistasis/
│   │   ├── __init__.py
│   │   ├── pairwise.py            # Pairwise tests
│   │   ├── mediation.py           # Mediation analysis
│   │   └── sign.py                # Sign epistasis
│   ├── validation/
│   │   ├── __init__.py
│   │   └── who.py                 # WHOValidator
│   ├── output/
│   │   ├── __init__.py
│   │   ├── visualizer.py
│   │   └── reporter.py
│   ├── pipeline/
│   │   ├── __init__.py
│   │   └── runner.py              # PipelineRunner
│   └── utils/
│       ├── __init__.py
│       ├── cache.py
│       ├── logging.py
│       └── reproducibility.py
│
├── scripts/
│   ├── run_pipeline.py            # Main entry point
│   └── validate_checkpoint.py     # Checkpoint validation
│
├── tests/
│   ├── conftest.py
│   ├── test_data/
│   ├── test_loader.py
│   ├── test_preprocessor.py
│   ├── test_gwas.py
│   ├── test_abess.py
│   └── test_pipeline.py
│
├── notebooks/
│   ├── 01_exploration.ipynb
│   └── 02_results_review.ipynb
│
└── data/
    ├── raw/                        # Original data (gitignored)
    ├── processed/                  # After QC (gitignored)
    └── external/                   # WHO catalogue, etc.
```

---

## Revision History

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | 2026-01-05 | Developer | Initial design |
