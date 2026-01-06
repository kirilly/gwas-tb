# Phase 1 Improvement Details

Method validation improvements for Phase 1 baseline GWAS, based on literature review.

**Related**: [03_tasks.md](03_tasks.md#task-17-validate-elastic-net-method)

---

## Summary

| Task | Description | Literature | Status |
|------|-------------|------------|--------|
| [1.7](#imp-1-alternative-gwas-methods) | Elastic net method validation | Saber & Shapiro 2020 | ⬜ |

---

## IMP-1: Alternative GWAS Methods

**Task**: 1.7
**Literature**: [Saber & Shapiro 2020](https://doi.org/10.1099/mgen.0.000337)

### Problem

Lambda GC = 0.139 in Phase 1 Rifampicin GWAS indicates **overcorrection** for population structure. This is a known issue in clonal bacterial GWAS.

### Literature Evidence

> "Q-Q plot deflation (lambda < 1) typically indicates overcorrection for population structure - you may be using too many MDS components or an overly stringent kinship matrix"
> — [pyseer documentation](https://pyseer.readthedocs.io/en/master/tutorial.html)

> "**All methods performed relatively poorly on highly clonal (low-recombining) genomes**... particularly for highly clonal populations, there may be a limit to what can be learned from GWAS approaches"
> — [Saber & Shapiro 2020](https://doi.org/10.1099/mgen.0.000337)

> "**The multi-locus elastic net approach was consistently amongst the highest-performing methods**, and had the highest power in detecting causal variants with both low and high effect sizes"
> — [Saber & Shapiro 2020](https://doi.org/10.1099/mgen.0.000337)

> "Elastic net has higher power than the linear mixed model and a lower false-positive rate than fixed-effect models in clonal populations"
> — [pyseer paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6289128/)

### Proposed Solutions

**Option 1: Elastic Net Method** (Recommended)

pyseer supports elastic net via `--wg enet`. Benefits:
- Higher power than LMM in clonal populations
- Lower FP rate than fixed-effect models
- Built-in feature selection (sparse coefficients)

**Option 2: Configurable PC Count**

Simpler approach - reduce number of principal components:

```python
# In _run_lmm():
if self.config.n_components is not None:
    n = self.config.n_components
    eigenvalues = eigenvalues[-n:]  # Top N (largest)
    eigenvectors = eigenvectors[:, -n:]
```

**Option 3: Ridge Regression**

Alternative regularization with adjustable strength.

### Proposed Code Changes

**File**: `src/config/config.py`

```python
@dataclass
class GWASConfig:
    """GWAS analysis configuration."""

    method: str = "lmm"  # Options: "lmm", "enet", "ridge", "simple"
    n_components: int | None = None  # None = use all, int = use top N eigencomponents
    alpha: float = 1.0  # Elastic net mixing (1.0 = LASSO, 0.0 = Ridge)
    l1_ratio: float = 0.5  # L1/L2 mixing for elastic net
    covariates: list[str] = field(default_factory=list)
    p_threshold: float = 5e-8
    fdr_threshold: float = 0.05
    min_maf: float = 0.01
```

**File**: `src/gwas/runner.py`

```python
def _run_elastic_net(
    self,
    snps: pd.DataFrame,
    phenotype: pd.Series,
    kinship: np.ndarray,
    covariates: pd.DataFrame | None,
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    """Run elastic net GWAS using pyseer.

    Elastic net offers:
    - Higher power than LMM in clonal populations
    - Lower FP rate than fixed-effect models
    - Built-in feature selection (sparse coefficients)
    """
    from pyseer.enet import fit_enet

    # Implementation using pyseer's elastic net API
    # ...
```

### Testing Plan

1. Run GWAS with current LMM (baseline)
2. Run GWAS with elastic net
3. Run GWAS with reduced PCs (10, 20, 50)
4. Compare: Lambda GC, top hits overlap, rpoB p-value

### Expected Outcomes

| Method | Expected Lambda GC | Top Hit Preserved |
|--------|-------------------|-------------------|
| LMM (current) | 0.139 | Yes |
| Elastic net | 0.7-1.0 | Yes |
| LMM (10 PCs) | 0.5-0.8 | Yes |

---

## Literature References

### Primary Source

1. **Saber & Shapiro 2020** - BacGWASim: Benchmarking bacterial GWAS
   - DOI: [10.1099/mgen.0.000337](https://doi.org/10.1099/mgen.0.000337)
   - PMID: 32100713
   - Key findings:
     - Elastic net consistently outperformed single-locus models
     - All methods struggled with highly clonal genomes
     - Sample size of 2000+ needed for strong effect detection

### Supporting Sources

2. **Lees et al. 2018** - pyseer paper
   - DOI: [10.1093/bioinformatics/bty539](https://doi.org/10.1093/bioinformatics/bty539)
   - Key: pyseer implementation details

---

## Status

⬜ Awaiting bioinformatician feedback at CP1

---

## Related Improvements (Other Phases)

The following improvements have been moved to their respective phases:

- **Gene Annotation Module** → [Phase 2, Task 2.0](phase2_improvements.md#gene-annotation-module)
- **WHO Catalogue Integration** → [Phase 2, Task 2.4](03_tasks.md#task-24-implement-whovalidator)
- **Lineage-Stratified GWAS** → [Phase 2, Task 2.6](phase2_improvements.md#lineage-stratified-gwas)
- **Convergent Evolution (phyC)** → [Phase 4, Task 4.5](phase4_improvements.md#convergent-evolution-analysis)
- **PFAM Domain Aggregation** → [Phase 4, Task 4.6](phase4_improvements.md#pfam-domain-aggregation)
- **3D Structure Validation** → [Phase 4, Task 4.7](phase4_improvements.md#3d-structure-validation)
