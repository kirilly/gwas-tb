# Phase 1 Code Improvements

Based on Phase 1 Rifampicin GWAS results (2026-01-05), these improvements have been identified for implementation after bioinformatician feedback.

**Related**: [03_tasks.md](03_tasks.md#potential-improvements-phase-1) | [Phase 1 Report](../results/reports/phase1/2026-01-05_rifampicin_baseline.md)

---

## Summary

| ID | Issue | Solution | Priority | Effort | When |
|----|-------|----------|----------|--------|------|
| [IMP-1](#imp-1-alternative-gwas-methods) | Lambda GC = 0.139 (overcorrection) | Add elastic net method | High | Medium | After CP1 |
| [IMP-2](#imp-2-gene-annotation-module) | Missing gene annotations | Create GeneMapper module | High | Low | Can do now |
| [IMP-3](#imp-3-who-catalogue-integration) | No WHO catalogue validation | Create WHOCatalogue module | Medium | Low | After IMP-2 |
| [IMP-4](#imp-4-effect-size-documentation) | Effect size attenuation | Documentation only | Low | Minimal | Done |

---

## IMP-1: Alternative GWAS Methods

### Problem

Lambda GC = 0.139 in Phase 1 Rifampicin GWAS indicates overcorrection for population structure. This is a known issue in clonal bacterial GWAS.

### Literature Evidence

> "Q-Q plot deflation (lambda < 1) typically indicates overcorrection for population structure - you may be using too many MDS components or an overly stringent kinship matrix"
> — [pyseer documentation](https://pyseer.readthedocs.io/en/master/tutorial.html)

> "All methods performed relatively poorly on highly clonal (low-recombining) genomes... particularly for highly clonal populations, there may be a limit to what can be learned from GWAS approaches"
> — [Benchmarking bacterial GWAS methods, PMC7200059](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7200059/)

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

### Status

⬜ Awaiting bioinformatician feedback at CP1

---

## IMP-2: Gene Annotation Module

### Problem

GWAS results show variant positions (e.g., `3650378_G_A`) without gene context. Manual annotation is required to identify genes like rpoB, fdxA, etc.

### Required Data

- **H37Rv GFF3**: NCBI RefSeq NC_000962.3 annotation
- **Download**: `ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.gff.gz`
- **Location**: `data/external/h37rv.gff3`

### Proposed Code

**File**: `src/annotation/__init__.py`

```python
from .gene_mapper import GeneMapper
from .who_catalogue import WHOCatalogue

__all__ = ["GeneMapper", "WHOCatalogue"]
```

**File**: `src/annotation/gene_mapper.py`

```python
"""Map genomic positions to M. tuberculosis H37Rv genes."""

from dataclasses import dataclass
from pathlib import Path

import pandas as pd


@dataclass
class GeneInfo:
    """Gene annotation information."""

    gene_name: str  # e.g., "rpoB"
    rv_number: str  # e.g., "Rv0667"
    start: int
    end: int
    strand: str  # "+" or "-"
    product: str  # e.g., "DNA-directed RNA polymerase subunit beta"
    locus_tag: str


class GeneMapper:
    """Map genomic positions to M. tuberculosis H37Rv genes."""

    def __init__(self, gff_path: str | Path) -> None:
        """Load H37Rv GFF3 file.

        Args:
            gff_path: Path to H37Rv GFF3 annotation file
        """
        self.genes = self._parse_gff(gff_path)
        self._build_index()

    def _parse_gff(self, path: Path) -> pd.DataFrame:
        """Parse GFF3 file into DataFrame."""
        # Parse GFF3 format:
        # seqid, source, type, start, end, score, strand, phase, attributes
        ...

    def _build_index(self) -> None:
        """Build interval index for fast position lookup."""
        from intervaltree import IntervalTree
        self._tree = IntervalTree()
        for _, gene in self.genes.iterrows():
            self._tree[gene.start:gene.end] = gene

    def annotate_position(self, position: int) -> GeneInfo | None:
        """Get gene info for a genomic position.

        Args:
            position: 1-based genomic coordinate

        Returns:
            GeneInfo if position is within a gene, None otherwise
        """
        overlaps = self._tree[position]
        if overlaps:
            gene = list(overlaps)[0].data
            return GeneInfo(
                gene_name=gene.gene_name,
                rv_number=gene.rv_number,
                start=gene.start,
                end=gene.end,
                strand=gene.strand,
                product=gene.product,
                locus_tag=gene.locus_tag,
            )
        return None

    def annotate_variant(self, variant_id: str) -> GeneInfo | None:
        """Get gene info for a variant ID (e.g., '3650378_G_A').

        Args:
            variant_id: Variant identifier with position prefix

        Returns:
            GeneInfo if position is within a gene, None otherwise
        """
        position = int(variant_id.split("_")[0])
        return self.annotate_position(position)

    def annotate_results(self, results: pd.DataFrame) -> pd.DataFrame:
        """Add gene annotations to GWAS results DataFrame.

        Args:
            results: DataFrame with 'variant' column

        Returns:
            DataFrame with added gene_name, rv_number, product columns
        """
        annotations = []
        for variant in results["variant"]:
            info = self.annotate_variant(variant)
            if info:
                annotations.append({
                    "gene_name": info.gene_name,
                    "rv_number": info.rv_number,
                    "product": info.product,
                })
            else:
                annotations.append({
                    "gene_name": "intergenic",
                    "rv_number": "",
                    "product": "",
                })

        anno_df = pd.DataFrame(annotations)
        return pd.concat([results.reset_index(drop=True), anno_df], axis=1)
```

### Dependencies

Add to `pixi.toml`:
```toml
[dependencies]
intervaltree = ">=3.1"
```

### Key H37Rv Gene Coordinates (Reference)

| Gene | Rv Number | Start | End | Function |
|------|-----------|-------|-----|----------|
| rpoB | Rv0667 | 759,807 | 763,325 | RNA polymerase β - RIF target |
| rpoC | Rv0668 | 763,370 | 767,320 | RNA polymerase β' - compensatory |
| rpoA | Rv3457c | 3,877,464 | 3,878,507 | RNA polymerase α - compensatory |
| katG | Rv1908c | 2,153,889 | 2,156,111 | Catalase-peroxidase - INH target |
| inhA | Rv1484 | 1,674,202 | 1,675,011 | Enoyl-ACP reductase - INH target |
| gyrA | Rv0006 | 7,302 | 9,818 | DNA gyrase A - FQ target |
| gyrB | Rv0005 | 5,240 | 7,267 | DNA gyrase B - FQ target |
| embB | Rv3795 | 4,246,514 | 4,249,810 | Arabinosyltransferase - EMB target |
| pncA | Rv2043c | 2,288,681 | 2,289,241 | Pyrazinamidase - PZA target |
| rrs | MTB000019 | 1,471,846 | 1,473,382 | 16S rRNA - STR/AMK target |

### Status

⬜ Can implement independently

---

## IMP-3: WHO Catalogue Integration

### Problem

Manual comparison of GWAS hits against known resistance mutations is tedious and error-prone.

### Required Data

- **WHO Catalogue 2024**: Excel format from WHO
- **Download**: [WHO TB mutation catalogue](https://www.who.int/publications/i/item/9789240082410)
- **Location**: `data/external/who_catalogue_2024.xlsx`

### Proposed Code

**File**: `src/annotation/who_catalogue.py`

```python
"""Compare variants against WHO 2024 mutation catalogue."""

from dataclasses import dataclass
from pathlib import Path

import pandas as pd


@dataclass
class WHOEntry:
    """WHO catalogue entry."""

    drug: str
    gene: str
    mutation: str  # e.g., "S450L"
    variant_type: str  # "SNP", "Deletion", etc.
    confidence: str  # "High", "Moderate", "Low", "Indeterminate"
    position: int | None  # Genomic position if available
    evidence: str  # Evidence description


class WHOCatalogue:
    """WHO 2024 TB mutation catalogue lookup."""

    def __init__(self, catalogue_path: str | Path) -> None:
        """Load WHO catalogue.

        Args:
            catalogue_path: Path to WHO catalogue Excel file
        """
        self.catalogue = self._load_catalogue(catalogue_path)
        self._build_lookup()

    def _load_catalogue(self, path: Path) -> pd.DataFrame:
        """Load and parse WHO Excel catalogue."""
        # Parse the WHO format (multiple sheets per drug)
        ...

    def _build_lookup(self) -> None:
        """Build lookup dictionaries for fast access."""
        self._by_gene = {}  # gene -> list of mutations
        self._by_position = {}  # position -> entry
        ...

    def lookup_position(self, position: int, drug: str | None = None) -> list[WHOEntry]:
        """Find WHO entries for a genomic position."""
        ...

    def lookup_gene_mutation(self, gene: str, mutation: str, drug: str | None = None) -> WHOEntry | None:
        """Find WHO entry for a specific gene/mutation."""
        ...

    def validate_gwas_results(
        self,
        results: pd.DataFrame,
        gene_mapper: "GeneMapper",
        drug: str,
    ) -> pd.DataFrame:
        """Add WHO validation columns to GWAS results."""
        ...

    def calculate_sensitivity_ppv(
        self,
        results: pd.DataFrame,
        p_threshold: float = 5e-8,
    ) -> dict:
        """Calculate sensitivity and PPV against WHO catalogue."""
        ...
```

### Status

⬜ Depends on IMP-2 (gene mapper)

---

## IMP-4: Effect Size Documentation

### Problem

OR values (1.8-1.9) are lower than clinical studies (4.8).

### Solution

This is expected behavior, not a bug. Phylogenetic correction removes lineage-associated variance, attenuating effect sizes.

Already documented in [Phase 1 Report](../results/reports/phase1/2026-01-05_rifampicin_baseline.md#effect-sizes-or-18-19-lower-than-clinical-studies).

### Status

✅ Complete (documentation only)

---

## Future: Scientific Literature MCP Integration

For automated fact-checking against published literature:

| API | MCP Server | Config |
|-----|------------|--------|
| PubMed | `@ncukondo/pubmed-mcp` | `NCBI_EMAIL: 2flashback@gmail.com` |
| Semantic Scholar | `semantic-scholar-mcp` | API key required |
| OpenAlex | `openalex-research-mcp` | No auth required |

**Installation (when ready)**:

```json
// ~/.claude/.mcp.json
{
  "mcpServers": {
    "pubmed": {
      "command": "npx",
      "args": ["-y", "@ncukondo/pubmed-mcp"],
      "env": {
        "NCBI_EMAIL": "2flashback@gmail.com"
      }
    }
  }
}
```

### Status

⬜ Not implemented - for future enhancement
