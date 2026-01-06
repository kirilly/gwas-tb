# Phase 2 Improvement Details

Detailed specifications for Phase 2 multi-drug analysis improvements, based on literature review (2026-01-06).

**Related**: [03_tasks.md](03_tasks.md#phase-2-multi-drug-analysis)

---

## Summary

| Task | Description | Literature | Status |
|------|-------------|------------|--------|
| [2.0](#gene-annotation-module) | Gene annotation module | Phelan 2019 | ⬜ |
| [2.6](#lineage-stratified-gwas) | Lineage-stratified GWAS | Oppong 2019 | ⬜ |

---

## Gene Annotation Module

**Task**: 2.0
**Literature**: [Phelan et al. 2019](https://doi.org/10.1038/s41598-019-45566-5)

### Problem

GWAS results contain genomic positions (e.g., `3650378_G_A`) without gene context. Manual annotation is error-prone and slow.

### Required Data

| File | Source | Location |
|------|--------|----------|
| H37Rv GFF3 | NCBI RefSeq NC_000962.3 | `data/external/h37rv.gff3` |

**Download**:
```bash
wget -O data/external/h37rv.gff3.gz \
  "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.gff.gz"
gunzip data/external/h37rv.gff3.gz
```

### Key Gene Coordinates (H37Rv Reference)

| Gene | Rv Number | Start | End | Drug Target |
|------|-----------|-------|-----|-------------|
| rpoB | Rv0667 | 759,807 | 763,325 | Rifampicin |
| rpoC | Rv0668 | 763,370 | 767,320 | Compensatory |
| rpoA | Rv3457c | 3,877,464 | 3,878,507 | Compensatory |
| katG | Rv1908c | 2,153,889 | 2,156,111 | Isoniazid |
| inhA | Rv1484 | 1,674,202 | 1,675,011 | Isoniazid |
| gyrA | Rv0006 | 7,302 | 9,818 | Fluoroquinolones |
| gyrB | Rv0005 | 5,240 | 7,267 | Fluoroquinolones |
| embB | Rv3795 | 4,246,514 | 4,249,810 | Ethambutol |
| pncA | Rv2043c | 2,288,681 | 2,289,241 | Pyrazinamide |
| rrs | MTB000019 | 1,471,846 | 1,473,382 | Aminoglycosides |

### Implementation

**File**: `src/annotation/gene_mapper.py`

```python
"""Map genomic positions to M. tuberculosis H37Rv genes."""

from dataclasses import dataclass
from pathlib import Path
import pandas as pd

@dataclass
class GeneInfo:
    """Gene annotation information."""
    gene_name: str      # e.g., "rpoB"
    rv_number: str      # e.g., "Rv0667"
    start: int
    end: int
    strand: str         # "+" or "-"
    product: str        # e.g., "DNA-directed RNA polymerase subunit beta"

class GeneMapper:
    """Map genomic positions to H37Rv genes using IntervalTree."""

    def __init__(self, gff_path: str | Path) -> None:
        self.genes = self._parse_gff(gff_path)
        self._build_index()

    def annotate_variant(self, variant_id: str) -> GeneInfo | None:
        """Get gene info for variant ID (e.g., '3650378_G_A')."""
        position = int(variant_id.split("_")[0])
        return self.annotate_position(position)

    def annotate_results(self, results: pd.DataFrame) -> pd.DataFrame:
        """Add gene_name, rv_number, product columns to GWAS results."""
        ...
```

### Dependencies

```toml
# pixi.toml
[dependencies]
intervaltree = ">=3.1"
```

---

## Lineage-Stratified GWAS

**Task**: 2.6
**Literature**: [Oppong et al. 2019](https://doi.org/10.1186/s12864-019-5615-3)

### Problem

Combined GWAS across all lineages may miss lineage-specific associations or be confounded by population structure differences between lineages.

### Key Finding from Literature

> "Unique associations with XDR in lineage-specific analyses provide evidence of **diverging evolutionary trajectories** between lineages 2 and 4 in response to antimicrobial drug therapy"

### M. tuberculosis Lineage Distribution

| Lineage | Name | Geographic Distribution | % in Dataset |
|---------|------|------------------------|--------------|
| L1 | Indo-Oceanic/Manila | Philippines, Southeast Asia | ~10% |
| L2 | East Asian/Beijing | China, Russia, Central Asia | ~40% |
| L3 | East African/Indian | East Africa, Central Asia | ~5% |
| L4 | Euro-American | Global | ~40% |

### Implementation

```python
def run_stratified(
    self,
    snps: pd.DataFrame,
    phenotype: pd.Series,
    lineages: pd.Series,
    kinship: np.ndarray,
) -> dict[str, GWASResult]:
    """Run GWAS separately for each lineage.

    Args:
        lineages: Series mapping sample_id to lineage (L1, L2, L3, L4)

    Returns:
        Dict of lineage -> GWASResult
    """
    results = {}
    for lineage in lineages.unique():
        mask = lineages == lineage
        if mask.sum() < 100:  # Skip small lineages
            continue
        results[lineage] = self.run(
            snps.loc[mask],
            phenotype.loc[mask],
            kinship[np.ix_(mask, mask)],
        )
    return results
```

### Expected Outcomes

1. **L2-specific**: Beijing family may have unique resistance pathways
2. **L4-specific**: Euro-American lineage compensatory mutations
3. **Shared**: Core resistance genes (rpoB, katG) should appear in both

---

## Literature References

### Primary Sources

1. **Phelan et al. 2019** - Manila strain WGS
   - DOI: [10.1038/s41598-019-45566-5](https://doi.org/10.1038/s41598-019-45566-5)
   - PMID: 31243306
   - Key: Molecular barcode, drug resistance mutations

2. **Oppong et al. 2019** - Lineage-specific GWAS
   - DOI: [10.1186/s12864-019-5615-3](https://doi.org/10.1186/s12864-019-5615-3)
   - PMID: 30922221
   - Key: L2 vs L4 diverging trajectories, PhyC method

### Supporting Sources

3. **Cohen et al. 2019** - WGS for TB drug resistance (Review)
   - DOI: [10.1186/s13073-019-0660-8](https://doi.org/10.1186/s13073-019-0660-8)
   - PMID: 31345251
   - Key: Comprehensive review of resistance mechanisms

4. **Worakitchanon et al. 2024** - Large indels and virulence
   - DOI: [10.1016/j.chom.2024.10.004](https://doi.org/10.1016/j.chom.2024.10.004)
   - PMID: 39471821
   - Key: Non-canonical resistance genes, indel detection
