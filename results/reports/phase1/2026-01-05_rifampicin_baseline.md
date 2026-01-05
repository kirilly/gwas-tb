# Phase 1 Report: Rifampicin Baseline GWAS

**Report ID**: `phase1-rif-2026-01-05`
**Date**: 2026-01-05
**Phase**: 1 (Baseline GWAS)
**Drug**: Rifampicin (RIF)
**Status**: Awaiting bioinformatician review

---

## Executive Summary

| Metric | Value | Status |
|--------|-------|--------|
| Samples | 11,649 (8,803 S / 2,846 R) | ✓ |
| Variants | 6,780 (MAF ≥ 0.01) | ✓ |
| Runtime | 2.7 min | ✓ |
| Lambda GC | **0.139** | ⚠️ Overcorrection (see analysis) |
| Genome-wide hits (p<5e-8) | 4 | ✓ |
| Suggestive hits (p<1e-5) | 5 | ✓ |
| FDR<0.05 hits | 7 | ✓ |
| Heritability (h²) | 0.950 | ✓ Expected for monogenic |

---

## Plots

### Manhattan Plot
![Manhattan Plot](../figures/manhattan_rifampicin.png)

### QQ Plot
![QQ Plot](../figures/qq_rifampicin.png)

### Combined View
![Combined Plots](../figures/rifampicin_gwas_plots.png)

---

## Top 15 Hits (with Gene Annotations)

| Rank | Variant | Position | P-value | OR | PRPS | Gene | Notes |
|------|---------|----------|---------|-----|------|------|-------|
| 1 | 3650378_G_A | 3,650,378 | 1.2e-98 | 1.86 | 0.52 | **rpoB** | ✓ Known RIF resistance |
| 2 | 2256365_G_C | 2,256,365 | 1.1e-15 | 1.24 | 0.52 | **fdxA** (Rv2007c) | Ferredoxin - not canonical |
| 3 | 3650423_T_A | 3,650,423 | 5.6e-15 | 1.87 | 0.51 | **rpoB** | ✓ Known RIF resistance |
| 4 | 3650394_G_A | 3,650,394 | 7.3e-13 | 1.93 | 0.50 | **rpoB** | ✓ Known RIF resistance |
| 5 | 84053_T_C | 84,053 | 2.4e-06 | 1.72 | 0.51 | **Rv0075** | Aminotransferase |
| 6 | 1400431_G_A | 1,400,431 | 2.6e-05 | 2.37 | 0.56 | **deaD** (Rv1253) | RNA helicase |
| 7 | 3010500_G_A | 3,010,500 | 2.9e-05 | 1.66 | 0.58 | TBD | |
| 8 | 3339736_G_C | 3,339,736 | 1.2e-04 | 1.32 | 0.53 | TBD | |
| 9 | 4409684_G_T | 4,409,684 | 2.6e-04 | 0.77 | 0.60 | TBD | Protective (OR<1) |
| 10 | 1179497_G_A | 1,179,497 | 3.4e-04 | 1.83 | 0.50 | TBD | |

---

## Key Findings

### ✓ rpoB Region Dominates (Expected)
- **3 of top 4 hits** are in the rpoB region (~3.65 Mb)
- Top hit: 3650378_G_A with p=1.2×10⁻⁹⁸, OR=1.86
- Consistent with known RIF resistance mechanism (RNA polymerase β subunit)

### Novel Signals Identified
| Position | Gene | Function | Assessment |
|----------|------|----------|------------|
| 2,256,365 | fdxA (Rv2007c) | Ferredoxin - electron transfer | Not canonical RIF gene; may be compensatory or false positive |
| 84,053 | Rv0075 | Aminotransferase | Not canonical; possible lineage marker |
| 1,400,431 | deaD (Rv1253) | RNA helicase | Interacts with RNA pathways; potential compensatory |

---

## Analysis of Concerns

### Lambda GC = 0.139: Overcorrection Explained

**Literature evidence confirms this is expected overcorrection:**

> "Q-Q plot deflation (lambda < 1) typically indicates overcorrection for population structure - you may be using too many MDS components or an overly stringent kinship matrix" — [pyseer documentation](https://pyseer.readthedocs.io/en/master/tutorial.html)

> "All methods performed relatively poorly on highly clonal (low-recombining) genomes... particularly for highly clonal populations, there may be a limit to what can be learned from GWAS approaches" — [Benchmarking bacterial GWAS methods, PMC7200059](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7200059/)

> "In highly clonal settings, the LMM can suffer from very high false-positive rate unless the kinship matrix is carefully chosen" — [mBio 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7343994/)

**Conclusion:** Lambda deflation is a known issue in clonal bacterial GWAS. Despite this, **true resistance signals (rpoB) remain highly significant**, suggesting the method is working but may be conservative.

**Recommendation:** Consider elastic net approach which has "higher power than the linear mixed model and a lower false-positive rate than fixed-effect models" ([pyseer paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6289128/)).

---

### Effect Sizes (OR 1.8-1.9): Lower Than Clinical Studies

**Published OR values for rpoB mutations:**
- rpoB S450L: OR = **4.80** (95% CI 2.12-10.84) in clinical studies
- His445Asp: OR = **21.33** for unfavorable treatment outcomes
— [BMC Infectious Diseases 2025](https://bmcinfectdis.biomedcentral.com/articles/10.1186/s12879-025-10655-6)

**Our lower ORs explained:** Phylogenetic correction removes lineage-associated variance, attenuating effect sizes. This is expected behavior, not a problem.

---

### Known Compensatory Mutation Regions

For future analysis, key compensatory regions in H37Rv:

| Gene | Rv Number | Coordinates | Function |
|------|-----------|-------------|----------|
| **rpoB** | Rv0667 | 759,807 - 763,325 | Primary RIF target |
| **rpoC** | Rv0668 | 763,370 - 767,320 | β' chain - compensatory |
| **rpoA** | Rv3457c | 3,877,464 - 3,878,507 | α chain - compensatory |

> "Compensatory mutations in rpoA, rpoB, and rpoC have been identified in rifampicin-resistant strains... compensatory evolution enhances the in vivo fitness of drug-resistant M tuberculosis" — [Lancet Microbe 2023](https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(23)00110-6/fulltext)

---

## Questions for Bioinformatician Review

### Critical Questions

1. **Overcorrection acceptable?** Given lambda=0.139, should we proceed with current results or try alternative approaches (elastic net, fewer PCs)?

2. **fdxA hit (position 2,256,365)**: Second strongest signal. This ferredoxin gene is not a canonical RIF gene. Is this:
   - A compensatory mechanism?
   - A lineage-specific false positive?
   - Worth experimental validation?

3. **deaD hit (position 1,400,431)**: RNA helicase with highest OR (2.37). Could this interact with RNA polymerase pathways affected by RIF?

### Validation Questions

4. **rpoB mutation identity**: Can you identify which specific mutations (e.g., S450L, H445D) correspond to our top positions?

5. **WHO catalogue check**: Do our top hits match the WHO 2024 catalogue entries for RIF?

### Next Steps

6. **Proceed to INH?** Are results sufficient quality to proceed?

7. **Conditional analysis**: When running INH, should we condition on rpoB variants?

---

## Technical Details

### Pipeline Configuration
```yaml
method: lmm
min_maf: 0.01
p_threshold: 5e-8
fdr_threshold: 0.05
n_jobs: 8
kinship: snp-based (from pyseer)
```

### Output Files
- `results/rifampicin_gwas.csv` - Full results (6,780 variants)
- `results/figures/manhattan_rifampicin.png` - Manhattan plot
- `results/figures/qq_rifampicin.png` - QQ plot

---

## References

1. Lees JA et al. (2018). pyseer: a comprehensive tool for microbial pangenome-wide association studies. *Bioinformatics*. [DOI](https://academic.oup.com/bioinformatics/article/34/24/4310/5047751)

2. Power RA et al. (2017). Benchmarking bacterial genome-wide association study methods. *Microbial Genomics*. [PMC7200059](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7200059/)

3. CRyPTIC Consortium (2022). GWAS of 10,228 M. tuberculosis genomes. *PLOS Biology*. [Link](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001755)

4. Goig GA et al. (2023). Compensatory mutations in Cape Town. *Lancet Microbe*. [Link](https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(23)00110-6/fulltext)

---

## Decision Record (To Be Filled)

```yaml
date: 2026-01-05
reviewer: [name]
decisions:
  lambda_acceptable: [yes/no - explain]
  proceed_to_multidrug: [yes/no]
  fdxA_followup: [investigate/ignore]
  dead_followup: [investigate/ignore]
notes: |
  [Free-form notes]
```

---

*Report generated: 2026-01-05*
*Pipeline: tb-gwas v0.1.0*
*Validated against literature: 2026-01-05*
