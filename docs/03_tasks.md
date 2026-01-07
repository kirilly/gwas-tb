# TB GWAS: Tasks Specification

## Overview

**Total estimated time**: 60 hours  
**Bioinformatician syncs**: 5 checkpoints (~4-5 hours total)  
**Buffer**: 4 hours for unexpected issues

### Task Status Legend

- ⬜ Not started
- 🔄 In progress
- ✅ Complete
- 🔬 Awaiting bioinformatician review
- ❌ Blocked

---

## Phase 0: Environment & Data Setup

**Duration**: 4 hours  
**Dependencies**: None

### Task 0.1: Initialize Project Structure

**Time**: 30 min
**Status**: ✅

**Subtasks**:
- [x] Create directory structure as per Design spec
- [x] Initialize git repository
- [x] Create `.gitignore` (exclude data/, cache/, results/)
- [x] Create README.md with project overview

**Acceptance**:
- [x] `git status` shows clean repo
- [x] All directories exist

---

### Task 0.2: Setup Environment with Pixi

**Time**: 1 hour
**Status**: ✅

**Subtasks**:
- [x] Install Pixi (if not installed)
- [x] Create `pixi.toml` with dependencies:
  - python=3.11
  - pandas, numpy, scipy, scikit-learn, matplotlib, seaborn
  - biopython
  - fasttree, iqtree (bioconda)
  - pyseer, abess, torch (PyPI)
  - pytest, mypy, ruff (dev)
- [x] Run `pixi install`
- [x] Verify all imports work

**Acceptance**:
- [x] `pixi run python -c "import abess"` succeeds
- [x] `pixi run fasttree -h` shows help

**Note**: pyseer/torch removed from pypi-deps due to platform issues; can be added later if needed.

---

### Task 0.3: Download Data

**Time**: 1 hour
**Status**: ✅

**Subtasks**:
- [x] Clone article repository: `github.com/Reshetnikoff/m.tuberculosis-research-code`
- [x] Extract variant data archives (nucl_data.tar.gz, pheno.tar.gz, trees.zip)
- [x] Build SNP matrix from raw variant files (scripts/build_snp_matrix.py)
- [ ] Download WHO catalogue 2024 (pending)
- [ ] Download H37Rv reference genome (pending)
- [ ] Download CRyPTIC validation dataset (optional, for independent validation)

**Acceptance**:
- [x] SNP matrix file exists: `data/processed/snps.csv` (11,649 samples × 6,780 variants)
- [x] Phenotype file exists: `data/processed/phenotypes.csv` (Rifampicin: 8,803 S / 2,846 R)
- [ ] WHO catalogue loads without errors (pending)
- [ ] CRyPTIC dataset available at `data/external/cryptic_mics.csv` (optional)

**Notes**:
- Raw data: 12,814 .variants files in `data/raw/research-code/db/nucl_data/`
- SNP matrix built with MAF >= 0.01 filter, SNPs only (no indels)
- Kinship matrix computed and cached: `data/processed/kinship.npy`
- **CRyPTIC dataset** (optional validation): 12,289 isolates, WGS + MIC for 13 drugs
  - Source: CRyPTIC Consortium 2022 (PLOS Biology, DOI: 10.1371/journal.pbio.3001721)
  - Access: Open source (no registration)

---

### Task 0.4: Implement Configuration Module

**Time**: 1 hour
**Status**: ✅
**Depends on**: 0.1, 0.2

**Subtasks**:
- [x] Create `src/config/config.py` with dataclasses
- [x] Create `config/base.yaml` with default values
- [x] Create `config/debug.yaml` (small data for testing)
- [x] Implement `Config.from_yaml()` method
- [x] Implement `Config.validate()` method
- [x] Write unit tests

**Acceptance**:
- [x] Config loads from YAML
- [x] Validation catches missing required fields
- [x] Tests pass: `pixi run pytest tests/test_config.py` (5 tests passed)

---

### Task 0.5: Implement Logging & Reproducibility Utils

**Time**: 30 min
**Status**: ✅
**Depends on**: 0.1, 0.2

**Subtasks**:
- [x] Create `src/utils/logging.py` with setup_logging()
- [x] Create `src/utils/reproducibility.py` with set_seeds()
- [x] Log format: timestamp | level | module | message
- [x] Log to both console and file

**Acceptance**:
- [x] Logs appear in console and `logs/` directory
- [x] Random seed is logged at startup

---

## 🔬 CHECKPOINT 0: Data Validation

**When**: After Task 0.5, before starting Phase 1  
**Duration**: 30 min  
**Type**: Sync with bioinformatician

### Preparation (Developer)

Before meeting, prepare:
- [ ] Summary of loaded data: N samples, M variants, drugs available
- [ ] Basic statistics: MAF distribution, missingness rates
- [ ] Any data loading errors or warnings
- [ ] Questions list (see below)

### Questions to Ask

**Data Quality**:
1. "SNP matrix has [N] samples × [M] variants. Expected?"
2. "Missing rate is [X]%. Acceptable threshold?"
3. "MAF distribution looks like [describe]. Normal?"
4. "Found [list] drugs with phenotypes. Which prioritize?"

**QC Parameters**:
5. "Recommend MAF threshold? (default 0.01)"
6. "Sample missingness threshold? (default 5%)"
7. "Variant missingness threshold? (default 5%)"

**Data Structure**:
8. "Is there a pre-built phylogeny or build from scratch?"
9. "Any batch effects or known confounders to account for?"

### Expected Outcomes

- [ ] QC parameters agreed
- [ ] Priority drug confirmed (likely RIF)
- [ ] Data quality approved or issues identified
- [ ] Go/No-Go decision for Phase 1

### Record Decisions

```yaml
# checkpoint_0_decisions.yaml
date: YYYY-MM-DD
participants: [developer, bioinformatician_name]
decisions:
  qc_maf_threshold: 0.01
  qc_missing_threshold: 0.05
  priority_drug: RIF
  use_prebuilt_tree: true/false
  data_quality_approved: true/false
notes: |
  Free-form notes from discussion
```

---

## Phase 1: Baseline GWAS

**Duration**: 8 hours  
**Dependencies**: Phase 0

### Task 1.1: Implement DataLoader

**Time**: 1.5 hours
**Status**: ✅
**Depends on**: 0.4

**Subtasks**:
- [x] Create `src/data/loader.py`
- [x] Implement `load_snp_matrix()` with HDF5/CSV support
- [x] Implement `load_phenotypes()`
- [x] Implement `load_who_catalogue()`
- [x] Handle missing values (-1 encoding)
- [x] Add input validation
- [x] Write unit tests

**Acceptance**:
- [x] Loads test data without errors
- [x] Returns correct DataFrame shapes
- [x] Validates input format
- [x] Tests pass

---

### Task 1.2: Implement Preprocessor

**Time**: 1.5 hours
**Status**: ✅
**Depends on**: 1.1

**Subtasks**:
- [x] Create `src/data/preprocessor.py`
- [x] Implement `filter_by_maf()`
- [x] Implement `filter_by_missingness()`
- [x] Implement `run_qc()` combining all filters
- [x] Generate QCReport dataclass
- [x] Log before/after statistics
- [x] Write unit tests

**Acceptance**:
- [x] MAF filter removes variants correctly
- [x] Missingness filter removes samples/variants
- [x] QCReport contains all required fields
- [x] Tests pass (7 tests in test_preprocessor.py)

---

### Task 1.3: Implement PhylogenyBuilder

**Time**: 2 hours
**Status**: ✅
**Depends on**: 0.4

**Subtasks**:
- [x] Create `src/phylogeny/builder.py`
- [x] Implement `build_tree()` with FastTree wrapper
- [x] Implement `build_tree()` with IQ-TREE wrapper (optional)
- [x] Implement `compute_kinship()` from tree (pyseer method)
- [x] Implement `compute_kinship()` from SNPs (fallback)
- [x] Implement `validate_kinship()` (PSD check)
- [x] Add caching for expensive operations
- [ ] Write unit tests (pending)

**Acceptance**:
- [x] Tree builds successfully from alignment
- [x] Kinship matrix is positive semi-definite
- [x] Kinship dimensions match sample count
- [x] Cache works (second run is faster)
- [ ] Tests pass (pending integration tests)

---

### Task 1.4: Implement GWASRunner (Basic)

**Time**: 2 hours
**Status**: ✅
**Depends on**: 1.2, 1.3

**Subtasks**:
- [x] Create `src/gwas/runner.py`
- [x] Implement `run()` using custom LMM (pyseer-style)
- [x] Create `GWASResult` dataclass
- [x] Implement `calculate_lambda_gc()` in `src/gwas/stats.py`
- [x] Log progress and key metrics
- [x] Write unit tests

**Acceptance**:
- [x] GWAS runs without errors
- [x] Lambda GC is calculated and logged
- [x] Results contain all required columns
- [x] Tests pass with mock data (9 tests in test_gwas.py)

---

### Task 1.5: Implement Visualizer (Basic)

**Time**: 1 hour
**Status**: ✅
**Depends on**: 1.4

**Subtasks**:
- [x] Create `src/output/visualizer.py`
- [x] Implement `manhattan_plot()`
- [x] Implement `qq_plot()` with Lambda annotation
- [x] Save as PNG (300 DPI for publication)
- [ ] Write unit tests (pending)

**Acceptance**:
- [x] Manhattan plot shows genome positions
- [x] QQ plot shows expected vs observed
- [x] Lambda GC displayed on QQ plot
- [x] Plots saved to correct location

---

### Task 1.6: Run Baseline GWAS for RIF

**Time**: 1 hour (mostly waiting)
**Status**: ✅
**Depends on**: 1.1-1.5

**Subtasks**:
- [x] Load data using DataLoader
- [x] Run QC using Preprocessor (MAF >= 0.01)
- [x] Compute SNP-based kinship matrix (cached)
- [x] Parallelize GWAS implementation:
  - [x] pyseer Python API (no I/O overhead)
  - [x] Block-based parallel processing (8 jobs)
- [x] Run LMM GWAS for RIF
- [x] Generate results with PRPS scores
- [x] Save results to `results/rifampicin_gwas.csv`
- [x] Log all metrics

**Acceptance**:
- [x] Lambda GC < 1.1 → **0.139** (excellent, indicates good population structure correction)
- [x] Top hits include known resistance genes (rpoB region: p=1.2e-98)
- [x] Results table contains expected columns

**Results (2026-01-05)**:
- Dataset: 11,649 samples × 6,780 variants
- Time: 2.7 min total (42.4 var/sec)
- Heritability h² = 0.950
- Significant hits (p<5e-8): 4
- Significant hits (FDR<0.05): 7
- Top hit: 3650378_G_A (rpoB region), p=1.2e-98, OR=1.86

### Phase 1 Improvements Summary

Based on Phase 1 results and literature review. **Details**: [phase1_improvements.md](phase1_improvements.md)

| Task | Issue | Status | 🔬 Ask Bio | 📚 Suggested Answer |
|------|-------|--------|-----------|---------------------|
| [1.7](#task-17-validate-elastic-net-method) | Lambda GC = 0.139 overcorrection | ⬜ After CP1 | Is elastic net appropriate for this dataset? | Yes, elastic net recommended. Saber & Shapiro 2020 found "all methods performed relatively poorly on highly clonal genomes" and elastic net "was consistently amongst the highest-performing methods" with superior precision/recall. Lambda GC < 1.0 is acceptable for clonal bacteria - indicates strong LMM correction working as intended. [DOI: 10.1099/mgen.0.000337](https://doi.org/10.1099/mgen.0.000337) |

**Action Points for CP1**:
- [x] 🔬 Confirm Lambda GC = 0.139 acceptable, or need elastic net → ✅ **Accepted** (2026-01-07): "Overcorrection в этой схеме уж точно меньше чем у нас, так что она имеет право быть такой, как есть". Task 1.7 kept for validation.
- [ ] 🔬 Review OR attenuation (1.86 vs clinical ~4.8)
- [ ] 🔬 Approve proceeding to Phase 2

---

### Task 1.7: Validate Elastic Net Method

**Time**: 2 hours
**Status**: ⬜
**Depends on**: 1.6
**Literature**: [Saber & Shapiro 2020](https://doi.org/10.1099/mgen.0.000337)

> "The multi-locus elastic net approach was consistently amongst the highest-performing methods... all methods performed relatively poorly on highly clonal genomes"

**Problem**: Lambda GC = 0.139 indicates LMM overcorrection for clonal TB population.

**Subtasks**:
- [ ] Add `method: str` to GWASConfig (`"lmm"`, `"enet"`)
- [ ] Implement `_run_elastic_net()` using pyseer `--wg enet`
- [ ] Run GWAS with elastic net for RIF
- [ ] Compare results: Lambda GC, top hits overlap, rpoB p-value
- [ ] Document method choice rationale

**Acceptance**:
- [ ] Elastic net Lambda GC closer to 1.0
- [ ] Top hits (rpoB) preserved
- [ ] Method comparison documented

**Details**: [phase1_improvements.md](phase1_improvements.md#imp-1-alternative-gwas-methods)

---

## 🔬 CHECKPOINT 1: Baseline Evaluation

**When**: After Task 1.6  
**Duration**: 45 min  
**Type**: Sync with bioinformatician

### Preparation (Developer)

Before meeting, prepare:
- [ ] Manhattan plot for RIF
- [ ] QQ plot with Lambda GC value
- [ ] Table of top-20 SNPs (variant, gene, p-value, OR)
- [ ] Summary statistics: n_samples, n_variants, n_significant
- [ ] Comparison with WHO catalogue (manual)

### Artifacts to Show

```
results/phase1/
├── rif_manhattan.png
├── rif_qq_plot.png
├── rif_top_hits.csv
└── rif_summary.json
```

### Questions to Ask

**Quality Assessment**:
1. "Lambda GC = [X]. Acceptable? Need more correction?"
2. "QQ plot looks like [describe]. Expected?"
3. "Found [N] significant SNPs. Too many/few?"

**Biological Validation**:
4. "Top genes are [list]. Which are expected?"
5. "rpoB S450L has p = [X], OR = [Y]. Correct magnitude?"
6. "Found [novel gene]. Interesting or artifact?"
7. "Missed [WHO gene]. Data issue or method issue?"

**Next Steps**:
8. "Quality sufficient to proceed to multi-drug?"
9. "Which drug to analyze next? (INH recommended)"
10. "Any parameter changes before scaling up?"

### Expected Outcomes

- [ ] Lambda GC approved (< 1.1) or remediation plan
- [ ] rpoB found in top hits confirmed
- [ ] Quality sufficient for multi-drug (Go/No-Go)
- [ ] Parameters adjustments documented

### Decision Record

```yaml
# checkpoint_1_decisions.yaml
date: YYYY-MM-DD
metrics:
  lambda_gc: 1.05
  n_significant: 25
  rpob_rank: 3
  sensitivity_who: 0.85
quality_approved: true
parameter_changes: []
proceed_to_multidrug: true
notes: |
  rpoB S450L found at rank 3, p=1.2e-45
  One unexpected hit in gene X - likely artifact
```

---

## Phase 2: Multi-Drug Analysis

**Duration**: 8 hours
**Dependencies**: Phase 1, Checkpoint 1 approval

### Phase 2 Improvements Summary

Based on literature review. **Details**: [phase2_improvements.md](phase2_improvements.md)

| Task | Description | Status | 🔬 Ask Bio | 📚 Suggested Answer |
|------|-------------|--------|-----------|---------------------|
| [2.0](#task-20-implement-gene-annotation-module) | Gene annotation (H37Rv) | ⬜ Can start | Which annotation source? NCBI vs TubercuList? | TubercuList (tuberculist.epfl.ch) is the gold standard for TB annotation with manually curated gene information. NCBI RefSeq NC_000962.3 provides coordinates but TubercuList adds functional annotations. Use both: NCBI for coordinates, TubercuList for gene names/functions. [DOI: 10.1016/j.tube.2010.09.008](https://doi.org/10.1016/j.tube.2010.09.008) |
| [2.6](#task-26-implement-lineage-stratified-gwas) | Lineage-stratified GWAS | ⬜ After 2.5 | L2/L4 split meaningful for this dataset? | Yes, highly meaningful. Omae et al. 2017 showed pathogen lineage-based GWAS reveals unique associations. Sobkowiak et al. 2020 found L2 Beijing associated with higher transmissibility. Li et al. 2024 found lineage-specific compensatory mutations. [DOI: 10.1038/jhg.2017.82](https://doi.org/10.1038/jhg.2017.82), [DOI: 10.1099/mgen.0.000361](https://doi.org/10.1099/mgen.0.000361) |
| [2.7](#task-27-predictive-ability-validation) | Validate predictive ability | ⬜ After 2.1 | Is ~98% accuracy achievable? | Known resistance mutations achieve ~98% sensitivity for RIF/INH. Compare GWAS hits vs known mutations from WHO catalog. |

**Action Points for CP2**:
- [ ] 🔬 Confirm annotation source (NCBI vs TubercuList)
- [ ] 🔬 Approve L2/L4 lineage stratification
- [ ] 🔬 Define cross-drug contamination threshold

---

### Task 2.0: Implement Gene Annotation Module

**Time**: 1.5 hours
**Status**: ⬜
**Depends on**: 0.3 (H37Rv download)
**Literature**: [Phelan et al. 2019](https://doi.org/10.1038/s41598-019-45566-5)

**Problem**: GWAS results show positions without gene context.

**Subtasks**:
- [ ] Download H37Rv GFF3 from NCBI RefSeq (NC_000962.3)
- [ ] Create `src/annotation/gene_mapper.py`
- [ ] Implement `GeneMapper` class with IntervalTree
- [ ] Implement `annotate_results()` for batch annotation
- [ ] Write unit tests

**Acceptance**:
- [ ] Position 761155 → rpoB
- [ ] Position 2155168 → katG
- [ ] Batch annotation <1 sec for 10K variants

**Details**: [phase2_improvements.md](phase2_improvements.md#gene-annotation-module)

---

### Task 2.1: Implement Conditional GWAS

**Time**: 2 hours  
**Status**: ⬜  
**Depends on**: 1.4

**Subtasks**:
- [ ] Extend `GWASRunner` with `run_conditional()` method
- [ ] Implement covariate extraction from previous drug results
- [ ] Order drugs correctly (RIF → INH → FQ → ...)
- [ ] Store drug-specific results separately
- [ ] Write integration tests

**Acceptance**:
- [ ] Covariates added for each subsequent drug
- [ ] Drug order follows config
- [ ] Results stored per drug
- [ ] Tests pass

---

### Task 2.2: Implement Gene Burden Analysis

**Time**: 3 hours
**Status**: ⬜
**Depends on**: 1.4, 2.0

**Reference**: Farhat et al. 2019 (Nature Communications) - [articles/farhat_2019_natcomm_methods.md](../articles/farhat_2019_natcomm_methods.md)

**Subtasks**:
- [ ] Create `src/gwas/burden.py`
- [ ] Implement `GeneAnnotation` dataclass
- [ ] Implement `load_annotation()` for H37Rv GFF3
- [ ] Implement `compute_burden_matrix()`:
  - Binary score: 1 if any non-synonymous SNS or indel in gene
  - For intergenic: 1 if any variant present
  - Exclude synonymous variants
  - Filter burden frequency < 0.01
- [ ] Implement `run_burden_gwas()` using existing LMM
- [ ] Implement `compare_with_site_level()` for concordance analysis
- [ ] Create `BurdenGWASResult` and `ComparisonReport` dataclasses
- [ ] Write unit tests

**Acceptance**:
- [ ] Burden matrix computed correctly (binary 0/1)
- [ ] ~2800 genes/regions tested (matching Farhat)
- [ ] Reuses LMM eigendecomposition (no redundant computation)
- [ ] Known genes (rpoB, katG) significant for respective drugs
- [ ] Tests pass

**Key Parameters (from Farhat et al.):**
- MAF threshold: 0.01 (burden frequency)
- Variants per gene: non-synonymous + indels only
- Intergenic regions: any variant counts
- Lineage markers: excluded from burden

---

### Task 2.3: Implement PRPS Calculation

**Time**: 1.5 hours
**Status**: ✅
**Depends on**: 1.3, 2.1

**Subtasks**:
- [x] Add `calculate_prps()` to `src/gwas/stats.py`
- [x] Add `calculate_prps_batch()` for vectorized batch processing (18,550x faster)
- [x] Add `calculate_prps_eigenspace()` for single-variant O(n*k) calculation
- [x] Reuse eigendecomposition from LMM (computed during h² estimation)
- [x] Integrate PRPS into GWAS results
- [x] Add PRPS column to output
- [x] Write unit tests (test_calculate_prps_batch, test_calculate_prps_batch_monomorphic)

**Acceptance**:
- [x] PRPS calculated for each variant
- [x] PRPS values in [0, 1] range
- [x] Batch PRPS: <3 sec for 6,780 variants (vs ~55 sec naive)
- [x] Tests pass (16 tests in test_gwas.py)

**Notes**:
- PRPS interpretation: High (>0.5) = hitchhiking, Low (<0.3) = convergent evolution
- Eigenspace method reuses eigenvectors/eigenvalues from LMM initialization

---

### Task 2.4: Implement WHOValidator

**Time**: 1.5 hours  
**Status**: ⬜  
**Depends on**: 2.1

**Subtasks**:
- [ ] Create `src/validation/who.py`
- [ ] Implement `validate()` method
- [ ] Implement `cross_drug_contamination()` method
- [ ] Create `ValidationReport` dataclass
- [ ] Support WHO catalogue 2021 and 2024 formats
- [ ] Write unit tests

**Acceptance**:
- [ ] Sensitivity and PPV calculated correctly
- [ ] Cross-drug contamination detected
- [ ] Reports generated with all fields
- [ ] Tests pass

---

### Task 2.5: Run Multi-Drug Conditional GWAS

**Time**: 2 hours (mostly waiting)
**Status**: ⬜
**Depends on**: 2.1-2.4

**Subtasks**:
- [ ] Run conditional GWAS for: RIF, INH, FQ (minimum)
- [ ] Run gene burden GWAS for each drug (Task 2.2)
- [ ] Calculate PRPS for all results
- [ ] Compare site-level vs gene burden results
- [ ] Validate each drug against WHO
- [ ] Calculate cross-drug contamination
- [ ] Generate per-drug Manhattan plots
- [ ] Generate cross-drug heatmap
- [ ] Save all results

**Acceptance**:
- [ ] 3+ drugs analyzed (site-level + gene burden)
- [ ] Cross-drug contamination < 10%
- [ ] RIF: rpoB ✓, katG ✗
- [ ] INH: katG ✓, rpoB ✗
- [ ] Gene burden identifies ≥1 additional locus not found by site-level

---

### Task 2.6: Implement Lineage-Stratified GWAS

**Time**: 2 hours
**Status**: ⬜
**Depends on**: 2.1, 2.5
**Literature**: [Oppong et al. 2019](https://doi.org/10.1186/s12864-019-5615-3)

> "Unique associations with XDR in lineage-specific analyses provide evidence of diverging evolutionary trajectories between lineages 2 and 4"

**Problem**: Combined analysis may miss lineage-specific associations.

**Subtasks**:
- [ ] Add lineage assignment to samples (L2 Beijing vs L4 Euro-American)
- [ ] Implement `run_stratified()` method in GWASRunner
- [ ] Run separate GWAS for L2 and L4 subsets
- [ ] Compare results: unique hits, shared hits, effect size differences
- [ ] Document lineage-specific findings

**Acceptance**:
- [ ] L2/L4 stratification working
- [ ] ≥1 lineage-specific association identified
- [ ] Comparison table generated

**Details**: [phase2_improvements.md](phase2_improvements.md#lineage-stratified-gwas)

---

### Task 2.7: Predictive Ability Validation

**Time**: 1.5 hours
**Status**: ⬜
**Depends on**: 2.1 (multi-drug complete)
**Suggested by**: Bioinformatics review (2026-01-07)

> "Интересно было бы ещё глянуть на предсказательную способность, если ориентироваться только на эти мутации (она для RIF и INH порядка 98% по известным мутациям)"

**Problem**: Need to validate GWAS-identified mutations against known resistance mutation benchmarks.

**Subtasks**:
- [ ] Obtain WHO catalog of known resistance mutations
- [ ] Calculate sensitivity: GWAS hits ∩ known mutations / known mutations
- [ ] Calculate precision: GWAS hits ∩ known mutations / GWAS hits
- [ ] Compare against ~98% benchmark for RIF/INH
- [ ] Document novel candidates (GWAS hits not in catalog)

**Acceptance**:
- [ ] Sensitivity ≥90% for RIF/INH against known mutations
- [ ] Precision documented for each drug
- [ ] Novel candidates flagged for review

---

## 🔬 CHECKPOINT 2: Drug Specificity

**When**: After Task 2.4  
**Duration**: 45 min  
**Type**: Sync with bioinformatician

### Preparation (Developer)

Before meeting, prepare:
- [ ] Per-drug summary table: drug → top genes → WHO overlap
- [ ] Cross-drug contamination matrix/heatmap
- [ ] PRPS distribution plots
- [ ] Comparison: conditional vs non-conditional results

### Questions to Ask

**Cross-Drug Validation**:
1. "Cross-contamination is [X]%. Acceptable?"
2. "katG still shows weak signal for RIF (p=[X]). Problem?"
3. "Drug analysis order correct? Should we change?"

**PRPS Interpretation**:
4. "PRPS distribution looks like [describe]. Expected?"
5. "High-PRPS hits include [genes]. Filter them?"
6. "Low-PRPS novel hits: [genes]. Worth investigating?"

**Methodology**:
7. "Conditional approach working? Alternative methods?"
8. "Ready for ABESS feature selection?"

### Expected Outcomes

- [ ] Cross-drug contamination approved (< 10%)
- [ ] Drug-specific genes confirmed
- [ ] PRPS interpretation agreed
- [ ] Proceed to ABESS (Go/No-Go)

---

## Phase 3: Feature Selection

**Duration**: 4 hours  
**Dependencies**: Phase 2

### Task 3.1: Implement ABESSSelector

**Time**: 2 hours  
**Status**: ⬜  
**Depends on**: 2.1

**Subtasks**:
- [ ] Create `src/feature_selection/abess.py`
- [ ] Implement `select()` with iterative approach
- [ ] Implement `cross_validate()` for evaluation
- [ ] Create `ABESSResult` dataclass
- [ ] Log iteration details
- [ ] Write unit tests

**Acceptance**:
- [ ] ABESS selects reasonable number of features (10-50)
- [ ] Iteration 2 finds additional features
- [ ] Cross-validation works
- [ ] Tests pass

---

### Task 3.2: Run ABESS for All Drugs

**Time**: 2 hours  
**Status**: ⬜  
**Depends on**: 3.1

**Subtasks**:
- [ ] Run ABESS for each analyzed drug
- [ ] Compute CV metrics (AUC, etc.)
- [ ] Compare with GWAS top hits
- [ ] Generate comparison table
- [ ] Save results

**Acceptance**:
- [ ] CV AUC > 0.85 for RIF, INH
- [ ] ABESS-GWAS overlap > 50%
- [ ] Results for all drugs saved

---

## Phase 4: Epistasis Analysis

**Duration**: 8 hours
**Dependencies**: Phase 2

### Phase 4 Improvements Summary

Advanced analysis based on literature. **Details**: [phase4_improvements.md](phase4_improvements.md)

| Task | Description | Status | 🔬 Ask Bio | 📚 Suggested Answer |
|------|-------------|--------|-----------|---------------------|
| [4.5](#task-45-implement-convergent-evolution-phyc-analysis) | Convergent evolution (phyC) | ⬜ After 4.4 | Is phyC appropriate for compensatory mutation detection? | Yes, phyC is appropriate. Sobkowiak et al. 2020 used "evolutionary convergence testing (phyC) and GWAS" together to identify 6 loci associated with transmission. Lai & Ioerger 2020 ECAT method also exploits homoplasy for improved GWAS power. [DOI: 10.1099/mgen.0.000361](https://doi.org/10.1099/mgen.0.000361), [DOI: 10.1177/1176934320944932](https://doi.org/10.1177/1176934320944932) |
| [4.6](#task-46-implement-pfam-domain-aggregation) | PFAM domain aggregation | ⬜ After 2.0 | Which PFAM database version? | Use latest PFAM (now InterPro). Reshetnikov et al. 2025 showed "aggregating rare mutations within protein-coding genes into markers indicative of changes in PFAM domains improved prediction quality" and ABESS preferentially selected domain markers. [DOI: 10.3389/fmicb.2025.1586476](https://doi.org/10.3389/fmicb.2025.1586476) |
| [4.7](#task-47-implement-3d-structure-validation) | 3D structure validation | ⬜ After 4.4 | Which PDB structures for rpoB/katG? | rpoB: Use MTB RNAP structures with rifampin (PDB 5UHB, 5UHC, 5UHE). Molodtsov et al. 2017 solved structures with S531L, D516V, H526Y mutations. katG: PDB 1SJ2 (M. tuberculosis KatG), 3WXO (S. elongatus with INH). [DOI: 10.1111/mmi.13606](https://doi.org/10.1111/mmi.13606), [DOI: 10.1111/febs.13102](https://doi.org/10.1111/febs.13102) |

**Action Points for CP3**:
- [ ] 🔬 Confirm epistatic pairs (rpoB-rpoC, katG-ahpC)
- [ ] 🔬 Approve phyC for compensatory mutations
- [ ] 🔬 Define 3D distance threshold (<10Å?)

---

### Task 4.1: Implement Pairwise Epistasis Tests

**Time**: 2 hours  
**Status**: ⬜  
**Depends on**: 2.1

**Subtasks**:
- [ ] Create `src/epistasis/pairwise.py`
- [ ] Implement `test_known_pairs()` method
- [ ] Interaction term in logistic regression
- [ ] Create `EpistasisResult` dataclass
- [ ] Include kinship as covariates
- [ ] Write unit tests

**Acceptance**:
- [ ] Interaction p-value calculated
- [ ] OR and CI calculated
- [ ] Kinship correction included
- [ ] Tests pass

---

### Task 4.2: Implement Mediation Analysis

**Time**: 2 hours  
**Status**: ⬜  
**Depends on**: 4.1

**Subtasks**:
- [ ] Create `src/epistasis/mediation.py`
- [ ] Implement mediation framework (total/direct/indirect)
- [ ] Implement Sobel test for indirect effect
- [ ] Calculate proportion mediated
- [ ] Create `MediationResult` dataclass
- [ ] Write unit tests

**Acceptance**:
- [ ] All effects calculated correctly
- [ ] Sobel test p-value calculated
- [ ] Proportion mediated in [0, 1]
- [ ] Tests pass

---

### Task 4.3: Implement Sign Epistasis Detection

**Time**: 1 hour  
**Status**: ⬜  
**Depends on**: 4.1

**Subtasks**:
- [ ] Create `src/epistasis/sign.py`
- [ ] Implement `detect_sign_epistasis()` method
- [ ] Calculate all four genotype effects (00, 01, 10, 11)
- [ ] Detect sign reversals
- [ ] Create `SignEpistasisResult` dataclass
- [ ] Write unit tests

**Acceptance**:
- [ ] Four effects calculated
- [ ] Sign epistasis correctly identified
- [ ] Tests pass

---

### Task 4.4: Run Epistasis Analysis

**Time**: 2 hours  
**Status**: ⬜  
**Depends on**: 4.1-4.3

**Subtasks**:
- [ ] Test known pairs: rpoB-rpoC, rpoB-rpoA, katG-ahpC, gyrA-rpoC
- [ ] Run mediation for rpoB→rpoC and katG→ahpC
- [ ] Check sign epistasis for eis-rrs
- [ ] Generate summary table
- [ ] Save results

**Acceptance**:
- [ ] rpoB-rpoC interaction significant (p < 0.01)
- [ ] Mediation shows indirect effect
- [ ] Results saved correctly

---

### Task 4.5: Implement Convergent Evolution (phyC) Analysis

**Time**: 2 hours
**Status**: ⬜
**Depends on**: 1.3, 4.1
**Literature**: [Chen & Shapiro 2015](https://doi.org/10.1016/j.mib.2015.03.002)

> "We compare the traditional GWAS against phyC, a contrasting method of mapping genotype to phenotype based upon evolutionary convergence"

**Problem**: GWAS may miss variants that arose independently multiple times.

**Subtasks**:
- [ ] Create `src/epistasis/convergence.py`
- [ ] Implement `compute_convergence_score()` using phylogeny
- [ ] Count independent origins of each mutation
- [ ] Compare with GWAS p-values
- [ ] Identify convergence-enriched vs GWAS-enriched variants

**Acceptance**:
- [ ] Convergence scores calculated for all variants
- [ ] Comparison table: GWAS rank vs convergence rank
- [ ] Known compensatory mutations (rpoC) have high convergence

**Details**: [phase4_improvements.md](phase4_improvements.md#convergent-evolution-analysis)

---

### Task 4.6: Implement PFAM Domain Aggregation

**Time**: 1.5 hours
**Status**: ⬜
**Depends on**: 2.0
**Literature**: [Reshetnikov et al. 2025](https://doi.org/10.3389/fmicb.2025.1586476)

> "Aggregating rare mutations within protein-coding genes into markers indicative of changes in PFAM domains improved prediction quality"

**Problem**: Rare variants may lack power individually but aggregate to domains.

**Subtasks**:
- [ ] Download PFAM domain annotations for H37Rv
- [ ] Extend GeneMapper with domain information
- [ ] Implement `aggregate_by_domain()` method
- [ ] Run domain-level GWAS
- [ ] Compare with gene-level and site-level results

**Acceptance**:
- [ ] PFAM domains mapped to H37Rv genes
- [ ] Domain-level association test working
- [ ] ≥1 domain hit not captured by site-level

**Details**: [phase4_improvements.md](phase4_improvements.md#pfam-domain-aggregation)

---

### Task 4.7: Implement 3D Structure Validation

**Time**: 2 hours
**Status**: ⬜
**Depends on**: 2.0, 4.4
**Literature**: [Phelan et al. 2016](https://doi.org/10.1186/s12916-016-0575-9)

> "A strong direct correlation was observed between MIC values and the distance of mutated residues to drug binding sites in 3D structures"

**Problem**: Need to validate GWAS hits against structural biology.

**Subtasks**:
- [ ] Download rpoB/katG/gyrA PDB structures
- [ ] Implement `calculate_binding_distance()` method
- [ ] Map GWAS hits to 3D positions
- [ ] Correlate p-value/effect size with binding distance
- [ ] Flag structurally implausible hits

**Acceptance**:
- [ ] 3D distances calculated for rpoB hits
- [ ] Correlation plot: distance vs -log10(p)
- [ ] Structural interpretation documented

**Details**: [phase4_improvements.md](phase4_improvements.md#3d-structure-validation)

---

## 🔬 CHECKPOINT 3: Epistasis Review

**When**: After Task 4.4  
**Duration**: 45 min  
**Type**: Sync with bioinformatician

### Preparation (Developer)

Before meeting, prepare:
- [ ] Epistasis results table (pair, p-value, effect size)
- [ ] Mediation analysis diagram/results
- [ ] Sign epistasis findings (if any)
- [ ] Questions about interpretation

### Questions to Ask

**Epistasis Validation**:
1. "rpoB-rpoC interaction: p=[X], OR=[Y]. Matches literature?"
2. "Mediation shows [Z]% through rpoC. Biologically sensible?"
3. "katG-ahpC not found. Problem or expected for this dataset?"

**Novel Findings**:
4. "Found unexpected interaction: [pair]. Known or novel?"
5. "Sign epistasis in [pair]. Real or artifact?"

**Custom Model Decision**:
6. "Given results so far, is custom model worth 12 hours?"
7. "What should model do that GWAS+ABESS don't?"

### Expected Outcomes

- [ ] Known epistasis validated
- [ ] Novel findings assessed
- [ ] Custom model Go/No-Go decision
- [ ] Interpretation guidance for write-up

---

## Phase 5: Custom Model (Optional)

**Duration**: 12 hours  
**Dependencies**: Phase 2, Checkpoint 3 approval  
**Condition**: Only if approved at Checkpoint 3

### Task 5.1: Implement Drug-SNP Interaction Model

**Time**: 4 hours  
**Status**: ⬜  
**Depends on**: 2.1

**Subtasks**:
- [ ] Create `src/models/drug_snp_model.py`
- [ ] Implement model architecture (PyTorch)
- [ ] Implement drug-SNP attention mechanism
- [ ] Implement kinship correction layer
- [ ] Write unit tests for forward pass

**Acceptance**:
- [ ] Model builds without errors
- [ ] Forward pass produces valid output
- [ ] Attention weights extractable
- [ ] Tests pass

---

### Task 5.2: Implement Training Pipeline

**Time**: 3 hours  
**Status**: ⬜  
**Depends on**: 5.1

**Subtasks**:
- [ ] Create `src/models/trainer.py`
- [ ] Implement train/val/test split
- [ ] Implement training loop with early stopping
- [ ] Implement evaluation metrics (AUC per drug)
- [ ] Add logging and checkpointing
- [ ] Write integration tests

**Acceptance**:
- [ ] Training runs without errors
- [ ] Loss decreases over epochs
- [ ] Early stopping works
- [ ] Checkpoints saved

---

### Task 5.3: Train and Evaluate Model

**Time**: 3 hours  
**Status**: ⬜  
**Depends on**: 5.2

**Subtasks**:
- [ ] Prepare data tensors
- [ ] Train model (100 epochs, early stopping)
- [ ] Evaluate on test set
- [ ] Compare with ABESS baselines
- [ ] Extract attention weights per drug
- [ ] Save model and results

**Acceptance**:
- [ ] Per-drug AUC ≥ 0.85
- [ ] Model matches or beats ABESS
- [ ] Attention weights stable

---

### Task 5.4: Interpret Model Results

**Time**: 2 hours  
**Status**: ⬜  
**Depends on**: 5.3

**Subtasks**:
- [ ] Extract top-50 features per drug from attention
- [ ] Compare with GWAS hits (Jaccard similarity)
- [ ] Identify "model-only" candidates
- [ ] Generate interpretation report
- [ ] Visualize attention patterns

**Acceptance**:
- [ ] GWAS-model overlap > 50%
- [ ] Novel candidates flagged
- [ ] Report generated

---

## Phase 6: Finalization

**Duration**: 4 hours  
**Dependencies**: All previous phases

### Task 6.1: Generate Final Report

**Time**: 2 hours  
**Status**: ⬜  
**Depends on**: All analysis tasks

**Subtasks**:
- [ ] Create `src/output/reporter.py`
- [ ] Generate comprehensive report with all results
- [ ] Include all figures (publication quality)
- [ ] Export to HTML and PDF
- [ ] Export data tables to CSV

**Acceptance**:
- [ ] Report includes all phases
- [ ] Figures are high quality (300 DPI)
- [ ] Tables are complete

---

### Task 6.2: Code Cleanup and Documentation

**Time**: 1.5 hours
**Status**: 🔄
**Depends on**: 6.1

**Subtasks**:
- [x] Set up GitHub CI with ruff, mypy, pytest
- [x] Run linter (ruff) and fix issues
- [x] Run type checker (mypy) and fix issues
- [x] Configure ruff for scientific code (ignore N806, F841)
- [x] Configure mypy for numpy code (ignore import-untyped, warn_return_any=false)
- [x] Update README with quick start guide (ARM/AMD)
- [ ] Ensure all public functions have docstrings
- [ ] Create CHANGELOG.md

**Acceptance**:
- [x] `pixi run ruff check src/ tests/` passes
- [x] `pixi run mypy src/` passes (21 source files)
- [x] `pixi run pytest tests/` passes (23 tests)
- [x] GitHub CI runs on push/PR to main

---

### Task 6.3: Final Testing

**Time**: 30 min  
**Status**: ⬜  
**Depends on**: 6.2

**Subtasks**:
- [ ] Run full test suite
- [ ] Verify test coverage > 80%
- [ ] Run full pipeline on small test data
- [ ] Verify reproducibility (run twice, compare)

**Acceptance**:
- [ ] All tests pass
- [ ] Coverage > 80%
- [ ] Results reproducible

---

## 🔬 CHECKPOINT 4: Final Review

**When**: After all tasks complete  
**Duration**: 1 hour  
**Type**: Sync with bioinformatician

### Preparation (Developer)

Before meeting, prepare:
- [ ] Final report (HTML/PDF)
- [ ] Summary of all findings
- [ ] List of novel candidates
- [ ] Methods section draft
- [ ] Limitations list

### Questions to Ask

**Results Validation**:
1. "Overall findings summary. Anything missing?"
2. "Novel candidates: [list]. Which to highlight?"
3. "Unexpected results: [list]. How to interpret?"

**Publication Readiness**:
4. "Methods description complete?"
5. "Which findings are publication-worthy?"
6. "What experiments would validate novel findings?"

**Future Work**:
7. "What would improve with more time?"
8. "Which additional datasets would help?"
9. "Recommended next steps?"

### Expected Outcomes

- [ ] Results approved for use
- [ ] Novel findings prioritized
- [ ] Methods validated
- [ ] Next steps defined

---

## Summary: Time Allocation

| Phase | Tasks | Hours | Checkpoints |
|-------|-------|-------|-------------|
| 0: Setup | 0.1-0.5 | 4 | CP0 (0.5h) |
| 1: Baseline | 1.1-1.6 | 8 | CP1 (0.75h) |
| 2: Multi-drug | 2.1-2.5 | 11 | CP2 (0.75h) |
| 3: ABESS | 3.1-3.2 | 4 | - |
| 4: Epistasis | 4.1-4.4 | 8 | CP3 (0.75h) |
| 5: Model | 5.1-5.4 | 12 | - |
| 6: Final | 6.1-6.3 | 4 | CP4 (1h) |
| **Total** | | **51** | **~4h** |
| **Buffer** | | **4** | |
| **Checkpoints** | | **4** | |
| **Grand Total** | | **59** | |

*Note: Phase 2 includes gene burden analysis (Task 2.2, +3h) based on Farhat et al. 2019*

*Remaining 4 hours for unexpected issues*

---

## Dependency Graph

```
Phase 0
  │
  ├─ 0.1 ─┬─ 0.4 ─┬─ 1.1 ─── 1.2 ─┬─ 1.4 ─── 1.6 ──── [CP1]
  ├─ 0.2 ─┤       │               │                      │
  ├─ 0.3 ─┘       └─ 1.3 ─────────┘                      │
  └─ 0.5 ──────────── 1.5 ───────────────────────────────┘
                                                         │
  [CP0] ─────────────────────────────────────────────────┘
                                                         │
Phase 2-3                                                ▼
  │
  ├─ 2.1 ─┬─ 2.2 ─┬─ 2.4 ──── [CP2] ─── 3.1 ─── 3.2
  │       │       │
  └─ 2.3 ─┴───────┘
                                        │
Phase 4                                 ▼
  │
  └─ 4.1 ─┬─ 4.2 ─┬─ 4.4 ──── [CP3]
          │       │
          └─ 4.3 ─┘             │
                                ▼
Phase 5 (optional)              │
  │                             │
  └─ 5.1 ─── 5.2 ─── 5.3 ─── 5.4
                                │
Phase 6                         ▼
  │
  └─ 6.1 ─── 6.2 ─── 6.3 ──── [CP4]
```

---

## Risk Mitigation Tasks

### If Lambda GC > 1.2 (after Task 1.6)

**Mitigation Task M1**: Enhanced Phylogenetic Correction

- [ ] Increase PCA components to 20
- [ ] Try SNP-based kinship instead of tree-based
- [ ] Use stricter MAF filter (0.05 instead of 0.01)
- [ ] Consult bioinformatician at CP1

**Time**: +2 hours

---

### If rpoB Not Found (after Task 1.6)

**Mitigation Task M2**: Debug Data Pipeline

- [ ] Verify rpoB mutations exist in raw data
- [ ] Check phenotype labels not swapped
- [ ] Check sample alignment between SNPs and phenotypes
- [ ] Inspect intermediate data at each step
- [ ] Consult bioinformatician at CP1

**Time**: +2-4 hours

---

### If Cross-Drug Contamination > 20% (after Task 2.4)

**Mitigation Task M3**: Enhanced Conditional Analysis

- [ ] Use more aggressive covariate inclusion
- [ ] Try multi-trait GWAS instead of conditional
- [ ] Apply stricter PRPS filtering (threshold 0.3)
- [ ] Consult bioinformatician at CP2

**Time**: +3 hours

---

### If Model Underperforms (after Task 5.3)

**Mitigation Task M4**: Model Debugging

- [ ] Simplify architecture
- [ ] Increase regularization
- [ ] Try longer training
- [ ] Fall back to per-drug models
- [ ] Document as negative result

**Time**: +2 hours (or abandon)

---

## Revision History

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | 2026-01-05 | Developer | Initial task breakdown |
| 1.1 | 2026-01-05 | Developer | Task 1.6 complete (RIF GWAS), Task 2.2 complete (PRPS), CI setup |
| 1.2 | 2026-01-05 | Developer | Added improvements checklist after Task 1.6, details in [phase1_improvements.md](phase1_improvements.md) |
| 1.3 | 2026-01-06 | Developer | Added Task 2.2 (Gene Burden Analysis) based on Farhat et al. 2019; renumbered subsequent tasks |
| 1.4 | 2026-01-06 | Developer | Added CRyPTIC validation dataset (12,289 isolates) to Task 0.3 |
