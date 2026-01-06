# TB GWAS: Requirements Specification

## Project Overview

**Project**: Genome-Wide Association Study для поиска генетических вариантов, ассоциированных с лекарственной устойчивостью Mycobacterium tuberculosis

**Domain**: Bioinformatics / Computational Genomics

**Scientific Context**: 
- WHO estimates 10 million TB cases annually, ~500,000 with drug resistance
- Understanding genetic basis of resistance is critical for diagnostics and treatment
- Standard GWAS methods fail for bacteria due to clonal population structure

---

## Stakeholders

| Role | Needs | Success Metric |
|------|-------|----------------|
| **Researcher** | Find resistance-associated variants | Validated discoveries matching WHO catalogue |
| **Bioinformatician** | Scientifically sound methodology | Lambda GC < 1.1, proper phylogenetic correction |
| **Developer** | Clear specs, testable components | Modular code, >80% test coverage |
| **Reviewer** | Reproducible results | Docker/Pixi environment, fixed seeds |

---

## Epic 1: Single-Drug GWAS with Phylogenetic Correction

### User Story 1.1: Basic Association Analysis

**As a** researcher  
**I want to** identify genetic variants associated with resistance to a single antibiotic  
**So that** I can understand the molecular mechanisms of drug resistance

**Acceptance Criteria:**

| ID | Criterion | Validation Method |
|----|-----------|-------------------|
| AC-1.1.1 | System processes SNP matrix with ≥500 samples and ≥10,000 variants | Load test dataset, verify dimensions |
| AC-1.1.2 | Results include p-value, odds ratio, confidence interval for each variant | Check output schema |
| AC-1.1.3 | Known resistance genes (rpoB for RIF, katG for INH) appear in top-20 hits | Compare with WHO catalogue |
| AC-1.1.4 | Processing completes within 2 hours on standard laptop (M1/M4 Mac, 16GB RAM) | Benchmark on reference dataset |

**Scientific References:**
- Farhat et al. 2013 (Nature Genetics) - phylogenetic convergence approach
- CRyPTIC Consortium 2022 (PLOS Biology) - large-scale TB GWAS

---

### User Story 1.2: Phylogenetic Correction

**As a** bioinformatician  
**I want to** account for population structure in the analysis  
**So that** I avoid false positive associations due to shared ancestry

**Acceptance Criteria:**

| ID | Criterion | Validation Method |
|----|-----------|-------------------|
| AC-1.2.1 | Genomic inflation factor (Lambda GC) < 1.1 after correction | Calculate from p-value distribution |
| AC-1.2.2 | QQ plot shows minimal deviation from expected distribution (except true signals) | Visual inspection + statistical test |
| AC-1.2.3 | Support both tree-based kinship and SNP-based kinship matrix | Test both methods |
| AC-1.2.4 | Kinship matrix is positive semi-definite (valid for LMM) | Eigenvalue check |

**Scientific References:**
- Earle et al. 2016 (Nature Microbiology) - bugwas, importance of LMM for bacteria
- Power et al. 2017 (Nature Communications) - benchmarking bacterial GWAS methods
- Lees et al. 2018 (Bioinformatics) - pyseer methodology

---

### User Story 1.3: Validation Against Known Mutations

**As a** researcher  
**I want to** validate results against WHO mutation catalogue  
**So that** I can assess the quality of my analysis

**Acceptance Criteria:**

| ID | Criterion | Validation Method |
|----|-----------|-------------------|
| AC-1.3.1 | Sensitivity ≥ 80% for Tier 1 WHO mutations | Compare significant hits with catalogue |
| AC-1.3.2 | Positive Predictive Value (PPV) ≥ 70% | Manual review of top-50 hits |
| AC-1.3.3 | Report clearly separates "confirmed", "novel", and "missed" variants | Structured output format |
| AC-1.3.4 | Support WHO catalogue versions 2021 and 2024 | Configurable catalogue path |

**Scientific References:**
- WHO Catalogue 2021/2024 - ground truth for known mutations
- Walker et al. 2015 (Lancet ID) - TBProfiler validation

---

## Epic 2: Multi-Drug Conditional Analysis

### User Story 2.1: Drug-Specific Association Discovery

**As a** researcher  
**I want to** find variants specific to each drug (not cross-drug correlations)  
**So that** I understand true causal relationships, not confounders

**Acceptance Criteria:**

| ID | Criterion | Validation Method |
|----|-----------|-------------------|
| AC-2.1.1 | rpoB mutations appear significant ONLY for RIF (not INH, EMB) | Cross-tabulate results |
| AC-2.1.2 | katG mutations appear significant ONLY for INH (not RIF, FQ) | Cross-tabulate results |
| AC-2.1.3 | Cross-drug contamination rate < 10% | Calculate overlap of top hits |
| AC-2.1.4 | Conditional analysis includes previous drugs' top hits as covariates | Inspect model specification |

**Scientific References:**
- Reshetnikoff et al. 2025 (Frontiers Microbiology) - cross-drug association problem
- CRyPTIC 2022 - multi-drug phenotyping

---

### User Story 2.2: Gene Burden Analysis

**As a** researcher
**I want to** aggregate rare variants per gene into burden scores
**So that** I can detect resistance genes where multiple different mutations contribute

**Acceptance Criteria:**

| ID | Criterion | Validation Method |
|----|-----------|-------------------|
| AC-2.2.1 | Binary burden score: 1 if any non-syn SNS or indel in gene, 0 otherwise | Verify scoring logic |
| AC-2.2.2 | Gene-level associations tested alongside site-level | Run both analyses |
| AC-2.2.3 | Non-coding regions (promoters) tested with any variant = 1 | Verify intergenic handling |
| AC-2.2.4 | Burden frequency ≥ 0.01 filter applied | Check MAF filtering |
| AC-2.2.5 | Known genes (rpoB, katG) rank in top-20 for respective drugs | Compare results |

**Scientific References:**
- Farhat et al. 2019 (Nature Communications) - Gene burden GWAS for MTB resistance
- GEMMA gene-based tests - Standard burden test methodology

**Rationale:**
- Site-level GWAS may miss genes with multiple rare causal variants
- Farhat et al. found 13 non-canonical loci using burden approach (including ubiA, whiB6)
- Aggregation increases power for rare variant effects (e.g., ccsA with AF 0.03)

---

### User Story 2.4: PRPS Filtering

**As a** bioinformatician
**I want to** filter variants by Phylogeny-Related Parallelism Score
**So that** I prioritize convergent evolution signals over hitchhiking mutations

**Acceptance Criteria:**

| ID | Criterion | Validation Method |
|----|-----------|-------------------|
| AC-2.4.1 | PRPS score calculated for each significant variant | Check output includes PRPS column |
| AC-2.4.2 | High-PRPS variants (>0.5) flagged as "potentially spurious" | Verify flagging logic |
| AC-2.4.3 | True positive rate higher in low-PRPS subset | Compare with WHO catalogue |
| AC-2.4.4 | PRPS calculation uses actual phylogenetic tree | Verify tree is input to function |

**Scientific References:**
- Farhat et al. 2013 - evolutionary convergence as signal
- Collins & Didelot 2018 - TreeWAS methodology

---

### User Story 2.5: Feature Selection with ABESS

**As a** researcher  
**I want to** use best-subset selection to find optimal SNP combinations  
**So that** I account for correlations between variants and find weak signals

**Acceptance Criteria:**

| ID | Criterion | Validation Method |
|----|-----------|-------------------|
| AC-2.3.1 | ABESS selects 10-50 features per drug (not thousands) | Check selected feature count |
| AC-2.3.2 | Iterative approach runs 2+ iterations to find weak associations | Verify iteration count in logs |
| AC-2.3.3 | Cross-validation AUC > 0.85 for well-studied drugs (RIF, INH) | 5-fold CV evaluation |
| AC-2.3.4 | ABESS results overlap >50% with GWAS top hits | Jaccard similarity |

**Scientific References:**
- Zhu et al. 2020 (PNAS) - ABESS algorithm
- Reshetnikoff et al. 2025 - ABESS vs HHS comparison for TB

---

## Epic 3: Epistasis Analysis

### User Story 3.1: Known Epistatic Pairs

**As a** researcher  
**I want to** detect interactions between known compensatory mutation pairs  
**So that** I validate the pipeline can find epistasis

**Acceptance Criteria:**

| ID | Criterion | Validation Method |
|----|-----------|-------------------|
| AC-3.1.1 | rpoB-rpoC interaction detected with p < 0.01 | Statistical test output |
| AC-3.1.2 | katG-ahpC interaction detected with p < 0.01 | Statistical test output |
| AC-3.1.3 | Interaction effect size (OR) reported with confidence interval | Check output schema |
| AC-3.1.4 | At least 4 known pairs tested from literature | Verify test coverage |

**Known pairs to test:**
- rpoB + rpoC (RIF fitness compensation)
- rpoB + rpoA (RIF fitness compensation)
- katG + ahpC (INH oxidative stress)
- gyrA + rpoC (FQ-RIF cross-resistance, Moura 2025)

**Scientific References:**
- Comas et al. 2012 (Nature Genetics) - compensatory mutations
- Moura et al. 2024 (Scientific Reports) - rpoC epistasis

---

### User Story 3.2: Mediation Analysis

**As a** bioinformatician  
**I want to** quantify how much compensatory mutations mediate primary mutation effects  
**So that** I understand the causal pathway

**Acceptance Criteria:**

| ID | Criterion | Validation Method |
|----|-----------|-------------------|
| AC-3.2.1 | Indirect effect (through mediator) calculated and tested | Check mediation output |
| AC-3.2.2 | Proportion mediated reported as percentage | Verify calculation |
| AC-3.2.3 | Sobel test p-value provided for indirect effect | Statistical significance |
| AC-3.2.4 | rpoC mediates >20% of rpoB effect on fitness/transmission | Compare with literature |

**Scientific References:**
- Casali et al. 2014 (Lancet ID) - compensatory evolution
- CompMut-TB 2024 - mediation framework

---

### User Story 3.3: Sign Epistasis Detection

**As a** researcher  
**I want to** detect cases where two beneficial mutations are deleterious together  
**So that** I understand evolutionary constraints on resistance

**Acceptance Criteria:**

| ID | Criterion | Validation Method |
|----|-----------|-------------------|
| AC-3.3.1 | Sign epistasis detected for eis-rrs pair (AMK/KAN) | Statistical test |
| AC-3.3.2 | All four genotype combinations (00, 01, 10, 11) have sufficient samples (n≥10) | Sample size check |
| AC-3.3.3 | Effect sizes reported for each genotype | Output includes all effects |
| AC-3.3.4 | Warning generated when sign epistasis conflicts with expectations | Logging verification |

**Scientific References:**
- WHO Catalogue 2023 - eis promoter interactions
- Borrell & Gagneux 2011 - epistasis in TB evolution

---

## Epic 4: Custom Model (Optional)

### User Story 4.1: Multi-Drug Joint Prediction

**As a** researcher  
**I want to** predict resistance to multiple drugs simultaneously  
**So that** I can leverage shared information across drugs

**Acceptance Criteria:**

| ID | Criterion | Validation Method |
|----|-----------|-------------------|
| AC-4.1.1 | Single model predicts all drugs (not N separate models) | Architecture inspection |
| AC-4.1.2 | Per-drug AUC ≥ 0.85 on held-out test set | Evaluation metrics |
| AC-4.1.3 | Training completes in < 4 hours on standard hardware | Benchmark |
| AC-4.1.4 | Model outperforms or matches ABESS per-drug models | Comparative evaluation |

---

### User Story 4.2: Interpretable Drug-SNP Importance

**As a** bioinformatician  
**I want to** extract which SNPs the model considers important for each drug  
**So that** I can validate biological plausibility

**Acceptance Criteria:**

| ID | Criterion | Validation Method |
|----|-----------|-------------------|
| AC-4.2.1 | Attention weights / feature importance extractable per drug | API method exists |
| AC-4.2.2 | Top-50 model features overlap >50% with GWAS hits | Jaccard calculation |
| AC-4.2.3 | "Model-only" candidates flagged for biological review | Output format |
| AC-4.2.4 | Importance scores are stable across training runs (CV std < 0.1) | Cross-validation |

---

## Non-Functional Requirements

### NFR-1: Reproducibility

| ID | Requirement | Validation |
|----|-------------|------------|
| NFR-1.1 | All random seeds fixed and logged | Grep logs for seed values |
| NFR-1.2 | Environment fully specified (Pixi/Docker) | Fresh install test |
| NFR-1.3 | Results identical across runs with same seed | Diff outputs |
| NFR-1.4 | Git commit hash recorded with each run | Check output metadata |

### NFR-2: Performance

| ID | Requirement | Validation | Achieved |
|----|-------------|------------|----------|
| NFR-2.1 | Process 500 samples × 50K variants in < 30 min | Benchmark | ✅ 11,649 × 6,780 in 2.7 min |
| NFR-2.2 | RAM usage < 16 GB for standard analysis | Memory profiling | ✅ ~4 GB peak |
| NFR-2.3 | Support parallel execution (configurable n_jobs) | Test with n_jobs=1,4,8 | ✅ 8 jobs default |
| NFR-2.4 | Intermediate results cached to avoid recomputation | Verify cache hits | ✅ Eigendecomp cached |

**Actual Performance Benchmarks (11,649 samples):**

| Operation | Target | Achieved | Notes |
|-----------|--------|----------|-------|
| Eigendecomposition | - | ~2 min | One-time, cached |
| GWAS (6,780 variants) | <30 min | 2.7 min | 143 var/sec (pyseer API) |
| PRPS calculation | - | <1 sec | 18,550x faster than naive |
| h² optimization | - | ~2 min | Per phenotype |
| Full pipeline | <2 hours | ~5 min | After eigendecomp |

### NFR-3: Usability

| ID | Requirement | Validation |
|----|-------------|------------|
| NFR-3.1 | Single command to run full pipeline | Test CLI |
| NFR-3.2 | Configuration via YAML files (not code changes) | Config-driven run |
| NFR-3.3 | Meaningful error messages for common failures | Test error scenarios |
| NFR-3.4 | Progress logging with ETA | Visual verification |

### NFR-4: Maintainability

| ID | Requirement | Validation |
|----|-------------|------------|
| NFR-4.1 | Test coverage > 80% for core modules | pytest-cov report |
| NFR-4.2 | Type hints on all public functions | mypy check |
| NFR-4.3 | Docstrings with examples for key functions | Documentation review |
| NFR-4.4 | Modular design: each component independently testable | Unit test structure |

---

## Data Requirements

### Input Data

| Data | Format | Source | Required |
|------|--------|--------|----------|
| SNP matrix | CSV/HDF5 | Pre-processed from VCF | Yes |
| Phenotypes | CSV | Binary R/S per drug | Yes |
| Core genome alignment | FASTA | For phylogeny building | Yes* |
| Pre-built phylogeny | Newick | Alternative to alignment | Yes* |
| WHO catalogue | Excel/CSV | WHO website | Yes |
| CRyPTIC dataset | CSV/Parquet | PLOS Biology 2022 | Optional** |
| Reference genome | FASTA | H37Rv (GenBank) | Optional |

*Either alignment OR pre-built tree required

**CRyPTIC Consortium Dataset (Validation):
- **Size**: 12,289 MTB clinical isolates with WGS + MIC for 13 drugs
- **Source**: CRyPTIC Consortium 2022 (PLOS Biology)
- **DOI**: [10.1371/journal.pbio.3001721](https://doi.org/10.1371/journal.pbio.3001721)
- **Access**: Open source (no registration required)
- **Use case**: Independent validation dataset with quantitative MIC phenotypes

### Output Data

| Output | Format | Description |
|--------|--------|-------------|
| GWAS results | CSV | Per-variant statistics |
| Validation report | JSON | WHO comparison metrics |
| Manhattan plot | PNG/PDF | Genome-wide visualization |
| QQ plot | PNG/PDF | P-value distribution |
| Epistasis table | CSV | Interaction statistics |
| Model weights | PyTorch | Trained model (if applicable) |
| Run metadata | JSON | Git hash, seeds, config |

---

## Constraints

### Technical Constraints

- Must run on Apple Silicon (M1/M4) and x86_64
- No GPU required (CPU-only acceptable)
- Network access not required after initial setup
- Compatible with Python 3.10+

### Scientific Constraints

- LMM required for phylogenetic correction (simple Fisher test not acceptable)
- Multiple testing correction required (Bonferroni or FDR)
- Validation against external catalogue required
- Methods must be published/peer-reviewed

### Project Constraints

- 60 hours total development time
- Regular checkpoints with bioinformatician
- Results must be publication-ready

---

## Glossary

| Term | Definition |
|------|------------|
| GWAS | Genome-Wide Association Study |
| SNP | Single Nucleotide Polymorphism |
| LMM | Linear Mixed Model |
| MAF | Minor Allele Frequency |
| Lambda GC | Genomic inflation factor |
| PRPS | Phylogeny-Related Parallelism Score |
| Epistasis | Non-additive interaction between genetic variants |
| Kinship matrix | Matrix of pairwise genetic relatedness |

---

## Revision History

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | 2026-01-05 | Developer | Initial specification |
| 1.1 | 2026-01-06 | Developer | Added User Story 2.2 (Gene Burden Analysis) based on Farhat et al. 2019 |
| 1.2 | 2026-01-06 | Developer | Updated NFR-2 with achieved performance benchmarks |
| 1.3 | 2026-01-06 | Developer | Added CRyPTIC Consortium dataset (12,289 isolates) as validation data |
