# GWAS Reports Index

This directory contains analysis reports for each phase of the TB GWAS pipeline.

## Report Structure

```
reports/
├── README.md              # This index
├── phase1/                # Baseline GWAS reports
│   └── YYYY-MM-DD_drug_baseline.md
├── phase2/                # Multi-drug conditional GWAS
│   └── YYYY-MM-DD_drug_conditional.md
├── phase3/                # ABESS feature selection
│   └── YYYY-MM-DD_drug_abess.md
├── phase4/                # Epistasis analysis
│   └── YYYY-MM-DD_epistasis.md
└── checkpoints/           # Bioinformatician review records
    └── checkpoint_N_decisions.yaml
```

---

## Reports

### Phase 1: Baseline GWAS

| Date | Drug | Report | Status | Key Findings |
|------|------|--------|--------|--------------|
| 2026-01-05 | Rifampicin | [phase1/2026-01-05_rifampicin_baseline.md](phase1/2026-01-05_rifampicin_baseline.md) | 🔬 Awaiting Review | rpoB top hit (p=1.2e-98), λGC=0.139 |

### Phase 2: Multi-Drug Conditional GWAS

| Date | Drug | Report | Status | Key Findings |
|------|------|--------|--------|--------------|
| - | INH | - | ⬜ Not started | - |
| - | FQ | - | ⬜ Not started | - |

### Phase 3: ABESS Feature Selection

| Date | Drug | Report | Status | Key Findings |
|------|------|--------|--------|--------------|
| - | All | - | ⬜ Not started | - |

### Phase 4: Epistasis Analysis

| Date | Analysis | Report | Status | Key Findings |
|------|----------|--------|--------|--------------|
| - | Pairwise | - | ⬜ Not started | - |

---

## Status Legend

- ⬜ Not started
- 🔄 In progress
- 🔬 Awaiting bioinformatician review
- ✅ Approved
- ❌ Needs revision

---

## Quick Links

### Latest Reports
- **Phase 1 (Current)**: [Rifampicin Baseline GWAS](phase1/2026-01-05_rifampicin_baseline.md)

### Figures
- [Manhattan Plot - Rifampicin](../figures/manhattan_rifampicin.png)
- [QQ Plot - Rifampicin](../figures/qq_rifampicin.png)
- [Combined Plots](../figures/rifampicin_gwas_plots.png)

### Raw Results
- [Rifampicin GWAS Results CSV](../rifampicin_gwas.csv)

---

## Review Checklist

For each report, the bioinformatician should verify:

- [ ] Lambda GC is acceptable (0.9-1.1 typical, discuss if outside range)
- [ ] Top hits include expected genes for the drug
- [ ] Effect sizes (OR) are biologically plausible
- [ ] No unexpected cross-drug contamination
- [ ] PRPS distribution is interpretable
- [ ] Gene annotations completed for significant hits
- [ ] Decision recorded in checkpoint file

---

*Last updated: 2026-01-05*
