# TB GWAS Project

GWAS pipeline for TB drug resistance with LMM phylogenetic correction.

## Key Files

- `src/gwas/runner.py` - GWAS runner (pyseer Python API)
- `src/gwas/stats.py` - PRPS, lambda GC, FDR
- `docs/02_design.md` - Architecture and conventions
- `docs/runbook.md` - Testing and troubleshooting

## Performance

| Operation | Time |
|-----------|------|
| Eigendecomposition | ~2 min |
| GWAS (1000 var) | ~7 sec |
| PRPS (1000 var) | <1 sec |

## Development Guidelines

**Bioinformatics issues**: Always use debug-specialist agent for research first:
- Research the specific bioinformatics issue before implementing
- Test solutions with increasing data sizes (100 → 1000 → 10000)
- Verify correctness with correlation tests against reference implementations
