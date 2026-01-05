# TB GWAS Runbook

## Testing

### Unit Tests
```bash
~/.pixi/bin/pixi run pytest tests/
```

### Quick Validation (100 samples)
```bash
~/.pixi/bin/pixi run python scripts/validate_gwas.py --n-samples 100
```

### Full Validation (11K samples)
```bash
~/.pixi/bin/pixi run python scripts/validate_gwas.py
```

### Performance Testing
Test with increasing data sizes to validate optimizations:
```bash
# Small scale
~/.pixi/bin/pixi run python -c "
from src.gwas.stats import calculate_prps_batch
import numpy as np
# 1000 samples × 100 variants
"

# Medium scale
# 5000 samples × 500 variants

# Full scale
# 11000 samples × 1000 variants
```

## Troubleshooting

### pyseer installation fails
- On ARM Macs: pixi automatically uses Rosetta (osx-64 platform)
- Check `pixi.toml` has `platforms = ["osx-64", "linux-64"]`
- Verify: `~/.pixi/bin/pixi run pyseer --version`

### GWAS too slow
- Ensure using pyseer Python API (default), not CLI
- Check eigendecomposition is cached
- Use batch PRPS calculation (`calculate_prps_batch`)
- Expected: ~143 var/sec for GWAS, <1 sec for PRPS batch

### High Lambda GC (>1.1)
- Population structure not fully corrected
- Check kinship matrix is positive semi-definite
- Consider adding more covariates
- Review sample alignment between SNPs and kinship

### Memory issues
- For 11K samples: ~1GB for kinship matrix
- Eigendecomposition: ~2GB peak
- Use chunked processing for >20K samples

## Common Operations

### Run Full GWAS
```python
from src.gwas.runner import GWASRunner
from src.config import GWASConfig
import pandas as pd
import numpy as np

snps = pd.read_csv('data/processed/snps.csv', index_col=0)
pheno = pd.read_csv('data/processed/phenotypes.csv', index_col=0)
kinship = np.load('data/processed/kinship.npy')

config = GWASConfig(method='lmm', min_maf=0.01)
runner = GWASRunner(config, n_jobs=8)
result = runner.run(snps, pheno['Rifampicin'], kinship, drug_name='Rifampicin')
```

### Check Results
```python
print(f'Lambda GC: {result.lambda_gc:.3f}')
print(f'Significant: {result.n_significant}')
print(result.get_top_hits(10))
```
