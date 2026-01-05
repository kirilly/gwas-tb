# TB GWAS

Genome-Wide Association Study pipeline for identifying genetic variants associated with drug resistance in *Mycobacterium tuberculosis*.

WHO estimates 10 million TB cases annually, with ~500,000 showing drug resistance. This pipeline implements phylogenetically-corrected GWAS to find resistance-associated mutations while accounting for bacterial population structure.

## Quick Start

### Installation

```bash
# Install pixi (package manager)
curl -fsSL https://pixi.sh/install.sh | bash

# Clone and install dependencies
git clone https://github.com/kirilly/gwas-tb.git
cd gwas-tb
pixi install
```

**Platform notes:**
- **ARM Mac (M1/M2/M3)**: Uses Rosetta emulation automatically for pyseer
- **AMD64 (Linux/Intel Mac)**: Native execution

### Run GWAS Analysis

```python
import pandas as pd
import numpy as np
from src.config import GWASConfig
from src.gwas.runner import GWASRunner

# Load data
snps = pd.read_csv('data/processed/snps.csv', index_col=0)
phenotype = pd.read_csv('data/processed/phenotypes.csv', index_col=0)['Rifampicin']
kinship = np.load('data/processed/kinship.npy')

# Run GWAS
config = GWASConfig(method='lmm', min_maf=0.01)
runner = GWASRunner(config, n_jobs=8)
result = runner.run(snps, phenotype, kinship, drug_name='Rifampicin')

# Results
print(f'Lambda GC: {result.lambda_gc:.3f}')
print(f'Significant (p<5e-8): {result.n_significant}')
print(result.get_top_hits(10))
```

### Performance

| Dataset | GWAS | PRPS | Total |
|---------|------|------|-------|
| 11K samples × 7K variants | ~2.5 min | <3 sec | ~3 min |

## Documentation

See `docs/` for detailed specifications:
- `01_requirements.md` - Requirements and acceptance criteria
- `02_design.md` - Architecture and component design
- `03_tasks.md` - Implementation tasks and checkpoints
- `runbook.md` - Testing and troubleshooting
