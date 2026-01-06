# Phase 4 Improvement Details

Advanced analysis improvements for epistasis and compensatory mutation detection, based on literature review (2026-01-06).

**Related**: [03_tasks.md](03_tasks.md#phase-4-epistasis-analysis)

---

## Summary

| Task | Description | Literature | Status |
|------|-------------|------------|--------|
| [4.5](#convergent-evolution-analysis) | Convergent evolution (phyC) | Chen & Shapiro 2015 | ⬜ |
| [4.6](#pfam-domain-aggregation) | PFAM domain aggregation | Reshetnikov 2025 | ⬜ |
| [4.7](#3d-structure-validation) | 3D structure validation | Phelan 2016 | ⬜ |

---

## Convergent Evolution Analysis

**Task**: 4.5
**Literature**: [Chen & Shapiro 2015](https://doi.org/10.1016/j.mib.2015.03.002)

### Problem

Traditional GWAS identifies variants correlated with phenotype, but may miss variants that arose **independently multiple times** on different lineages (convergent evolution).

### Key Finding from Literature

> "We compare the traditional GWAS against **phyC**, a contrasting method of mapping genotype to phenotype based upon evolutionary convergence"

### Why This Matters for TB

Compensatory mutations (e.g., rpoC mutations after rpoB resistance) often:
- Arise independently on multiple branches
- Have lower GWAS p-values due to lower frequency
- But have HIGH convergence scores (multiple independent origins)

### Implementation

```python
def compute_convergence_score(
    tree: TreeNode,
    variant_presence: np.ndarray,
    phenotype: np.ndarray,
) -> float:
    """Count independent origins of a mutation associated with phenotype.

    Algorithm:
    1. Map mutations and phenotypes to tree leaves
    2. Reconstruct ancestral states (parsimony)
    3. Count number of independent mutation events
    4. Weight by phenotype association at each origin

    Returns:
        Convergence score (higher = more independent origins)
    """
    # Use Fitch parsimony for ancestral reconstruction
    ancestral = reconstruct_ancestral_states(tree, variant_presence)

    # Count transitions from 0→1 that occur on resistant branches
    n_origins = 0
    for node in tree.traverse():
        if node.is_leaf():
            continue
        parent_state = ancestral[node.up.name] if node.up else 0
        node_state = ancestral[node.name]
        if parent_state == 0 and node_state == 1:
            # Check if associated with resistance
            if is_resistance_associated(node, phenotype):
                n_origins += 1

    return n_origins
```

### Expected Results

| Variant Type | GWAS Rank | Convergence Rank | Interpretation |
|--------------|-----------|------------------|----------------|
| rpoB S450L | 1 | 1 | High both |
| rpoC compensatory | 50-100 | 5-10 | Low GWAS, high convergence |
| Lineage marker | 10-20 | 200+ | High GWAS, low convergence (hitchhiking) |

---

## PFAM Domain Aggregation

**Task**: 4.6
**Literature**: [Reshetnikov et al. 2025](https://doi.org/10.3389/fmicb.2025.1586476)

### Problem

Rare variants lack statistical power individually, but may aggregate meaningfully at the **protein domain level**.

### Key Finding from Literature

> "Aggregating rare mutations within protein-coding genes into markers indicative of changes in **PFAM domains improved prediction quality**"

### PFAM Domains in Key TB Genes

| Gene | PFAM Domain | Domain Function |
|------|-------------|-----------------|
| rpoB | RNA_pol_Rpb1_2 | RNA polymerase active site |
| rpoB | RNA_pol_Rpb1_3 | Rifampicin binding pocket |
| katG | Peroxidase | Catalytic site (INH activation) |
| gyrA | DNA_gyraseA | Quinolone resistance determining region |
| embB | Glyco_transf_4 | Ethambutol binding site |

### Implementation

```python
def aggregate_by_domain(
    variants: pd.DataFrame,
    domain_map: dict[str, list[tuple[int, int, str]]],  # gene -> [(start, end, domain_name)]
) -> pd.DataFrame:
    """Aggregate variants by PFAM domain.

    Args:
        variants: DataFrame with position, gene columns
        domain_map: Gene to list of (start, end, domain_name) tuples

    Returns:
        DataFrame with domain-level binary indicators
    """
    domain_matrix = {}

    for gene, domains in domain_map.items():
        gene_variants = variants[variants["gene"] == gene]
        for start, end, domain_name in domains:
            key = f"{gene}:{domain_name}"
            # Binary: any non-synonymous variant in domain
            has_variant = gene_variants["position"].between(start, end).any()
            domain_matrix[key] = has_variant

    return pd.DataFrame(domain_matrix)
```

### Required Data

```bash
# Download PFAM annotations for H37Rv proteins
wget -O data/external/h37rv_pfam.tsv \
  "https://www.uniprot.org/uniprotkb?query=proteome:UP000001584&format=tsv&fields=accession,gene_names,ft_domain"
```

---

## 3D Structure Validation

**Task**: 4.7
**Literature**: [Phelan et al. 2016](https://doi.org/10.1186/s12916-016-0575-9)

### Problem

GWAS identifies statistical associations, but doesn't validate **mechanistic plausibility**. Mutations far from drug binding sites are less likely to confer resistance.

### Key Finding from Literature

> "A strong direct correlation was observed between the minimum inhibitory concentration values and the **distance of the mutated residues in 3D structures** to their respective drug binding sites"

### Available PDB Structures

| Gene | PDB ID | Bound Drug | Resolution |
|------|--------|------------|------------|
| rpoB | 5UHB | Rifampicin | 3.4 Å |
| katG | 1SJ2 | Isoniazid | 2.0 Å |
| gyrA | 5BS8 | Moxifloxacin | 2.8 Å |
| InhA | 4TZK | Isoniazid | 1.9 Å |
| EmbB | 7BVF | Ethambutol | 2.3 Å |

### Implementation

```python
from Bio.PDB import PDBParser, NeighborSearch

def calculate_binding_distance(
    pdb_path: str,
    mutation_residue: int,
    ligand_id: str,
) -> float:
    """Calculate distance from mutation site to drug binding site.

    Args:
        pdb_path: Path to PDB file
        mutation_residue: Residue number of mutation
        ligand_id: Three-letter code of bound drug (e.g., "RFP" for rifampicin)

    Returns:
        Minimum distance in Angstroms to any ligand atom
    """
    parser = PDBParser()
    structure = parser.get_structure("protein", pdb_path)

    # Find ligand atoms
    ligand_atoms = []
    for residue in structure.get_residues():
        if residue.get_resname() == ligand_id:
            ligand_atoms.extend(residue.get_atoms())

    # Find mutation residue
    mutation_atoms = []
    for residue in structure.get_residues():
        if residue.get_id()[1] == mutation_residue:
            mutation_atoms.extend(residue.get_atoms())

    # Calculate minimum distance
    min_dist = float("inf")
    for m_atom in mutation_atoms:
        for l_atom in ligand_atoms:
            dist = m_atom - l_atom  # Biopython distance operator
            min_dist = min(min_dist, dist)

    return min_dist
```

### Expected Correlation

```
High-effect mutations: rpoB S450L → 4.2 Å from rifampicin
Low-effect mutations: rpoB far-site → 15+ Å from rifampicin
```

### Validation Threshold

| Distance | Interpretation |
|----------|----------------|
| < 5 Å | Direct drug contact - high confidence |
| 5-10 Å | Binding pocket - plausible |
| 10-15 Å | Secondary effects - needs evidence |
| > 15 Å | Unlikely direct effect - flag for review |

---

## Literature References

### Primary Sources

1. **Chen & Shapiro 2015** - Bacterial GWAS advent
   - DOI: [10.1016/j.mib.2015.03.002](https://doi.org/10.1016/j.mib.2015.03.002)
   - PMID: 25835153
   - Key: phyC convergence method, GWAS comparison

2. **Reshetnikov et al. 2025** - Feature selection for TB GWAS
   - DOI: [10.3389/fmicb.2025.1586476](https://doi.org/10.3389/fmicb.2025.1586476)
   - PMID: 40606161
   - Key: ABESS, PFAM domain aggregation

3. **Phelan et al. 2016** - WGS + protein structure modeling
   - DOI: [10.1186/s12916-016-0575-9](https://doi.org/10.1186/s12916-016-0575-9)
   - PMID: 27005572
   - Key: 3D distance correlation with MIC

### Supporting Sources

4. **Saber & Shapiro 2020** - Benchmarking bacterial GWAS
   - DOI: [10.1099/mgen.0.000337](https://doi.org/10.1099/mgen.0.000337)
   - PMID: 32100713
   - Key: BacGWASim, elastic net superiority

5. **Fang et al. 2025** - WGS + ML for RR-TB
   - DOI: [10.3389/fcimb.2025.1641385](https://doi.org/10.3389/fcimb.2025.1641385)
   - PMID: 40851801
   - Key: Random Forest, embB M306I, katG S315T
