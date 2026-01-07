# Open-Source GWAS Clinical Service for Russia

Research compiled: 2026-01-07

---

## Executive Summary

This document consolidates research on establishing a non-profit WGS/NGS clinical genomics service for Russian clinics, with initial focus on TB drug resistance and potential expansion to other high-value applications.

**Verdict**: Viable with hybrid funding model, but requires careful navigation of sanctions, equipment sourcing (MGI/BGI), and data localization requirements.

---

## Part 1: Global Context

### Clinical GWAS Landscape

| Aspect | Status |
|--------|--------|
| Primary application | Polygenic Risk Scores (PRS) |
| Key gap | 95%+ GWAS variants are non-coding, mechanisms unknown |
| Population bias | Only 20% non-European samples in GWAS studies |
| Regulatory | FDA/CLIA (US), CE-IVDR (EU), GOST R 72067 (Russia from 2026) |

### Existing Tools Comparison

| Tool | Type | Clinical Readiness | Link |
|------|------|-------------------|------|
| **TBProfiler** | TB resistance | High (WHO catalog, 94% MDR-TB sensitivity) | [tbdr.lshtm.ac.uk](https://tbdr.lshtm.ac.uk/) |
| **pyseer** | Bacterial GWAS | Research-grade | [pyseer.readthedocs.io](https://pyseer.readthedocs.io/) |
| **BOLT-LMM** | Large cohort GWAS | Research (needs 5000+ samples) | [alkesgroup.broadinstitute.org](https://alkesgroup.broadinstitute.org/BOLT-LMM/) |
| **Fabric Genomics** | Clinical interpretation | FDA-recognized | [fabricgenomics.com](https://www.fabricgenomics.com/) |
| **SOPHiA DDM** | Clinical genomics | High (750+ institutions) | [sophiagenetics.com](https://www.sophiagenetics.com/) |
| **WHO Mutation Catalogue** | TB resistance reference | Standard | [who.int](https://www.who.int/publications/i/item/9789240082410) |

**Gap**: No open-source clinical-grade GWAS with audit trails, EHR integration, and regulatory compliance.

### TB Clinical Significance

| Metric | Value |
|--------|-------|
| Accuracy (RIF/INH) | 98% sensitivity/specificity |
| Speed advantage | WGS 36 days faster than phenotypic DST |
| WHO recommendation | NGS for TB drug resistance (2024 guidance) |
| Reference implementation | UK universal WGS for all TB (2017) |

### Market Size

| Segment | 2024-2025 | 2030-2034 | CAGR |
|---------|-----------|-----------|------|
| Clinical genomics | $1.06B | $5.3B | 17.5% |
| TB diagnostics | $2.5B | $3.4B | 6% |
| WGS market | $2.2B | $16.5B | 22% |

---

## Part 2: Russia-Specific Analysis

### TB Burden in Russia

| Indicator | Value |
|-----------|-------|
| TB incidence | ~48-83 per 100,000 |
| Primary MDR-TB rate | 26.4% (Omsk, 2021) |
| Pre-XDR among MDR | 40.4% |
| XDR among MDR | 19.2% |
| XDR trend | Increased 16x since 2012 |
| Dominant genotype | 51% Beijing lineage |

Russia is among WHO's 30 high-burden countries for MDR/RR-TB.

### Current TB Diagnostics (Not WGS)

| Method | Platform | Use |
|--------|----------|-----|
| GeneXpert MTB/RIF | Cepheid | Primary diagnosis, RIF resistance |
| Multiplex PCR | Various | Drug resistance |
| DNA strips | Hain Lifescience | Line probe assays |
| Biological microchips | Domestic | Spoligotyping |

**WGS Status**: Research only, not clinical routine.

### Key Institutions

| Institution (Russian) | English | Focus | Link |
|----------------------|---------|-------|------|
| ЦНИИТ | Central TB Research Institute | WHO Collaborating Center, molecular diagnostics | - |
| ЦНИИ эпидемиологии | Central Research Institute of Epidemiology | AMR genomic surveillance | [crie.ru](https://www.crie.ru/en/) |
| НИЦ им. Гамалеи | Gamaleya Institute | GamTBvac, molecular biology | [gamaleya.org](https://gamaleya.org/en/) |
| ГНЦ ВБ "Вектор" | Vector State Research Center | Biosafety, virus genomics | [vector.nsc.ru](https://vector.nsc.ru/en/) |
| РЦМГ | Research Centre for Medical Genetics | Rare diseases, clinical genetics | [med-gen.ru](https://med-gen.ru/en/) |

### Key Researchers

- **A.E. Ergeshov** (А.Э. Эргешов) - CTRI, drug-resistant TB diagnostics
- **L.N. Chernousova** (Л.Н. Черноусова) - CTRI, molecular diagnostics
- **S.N. Andreevskaya** (С.Н. Андреевская) - CTRI, molecular epidemiology

### Regulatory Environment

| Regulation | Scope |
|------------|-------|
| Federal Law 152-FZ | Personal data protection |
| GOST R 72067-2025 | Medical lab standards (effective Jan 2026) |
| Data localization | All personal data must stay in Russia (July 2025) |

**Implication**: Cloud-based international analysis tools may face compliance issues. Local infrastructure required.

### Sanctions Impact (Critical)

| Challenge | Status |
|-----------|--------|
| Equipment dependency | 80% foreign pre-2022 |
| Reagent dependency | 90% Western |
| Illumina/Thermo | Paused sales in Russia |
| Price impact | DNA sequencing consumables tripled |

**Solution**: MGI/BGI (Chinese) equipment via OOO Helicon distributor. MGISEQ-2000 comparable to Illumina HiSeq 2500.

---

## Part 3: Other WGS Opportunities

### Priority Applications for Russian Clinics

| Application | Demand | Revenue | Complexity |
|-------------|--------|---------|------------|
| Oncology tumor profiling | Highest | High | High |
| Rare disease diagnostics | High (govt priority) | Medium | Medium |
| NIPT | Growing | Medium | Low |
| AMR surveillance | Medium | Low-Medium | Low |
| Pharmacogenomics | Emerging | Medium | Medium |

### Government Initiatives

- **"100,000 + Me"**: National genome initiative, 100K genomes by end 2025
- **1 million genomes by 2030**: Long-term target
- **EXAMEN Project**: WES-based newborn screening pilot (found 0.9% actionable variants)
- **mRNA cancer vaccines**: Clinical launch October 2025

---

## Part 4: Business Model

### Successful Russian Model: RUSeq Consortium

- First integration of genetic data across major Russian labs
- Partners: Kulakov Centre, City Hospital No. 40, GENETICO Ltd., CerbaLab Ltd.
- Output: 7,452 exome samples, open Russian allele frequency database
- Funding: Joint commercial partner funding

### Recommended Hybrid Model

1. **RSF grants** for infrastructure (up to 30M rubles/year)
2. **Sample-for-service** with clinics
3. **Bioinformatics consulting** as revenue stream
4. **Build Russian variant database** as research output

### Funding Sources

| Source | Amount | Notes |
|--------|--------|-------|
| RSF small grants | Up to 1.5M rubles/year | |
| RSF research teams | 20-30M rubles (2023-2026) | |
| RSF world-class labs | Up to 30M rubles ($357K) | 50%+ team under 39 years |
| International | Limited | China, Vietnam, Iran, India only |

### Competitive Landscape

| Company | Focus | Notes |
|---------|-------|-------|
| Genotek | 5000+ gene tests | Largest in Eastern Europe, $7.6M funding |
| Atlas Biomed | DTC genomics | UK-based, Russian origins |
| INVITRO | Clinical trials | 300+ trials since 2003 |

---

## Part 5: Implementation Roadmap

### Short-Term (6-12 months)

1. Contact RUSeq consortium participants
2. Attend GENETICS 2025 conference (November 5-7)
3. Assess MGI/BGI equipment via OOO Helicon
4. Partner with Kulakov Centre (EXAMEN experience)

### Medium-Term (1-2 years)

1. Start with rare disease + TB diagnostics
2. Build Russian variant frequency database
3. Develop bioinformatics training component
4. Seek RSF funding for infrastructure

### Critical Success Factors

1. Secure MGI/BGI equipment and reagent supply
2. Build multidisciplinary interpretation team
3. Partner with established institution for compliance
4. Focus on government-funded applications

---

## Part 6: Key Contacts

### Institutions

| Organization | Contact Purpose | Link |
|--------------|-----------------|------|
| ЦНИИТ (CTRI) | TB genomics partnership | - |
| ЦНИИ эпидемиологии | General genomics, AMR | [crie.ru](https://www.crie.ru/en/) |
| Kulakov Centre | RUSeq model, rare diseases | [ncagp.ru](https://www.ncagp.ru/en/) |
| OOO Helicon | MGI equipment/reagents | [helicon.ru](https://helicon.ru/) |
| RSF | Grant funding | [rscf.ru](https://rscf.ru/en/) |
| РЦМГ | Rare disease diagnostics | [med-gen.ru](https://med-gen.ru/en/) |
| Genotek | Commercial genomics (competitor) | [genotek.ru](https://www.genotek.ru/) |

### Conferences

| Event | Date | Focus | Link |
|-------|------|-------|------|
| GENETICS 2025 | Nov 5-7, 2025 | Vavilov Institute | [vavilovgenetics.org](http://vavilovgenetics.org/en) |
| XIX Congress Reproductive Medicine | Jan 21-24, 2025 | Reproductive genetics | - |
| TB-Profiler Workshop | Various | WHO-backed training | [tbdr.lshtm.ac.uk](https://tbdr.lshtm.ac.uk/) |

### Questions Requiring Direct Contact

**To CTRI:**
- Current WGS capabilities
- Willingness for sample-sharing partnership

**To MGI/Helicon:**
- Equipment pricing in Russia
- Reagent supply reliability

**To RSF:**
- Funding for clinical genomics services (not just research)

**To Regional TB dispensaries:**
- What tests are currently sent abroad?
- What turnaround time/cost is acceptable?

---

## Part 7: Risk Assessment

| Risk | Severity | Mitigation |
|------|----------|------------|
| Sanctions (equipment) | High | Use MGI/BGI via Helicon |
| Sanctions (reagents) | High | Chinese supply chain |
| Data localization | Medium | Local infrastructure from start |
| Regulatory changes | Medium | Partner with established institution |
| Funding sustainability | Medium | Hybrid grant + service model |

---

## Sources

### Global GWAS/Clinical Genomics
- [WHO Global TB Report 2024](https://www.who.int/teams/global-programme-on-tuberculosis-and-lung-health/tb-reports/global-tuberculosis-report-2024)
- [Saber & Shapiro 2020 - BacGWASim](https://doi.org/10.1099/mgen.0.000337) - elastic net for clonal bacteria
- [Phelan et al. 2016](https://doi.org/10.1186/s12916-016-0575-9) - 3D structure validation
- [CRyPTIC Consortium 2022](https://doi.org/10.1371/journal.pbio.3001755) - 10,228 TB genomes

### Russia-Specific
- [Western Siberia TB Epidemiology (MDPI 2023)](https://www.mdpi.com/2076-2607/11/2/425)
- [ЦНИИТ publications - Ergeshov et al.](https://vestnikramn.spr-journal.ru/jour/article/view/1163)
- [RUSeq Consortium (Oxford Academic)](https://academic.oup.com/nsr/article/11/10/nwae326/7758245)
- [Russian Science Foundation](https://rscf.ru/en/)
- [Russia 1M Genomes Initiative (Bio-IT World)](https://www.bio-itworld.com/news/2025/08/07/russia-targets-national-database-of-1-million-genomes-by-2030)
- [EXAMEN WES Newborn Screening Pilot](https://pubmed.ncbi.nlm.nih.gov/39033325/)

### Business Models
- [Seqera/Nextflow - $26M Series B](https://seqera.io/blog/seqera-raises-26m-series-b/)
- [Galaxy Project Funding](https://galaxyproject.org/galaxy-project/)
- [MongoDB Atlas Transformation](https://www.sandeepmvp.com/how-mongodb-transitioned-to-a-saas-business-with-atlas-svp-of-product-at-mongodb/)

### Regulatory & Sanctions
- [Russia Data Protection Law 152-FZ](https://securiti.ai/russian-federal-law-no-152-fz/)
- [New Medical Laboratory Standards GOST R 72067-2025](https://certru.ru/en/new-standard-for-medical-laboratories-in-russia-from-january-1-2026/)
- [Sanctions Impact on Russian Labs](https://sciencebusiness.net/news/russian-labs-run-out-equipment-sanctions-begin-bite)

---

## Appendix: Acronyms

| Acronym | Meaning |
|---------|---------|
| CTRI/ЦНИИТ | Central TB Research Institute |
| RSF | Russian Science Foundation |
| MDR-TB | Multi-drug resistant TB |
| XDR-TB | Extensively drug-resistant TB |
| WGS | Whole Genome Sequencing |
| NGS | Next-Generation Sequencing |
| NIPT | Non-Invasive Prenatal Testing |
| AMR | Antimicrobial Resistance |
