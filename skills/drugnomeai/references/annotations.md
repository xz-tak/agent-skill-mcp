# DrugnomeAI Biological Annotation Reference

Definitions for the classification systems behind the 11 DrugnomeAI runs, based on:
- Vitsios & Petrovski, *Commun Biol* 2022; doi:10.1038/s42003-022-04245-4
- NIH/NCATS PHAROS (Illuminating the Druggable Genome)
- Open Targets tractability assessments
- Codebase data files

---

## 1. PHAROS Target Development Level (TDL)

| Category | Full Name | Definition | Criteria | Genes in Data |
|----------|-----------|------------|----------|---------------|
| **Tclin** | Target Clinical | Targets of approved drugs with known mechanism of action | Has >=1 approved drug (DrugCentral) with a known MoA | 703 |
| **Tchem** | Target Chemical | Targets with potent chemical tools but no approved drug | Has bioactive compound with activity <=30nM (ChEMBL) or selective compound (<=100nM) but no approved drug | 1,900 |
| **Tbio** | Target Biological | Targets with functional annotation but no chemical tools at Tchem threshold | Has OMIM confirmed phenotype or GO experimental annotation, but doesn't meet Tchem criteria | 12,242 |
| **Tdark** | Target Dark | Understudied/uncharacterized proteins | Fails all above criteria; limited literature, no qualifying GO, no OMIM phenotype | 5,236 |

**Source:** NIH/NCATS PHAROS database (Illuminating the Druggable Genome project)
**Data file:** `drugnome_ai/data/PHAROS/pharos_merged_03112026.csv` (n=20,081)

---

## 2. Open Targets Druggability Tiers

| Tier | Definition | Criteria | Genes in Data |
|------|------------|----------|---------------|
| **Tier 1** | Clinical precedence | Target of an approved drug or clinical candidate (small molecule or biologic); highest tractability evidence | 1,427 |
| **Tier 2** | Discovery precedence | Active compound in biochemical/cell assays OR predicted highly druggable structure, but no clinical candidate yet | 682 |
| **Tier 3A** | Predicted tractable (structure) | Predicted druggable by protein structure: has modeled binding pocket, belongs to a druggable protein family, or structural homolog is a drug target | 870 |
| **Tier 3B** | Predicted tractable (other) | Weaker tractability evidence: gene is in a druggable pathway, has favorable expression profile, or other indirect druggability indicators | 1,500 |

**Source:** Open Targets tractability assessments
**Data file:** `drugnome_ai/data/labels/gene_druggable_labels.csv` (n=4,479)
**Additional columns in label data:** `small_mol_druggable` (Y/N), `bio_druggable` (Y/N), `adme_gene` (Y/N)

---

## 3. Therapeutic Modalities

| Modality | Definition | Target Class | Gene List Source | Genes |
|----------|------------|-------------|-----------------|-------|
| **Small Molecule** | Traditional small chemical compounds (MW <900 Da) that bind intracellular/extracellular targets | Enzymes, GPCRs, ion channels, nuclear receptors, kinases | Open Targets B1-B3 clinical buckets (approved to Phase 1) | 1,122 |
| **Antibody** | Monoclonal antibodies or antibody-derived biologics targeting cell-surface or secreted proteins | Extracellular/membrane proteins, secreted ligands, receptors | Open Targets antibody tractability annotations | 319 |
| **PROTAC** | Proteolysis-targeting chimeras that recruit E3 ubiquitin ligases to degrade intracellular targets | Intracellular proteins (transcription factors, scaffolds, oncoproteins) | Curated list of PROTAC-validated and predicted targets | 265 |

**Gene list files:** `misc/gene_lists/small_moelcules_genes.txt`, `antibody_genes.txt`, `protac.txt`

---

## 4. Full 11-Run Annotation

| # | Run Name | Dimension | CLI Flag | Positive Label Definition | Positive Genes | Biological Question Answered |
|---|----------|-----------|----------|---------------------------|----------------|------------------------------|
| 1 | `pharos_tclin` | PHAROS | `-p tclin` | Genes classified Tclin (approved drug targets) | ~613* | "Does this gene look like a known clinical drug target?" |
| 2 | `pharos_tchem` | PHAROS | `-p tchem` | Genes classified Tchem (chemical tool targets) | ~1,598* | "Does this gene have properties similar to targets with potent chemical tools?" |
| 3 | `pharos_tclin_tchem` | PHAROS | `-p tclin tchem` | Genes classified Tclin OR Tchem | ~2,211* | "Does this gene resemble targets with any drug/chemical engagement?" |
| 4 | `tier_1` | Tier | `-t 1` | Tier 1 druggable genes (clinical precedence) | 1,427 | "Does this gene look like a target with clinical drug precedence?" |
| 5 | `tier_1_2` | Tier | `-t 1 2` | Tier 1 OR Tier 2 genes | 2,109 | "Does this gene resemble a target with clinical or discovery-level tractability?" |
| 6 | `tier_1_2_3A` | Tier | `-t 1 2 3A` | Tier 1, 2, OR 3A genes | 2,979 | "Does this gene show structural or clinical druggability features?" |
| 7 | `tier_1_2_3A_3B` | Tier | `-t 1 2 3A 3B` | All 4 tiers (any tractability evidence) | 4,479 | "Does this gene have any druggability indicator at all?" |
| 8 | `tclin_tier1_intersect` | Intersection | `-k tclin_tier1_intersect.txt` | Genes that are BOTH Tclin AND Tier 1 (strictest set) | 524 | "Does this gene match the gold-standard: approved drug + clinical precedence?" |
| 9 | `modality_small_mol` | Modality | `-k small_moelcules_genes.txt` | Small molecule tractable targets | ~885* | "Is this gene druggable by a small molecule?" |
| 10 | `modality_antibody` | Modality | `-k antibody_genes.txt` | Antibody tractable targets | ~251* | "Is this gene targetable by an antibody/biologic?" |
| 11 | `modality_protac` | Modality | `-k protac.txt` | PROTAC-degradable targets | ~266* | "Is this gene a candidate for PROTAC-mediated degradation?" |

*Asterisk: actual positive gene count after merging with HGNC gene list (may differ slightly from raw list count)

---

## 5. PU Learning Framework Interpretation

| Term | Meaning |
|------|---------|
| **Positive (P)** | Genes with known druggability evidence for that run's definition (known_gene=1) |
| **Unlabelled (U)** | All remaining genes -- NOT assumed negative, just lacking evidence. The key insight of PU learning. |
| **Probability (proba)** | Mean prediction probability across all stochastic iterations; higher = more likely to be druggable |
| **Percentile (perc)** | Rank percentile (0-100) among all 20,080 genes; 95th = top 5% |
| **Novel gene** | A gene that was in the Unlabelled set (not Positive); high-scoring novel genes are predictions |
| **Known gene** | A gene that was in the Positive training set; high scores validate the model, low scores flag potential mislabels |

---

## 6. Interpreting Per-Run Results for a Gene

When reviewing a gene's profile across the 11 runs:

| Pattern | Interpretation |
|---------|----------------|
| High proba across all PHAROS runs | Gene has features consistent with established drug targets (even if not yet one) |
| High in Tclin but low in Tchem | Gene looks like a clinical target but lacks features of chemical-tool targets (may be biologic-only) |
| High across all tier runs | Strong structural/tractability signal |
| High in Tier 1 but drops off in Tier 1+2+3A+3B | Gene strongly resembles top-tier targets specifically, not generalized tractability |
| High in Tier 3A/3B runs but low in Tier 1 | Gene has predicted structural druggability but doesn't match established drug target features |
| Highest modality = Antibody | Gene likely encodes an extracellular/secreted protein -- classical antibody target |
| Highest modality = Small Molecule | Gene likely has an intracellular druggable pocket -- classical enzyme/receptor/kinase |
| Highest modality = PROTAC | Gene may lack a traditional binding pocket but is degradable -- suits undruggable targets |
| High modality_specificity (>0.15) | Strong evidence for one specific modality over others |
| Low modality_specificity (<0.05) | Multi-modal candidate -- viable via multiple approaches |
| Novel everywhere + high composite | **Discovery candidate** -- model predicts druggability despite zero prior evidence |
| Known in some runs, novel in others | Gene has partial evidence -- known for some criteria but predicted novel for others |
