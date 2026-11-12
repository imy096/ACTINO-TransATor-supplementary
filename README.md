# ACTINO–TransATor | Supporting data and reproducibility materials

**Dataset for the Option B comparative genomics manuscript (working title):**
*Biosynthetic Protocluster Composition in a Twenty-Assembly Actinobacterial Collection: A Comparative Study of the CCDA1 trans-AT PKS Region*.

## What is in this repository?

This repository contains curated data and scripts supporting the **computational results** reported in the manuscript. It is intended to be linked from the article's **Availability of data and materials** section once the authors have reviewed the files and publicly released a stable version.

| Folder/file | Description |
|---|---|
| `supplementary/` | Tables S1–S12 and detailed supplementary methods explaining samples, results and definitions |
| `source_data/` | Row-level antiSMASH summaries; selected antiSMASH GenBank record; 18 comparative GenBank genomic neighborhoods; IPG mappings; BLASTP, Pfam and motif outputs; FASTA sequences |
| `figures/` | Six publication figures, each in vector SVG, vector PDF and PNG formats |
| `reproduce_figures.py` | Source code for regenerating the six figures from the included data |
| `build_supplementary.py` | Script for assembling supplementary material from curated tables |
| `requirements.txt` | Python package requirements for these scripts |
| `references/` | Bibliography records and reference metadata used in the manuscript package |
| `documentation/` | Claim-to-evidence map, provenance notes and independent numerical checks |
| `FILE_INVENTORY.tsv`, `SHA256SUMS.txt` | Inventory and cryptographic checksums of distributed files |

## Principal analysis scopes

**Main inventory:** 20 actinobacterial assemblies; 674 predicted antiSMASH regions and 878 antiSMASH *protoclusters*. These are annotation objects and must not be interpreted as 878 independent chemical products or nonredundant biosynthetic gene-cluster families.

**CCDA1 comparative investigation:** 18 assembly-specific reference loci from the external protein-comparison investigation. The 18 loci are **not** an additional 18 independent samples in the 20-genome collection. Gene-neighborhood annotations, large-protein alignments and ACP-like sequences have distinct, documented supporting files.

**Important:** The cited terminal ACP-like serines (S5766, S5909 and S6221) are candidate positions derived from sequence/profile alignment. No biochemical assay or metabolite production was performed.

## Reproduce the six publication figures

From the repository root:

```bash
python -m pip install -r requirements.txt
python reproduce_figures.py
```

This command regenerates visualization output from existing curated results. It does **not** rerun antiSMASH, original assembly calling, BLASTP or Pfam from raw genomes.

## Reuse of public reference data

Some sequence records, annotations and Pfam profiles originate from public external resources, including NCBI and InterPro/Pfam. Keep original accession identifiers and cite the underlying resources. InterPro's distributable Pfam data are provided under CC0; that does **not** automatically determine the licence chosen for original project scripts and tables.

## Data gaps and scientific limits

- The 674 antiSMASH region annotations and 878 protocluster objects are available as row-level machine-readable tables, **but not all 674 original antiSMASH GenBank region files or complete upstream antiSMASH output directories are included**.
- Original antiSMASH invocation parameters and input assembly provenance require final verification against primary records.
- No BiG-SCAPE gene-cluster-family dereplication has been conducted; do not claim a count of distinct gene-cluster families.
- Genomic predictions do not establish expressed enzyme activity, final metabolite identity or experimentally demonstrated novelty.
- Author metadata, permissions, release licence and any persistent DOI must be completed/approved by the authors before a final scholarly deposit.

## Versioning and article citation

Once the authors have validated this repository, create a **GitHub Release** (e.g. `v1.0.0`). For durable dataset citation, archive that public release with **Zenodo**, then cite its actual minted DOI in the article. Do not cite a placeholder as if it were a real DOI.

See [`GITHUB_ZENODO_RELEASE_GUIDE.md`](GITHUB_ZENODO_RELEASE_GUIDE.md) and [`DATA_AVAILABILITY_TEXT_TEMPLATE.md`](DATA_AVAILABILITY_TEXT_TEMPLATE.md).
