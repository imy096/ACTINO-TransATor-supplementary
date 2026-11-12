# Source provenance and reproducibility log

**Primary archive:** `ACTINO_TRANSATOR_MASTER_v050_STEP030J_AUDITED_2026-11-12.zip`, available separately to the researcher. The curated GitHub-ready package is derived from independently audited tables and small, relevant genomic-sequence excerpts from this archive. All derivative figures were created from these source files; none contains invented numerical data.

## Provenance levels

- **P1 — Original computational output extracted from cumulative archive:** antiSMASH STEP014 row tables, STEP001 assembly-QC data, STEP030H Pfam/HMMER domain files, STEP030I genomic IPG and GenBank outputs and BLASTP alignments, STEP030J published-reference comparison.
- **P2 — Derived, recalculated outputs:** 20-genome merged Table S1, numeric-QC JSON, six figures, supplementary table copies, and editorial summaries.
- **P3 — Scholarly contextual sources:** DOI-verified methodological and biochemical articles; bibliography and DOI links are documented in `references/verified_reference_metadata.tsv`. This level supports interpretation, not the numerical values in Tables S1–S12.
- **P4 — Unresolved external evidence:** all BioSample country/locality metadata, raw complete antiSMASH run archives from every assembly, prepared-but-unrun BiG-SCAPE family clustering, functional ACP assays and metabolite identification. These are not fabricated or presented as completed.

## Running the reproduction scripts

1. Install dependencies from `requirements.txt` in an isolated environment.
2. `python build_supplementary.py` reproduces the twelve supplementary tables and numeric-check JSON from the curated archived source tables.
3. `python reproduce_figures.py` rebuilds figures 1–6 in SVG/PDF/PNG formats; output is placed in `figures/`.
4. Scientific-document generation of the formatted Word file used `pandoc` with a local reference DOCX template; the canonical editable manuscript is `manuscript/ACTINO_Option_B_Full_Manuscript.md`. Two large tables remain editable in the Markdown source and original CSVs. Designed table plates were used in the review-export DOCX/PDF to avoid table-layout incompatibilities.
5. `docs/independent_numeric_checks.json` records the core audited values and should be regenerated after any new raw data are introduced.

## Data unit caveats

Genome assemblies and BioSamples are different concepts from `WP_` nonredundant protein accessions and assembly-specific CDS loci. The 18 comparator loci are not a random sample of genomes or necessarily independent evolutionary events. The 202/202 non-self positional HSP count excludes 13 self-comparisons and six unobservable positions. The PF00550 score threshold for exploratory candidates is not equivalent to a formal experimental assay.

## File checksums

See `docs/SHA256_MANIFEST.tsv` for each distributed file's SHA-256, length and relative path. Hashes are calculated over the exact bytes packaged, not inferred from earlier file names. The separately preserved MASTER v050 archive has its own independent SHA-256 checksum in the original audit.
