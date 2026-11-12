# Supplementary Information — ACTINO–TransATor Option B

**Companion to:** *Genome Mining of Twenty Actinobacterial Assemblies Reveals a Heterogeneous Predicted Biosynthetic Repertoire and a Conserved Leinamycin-Related PKS Locus*  
**Status:** Research draft, source-grounded, version 1.0, 12 November 2026.

## Supplementary Methods S1 — Data provenance and counting units

The primary source was the independently audited cumulative archive `ACTINO_TRANSATOR_MASTER_v050_STEP030J_AUDITED_2026-11-12.zip`. The present repository is a curated, reproducible subset of that archive; it is **not** a replacement for the complete upstream antiSMASH output archives. Table S1 combines the 20 genome FASTA-QC records with their antiSMASH per-genome classification counts. Tables S2 and S3 retain all 878 protocluster and 674 region records, respectively, including genomic coordinates, predicted product labels, contig-edge flags, and region-associated KnownClusterBlast reference fields. Table S4 contains the 911 candidate-cluster objects grouped by antiSMASH relationship type.

Three hierarchical annotation concepts must remain distinct: (a) region area (the broader genomic interval), (b) candidate-cluster object (a representation of one or more interacting or overlapping protoclusters), and (c) protocluster (a prediction based on a biosynthetic rule and its associated core/border). Their counts are not additive and none of the units is a measured metabolite. The archived product-class mapping was applied to the detailed antiSMASH product labels to generate PKS, NRPS, RiPP, terpene, other and saccharide groupings. **DO NOT** use “878 independent BGC families,” “878 predicted unique products,” or “878 experimentally verified pathways.”

The region table's edge flag is `region_touches_contig_edge`; the protocluster equivalent is `protocluster_touches_contig_edge`. These were summed from Boolean data fields, not inferred by visual inspection. The corresponding figures and checks are generated without making population-level assumptions.

## Supplementary Methods S2 — Accession mapping and deduplication

The selected case-study protein `MIK3493667.1` is from *Streptomyces* sp. CCDA1, GenBank assembly `GCA_060303205.1`. Fourteen additional `WP_` reference protein accessions were profiled. The NCBI IPG response can associate one nonredundant `WP_` sequence with multiple strain genomes; these should not be counted as independent protein sequences. Eighteen genome-assembly-specific reference loci were retained after reconciling `GCA_` and `GCF_` paired accession representations. Table S5 retains the selected protein, assembly and contig identifiers, CDS coordinates and strand. The original XML records from fourteen requests and all eighteen GenBank regions are provided in `source_data/IPG_XML_14` and `source_data/genomic_neighborhoods_18`.

The `GCA_040299725.1` GenBank-only locus carries an INSDC translation accession `MER6801813.1` rather than the RefSeq `WP_344643858.1` used to initiate IPG discovery. GenBank localization was checked by the accession **actually present in that CDS annotation**, plus exact coordinates and strand. An earlier attempt that required the `WP_` accession alone incorrectly flagged this valid GenBank-only record; the corrected mapping is the one used here.

## Supplementary Methods S3 — Synteny and missingness

Gene position was normalized relative to the central target protein: indices `−10 … −1` precede the target in its transcription direction; index `0` is the target; `+1 … +3` follow it. On negative-strand targets, genomic-coordinate order was reversed before comparison. Products were assigned to fourteen broad annotation families: DUF6039, ACP S-malonyltransferase, DUF6039, VOC, DUF3224, PIG-L, SAM-dependent methyltransferase, alpha/beta hydrolase, ACP, AMP/adenylation, target protein, HMG-CoA synthase, carrier protein, and KS-like protein. These positions are an analysis window, not independently proven orthologous operon boundaries. Table S6 retains all eighteen patterns and differences. Fifteen matched completely (14/14), two matched all eleven positions present in their contig intervals, and one displayed a partial disrupted target with shifted annotations.

The partial locus in `GCF_031181555.1` included three translated CDS products `WP_309061376.1` (330 aa), `WP_309061384.1` (5,494 aa) and `WP_309061375.1` (761 aa). Two exactly 100-base stretches of the contig were composed entirely of `N`. No literal `assembly_gap` feature appeared in the retrieved GenBank interval. This is evidence of unresolved assembly sequence but not of known physical lengths nor a demonstrable historical gene fragmentation event. Supplementary domain evidence is in `source_data/raw_alignments/PARTIAL_3CDS_Pfam17.domtblout`.

## Supplementary Methods S4 — BLASTP, coverage and self-comparisons

The complete study included eighteen reference neighborhoods with 610 annotated CDS features and 598 available translations. The BLASTP input files and archived TSV are present in `source_data/primary_sequences` and `source_data/raw_alignments`. Search thresholds in the retained comparative scripts used an E-value threshold of `1e−5`. The per-position exploratory criterion was at least 70% query and subject HSP coverage for the top-scoring positional hit; sequence identity, E-value and aligned span should be distinguished. Table S7 includes all 221 positional comparisons across seventeen loci: 215 hits and six missing/outside positions. **Thirteen** hits are self-comparisons of the CL1-1 reference and are excluded from the main between-locus statistics, yielding 202 evaluable non-self comparisons among 208 possible positions; the other six are missing due to contig boundaries. The 202 evaluable non-self comparisons all met the stated coverage criterion. This does not establish 202 formal one-to-one orthology relationships.

For large modular proteins, many BLASTP HSPs can align repeated catalytic modules at alternative positions. A strictly nonoverlapping and order-preserving HSP chain was selected to produce reproducible *colinear coverage* values. Coverage was calculated as the sum of selected HSP query residue spans divided by full query length, with the analogous calculation for subject length. These are neither overall sequence identities nor global alignment percentages. The summary values appear in Table S8; selected GnmT–CCDA1 and GnmT–Microbispora HSPs appear in Tables S10–S11. A global colinear alignment may yield different percentages and has not been claimed.

## Supplementary Methods S5 — Profile-HMM hit interpretation

The targeted HMMER experiment used a seventeen-family Pfam profile panel, saved in `source_data/PFAM_PANEL17.hmm`, with strict `--cut_ga` calls and a separate exploratory scan without `--cut_ga`. PF00550 terminal carrier-like fragments were investigated when the aligned segment was at least 45 amino acids and the domain score was at least 15 bits. This post hoc exploratory filter is not the official Pfam gathering threshold; for the PF00550 model used, the domain gathering score is 20.8 bits. The CCDA1 terminal fragment scored 27.2 bits, the *Microbispora* representative 20.3 and the published GnmT 19.9. GnmT and *Microbispora* candidates are below the chosen Pfam domain gathering score; reporting all three as experimentally established ACPs would overstate the evidence.

The serine aligned to Stockholm profile column 31 was mapped to the original untrimmed protein fragment, accounting for N-terminal residue trimming by `hmmalign --trim`. Confirmed coordinates were S5766 (CCDA1, `GLTSF`), S5909 (`WP_168066089.1`, `ALTSF`) and S6221 (GnmT, `GFTSF`). Only *sequence-level candidates* are established. Tandem ACP function, phosphopantetheinylation, carrier binding, β-branching activities and chemical product formation have not been measured for CCDA1.

## Supplementary Methods S6 — Reproducible scientific figures

Run `python reproduce_figures.py` from the repository root. The script uses only archived CSV/TSV/GenBank/Stockholm inputs and emits six scientific figures as editable SVG and print-ready PDF plus 360-dpi PNG. No data are simulated or generated to fill gaps. Figure 3 uses GenBank CDS starts, ends, strands, and relative target orientation. Figure 4 flags the partial target; Figure 5 distinguishes Pfam below-gathering exploratory calls by dashed outlines. Figure 6 is a sequence *alignment-context* graphic and does not imply successful biochemical loading.

### Figure captions and source traceability

**Figure 1:** Per-genome predicted protocluster composition. Source: `antismash_genome_category_totals.csv`; 20 genome rows.

**Figure 2:** AntiSMASH annotation hierarchy and edge exposure. Sources: `antismash_regions_674.csv`, `antismash_all_protoclusters_878.csv`, `candidate_cluster_kind_counts.csv`.

**Figure 3:** CDS genomic synteny of CCDA1, *Microbispora* CL1-1 and *Streptomyces fimbriatus*. Sources: archived GenBank files in `source_data`, including `CCDA1_antismash_region002.gbk`.

**Figure 4:** Annotation-level synteny and query colinear coverage across eighteen assemblies. Sources: `comparative_synteny_18.tsv` and `target_colinear_18.tsv`. Relative positions not observed are explicitly represented as unavailable.

**Figure 5:** Selected domain architectures of CCDA1, *Microbispora* and GnmT. Sources: strict and relaxed Pfam17 HMMER `domtblout` files.

**Figure 6:** Terminal ACP-like motif contexts in fifteen initial proteins plus GnmT. Sources: original Stockholm alignments; all sequences retain the profile-aligned serine.

## Complete Supplementary Table Index

| Table | File | Unit / description |
|---|---|---|
| S1 | `Table_S1_20_genomes_assemblies_and_inventory.csv` | All 20 accessions, BioSamples, assembly metrics, per-category protocluster counts |
| S2 | `Table_S2_all_878_protoclusters.csv` | Full row-level antiSMASH protocluster annotation with coordinates and category labels |
| S3 | `Table_S3_all_674_regions.csv` | All annotated region areas and their contig-edge, product and KnownClusterBlast fields |
| S4 | `Table_S4_candidate_cluster_types.csv` | Counts of `single`, `neighbouring`, `interleaved` and `chemical_hybrid` objects |
| S5 | `Table_S5_IPG_18_reference_loci.tsv` | Eighteen selected genomic loci mapped through IPG records |
| S6 | `Table_S6_synteny_18_loci.tsv` | Exact 14-position patterns and annotation deviations |
| S7 | `Table_S7_positional_BLASTP_17_loci.tsv` | 221 positional BLASTP comparisons, including self and missing positions |
| S8 | `Table_S8_large_target_colinear_coverage.tsv` | HSP-chain query and subject coverage across 18 genomic loci |
| S9 | `Table_S9_terminal_ACP_serines.tsv` | Three verified serine coordinates and their profile hit statistics |
| S10 | `Table_S10_GnmT_colinear_HSPs.tsv` | GnmT colinear selected HSP chain |
| S11 | `Table_S11_GnmT_BLASTP_summary.tsv` | GnmT-to-target protein BLASTP summaries |
| S12 | `Table_S12_annotated_product_type_counts.csv` | antiSMASH fine product labels and category counts |

## Data and reporting caveats

This supplementary package includes all available *relevant audit tables and selected original record excerpts* but is **not** a promise that every source dataset ever considered is present. Specifically, original complete antiSMASH output ZIPs for all 20 assemblies, all 674 region GenBank input records for BiG-SCAPE, experimentally measured products, verified collection localities, and author/funding declarations are not included. The original cumulative MASTER v050 ZIP is available as a separate offline archive and should be preserved unchanged. Repository deposit, authorship, licensing for upstream database files and journal disclosure policies require owner review.
