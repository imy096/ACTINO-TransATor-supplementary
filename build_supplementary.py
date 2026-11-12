#!/usr/bin/env python3
"""Build audited supplementary tables, reference records, and quality-control summaries."""
from pathlib import Path
from datetime import date
import pandas as pd, csv, json, shutil, hashlib
ROOT=Path(__file__).resolve().parent;S=ROOT/'source_data';P=ROOT/'supplementary';R=ROOT/'references';D=ROOT/'docs';
for x in [P,R,D]:x.mkdir(exist_ok=True)
q=pd.read_csv(S/'genome_assembly_qc_20.csv');g=pd.read_csv(S/'antismash_genome_category_totals.csv')
r=q.merge(g,left_index=True,right_index=True,validate='one_to_one')
assert (r.organism==r.organism_name).all()
cols=['genome_number','strain_y','organism','gca','biosample','bioproject','assembly_level','observed_size_bp','observed_n50_bp','contig_count','gc_percent','total_protoclusters','PKS','NRPS','RiPP','terpene','other','saccharide']
r[cols].rename(columns={'strain_y':'strain','gca':'assembly_accession'}).to_csv(P/'Table_S1_20_genomes_assemblies_and_inventory.csv',index=False)
cp={
'antismash_all_protoclusters_878.csv':'Table_S2_all_878_protoclusters.csv',
'antismash_regions_674.csv':'Table_S3_all_674_regions.csv',
'candidate_cluster_kind_counts.csv':'Table_S4_candidate_cluster_types.csv',
'comparative_assembly_loci_18.tsv':'Table_S5_IPG_18_reference_loci.tsv',
'comparative_synteny_18.tsv':'Table_S6_synteny_18_loci.tsv',
'neighbor_positional_BLAST17.tsv':'Table_S7_positional_BLASTP_17_loci.tsv',
 'target_colinear_18.tsv':'Table_S8_large_target_colinear_coverage.tsv',
 'terminal_acp_serine_3refs.tsv':'Table_S9_terminal_ACP_serines.tsv',
 'gnmt_colinear_chains.tsv':'Table_S10_GnmT_colinear_HSPs.tsv',
 'gnmt_blast_summary.tsv':'Table_S11_GnmT_BLASTP_summary.tsv',
 'protocluster_product_class_counts.csv':'Table_S12_annotated_product_type_counts.csv',
}
for a,b in cp.items():shutil.copy2(S/a,P/b)
# Reference metadata limited to verified primary/method sources.
refs=[
('Blin2025','Blin K et al.','2025','antiSMASH 8.0: extended gene cluster detection capabilities and analyses of chemistry, enzymology, and regulation','Nucleic Acids Research','53(W1)','W32-W38','10.1093/nar/gkaf334','PubMed PMID 40276974','Method'),
('Zdouc2025','Zdouc MM et al.','2025','MIBiG 4.0: advancing biosynthetic gene cluster curation through global collaboration','Nucleic Acids Research','53(D1)','D678-D690','10.1093/nar/gkae1115','PubMed PMID 39657789','Database'),
('Pan2017','Pan G et al.','2017','Discovery of the leinamycin family of natural products by mining actinobacterial genomes','PNAS','114(52)','E11131-E11140','10.1073/pnas.1716245115','PubMed PMID 29229819','Primary precedent'),
('Phillips2025','Phillips AP et al.','2025','Divergent tandem acyl carrier proteins necessitate in-series polyketide processing in the leinamycin family','Angewandte Chemie International Edition','64(2)','e202414165','10.1002/anie.202414165','PubMed PMID 39324406','Primary biochemical precedent'),
('Eddy2011','Eddy SR','2011','Accelerated Profile HMM Searches','PLoS Computational Biology','7(10)','e1002195','10.1371/journal.pcbi.1002195','PubMed PMID 22039361','HMMER method'),
('Camacho2009','Camacho C et al.','2009','BLAST+: architecture and applications','BMC Bioinformatics','10','421','10.1186/1471-2105-10-421','PubMed PMID 20003500','BLAST+ method'),
('Medema2011','Medema MH et al.','2011','antiSMASH: rapid identification, annotation and analysis of secondary metabolite biosynthesis gene clusters in bacterial and fungal genome sequences','Nucleic Acids Research','39','W339-W346','10.1093/nar/gkr466','PubMed PMID 21672958','Background/method'),
('NavarroMunoz2020','Navarro-Muñoz JC et al.','2020','A computational framework to explore large-scale biosynthetic diversity','Nature Chemical Biology','16','60-68','10.1038/s41589-019-0400-9','Nature publisher indexed page','GCF method not run'),
('Draisma2026','Draisma A et al.','2026','BiG-SCAPE 2.0 and BiG-SLiCE 2.0: scalable, accurate and interactive sequence clustering of metabolic gene clusters','Nature Communications','17','2000','10.1038/s41467-026-68733-5','Nature publisher indexed page','GCF method not run'),
]
with (R/'verified_reference_metadata.tsv').open('w',newline='') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['key','authors_display','year','title','journal','volume_issue','pages_or_article','doi','verification_source','role']);w.writerows(refs)
with (R/'references.bib').open('w') as f:
 for key,auth,yr,title,journal,vol,pages,doi,source,role in refs:
  f.write('@article{%s,\n  author = {%s},\n  title = {%s},\n  year = {%s},\n  journal = {%s},\n  volume = {%s},\n  pages = {%s},\n  doi = {%s}\n}\n\n' % (key,auth.replace(' et al.',' and others'),title,yr,journal,vol,pages,doi))
q1=pd.read_csv(P/'Table_S7_positional_BLASTP_17_loci.tsv',sep='\t')
non=q1[q1.Assembly!='GCF_012034195.1']
checks={
 'genome_n':len(q),
 'distinct_biosample_n':int(q.biosample.nunique()),
 'region_n':len(pd.read_csv(S/'antismash_regions_674.csv')),
 'protocluster_n':len(pd.read_csv(S/'antismash_all_protoclusters_878.csv')),
 'candidate_cluster_objects':int(pd.read_csv(S/'candidate_cluster_kind_counts.csv')['count'].sum()),
 'category_totals':{k:int(g[k].sum()) for k in ['PKS','NRPS','RiPP','terpene','other','saccharide']},
 'region_contig_edge_n':int((pd.read_csv(S/'antismash_regions_674.csv').region_touches_contig_edge.astype(str)=='True').sum()),
 'protocluster_contig_edge_n':int((pd.read_csv(S/'antismash_all_protoclusters_878.csv').protocluster_touches_contig_edge.astype(str)=='True').sum()),
 'reference_loci_n':len(pd.read_csv(S/'comparative_assembly_loci_18.tsv',sep='\t')),
 'neighbor_tests_with_self':len(q1),
 'neighbor_hits_with_self':int((q1.Status=='HIT').sum()),
 'nonself_possible_pairs':len(non),
 'nonself_hit_pairs':int((non.Status=='HIT').sum()),
 'nonself_outside_pairs':int((non.Status=='OUTSIDE').sum()),
 'nonself_high_coverage_pairs':int(((non.Status=='HIT')&(non.QueryHSPCoveragePct>=70)&(non.SubjectHSPCoveragePct>=70)).sum()),
 'biG_SCAPE_executed':False,
 'individual_geographical_metadata_verified':False,
 'metabolite_identity_experimentally_validated':False
}
assert checks['genome_n']==20 and checks['distinct_biosample_n']==20
assert checks['region_n']==674 and checks['protocluster_n']==878 and checks['candidate_cluster_objects']==911
assert sum(checks['category_totals'].values())==878
assert checks['nonself_possible_pairs']==208 and checks['nonself_hit_pairs']==202 and checks['nonself_outside_pairs']==6 and checks['nonself_high_coverage_pairs']==202
(D/'independent_numeric_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
print(json.dumps(checks,indent=2))
