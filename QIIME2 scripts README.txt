README

This readme covers the sequence file processing from duckweed cultured microbial communities.

Use care with calling of qiime, directories, script headers. This will likely have to be edited.
Not all files are available, please see figshare, 10.6084/m9.figshare.20311758. 
Raw input files to sequence analysis and selected intermediate files are provided at that figshare DOI.
gg_13_8_otus.tar.gz in the inputs folder must be unzipped before running script 03_fitNB_tax_filtertaxTabAndSeq_LmbBZ17.sh


00_ReadIn_LmbBZ17.sh
input files
	MANIFEST.csv
	fasta files listed in MANIFEST.csv
output
	Lmb_BZ17_demux.qza
	
01_trim_join_filter_Lmb_BZ17.sh
input
	Lmb_BZ17_demux.qza
output
	Lmb_BZ17_demux_view.qzv 
	Lmb_BZ17_trim.qza
	Lmb_BZ17_join.qza
	Lmb_BZ17_filter.qza
	Lmb_BZ17_filter-stats.qza
	Lmb_BZ17_filter.qzv

02_deblur_table_repseqs_stats_LmbBZ17.sh
input
	Lmb_BZ17_filter.qza
output
	Lmb_BZ17_deblurtable.qza
	Lmb_BZ17_deblurrep-seqs.qza
	Lmb_BZ17_deblur-stats.qza
	Lmb_BZ17_deblurtable.qzv
	Lmb_BZ17_deblurrep-seqs.qzv
	Lmb_BZ17_deblur-stats.qzv

03_fitNB_tax_filtertaxTabAndSeq_LmbBZ17.sh
input
	99_otus.fasta
	99_otu_taxonomy.txt
	Lmb_BZ17_deblurrep-seqs.qza
	Lmb_BZ17_deblurtable.qza
	Lemna_microbiome_barebones_meta_inbzt2_wkmcC_only.tsv
output
	otus99.qza
	ref-taxonomy99.qza
	ref-seqs99_402.qza
	classifier_402.qza
	Lmb_BZ17_taxonomy.qza
	Lmb_BZ17_taxonomy.qzv
	Lmb_BZ17_taxabarplot_before_filteraxa.qzv
	Lmb_BZ17_filteredtaxa-seqs.qza
	Lmb_BZ17_filteredtaxa-table.qza
	Lmb_BZ17_filteredtaxa-barplot.qzv

04_tree_and_diversity_LmbBZ17.sh
input
	sepp-refs-gg-13-8.qza
	Lmb_BZ17_filteredtaxa-seqs.qza
	Lmb_BZ17_filteredtaxa-table.qza
	Lemna_microbiome_barebones_meta_inbzt2_wkmcC_only.tsv
output
	Lmb_BZ17_sepptree.qza
	Lmb_BZ17_seppplace.qza
	Lmb_BZ17_sepp_keptfeatures_table.qza
	Lmb_BZ17_sepp_rmdfeatures_table.qza
	Lmb_BZ17_coremetrics/ [entire directory]
	Lmb_BZ17_alpha-rarefaction.qzv
	
05_balances_LmbBZ17.sh
input
	Lmb_BZ17_sepp_keptfeatures_table.qza
	Lmb_BZ17_sepptree.qza
	Lemna_microbiome_barebones_meta_inbzt2_wkmcC_only.tsv
output	
	Lmb_BZ17_balances_tree.qza
	Lmb_BZ17_hierarchy_tree.qza
	Lmb_BZ17_balances_heatmap.qzv
	
	