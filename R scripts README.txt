
This readme file documents the inputs and outputs of each R analysis script
Analyses include datafiles from experimental observations, as well as output files from sequence data processing
Not files are available, please see 10.6084/m9.figshare.20311758 for raw input files to sequence analysis and selected intermediate files.

bzt2growth.R
inputs
	locations_info_bzs2.csv
	AO BZS.2.ODall.csv
	BZS2 Frond Clara.csv
	end_map_errorchecked.csv
	end_data_errorchecked.csv
	Feb18-19_map_errorchecked.csv
	Feb18-19_data_errorchecked.csv
	Feb20-21map_errorchecked.csv
	Feb20-21dat.csv
	start map ALL.csv
	start data ALL.csv
	wellsizeBZS2photos.csv
outputs
	biodat.csv
	BZT2_no_inoc_effect.pdf *** FIGURE A.1
	dpixandOD_sigcats_rddens.pdf *** FIGURE 2
	frondsperday_nogenoRd.pdf *** FIGURE A.4
	frondsperday_nogenoRd_append.pdf *** FIGURE A.3
	
bzt3_stats_NLS.R
inputs
	bzt3_concs_original_full.xlsx
outputs
	kinetics_meansSEs_preds.pdf *** FIGURE A.2

sequencingstats.R
inputs
	BZTinocmicrobiome_sequencing_statistics.csv
	
prevExp products.R
inputs
	experiment A products.csv
	BZTinBZS1_trts_within_geno_plusgrowth.csv
	
microbiome bzt2.R
inputs
	Lmb_BZ17_sepp_keptfeatures_table.qza
	Lmb_BZ17_sepp_rmdfeatures_table.qza
	Lmb_BZ17_taxonomy.qza
	Lmb_BZ17_balances_tree.qza
	Lmb_BZ17_hierarchy_tree.qza
	Lmb_BZ17_coremetrics/shannon_vector.qza
	Lmb_BZ17_coremetrics/observed_features_vector.qza
	Lmb_BZ17_coremetrics/faith_pd_vector.qza
	Lemna_microbiome_barebones_meta_inbzt2_wkmcC_only.tsv 
	BZS2_transformation.csv
	biodat.csv	
outputs
	family relative abund.pdf *** FIGURE A.6
	products_plot_grid_nmperL.pdf *** used in FIGURE 1
	check_nodes.pdf  *** for code object inspection only
	check_nodes_genera.pdf *** for code object inspection only
	check_nodes_sp.pdf *** for code object inspection only
	ggtreeheat_interval_shortphylo_TPS.pdf *** used in FIGURE 3
	ggtreeheat_interval_shortphylo_bio.pdf *** FIGURE A.8
	check_cor_nodes1.pdf *** or others, see script text; for code object inspection only
	exploratorymicrobegroups_short.pdf *** FIGURE A.7
	diversity_bzt.pdf *** used in FIGURE 3
	diversity_growth.pdf *** FIGURE A.5
	
	