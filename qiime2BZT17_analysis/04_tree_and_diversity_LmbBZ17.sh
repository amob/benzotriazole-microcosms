#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --job-name=treeandiv
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --output=outandlog/%x_%j.txt
#SBATCH --error=errors/%x_%j.err


##GOAL of this script is to build a phylogenetic tree of ASVs in the dataset and calculate standard diversity metrics within and between samples
## the first command uses sepp-refs-gg-13-8.qza, which is a database downloaded from https://docs.qiime2.org/2021.4/data-resources/ November, 2020 

###get phylogeny info; note this depends on the header to set threads, a user may need to alter this.
qiime fragment-insertion sepp --i-reference-database /inputs/sepp-refs-gg-13-8.qza \
     --i-representative-sequences /outputs/Lmb_BZ17_filteredtaxa-seqs.qza \
     --o-tree /outputs/Lmb_BZ17_sepptree.qza --o-placements /outputs/Lmb_BZ17_seppplace.qza --p-threads $SLURM_CPUS_PER_TASK
###clean features from table that can't be placed in phylogeny
qiime fragment-insertion filter-features --i-table /outputs/Lmb_BZ17_filteredtaxa-table.qza \
	--i-tree /outputs/Lmb_BZ17_sepptree.qza --o-filtered-table /outputs/Lmb_BZ17_sepp_keptfeatures_table.qza \
	--o-removed-table /outputs/Lmb_BZ17_sepp_rmdfeatures_table.qza --verbose

#calculate basic set of diversity metrics
qiime diversity core-metrics-phylogenetic --i-table /outputs/Lmb_BZ17_sepp_keptfeatures_table.qza \
	--i-phylogeny /outputs/Lmb_BZ17_sepptree.qza --p-sampling-depth 1000 \
	--m-metadata-file /inputs/Lemna_microbiome_barebones_meta_inbzt2_wkmcC_only.tsv \
	--output-dir /outputs/Lmb_BZ17_coremetrics
#perform alpha rarefaction
qiime diversity alpha-rarefaction  --i-table /outputs/Lmb_BZ17_sepp_keptfeatures_table.qza \
	--i-phylogeny /outputs/Lmb_BZ17_sepptree.qza  --p-max-depth 5000 \
	--m-metadata-file /inputs/Lemna_microbiome_barebones_meta_inbzt2_wkmcC_only.tsv \
	--o-visualization /outputs/Lmb_BZ17_alpha-rarefaction.qzv
