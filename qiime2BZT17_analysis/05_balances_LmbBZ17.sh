#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --job-name=balances
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --output=outandlog/%x_%j.txt
#SBATCH --error=errors/%x_%j.err
#SBATCH --partition=debug

#GOAL: compute balances of relative abundance along the phylogeny for BZ17 dataset

qiime gneiss ilr-phylogenetic --i-table /outputs/Lmb_BZ17_sepp_keptfeatures_table.qza \
	--i-tree /outputs/Lmb_BZ17_sepptree.qza --o-balances /outputs/Lmb_BZ17_balances_tree.qza \
	--o-hierarchy /outputs/Lmb_BZ17_hierarchy_tree.qza  --verbose

#visualization of the balances
qiime gneiss dendrogram-heatmap  --i-table /outputs/Lmb_BZ17_sepp_keptfeatures_table.qza \
	--i-tree /outputs/Lmb_BZ17_hierarchy_tree.qza  \
	--m-metadata-file /inputs/Lemna_microbiome_barebones_meta_inbzt2_wkmcC_only.tsv  \
	--m-metadata-column fieldormaster  --p-color-map seismic  --o-visualization /outputs/Lmb_BZ17_balances_heatmap.qzv


