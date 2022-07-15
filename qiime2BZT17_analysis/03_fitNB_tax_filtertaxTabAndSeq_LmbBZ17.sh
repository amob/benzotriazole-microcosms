#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --job-name=taxonomy
#SBATCH --nodes=2
#SBATCH --ntasks=40
#SBATCH --output=outandlog/%x_%j.txt
#SBATCH --mail-type=END
#SBATCH --partition=debug

##the whole script should take only a few min to run


####GET CLASSIFIER
##import greengenes data
qiime tools import --type 'FeatureData[Sequence]' \
        --input-path /inputs/gg_13_8_otus/rep_set/99_otus.fasta --output-path /outputs/otus99.qza
## the 99_otu_taxonomy file from greengenes was modified to have a header row 'Feature ID', 'Taxon', tab separated, as this was required by qiime2
# input greengenes are from https://docs.qiime2.org/2021.4/data-resources/  downloaded november 2020

# the zipped greengenes refence database include both of the inputs used in the above and below command:  99_otus.fasta &  99_otu_taxonomy.txt 
qiime tools import --type 'FeatureData[Taxonomy]' \
	--input-path /inputs/gg_13_8_otus/taxonomy/99_otu_taxonomy.txt --output-path /outputs/ref-taxonomy99.qza

###trim it based on our reads, either with our without the truncation performed in deblur
qiime feature-classifier extract-reads --i-sequences /outputs/otus99.qza  \
	--p-f-primer  CCTACGGGNGGCWGCAG --p-r-primer  GACTACHVGGGTATCTAATCC \
	--p-trunc-len 402  --p-min-length 402 --o-reads /outputs/ref-seqs99_402.qza

##fit the classifiers
qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads /outputs/ref-seqs99_402.qza \
	--i-reference-taxonomy /outputs/ref-taxonomy99.qza --o-classifier /outputs/classifier_402.qza

####GET TAXONOMY CALLS IN DATASET
##identify the taxonomy in the dataset
qiime feature-classifier classify-sklearn --i-classifier /outputs/classifier_402.qza \
	--i-reads /outputs/Lmb_BZ17_deblurrep-seqs.qza --p-confidence 0.7 \
	--o-classification /outputs/Lmb_BZ17_taxonomy.qza

##get taxonomy visual
qiime metadata tabulate --m-input-file /outputs/Lmb_BZ17_taxonomy.qza \
	--o-visualization /outputs/Lmb_BZ17_taxonomy.qzv

##get barplots before filtering
qiime taxa barplot --i-table /outputs/Lmb_BZ17_deblurtable.qza \
	--i-taxonomy /outputs/Lmb_BZ17_taxonomy.qza --m-metadata-file /inputs/Lemna_microbiome_barebones_meta_inbzt2_wkmcC_only.tsv \
	--o-visualization /outputs/Lmb_BZ17_taxabarplot_before_filteraxa.qzv

##remove land plant chloroplast and mitochondria -- some mitochondria may be insect or lab technician contamination; most mitochondria will likely be duckweed, some small fraction may be interesting in another analysis
qiime taxa filter-seqs --i-sequences /outputs/Lmb_BZ17_deblurrep-seqs.qza --i-taxonomy /outputs/Lmb_BZ17_taxonomy.qza \
	--p-exclude Streptophyta,mitochondria --o-filtered-sequences /outputs/Lmb_BZ17_filteredtaxa-seqs.qza

qiime taxa filter-table --i-table /outputs/Lmb_BZ17_deblurtable.qza --i-taxonomy /outputs/Lmb_BZ17_taxonomy.qza \
	--p-exclude Streptophyta,mitochrondria --o-filtered-table /outputs/Lmb_BZ17_filteredtaxa-table.qza

##get barplots after filtering out chloroplast and mitochondria
qiime taxa barplot --i-table /outputs/Lmb_BZ17_filteredtaxa-table.qza --i-taxonomy /outputs/Lmb_BZ17_taxonomy.qza \
	--m-metadata-file /inputs//Lemna_microbiome_barebones_meta_inbzt2_wkmcC_only.tsv \
	--o-visualization /outputs/Lmb_BZ17_filteredtaxa-barplot.qzv

