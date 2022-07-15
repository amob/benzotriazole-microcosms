#!/bin/bash
#SBATCH --time=05:00:00
#SBATCH --job-name=deblur
#SBATCH --nodes=2
#SBATCH --ntasks=40
#SBATCH --output=outandlog/%x_%j.txt
#SBATCH --error=errors/%x_%j.err

#NOTE: deblur takes more time for more samples. For this dataset, on the [cluster] machine it ran on, only 1 hour, maybe less.

#GOAL: run deblur on trimmed, merged, and quality filtered sequences

#if re-running analysis, double check results of previous script in qiime2view  before setting truncate length in the below

#run deblur
qiime deblur denoise-16S --i-demultiplexed-seqs /outputs/Lmb_BZ17_filter.qza \
        --p-trim-length 402 --p-sample-stats --o-table /outputs/Lmb_BZ17_deblurtable.qza \
        --o-representative-sequences /outputs/Lmb_BZ17_deblurrep-seqs.qza \
        --o-stats /outputs/Lmb_BZ17_deblur-stats.qza

#summarize and visualize each output file 
qiime feature-table summarize --i-table /outputs/Lmb_BZ17_deblurtable.qza \
        --o-visualization /outputs/Lmb_BZ17_deblurtable.qzv

qiime feature-table tabulate-seqs --i-data /outputs/Lmb_BZ17_deblurrep-seqs.qza \
        --o-visualization /outputs/Lmb_BZ17_deblurrep-seqs.qzv

qiime deblur visualize-stats --i-deblur-stats /outputs/Lmb_BZ17_deblur-stats.qza \
        --o-visualization /outputs/Lmb_BZ17_deblur-stats.qzv

