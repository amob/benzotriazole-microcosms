#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --job-name=trimandfilter
#SBATCH --nodes=2
#SBATCH --ntasks=40
#SBATCH --output=outandlog/%x_%j.txt
#SBATCH --error=errors/%x_%j.err
#SBATCH --partition=debug

##GOAL: summarize demultiplexed reads in the BZS project, trim adapters, join pairs, quality filter and re-summarize
#quality filtering should be run before deblur, which is in the next step.

qiime demux summarize  --i-data /outputs/Lmb_BZ17_demux.qza \
   --o-visualization /outputs/Lmb_BZ17_demux_view.qzv
#outputs file Lmb_B17_demux.qzv, a qiime readable object, all .qzv can be visualizd on the https://view.qiime2.org website
#.qza objects can be uploaded to the site as well to view information on provenance.

qiime cutadapt trim-paired --i-demultiplexed-sequences /outputs/Lmb_BZ17_demux.qza \
   --p-cores 40 --p-front-f CCTACGGGNGGCWGCAG --p-front-r GACTACHVGGGTATCTAATCC \
   --o-trimmed-sequences /outputs/Lmb_BZ17_trim.qza 
#this is common region of the primer used -- the variable region comes "before" (in the reading direction) these start

#merging paired ends first is required, quality filtering is reccommended by deblur authors
##merge
qiime vsearch join-pairs --i-demultiplexed-seqs /outputs/Lmb_BZ17_trim.qza \
          --p-minovlen 15 --p-maxdiffs 10 --o-joined-sequences /outputs/Lmb_BZ17_join.qza
##quality filter
qiime quality-filter q-score  --i-demux /outputs/Lmb_BZ17_join.qza \
        --o-filtered-sequences /outputs/Lmb_BZ17_filter.qza \
        --o-filter-stats /outputs/Lmb_BZ17_filter-stats.qza
##demultiplex summmarize to see what truncation should be
qiime demux summarize --i-data /outputs/Lmb_BZ17_filter.qza \
           --o-visualization /outputs/Lmb_BZ17_filter.qzv
###based on this visualization, truncation was set at 402 (see summary table, can retain 98% of sequences with this length)
##read quality is good for longer, but deblur and taxonomy assignment require a fixed length of sequence across all reads
###TRUNCATION HAPPENS DURING DEBLUR, this number is determined here, but used in the next script.


