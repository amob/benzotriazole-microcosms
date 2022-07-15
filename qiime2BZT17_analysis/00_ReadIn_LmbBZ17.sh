#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --job-name=importandcut
#SBATCH --nodes=2
#SBATCH --ntasks=40
#SBATCH --output=outandlog/%x_%j.txt
#SBATCH --error=errors/%x_%j.err
#SBATCH --partition=debug


###FOR THIS AND ALL SCRIPTS
###Depending on your system, adjust header (above is for a common cluster queue system, you may want to delete), adjust call to qiime, adjust directories and paths.

##how to run qiime2 help; commented out
#qiime --help

##import data using sequencing files and a manifest file
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
 --input-path /inputs/Lemna_microbiome_raw_files/MANIFEST.csv \
 --input-format PairedEndFastqManifestPhred33 \
 --output-path /outputs/Lmb_BZ17_demux.qza

