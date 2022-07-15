
 
rawseqstat <- read.csv("BZTinocmicrobiome_sequencing_statistics.csv")
sum(rawseqstat$Number.of.Reads) # 1,541,469
sum(rawseqstat$Number.of.Bases) # 770,734,500
range(rawseqstat$Number.of.Reads) #69,234 117,565
mean(rawseqstat$Number.of.Reads) #90,674.65
range(rawseqstat$Number.of.Bases) #34617000 58782500

#stats after joining, trimming and quality filtering can be viewed on https://view.qiime2.org from file Lmb_BZ17_filter.qzv
# similar stats after 'deblur' can be viewed at the same site with file Lmb_BZ17_deblurtable.qzv
