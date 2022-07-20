# library(tidyverse)
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(SDMTools) #this is just for legend.gradient used in family abundance figure; may be difficult to install, if so simply comment out the library and legend.gradient()
# install.packages("remotes")
# remotes::install_version("SDMTools", "1.1-221")
library(MCMCglmm)
library(phangorn)#descendants 
library(ggtree)
library(ggnewscale) #for ggtree

#generic, useful functions

range01 <- function(x) {
	newnums=(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))
	return(newnums)
} 

HPDi <- function(vect,prob) {
	int <- HPDinterval(as.mcmc(vect),prob=prob)
	return(int)
} #shortcut for coda HPDinterval 

bufferX <- function(x,p) { 
	r<- range(x,na.rm=T)
	add <- c(-1,1)*p*(r[2]-r[1])
	return(r+add)
	}	

std.error <- function(dat, na.rm=TRUE) {sd(dat,na.rm=na.rm)/sqrt(length(dat))}#defaults to na.rm=T


################################################
####READ IN MICROBIOME DATA, BIOOLOGICAL DATA, COLLECTION INFO, PROCESS
################################################

feat.tab.qz 	<- read_qza("../qiime2BZT17_sequencefiles_and_analysis/outputs/Lmb_BZ17_sepp_keptfeatures_table.qza")#filtered taxa in the phylogeny (not sure how much this last filter tosses)	
feat.rmtab.qz 	<- read_qza("../qiime2BZT17_sequencefiles_and_analysis/outputs/Lmb_BZ17_sepp_rmdfeatures_table.qza")#
feat.tax.qz 	<- read_qza("../qiime2BZT17_sequencefiles_and_analysis/outputs/Lmb_BZ17_taxonomy.qza")##NOT FILTERED TO kept features, or even to rm streptophyta (chloroplast, which are removed from the dataset, along with mitochondrial sequences)
bal.phy.qz	<- read_qza("../qiime2BZT17_sequencefiles_and_analysis/outputs/Lmb_BZ17_balances_tree.qza")$data
hier.phy 	<- read_qza("../qiime2BZT17_sequencefiles_and_analysis/outputs/Lmb_BZ17_hierarchy_tree.qza")$data 
bal.phy 	<- bal.phy.qz[,order(colnames(bal.phy.qz))] 
feat.tab 	<- feat.tab.qz$data[,order(colnames(feat.tab.qz$data))] 
feat.rmtab 	<- feat.rmtab.qz$data[,order(colnames(feat.rmtab.qz$data))] 
dim(feat.rmtab) # empty. no taxa had to be removed when making phylogeny
feat.tax 	<- feat.tax.qz$data

#reporting stats
nrow(feat.rmtab) #no reads or ASVs removed during placement into tree
mean(colSums((feat.tab))) #[1] 30632.06
range(colSums((feat.tab))) #[1] 19544 47017
sum(colSums((feat.tab))) #520745
range(colSums(sign(feat.tab)))# 13 44
mean(colSums(sign(feat.tab))) # 25.23529
range(rowSums(sign(feat.tab)))
table(rowSums(sign(feat.tab)))
feat.tax[feat.tax$Feature.ID%in%names(which(rowSums(sign(feat.tab))>=5)),]


#alpha diversity metrics based on rarefied sampling in QIIME2; 
div.shannon.qz 	<- read_qza("../qiime2BZT17_sequencefiles_and_analysis/outputs/Lmb_BZ17_coremetrics/shannon_vector.qza")$data 
transfdivS <- div.shannon.qz[order(rownames(div.shannon.qz)),]
div.feats.qz 	<- read_qza("../qiime2BZT17_sequencefiles_and_analysis/outputs/Lmb_BZ17_coremetrics/observed_features_vector.qza")$data 
transfdivF <- div.feats.qz[order(rownames(div.feats.qz)),]
div.faith.qz 	<- read_qza("../qiime2BZT17_sequencefiles_and_analysis/outputs/Lmb_BZ17_coremetrics/faith_pd_vector.qza")$data 
transfdivPD <- div.faith.qz[order(rownames(div.faith.qz)),]
#reordering results:
rownames(div.shannon.qz)[order(rownames(div.shannon.qz))] == colnames(feat.tab)
colnames(bal.phy) == colnames(feat.tab)
#rarefied richness in # features or faith's pd is tightly correlated to raw richness: colSums(sign(feat.tab)))


#make proportion:
feat.prop <- feat.tab
ftsums <- colSums(feat.tab)
for(i in 1:ncol(feat.tab)){feat.prop[,i] <- feat.prop[,i]/ftsums[i]}
#split tax info into columns
taxlevs <- sapply(1:nrow(feat.tax), function(z) strsplit(as.character(feat.tax$Taxon[z]),split="; "))
tax <- matrix(NA,ncol=7,nrow=nrow(feat.tax))
for(i in 1:nrow(tax)){ tax[i,1:length(taxlevs[[i]])]<- taxlevs[[i]]}
rownames(tax) <- as.character(feat.tax$Feature.ID)
mastertab.s <- feat.prop
mastertax.s <- tax[rownames(mastertab.s),]


##read sample info, and sort accordingly
sampleinfo 		<- read.csv("Lemna_microbiome_barebones_meta_inbzt2_wkmcC_only.tsv",sep="\t",header=T, stringsAsFactors=F) 
newcolnamesM <- sapply(colnames(mastertab.s), function(x) sampleinfo$bzs2[sampleinfo$sample.id==x])
colnames(mastertab.s)<- newcolnamesM
colnames(bal.phy) <- newcolnamesM
famsr <- unique(mastertax.s[,5])
fams <- famsr[-c(which(famsr=="f__"),which(is.na(famsr)))]
famdf <- matrix(NA,ncol=length(fams)+1,nrow=17 )
for(i in 1:length(fams)){
	if(length(which(mastertax.s[,5]==fams[i]))>1) {
	famdf[,i] <- colSums(mastertab.s[which(mastertax.s[,5]==fams[i]),])
	} else { famdf[,i] <- mastertab.s[which(mastertax.s[,5]==fams[i]),] }
}
famdf[,length(fams)+1] <- colSums(mastertab.s[which(mastertax.s[,5]=="f__" | is.na(mastertax.s[,5]) ),])
colnames(famdf) <- c(fams,"Unknown")
rownames(famdf) <- colnames(mastertab.s)
famsim <- famdf[,names(sort(colSums(famdf),decreasing=T))[1:10]]
famsim <-cbind(famsim, rowSums(famdf[,names(sort(colSums(famdf),decreasing=T))[11:(length(fams))]]))
#op ten family names
names(sort(colSums(famdf),decreasing=T))[1:10]
#write out names
famnames <- c( "Pseudomonadaceae" , "Aeromonadaceae" ,"Enterobacteriaceae" , "Rhizobiaceae"  ,
				 "Xanthomonadaceae"  ,  "Sphingobacteriaceae",  "Bacillaceae"  , 
				 "Comamonadaceae"  , "Moraxellaceae",  "Chromatiaceae",
				   "Other"   )          
wb <- colorRampPalette(c(rgb(1,1,1),rgb(0,0,0)))
pdf("family relative abund.pdf",height=3,width=7)
layout(matrix(1:2,ncol=2),widths=c(1.5,0.5))
par(mar=c(4,11,1,1))
image(famsim,col=wb(50),xaxt="n",yaxt="n",ylab="",xlab="",zlim=c(0,1))
	axis(side = 1,at = seq(from=0,to=1,length.out=17), labels=rownames(famdf),las=2)
	axis(side=2, at=seq(from=0,to=1,length.out=11), labels=famnames,las=2)
par(mar=c(1,2,1,1))
plot(c(1,10)~c(0,15), xaxt="n", yaxt="n",bty="n",pch=NA,ylab="",xlab="")
	legend.gradient(cbind(c(1,1.75,1.75,1),c(2,2,9,9)),cols=wb(15),title="",limits=c("0%","100%"),cex=1.2)
	mtext("Relative abundance",side=2) #srt=270
dev.off()


sum.aeromonad <- colSums(mastertab.s[which(mastertax.s[,5]=="f__Aeromonadaceae"),])
sum.pseudomonad <- colSums(mastertab.s[which(mastertax.s[,5]=="f__Pseudomonadaceae"),])
sum.enterobact <- colSums(mastertab.s[which(mastertax.s[,5]=="f__Enterobacteriaceae"),])
sum.rhizob <- colSums(mastertab.s[which(mastertax.s[,5]=="f__Rhizobiaceae"),])
sum.anythingelse <- 1 - sum.aeromonad - sum.pseudomonad - sum.enterobact - sum.rhizob
 
#biological measurements
biodat <- read.csv("biodat.csv",header=T) # deltasqmm is change in sq mm of duckweed as measured by pixel area in images
biodat$lnOD600 <- log(biodat$od_dat.od600+ min(biodat$od_dat.od600[biodat$od_dat.od600>0],na.rm=T) )

#bzt transf data
transf <- read.csv("BZS2_transformation.csv",stringsAsFactors=F,header=T)
1-range(transf$conc_corrected) #report range
transf$Genotype[transf$Genotype=="Cnp"] <- "CnP" #to match other files
table(transf$Genotype)

#molar mass benzotriazole: 119.12
#molar masses: bzta, glybzt, bztaa, anilin, a3p, phen, mbzt, moxybz
molmassTPs <- c(207,282,249 ,94,110,181, 134,150) #daltons 
mgperL <- sapply(1:8, function(z) molmassTPs[z]*transf[,z+15]*(1000/2.5))  #molarmass * mols per 2.5 mL * 1000/2.5 = grams /L
#ecosar units: mg/L bzt is 40.7
fishLC50 <- c(3.66e6,4.21e4,1.61e5,40.3,11.4,80.9,652,45.6)#mg/L
toxunits <- rowSums(sapply(1:8, function(z) mgperL[,z]/fishLC50[z]))
toxuBZT <- toxunits + 40.7*transf$conc_corrected #proportion remaining * 1, starting value was 1 mg/L
cor(transf$conc_corrected,toxuBZT)#rho is 1

#create identifiers to merge bio/transformation/microbiome information
arabicplate <- as.numeric(sapply(1:nrow(biodat), function(z) strsplit(biodat$plate[z],"[.]")[[1]][1]))
splitpair <- sapply(1:nrow(biodat), function(z) strsplit(biodat$plate[z],"[.]")[[1]][2]) 
arabicpair <- ifelse(splitpair=="I",1,2)
same_sampleid <- paste(arabicplate,arabicpair,paste(biodat$row,biodat$column,sep=""),sep=".")
# subset
subbio <- biodat[biodat$genotype%in%transf$Genotype,]
subbioexact <- t(sapply(1:length(transf$id), function(id) biodat[which(same_sampleid == transf$id[id]),]))
#sorted idenically to transf object
subbioexactdf <- data.frame(matrix(ncol=ncol(biodat),nrow=nrow(transf)))
for(i in 1:ncol(biodat)){
	subbioexactdf[,i] <- unlist(subbioexact[,i])
}
colnames(subbioexactdf) <- colnames(biodat)

#### Subset, sort, merge.
#now subset inocula abundance data to only the genotypes for which there are transformation products
submastertab <- mastertab.s[,colnames(mastertab.s)%in%transf$Genotype]
	submastertab.s <- submastertab[-(which(rowSums(submastertab) ==0)),]
	submastertax <- mastertax.s[rownames(submastertab.s),]
subbalphy <- bal.phy[,colnames(bal.phy)%in%transf$Genotype]
#get transformation products matrix by genotype, transform via log when appropriate
products <- transf[,c(11,16:23)]
productsAA <- cbind(transf[,c(11,17)], rowSums(transf[,c(16,18)]), transf[,19:23])
productsAAnm <- productsAA
productsAAnm[,2:8] <- (productsAA[,2:8]/2.5)*10e09*1000 # from moles per well to nanomoles per liter
productsAAnm[,1] <- 100*(1-productsAA[,1])
colnames(productsAAnm) <- c( colnames(transf)[c(11,17)],"BZTa+BZTaa" , colnames(transf)[19:23] )
logprodsAA <- sapply(1:ncol(productsAAnm), function(x) log(productsAAnm[,x] + min(productsAAnm[productsAA[,x]>0,x])))
colnames(logprodsAA) <- colnames(productsAAnm)
chooselogprodsAA <- sapply(1:ncol(productsAAnm), function(z) shapiro.test(productsAAnm[,z])$statistic < 0.7 & shapiro.test(logprodsAA[,z])$statistic > shapiro.test(productsAAnm[,z])$statistic )
mixAA <- productsAAnm
mixAA[,which(chooselogprodsAA)] <- logprodsAA[,chooselogprodsAA]
mixAAtrt <- cbind(Genotype=transf$Genotype, Salt=transf$Salt, Microbes = transf$Microbes, rddens =transf$density/1000, km= transf$km, mixAA)
colnames(mixAAtrt)[colnames(mixAAtrt)=="BZTa+BZTaa"] <- "BZTaplusBZTaa"
mixAAtrt$OD600 <- subbioexactdf$od_dat.od600
mixAAtrt$lnOD600 <- subbioexactdf$lnOD600
mixAAtrt$deltasqmm <- subbioexactdf$deltasqmm
shapiro.test(toxunits); shapiro.test(log(toxunits))
mixAAtrt$logtoxu <- log(toxunits)
mixAAtrtbio <- mixAAtrt[!is.na(mixAAtrt$OD600),]
mixAAtrt$micrnum <- as.numeric(as.factor(mixAAtrt$Microbes))-1 
#genotyp biomeans by trts
biogenotrt <- sapply(sort(unique(biodat$genotype)), function(g) sapply(c(9:16,18:28,68), function(p) tapply(biodat[biodat$genotype==g,p],paste(biodat$BZT,biodat$Salt,biodat$microbe)[biodat$genotype==g],mean,na.rm=T) )) 
#trtmnt order is all 0 bzt before all + bzt; then all 0 s before + s, then all - before +, so 0 0 -; 0 0 +, 0S - ....
#traits are grouped in rows first 8 rows are x0.FN, order is x0.FN through x10.FN, od600 od750, then all imageJ columns (area, perim, mean...) 


################################################
###models of product responses to treatments, correlation with biological growth
################################################

rddens.s <- seq(from=min(mixAAtrt$rddens),to=max(mixAAtrt$rddens),length.out=1000)
dsqm.s <- seq(from=min(mixAAtrtbio$deltasqmm),to=max(mixAAtrtbio$deltasqmm),length.out=1000)
od.s <- seq(from=min(mixAAtrtbio$lnOD600),to=max(mixAAtrtbio$lnOD600),length.out=1000)
summary(MCMCglmm(conc_corrected~Microbes + Salt +rddens + rddens:Salt + rddens:Microbes + Microbes:Salt+ rddens:Microbes:Salt,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(conc_corrected~Microbes + Salt +rddens + rddens:Salt + rddens:Microbes + Microbes:Salt,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(conc_corrected~Microbes + Salt +rddens + rddens:Salt  + Microbes:Salt,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(conc_corrected~Microbes + Salt +rddens   + Microbes:Salt,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(conc_corrected~Microbes + Salt +rddens   ,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(conc_corrected~Salt +rddens   ,data=mixAAtrt,nitt=50000,verbose=F))
#refit
BZTremtrt <- MCMCglmm(conc_corrected~ Salt +rddens   ,data=mixAAtrt,nitt=100000,verbose=F) #salt neg and rddens pos
solbrt <- BZTremtrt$Sol
predBRTm <- sapply(rddens.s, function(z) mean(solbrt[,1] + solbrt[,2]*0.5 + solbrt[,3]*z) )
predBRTi <- sapply(rddens.s, function(z) HPDi(solbrt[,1] + solbrt[,2]*0.5 + solbrt[,3]*z, 0.95) )
summary(MCMCglmm(conc_corrected~deltasqmm + lnOD600 + lnOD600:deltasqmm,data=mixAAtrtbio,nitt=50000,verbose=F))
summary(MCMCglmm(conc_corrected~deltasqmm + lnOD600 ,data=mixAAtrtbio,nitt=50000,verbose=F))
summary(MCMCglmm(conc_corrected~deltasqmm  ,data=mixAAtrtbio,nitt=50000,verbose=F))
#refit
BZTdp <- (MCMCglmm(conc_corrected~deltasqmm ,data=mixAAtrtbio,nitt=100000,verbose=F))#pos
solbpx <- BZTdp$Sol
predBpxm <- sapply(dsqm.s, function(z) mean(solbpx[,1] + solbpx[,2]*z) )
predBpxi <- sapply(dsqm.s, function(z) HPDi(solbpx[,1] + solbpx[,2]*z, 0.95))
#no fit issues
#
summary(MCMCglmm(BZTaplusBZTaa~Microbes + Salt +rddens + rddens:Salt + rddens:Microbes + Microbes:Salt+ rddens:Microbes:Salt,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(BZTaplusBZTaa~Microbes + Salt +rddens + rddens:Salt + rddens:Microbes + Microbes:Salt,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(BZTaplusBZTaa~Microbes + Salt +rddens + rddens:Salt  + Microbes:Salt,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(BZTaplusBZTaa~Microbes + Salt +rddens  + Microbes:Salt,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(BZTaplusBZTaa~Microbes + Salt +rddens  ,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(BZTaplusBZTaa~ Salt +rddens  ,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(BZTaplusBZTaa~  rddens  ,data=mixAAtrt,nitt=50000,verbose=F))
#refit
BZTaatrt <- (MCMCglmm(BZTaplusBZTaa~ rddens  ,data=mixAAtrt,nitt=100000,verbose=F)) #rddens pos
solbaat <- BZTaatrt$Sol
predBATm <- sapply(rddens.s, function(z) mean(solbaat[,1] + solbaat[,2]*z) )
predBATi <- sapply(rddens.s, function(z) HPDi(solbaat[,1] + solbaat[,2]*z, 0.95) )
summary(MCMCglmm(BZTaplusBZTaa~deltasqmm + lnOD600 + lnOD600:deltasqmm,data=mixAAtrtbio,nitt=50000,verbose=F))
summary(MCMCglmm(BZTaplusBZTaa~deltasqmm + lnOD600 ,data=mixAAtrtbio,nitt=50000,verbose=F))
summary(MCMCglmm(BZTaplusBZTaa~ lnOD600 ,data=mixAAtrtbio,nitt=50000,verbose=F))
baaod <-(MCMCglmm(BZTaplusBZTaa~ lnOD600 ,data=mixAAtrtbio,nitt=100000,verbose=F)) #pos
solbaaod <- baaod$Sol
predBaom <- sapply(od.s, function(z) mean(solbaaod[,1] + solbaaod[,2]*z) )
predBaoi <- sapply(od.s, function(z) HPDi(solbaaod[,1] + solbaaod[,2]*z, 0.95))
#no fit issues
#
summary(MCMCglmm(scale(glycosylatedBZT)~Microbes + Salt +rddens + rddens:Salt + rddens:Microbes + Microbes:Salt+ rddens:Microbes:Salt,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(scale(glycosylatedBZT)~Microbes + Salt +rddens + rddens:Salt + rddens:Microbes + Microbes:Salt,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(scale(glycosylatedBZT)~Microbes + Salt +rddens + rddens:Salt + Microbes:Salt,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(scale(glycosylatedBZT)~Microbes + Salt +rddens + Salt:Microbes,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(scale(glycosylatedBZT)~Microbes + Salt +rddens  ,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(scale(glycosylatedBZT)~Salt + rddens   ,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(scale(glycosylatedBZT)~Salt    ,data=mixAAtrt,nitt=50000,verbose=F))#gly bzt is intercept
summary(MCMCglmm(scale(glycosylatedBZT)~1    ,data=mixAAtrt,nitt=50000,verbose=F))#gly bzt is intercept
summary(MCMCglmm(scale(glycosylatedBZT)~deltasqmm + lnOD600 + lnOD600:deltasqmm,data=mixAAtrtbio,nitt=50000,verbose=F))
summary(MCMCglmm(scale(glycosylatedBZT)~deltasqmm + lnOD600 ,data=mixAAtrtbio,nitt=50000,verbose=F)) #both sig! dsqm pos, od neg
#refit
glypo <- (MCMCglmm(scale(glycosylatedBZT)~deltasqmm + lnOD600,data=mixAAtrtbio,nitt=100000,verbose=F))#pos
solgpo <- glypo$Sol
predGpxm <- sapply(dsqm.s, function(z) mean(solgpo[,1] + solgpo[,2]*z + solgpo[,3]*mean(mixAAtrtbio$lnOD600)) )*sd(mixAAtrtbio$glycosylatedBZT) + mean(mixAAtrtbio$glycosylatedBZT)
predGpxi <- sapply(dsqm.s, function(z) HPDi(solgpo[,1] + solgpo[,2]*z + solgpo[,3]*mean(mixAAtrtbio$lnOD600), 0.95))*sd(mixAAtrtbio$glycosylatedBZT) + mean(mixAAtrtbio$glycosylatedBZT)
predGodm <- sapply(od.s, function(z) mean(solgpo[,1] + solgpo[,2]*mean(mixAAtrtbio$deltasqmm) + solgpo[,3]*z) )*sd(mixAAtrtbio$glycosylatedBZT) + mean(mixAAtrtbio$glycosylatedBZT)
predGodi <- sapply(od.s, function(z) HPDi(solgpo[,1] + solgpo[,2]* mean(mixAAtrtbio$deltasqmm)+ solgpo[,3]*z, 0.95))*sd(mixAAtrtbio$glycosylatedBZT) + mean(mixAAtrtbio$glycosylatedBZT)
#fit issues if unscaled
#
summary(MCMCglmm(scale(methylBZT)~Microbes + Salt +rddens + rddens:Salt + rddens:Microbes + Microbes:Salt+ rddens:Microbes:Salt,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(scale(methylBZT)~Microbes + Salt +rddens + rddens:Salt + rddens:Microbes + Microbes:Salt,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(scale(methylBZT)~Microbes + Salt +rddens + rddens:Salt + rddens:Microbes ,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(scale(methylBZT)~Microbes + Salt +rddens +  rddens:Microbes,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(scale(methylBZT)~Microbes + Salt +rddens ,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(scale(methylBZT)~Microbes + rddens  ,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(scale(methylBZT)~rddens   ,data=mixAAtrt,nitt=50000,verbose=F))#int
summary(MCMCglmm(scale(methylBZT)~1   ,data=mixAAtrt,nitt=50000,verbose=F))#int
summary(MCMCglmm(scale(methylBZT)~deltasqmm + lnOD600 + lnOD600:deltasqmm,data=mixAAtrtbio,nitt=50000,verbose=F))
summary(MCMCglmm(scale(methylBZT)~deltasqmm + lnOD600 ,data=mixAAtrtbio,nitt=50000,verbose=F))
summary(MCMCglmm(scale(methylBZT)~lnOD600 ,data=mixAAtrtbio,nitt=50000,verbose=F))#sig neg
#refit
mbod <-(MCMCglmm(scale(methylBZT)~ lnOD600 ,data=mixAAtrtbio,nitt=100000,verbose=F)) #pos
solmbod <- mbod$Sol
predMBom <- sapply(od.s, function(z) mean(solmbod[,1] + solmbod[,2]*z) )*sd(mixAAtrtbio$methylBZT) + mean(mixAAtrtbio$methylBZT)
predMBoi <- sapply(od.s, function(z) HPDi(solmbod[,1] + solmbod[,2]*z, 0.95))*sd(mixAAtrtbio$methylBZT) + mean(mixAAtrtbio$methylBZT)
#zero effective samples if unscaled in most models for most effects
#
summary(MCMCglmm(methoxyBZT~Microbes + Salt +rddens + rddens:Salt + rddens:Microbes + Microbes:Salt+ rddens:Microbes:Salt,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(methoxyBZT~Microbes + Salt +rddens + rddens:Salt + rddens:Microbes + Microbes:Salt,data=mixAAtrt,nitt=50000,verbose=F))#worse, but continue to see if better fit
summary(MCMCglmm(methoxyBZT~Microbes + Salt +rddens + rddens:Salt  + Microbes:Salt,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(methoxyBZT~Microbes + Salt +rddens + rddens:Salt  ,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(methoxyBZT~Microbes + Salt +rddens  ,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(methoxyBZT~Salt +rddens  ,data=mixAAtrt,nitt=50000,verbose=F))
#refit
moxBtrt <- (MCMCglmm(methoxyBZT~ Salt +rddens  ,data=mixAAtrt,nitt=100000,verbose=F))# salt neg and rddens neg 
solmBt <- moxBtrt$Sol
predMBTm <- sapply(rddens.s, function(z) mean(solmBt[,1] + solmBt[,2]*0.5 + solmBt[,3]*z) ) 
predMBTi <- sapply(rddens.s, function(z) HPDi(solmBt[,1] + solmBt[,2]*0.5 + solmBt[,3]*z, 0.95) )
summary(MCMCglmm(methoxyBZT~deltasqmm + lnOD600 + lnOD600:deltasqmm,data=mixAAtrtbio,nitt=50000,verbose=F))
summary(MCMCglmm(methoxyBZT~deltasqmm + lnOD600 ,data=mixAAtrtbio,nitt=50000,verbose=F))
summary(MCMCglmm(methoxyBZT~deltasqmm ,data=mixAAtrtbio,nitt=50000,verbose=F))#ns, int
summary(MCMCglmm(methoxyBZT~1 ,data=mixAAtrtbio,nitt=50000,verbose=F))#ns, int
#no fit issues
#
summary(MCMCglmm(amino_3_phenol~Microbes + Salt +rddens + rddens:Salt + rddens:Microbes + Microbes:Salt+ rddens:Microbes:Salt,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(amino_3_phenol~Microbes + Salt +rddens + rddens:Salt + rddens:Microbes + Microbes:Salt,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(amino_3_phenol~Microbes + Salt +rddens + rddens:Salt  + Microbes:Salt,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(amino_3_phenol~Microbes + Salt +rddens   + Microbes:Salt,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(amino_3_phenol~Microbes + Salt +rddens  ,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(amino_3_phenol~Salt +rddens  ,data=mixAAtrt,nitt=50000,verbose=F))#DIC not different from prev.
#refit
a3ptrt <- (MCMCglmm(amino_3_phenol~  Salt +rddens  ,data=mixAAtrt,nitt=100000,verbose=F)) #salt neg and rddens neg
sola3pt <- a3ptrt$Sol
predAPTm <- sapply(rddens.s, function(z) mean(sola3pt[,1] + sola3pt[,2]*0.5 + sola3pt[,3]*z) )
predAPTi <- sapply(rddens.s, function(z) HPDi(sola3pt[,1] + sola3pt[,2]*0.5 + sola3pt[,3]*z, 0.95) )
summary(MCMCglmm(amino_3_phenol~deltasqmm + lnOD600 + lnOD600:deltasqmm,data=mixAAtrtbio,nitt=50000,verbose=F))
summary(MCMCglmm(amino_3_phenol~deltasqmm + lnOD600 ,data=mixAAtrtbio,nitt=50000,verbose=F))
summary(MCMCglmm(amino_3_phenol~deltasqmm ,data=mixAAtrtbio,nitt=50000,verbose=F)) #n.s.
summary(MCMCglmm(amino_3_phenol~1 ,data=mixAAtrtbio,nitt=50000,verbose=F)) #
#no fit issues
#
summary(MCMCglmm(phenazine~Microbes + Salt +rddens + rddens:Salt + rddens:Microbes + Microbes:Salt+ rddens:Microbes:Salt,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(phenazine~Microbes + Salt +rddens + rddens:Salt + rddens:Microbes + Microbes:Salt,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(phenazine~Microbes + Salt +rddens  + rddens:Microbes + Microbes:Salt,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(phenazine~Microbes + Salt +rddens  + rddens:Microbes ,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(phenazine~Microbes + Salt +rddens   ,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(phenazine~ Salt +rddens   ,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(phenazine~ rddens   ,data=mixAAtrt,nitt=50000,verbose=F))
#refit
pztrt <- (MCMCglmm(phenazine ~  rddens ,data=mixAAtrt,nitt=100000,verbose=F)) #rddens NEG
solpzt <- pztrt$Sol
predPZTm <- sapply(rddens.s, function(z) mean(solpzt[,1] + solpzt[,2]*z) )
predPZTi <- sapply(rddens.s, function(z) HPDi(solpzt[,1] + solpzt[,2]*z, 0.95) )
summary(MCMCglmm(phenazine~deltasqmm + lnOD600 + lnOD600:deltasqmm,data=mixAAtrtbio,nitt=50000,verbose=F))
summary(MCMCglmm(phenazine~deltasqmm + lnOD600 ,data=mixAAtrtbio,nitt=50000,verbose=F))
summary(MCMCglmm(phenazine~ lnOD600 ,data=mixAAtrtbio,nitt=50000,verbose=F))
#refit
pzod <- (MCMCglmm(phenazine~ lnOD600 ,data=mixAAtrtbio,nitt=100000,verbose=F))#
solpzod <- pzod$Sol
predpzodm <- sapply(od.s, function(z) mean(solpzod[,1] + solpzod[,2]*z) )
predpzodi <- sapply(od.s, function(z) HPDi(solpzod[,1] + solpzod[,2]*z, 0.95))
#no fit issues
#
summary(MCMCglmm(scale(aniline)~Microbes + Salt +rddens + rddens:Salt + rddens:Microbes + Microbes:Salt+ rddens:Microbes:Salt,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(scale(aniline)~Microbes + Salt +rddens + rddens:Salt + rddens:Microbes + Microbes:Salt,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(scale(aniline)~Microbes + Salt +rddens  + rddens:Microbes + Microbes:Salt,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(scale(aniline)~Microbes + Salt +rddens  + rddens:Microbes ,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(scale(aniline)~Microbes + Salt +rddens ,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(scale(aniline)~ Salt +Microbes ,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(scale(aniline)~ Salt ,data=mixAAtrt,nitt=50000,verbose=F)) # without scale, salt not sig.
antrt <- (MCMCglmm(scale(aniline) ~  Salt ,data=mixAAtrt,nitt=100000,verbose=F)) #rddens NEG
summary(MCMCglmm(scale(aniline)~deltasqmm + lnOD600 + lnOD600:deltasqmm,data=mixAAtrtbio,nitt=400000,verbose=F))#
summary(MCMCglmm(scale(aniline)~deltasqmm + lnOD600,data=mixAAtrtbio,nitt=400000,verbose=F))
summary(MCMCglmm(scale(aniline)~deltasqmm,data=mixAAtrtbio,nitt=400000,verbose=F))
summary(MCMCglmm(scale(aniline)~1,data=mixAAtrtbio,nitt=400000,verbose=F))


summary(MCMCglmm(logtoxu ~Microbes + Salt +rddens + rddens:Salt + rddens:Microbes + Microbes:Salt+ rddens:Microbes:Salt,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(logtoxu ~Microbes + Salt +rddens + rddens:Salt + rddens:Microbes + Microbes:Salt,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(logtoxu ~Microbes + Salt +rddens + rddens:Microbes + Microbes:Salt,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(logtoxu ~Microbes + Salt +rddens + Microbes:Salt,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(logtoxu ~Microbes + Salt +rddens,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(logtoxu ~Microbes +rddens,data=mixAAtrt,nitt=50000,verbose=F))
summary(MCMCglmm(logtoxu ~ rddens,data=mixAAtrt,nitt=50000,verbose=F)) #sig neg
toxutrt <- MCMCglmm(logtoxu ~ rddens,data=mixAAtrt,nitt=100000,verbose=F) #sig neg
summary(MCMCglmm(logtoxu~deltasqmm + lnOD600 + lnOD600:deltasqmm,data=mixAAtrtbio,nitt=50000,verbose=F))#theoretically, stop there
summary(MCMCglmm(logtoxu~deltasqmm + lnOD600 ,data=mixAAtrtbio,nitt=50000,verbose=F))#theoretically, stop there
summary(MCMCglmm(logtoxu~deltasqmm  ,data=mixAAtrtbio,nitt=50000,verbose=F))#theoretically, stop there
summary(MCMCglmm(logtoxu~1 ,data=mixAAtrtbio,nitt=50000,verbose=F))#theoretically, stop there
#no fit issues
cor(mixAAtrt$logtoxu,mixAAtrt[,7:13])
soltoxu <- toxutrt$Sol
predTXm <- sapply(rddens.s, function(z) mean(soltoxu[,1] + soltoxu[,2]*z) ) 
predTXi <- sapply(rddens.s, function(z) HPDi(soltoxu[,1] + soltoxu[,2]*z, 0.95) )
(exp(predTXm[1])-exp(predTXm[1000]))/exp(predTXm[1])


################################################
####cockerham and muir partitioning variance across environment into rank vs scale changes
################################################

mixAAtrt2 <- mixAAtrt
mixAAtrt2$glycosylatedBZT <- scale(mixAAtrt$glycosylatedBZT)
mixAAtrt2$methylBZT <- scale(mixAAtrt$methylBZT)
mixAAtrt2$aniline <- scale(mixAAtrt$aniline)
mixAAtrtsm <- mixAAtrt2[mixAAtrt$Salt==.8 & mixAAtrt$micrnum==1,]
mixAAtrtsn <- mixAAtrt2[mixAAtrt$Salt==.8 & mixAAtrt$micrnum<1,]
mixAAtrtnm <- mixAAtrt2[mixAAtrt$Salt<.8 & mixAAtrt$micrnum==1,]
mixAAtrtnn <- mixAAtrt2[mixAAtrt$Salt<.8 & mixAAtrt$micrnum<1,]
percentrankvar <- c()
for(p in 6:13){
	msm <- mixAAtrtsm
	msn <- mixAAtrtsn
	mnm <- mixAAtrtnm
	mnn <- mixAAtrtnn
	colnames(msm)[p] <- "resp"; colnames(msn)[p] <- "resp";	colnames(mnm)[p] <- "resp"; colnames(mnn)[p] <- "resp"
	vgs <- unlist(lapply(list(msm,msn,mnm,mnn),function(dat)   
		mean(MCMCglmm(resp~1,random=~Genotype,data=dat,verbose=F,nitt=1000000,thin=10,burnin=10000)$VCV[,"Genotype"])
	))
	gmns <- sapply(list(msm,msn,mnm,mnn), function(dat) tapply(dat$resp,dat$Genotype,mean))
	vic <- sum(unlist( sapply(4:2, function(j) sapply( (j-1):1, function(i) 
				2*sqrt(vgs[i])*sqrt(vgs[j])*(1-cor(gmns[,i],gmns[,j] ) )
			  )   ) ) )/(4*(4-1))
	vhv <- sum(unlist( sapply(4:2, function(j) sapply( (j-1):1, function(i) 
				(sqrt(vgs[i])-sqrt(vgs[j]) )^2
			  )   ) ) )/(4*(4-1))
	percentrankvar[p-5] <- (vic/(vic+vhv))
}
####with 1000000 iterations in mcmc, varies by less than +/-3 percent across runs
# percentrankvar in various runs printed here for comparison:
# 0.3012633 0.1693064 0.8095921 0.4458892 0.4403463 0.9017446 0.3662411 0.3844023
# 0.3070844 0.1708263 0.8008258 0.4193345 0.4675983 0.9227594 0.3701301 0.3874386
# 0.2936997 0.1639262 0.8156199 0.4359500 0.4424861 0.9336793 0.3701622 0.3782598
# 0.3082690 0.1701508 0.7972499 0.4363310 0.4001421 0.8882549 0.3610946 0.3822920


################################################
####FIGURE FOR PRODUCT TRANSFORMATION AND VARIANCE PARTITIONING
################################################

#these are all only BZT+, by definition. for microbe analyses, we therefore only care about + microbes
mixAAdens		<- sapply(sort(unique(transf$Genotype)), function(g) mean(transf$density[transf$Genotype==g],na.rm=T)/1000)
mixAAkmC	<- sapply(sort(unique(transf$Genotype)), function(g) mean(transf$km[transf$Genotype==g],na.rm=T))
NLmixAAsaltm <- sapply(sort(unique(transf$Genotype)), function(g) sapply(1:ncol(mixAA), function(p) mean(productsAAnm[transf$Genotype==g & transf$Salt=="0.8" &  transf$Microbes=="Yes",p],na.rm=T) ))
NLmixAAnosaltm <- sapply(sort(unique(transf$Genotype)), function(g) sapply(1:ncol(mixAA), function(p) mean(productsAAnm[transf$Genotype==g & transf$Salt=="0" &  transf$Microbes=="Yes",p],na.rm=T)))
NLmixAAsaltn <- sapply(sort(unique(transf$Genotype)), function(g) sapply(1:ncol(mixAA), function(p) mean(productsAAnm[transf$Genotype==g & transf$Salt=="0.8" &  transf$Microbes=="No",p],na.rm=T) ))
NLmixAAnosaltn <- sapply(sort(unique(transf$Genotype)), function(g) sapply(1:ncol(mixAA), function(p) mean(productsAAnm[transf$Genotype==g & transf$Salt=="0" &  transf$Microbes=="No",p],na.rm=T)))
productsbioAAnm <- productsAAnm[!is.na(mixAAtrt$OD600),]
productsbioAAnm[,2:8] <- productsAAnm[!is.na(mixAAtrt$OD600),2:8]
transfbio <- transf[!is.na(mixAAtrt$OD600),]
NLbioAAgensaltm <- sapply(sort(unique(transf$Genotype)), function(g) sapply(1:ncol(mixAA), function(p) mean(productsbioAAnm[ transfbio$Genotype==g & transfbio$Salt=="0.8" &  transfbio$Microbes=="Yes",p],na.rm=T)))
NLbioAAgennosaltm <- sapply(sort(unique(transf$Genotype)), function(g) sapply(1:ncol(mixAA), function(p) mean(productsbioAAnm[ transfbio$Genotype==g & transfbio$Salt=="0" &  transfbio$Microbes=="Yes",p],na.rm=T)))
NLbioAAgensaltn <- sapply(sort(unique(transf$Genotype)), function(g) sapply(1:ncol(mixAA), function(p) mean(productsbioAAnm[ transfbio$Genotype==g & transfbio$Salt=="0.8" &  transfbio$Microbes=="No",p],na.rm=T)))
NLbioAAgennosaltn <- sapply(sort(unique(transf$Genotype)), function(g) sapply(1:ncol(mixAA), function(p) mean(productsbioAAnm[ transfbio$Genotype==g & transfbio$Salt=="0" &  transfbio$Microbes=="Yes",p],na.rm=T)))
NLbioodgentrt <- sapply(sort(unique(transf$Genotype)), function(g)  sapply(c(0,0.8), function(salt)
						sapply(c("Yes","No"), function(micr) mean(mixAAtrtbio$lnOD600[ mixAAtrtbio$Genotype==g & mixAAtrtbio$Salt==salt &mixAAtrtbio$Microbes==micr],na.rm=T)) )) 
			#genotypes in columns, treatments rows: first two at zero salt, and very first with microbes
NLbiopxgentrt <- sapply(sort(unique(transf$Genotype)), function(g)  sapply(c(0,0.8), function(salt)
						sapply(c("Yes","No"), function(micr) mean(mixAAtrtbio$deltasqmm[ mixAAtrtbio$Genotype==g & mixAAtrtbio$Salt==salt &mixAAtrtbio$Microbes==micr],na.rm=T)) )) 
			#genotypes in columns, treatments rows: first two at zero salt, and very first with microbes

saltBmns <- sapply(1:8, function(z) tapply(productsAAnm[,z], mixAAtrt$Salt, mean))
saltBses <- sapply(1:8, function(z) tapply(productsAAnm[,z], mixAAtrt$Salt, std.error))

ylimsB <- list( c(78,100), c(0,50), c(0,250), 		c(0,220) ,c(0,100), c(0,35),  			c(0,100), c(0,35), c(0,500), c(0,500))
atv <- list(c(80,90,100), c(0,25,50), c(0,100,200), c(0,100,200), c(0,40,80), c(0,15,30), c(0,40,80),c(0,15,30), c(0,250,500), c(0,250,500)   )

#reaction norm plots
pdf("products_plot_grid_nmperL.pdf",height=6.25,width=5.25)
#reaction norm plots
layout(matrix(1:40,ncol=5,byrow=T) ,widths=c(.6,.8,.8,.8,.8),heights=c(1,.25,1,1,1,1,1,1))
par(oma = c(5.25,1,3,1))
for(i in c(1,11,2,3,7,8,5,6)){
	if(i==11){
	plot(1:2~c(1:2),pch=NA,bty="n",yaxt="n",xaxt="n",xlab="",ylab="")
	plot(1:2~c(1:2),pch=NA,bty="n",yaxt="n",xaxt="n",xlab="",ylab="")
	plot(1:2~c(1:2),pch=NA,bty="n",yaxt="n",xaxt="n",xlab="",ylab="")
	plot(1:2~c(1:2),pch=NA,bty="n",yaxt="n",xaxt="n",xlab="",ylab="")
	plot(1:2~c(1:2),pch=NA,bty="n",yaxt="n",xaxt="n",xlab="",ylab="")
	} else{
	par(mar=c(0,2,0,0))
	plot(saltBmns[,i]~c(1:2),col=c(rgb(0.4,0,0),rgb(1,0,0)), xaxt="n",yaxt="n",pch=1,cex=1.5,xlim=c(0.75,2.25),
			ylim=ylimsB[[i]] )
		arrows(1:2, y0=saltBmns[,i] - saltBses[,i], y1=saltBmns[,i] + saltBses[,i],length=0,col=c(rgb(0.4,0,0),rgb(1,0,0)))
		if(i==6){axis(side=1,at=c(1,2),las=2,labels=c("0","0.8"))}
		if(i==6){mtext(side=1,"salt g/L",line=3,cex=0.85)}
		axis(side=2,at=atv[[i]],las=2)
		#bztr, methoxyBZT, and amino3phenol have sig salt effs
		if(i==1){text(1.5,ylimsB[[i]][2]*0.95,"*",cex=2)}
		if(i%in%c(5,8)){text(1.5,ylimsB[[i]][2]*0.75,"*",cex=2)}
	par(mar=c(0,0,0,0))
	plot(NLmixAAnosaltn[i,]~rep(1,times=ncol(NLmixAAnosaltn)), xlim=c(0.5,4.5),
		ylim=ylimsB[[i]], xlab="",xaxt="n",yaxt="n",ylab="" ,
		col=rgb(0,(1-range01(mixAAdens))*0.6,0,alpha=0.5),pch=16)
		points(NLmixAAsaltn[i,]~rep(2,times=ncol(NLmixAAnosaltn)), col=rgb(0,(1-range01(mixAAdens))*0.6,0,alpha=0.5),pch=16)
		points(NLmixAAnosaltm[i,]~rep(3,times=ncol(NLmixAAsaltn)), col=rgb(0,(1-range01(mixAAdens))*0.6,0,alpha=0.5),pch=16)
		points(NLmixAAsaltm[i,]~rep(4,times=ncol(NLmixAAsaltn)), col=rgb(0,(1-range01(mixAAdens))*0.6,0,alpha=0.5),pch=16)
		sapply(1:ncol(NLmixAAnosaltn), function(z) 
			lines(c(NLmixAAnosaltn[i,z],NLmixAAsaltn[i,z] ) ~ c(1,2),col= rgb(0,(1-range01(mixAAdens)[z])*0.6,0)) ) 
		sapply(1:ncol(NLmixAAsaltn), function(z) 
			lines(c(NLmixAAsaltn[i,z],NLmixAAnosaltm[i,z] ) ~ c(2,3),col= rgb(0,(1-range01(mixAAdens)[z])*0.6,0)) ) 
		sapply(1:ncol(NLmixAAsaltn), function(z) 
			lines(c(NLmixAAnosaltm[i,z],NLmixAAsaltm[i,z] ) ~ c(3,4),col= rgb(0,(1-range01(mixAAdens)[z])*0.6,0)) ) 
		abline(v=2.5,lty=2)
		if(i==6){axis(side=1,at=c(1,2,3,4),labels=c("0","0.8","0","0.8"),las=2)}
		if(i==6){mtext(side=1,"salt g/L",line=3,cex=0.85)}
		if(i==1){mtext("disrupted",side=3, line=0.75,cex=0.75,adj=0)}
		if(i==1){mtext("local",side=3, line=0.75,cex=0.75,adj=1)}
		if(i==1){mtext("micr.",side=3, line=0.1,cex=0.75,adj=0)}
		if(i==1){mtext("micr.",side=3, line=0.1,cex=0.75,adj=1)}
		if(i!=1){text(x=4,y=ylimsB[[i]][2]*.9,paste(round(percentrankvar[i]*100),"%"))
			} else{		text(x=4,y=ylimsB[[i]][2]*0.8,paste(round(percentrankvar[i]*100),"%"))}
	plot(c(NLmixAAnosaltn[i,],NLmixAAnosaltm[i,],NLmixAAsaltn[i,],NLmixAAsaltm[i,])~rep(mixAAdens,times=4), 
		ylim=ylimsB[[i]],
		 xlab="",xaxt="n",yaxt="n",ylab="" ,  pch=1, col=rgb(0,0,0,alpha=0.5))
		if(i==1){lines(predBRTm~rddens.s)}
		if(i==1){polygon(c(rddens.s,rev(rddens.s)), c(predBRTi[1,], rev(predBRTi[2,])),col=rgb(0,0,0,alpha=0.5),border=NA )  }
		if(i==3){lines(exp(predBATm)~rddens.s)}
		if(i==3){polygon(c(rddens.s,rev(rddens.s)), exp(c(predBATi[1,], rev(predBATi[2,]))),col=rgb(0,0,0,alpha=0.5),border=NA )  }
		if(i==8){lines(exp(predMBTm)~rddens.s)}
		if(i==8){polygon(c(rddens.s,rev(rddens.s)), exp(c(predMBTi[1,], rev(predMBTi[2,]))),col=rgb(0,0,0,alpha=0.5),border=NA )  }
		if(i==5){lines(exp(predAPTm)~rddens.s)}
		if(i==5){polygon(c(rddens.s,rev(rddens.s)), exp(c(predAPTi[1,], rev(predAPTi[2,]))),col=rgb(0,0,0,alpha=0.5),border=NA )  }
		if(i==6){lines(exp(predPZTm)~rddens.s)}
		if(i==6){polygon(c(rddens.s,rev(rddens.s)), exp(c(predPZTi[1,], rev(predPZTi[2,]))),col=rgb(0,0,0,alpha=0.5),border=NA )  }
		if(i==6){mtext(side=1,"road density",line=3,cex=0.85)}
		if(i==6){axis(side=1,las=2,gap.axis=1.5)}
		if(i==2){mtext(side=3,"transformation product, nM",line=0.1, adj=0.3)}
		if(i==1){mtext(side=3,"% benzotriazole removed",line=1.7, adj = 0.25)}
	plot( c(NLbioAAgennosaltm[i,],NLbioAAgennosaltn[i,],NLbioAAgensaltm[i,],NLbioAAgensaltn[i,] ) ~c(NLbiopxgentrt[1,] ,NLbiopxgentrt[2,] , NLbiopxgentrt[3,] ,NLbiopxgentrt[4,] ), ylim=ylimsB[[i]],
		 xlab="",xaxt="n",yaxt="n",ylab="" ,  pch=1, col=rgb(0,0,0,alpha=0.5))
	if(i==1){lines(predBpxm~dsqm.s)}
		if(i==1){polygon(c(dsqm.s,rev(dsqm.s)), c(predBpxi[1,], rev(predBpxi[2,])),col=rgb(0,0,0,alpha=0.5),border=NA )  }
		if(i==2){lines(predGpxm~dsqm.s)}
		if(i==2){polygon(c(dsqm.s,rev(dsqm.s)), c(predGpxi[1,], rev(predGpxi[2,])),col=rgb(0,0,0,alpha=0.5),border=NA )  }
		if(i==6){mtext(side=1,"plant growth",line=3,cex=0.85)}
		if(i==6){mtext(side=1,bquote("mm"^2),line=4.25,cex=0.85)}
		if(i==6){axis(side=1,las=2,gap.axis=1.5)}
	plot( c(NLbioAAgennosaltm[i,],NLbioAAgennosaltn[i,],NLbioAAgensaltm[i,],NLbioAAgensaltn[i,] ) ~c(NLbioodgentrt[1,] ,NLbioodgentrt[2,] , NLbioodgentrt[3,] ,NLbioodgentrt[4,] ), ylim=ylimsB[[i]],
		 xlab="",xaxt="n",yaxt="n",ylab="" ,  pch=1, col=rgb(0,0,0,alpha=0.5))
		if(i==3){lines(exp(predBaom)~od.s)}
		if(i==3){polygon(c(od.s,rev(od.s)), exp(c(predBaoi[1,], rev(predBaoi[2,]))),col=rgb(0,0,0,alpha=0.5),border=NA )  }
		if(i==2){lines(predGodm~od.s)}
		if(i==2){polygon(c(od.s,rev(od.s)), c(predGodi[1,], rev(predGodi[2,])),col=rgb(0,0,0,alpha=0.5),border=NA )  }
		if(i==7){lines(predMBom~od.s)}
		if(i==7){polygon(c(od.s,rev(od.s)), c(predMBoi[1,], rev(predMBoi[2,])),col=rgb(0,0,0,alpha=0.5),border=NA )  }
		if(i==6){lines(exp(predpzodm)~od.s)}
		if(i==6){polygon(c(od.s,rev(od.s)), exp(c(predpzodi[1,], rev(predpzodi[2,])) ),col=rgb(0,0,0,alpha=0.5),border=NA )  }
		if(i==6){mtext(side=1,"microbe growth",line=3,cex=0.85)}
		if(i==6){mtext(side=1,"ln(OD)",line=4.25,cex=0.85)}
		if(i==6){axis(side=1,las=2,gap.axis=1.5)}
	}
	}
dev.off()



################################################
####IS PRODUCT TRANSFORMATION RELATED TO MICROBIOME COMPOSISTION?
################################################

##prepare data
mixAAsaltm <- sapply(sort(unique(transf$Genotype)), function(g) sapply(1:ncol(mixAA), function(p) mean(mixAA[transf$Genotype==g & transf$Salt=="0.8" &  transf$Microbes=="Yes",p],na.rm=T) ))
mixAAnosaltm <- sapply(sort(unique(transf$Genotype)), function(g) sapply(1:ncol(mixAA), function(p) mean(mixAA[transf$Genotype==g & transf$Salt=="0" &  transf$Microbes=="Yes",p],na.rm=T)))

##prepare tree
hier.phy3 <- hier.phy
	treegenus <- sapply(1:250,function(T) tail(unlist(strsplit(mastertax.s[ hier.phy$tip.label[T], 6 ],"g__",fixed=T)),1))
	treesp <- sapply(1:250,function(T) tail(unlist(strsplit(mastertax.s[ hier.phy$tip.label[T], 7],"s__",fixed=T)),1))
	treefam <- sapply(1:250,function(T) tail(unlist(strsplit(mastertax.s[ hier.phy$tip.label[T], 5],"f__",fixed=T)),1))
	treefam[treefam=="[Weeksellaceae]"] <- "Weeksellaceae"
	treefam[treefam=="[Exiguobacteraceae]"] <- "Exiguobacteraceae"
	treefam[treefam=="[Chromatiaceae]"] <- "Chromatiaceae"
	tiplabeldat <- data.frame(label = hier.phy3$tip.label, fam=treefam, genus = treegenus,
                species = treesp)
hier.phy3$edge.length[is.na(hier.phy3$edge.length)] <- 0

#function for microbe phylogeny- products analysis
tipxout <- function(prdtrtmat, balances, tree ){ #prdtrtmat has rownames for products and colnames for genos
	ntip <- tree$Nnode + 1
	nnode <- tree$Nnode 
	cormat <- sapply(1:nrow(prdtrtmat), function(z) cor(prdtrtmat[z,],t(balances)) )
	colnames(cormat) <- rownames(prdtrtmat)
	rownames(cormat) <- rownames(balances)
	nodesout <- sapply(1:ncol(cormat), function(z)  abs(cormat[,z]) > quantile(abs(cormat[,z]),0.95)) 
	tipsout <- lapply(1:ncol(nodesout), function(z) sapply(row.names(cormat)[which(nodesout[,z])],
										 function(name)  Descendants(tree, node= ntip + which(tree$node.label==name),type="tips" ) )  )
	tipisoutR <- sapply(1:length(tipsout), function(l) sapply(1:ntip, function(T) sum(unlist(tipsout[[l]]) ==T) ) )
	return(tipisoutR)
}

#run analysis for transformation products
maanss <-rbind(mixAAsaltm,mixAAnosaltm)
	rownames(maanss) <- c( colnames(mixAA), paste(colnames(mixAA),"_S",sep="")  )
tipisout <- tipxout(maanss,subbalphy,hier.phy)
ntipisout <- list()
for(i in 1:1000){
	nmaanss <- maanss[,sample(1:ncol(maanss),repl=F)]
	ntipisout[[i]] <-  tipxout(nmaanss,subbalphy,hier.phy)
}
reord <- lapply(1:ncol(tipisout), function(prd) sapply(1:length(ntipisout), function(n) ntipisout[[n]][,prd]) )
int <- sapply(1:ncol(tipisout), function(p) sapply(1:nrow(tipisout),function(T) findInterval(tipisout[T,p],sort(reord[[p]][T,]) ,left.open=T )/length(reord[[p]][T,]) ) )
tipisoutpernode <- int
tipisoutpernode [tipisoutpernode<0.95] <- 0
reorderTOpn <- tipisoutpernode[,c(9,1,10,2,11,3,12,4,13,5,14,6,15,7,16,8)]#no salt first #(1,9,2,10,3,11,4,12,5,13,6,14,7,15,8,16)]
row.names(reorderTOpn)<- hier.phy3$tip.label
rotodfpn <- as.data.frame(reorderTOpn)
colnames(rotodfpn) <- c(" ","  "," ","  "," ","  "," ","  "," ","  "," ","  "," ","  "," ","  ")

###run analysis for microbe and plant growth
biogeno4corM <- (biogenotrt[c(42,44,46,48,154,156,158,160,50,52,54,56),colnames(biogenotrt)%in%transf$Genotype]) #alternating beacuse selecting only inoculated treatments
#the rows selected correspond to the treatments below - all are inoculated treatments
row.names(biogeno4corM) <- c("F10.00","F10.0S","F10.B0","F10.BS","dsqm.00","dsqm.0S","dsqm.B0","dsqm.BS","od600.00","od600.0S","od600.B0","od600.BS")
colnames(bal.phy)==colnames(biogeno4corM) #check names

tipisoutbio <- tipxout(biogeno4corM,bal.phy,hier.phy)
nBtipisout <- list()
for(i in 1:1000){
	nbiogeno4corM <- biogeno4corM[,sample(1:ncol(biogeno4corM),repl=F)]
	nBtipisout[[i]] <-  tipxout(nbiogeno4corM,bal.phy,hier.phy)
}
Breord <- lapply(1:ncol(tipisoutbio), function(prd) sapply(1:length(nBtipisout), function(n) nBtipisout[[n]][,prd]) )
Bint <- sapply(1:ncol(tipisoutbio), function(p) sapply(1:nrow(tipisoutbio),function(T) findInterval(tipisoutbio[T,p],sort(Breord[[p]][T,]) ,left.open=T )/length(Breord[[p]][T,]) ) )
Btipisoutpernode <- Bint
Btipisoutpernode [Btipisoutpernode<0.95] <- 0
row.names(Btipisoutpernode)<- hier.phy3$tip.label
biotodfpn <- as.data.frame(Btipisoutpernode)
colnames(biotodfpn) <- c(" ","  "," ","  "," ","  "," ","  "," ","  "," ","  ")


#collect some node numbers, names for pieces of the tree
hier.phy2 <- hier.phy
hier.phy2$tip.label <- as.character(treefam)
pdf("check_nodes.pdf", width = 9, height=24)
gtree <- ggtree(hier.phy2,ladderize=F,color="gray") + geom_tiplab()  + geom_nodelab() + scale_fill_manual(values =rep(c(rgb(0,0,0.4),rgb(0,0,0.55),rgb(0,0,0.7),rgb(0,0,0.95)), times=9) )
gtree + xlim(0,25)
dev.off()
hier.phyg <- hier.phy2
hier.phyg$tip.label <- as.character(treegenus)
pdf("check_nodes_genera.pdf", width = 9, height=24)
gtree <- ggtree(hier.phyg,ladderize=F,color="gray") + geom_tiplab()  + geom_nodelab() + scale_fill_manual(values =rep(c(rgb(0,0,0.4),rgb(0,0,0.55),rgb(0,0,0.7),rgb(0,0,0.95)), times=9) )
gtree + xlim(0,25)
dev.off()
hier.physp <- hier.phy2
hier.physp$tip.label <- as.character(treesp)
pdf("check_nodes_sp.pdf", width = 9, height=24)
gtree <- ggtree(hier.physp,ladderize=F,color="gray") + geom_tiplab()  + geom_nodelab() + scale_fill_manual(values =rep(c(rgb(0,0,0.4),rgb(0,0,0.55),rgb(0,0,0.7),rgb(0,0,0.95)), times=9) )
gtree + xlim(0,25)
dev.off()


#inspecting tree, visually determine match between nodes and taxonomy assigned family names
famnodes <-  c("y97","y96","y32","y22","y16", 
								"y21","y40",
								"y43","y18",
								"y78","y15")
famnodeNames <- c("Aeromonadaceae","Enterobacteriaceae","Pseudomonadaceae","Rhizobiaceae +","Comamonadaceae +", #Hyphomicrobiales is primarily rhizobiacae; burkholderiales is primarily comamonadaceae
							"Caulobacteriaceae","Bacillaceae +", #bacillaceae contains exiguobacteriaceae
						"Flavobacteriaceae","Moraxellaceae",
						"Chromatiaceae","Xanthomonadaceae")
		 # note that most recent greengenes taxonomy information does not fully reflect the most recent bacterial phylogeny
		 # Thus many (most) of these families are paraphyletic
famnodeLs <- sapply(famnodes, function(z) which(hier.phy3$node.label == z) )
hier.phynodenams <- hier.phy3
hier.phynodenams$node.label[!(hier.phynodenams$node.label%in%famnodes)] <- "  "
hier.phynodenams$node.label[famnodeLs] <- famnodeNames
tiplabeldat$fam2 <- tiplabeldat$fam
tiplabeldat$fam2[is.na(tiplabeldat$fam2)] <- "NA"
tiplabeldat$fam2[tiplabeldat$fam2==""] <- "unk"
filldat <- data.frame(id= unlist(famnodeLs)+250, type=famnodeNames)


###PLOT FIGURE
pdf("ggtreeheat_interval_shortphylo_TPS.pdf",width=7,height=5)
gtree <- ggtree(hier.phynodenams,ladderize=F,color="gray") + geom_hilight(data=filldat, mapping=aes(node=id, fill=type,extend=9),show.legend=FALSE)  + geom_nodelab() + scale_fill_manual(values =
			c(rgb(0,0.1,0.2),rgb(0.1,0,0.26),rgb(0,0,0.32),rgb(0,0.1,0.38),rgb(0.1,0,0.44), rgb(0,0,0.5),rgb(0,0.1,0.56),rgb(0.1,0,0.72),rgb(0,0,0.8),rgb(0,0.1,0.86), rgb(0.1,0,0.92)) )
p <- gtree + new_scale_fill()
p<-gheatmap(p, rotodfpn[,c(1:2)], offset=0,  width=.1,  high=rgb(0,0,0),low="gray")#colnames_angle=90, colnames_offset_y = 20,
p1 <- p + new_scale_fill()
p1<-gheatmap(p1, rotodfpn[,c(3:4)], offset=0.075,  width=.1 )#, colnames_angle=90, colnames_offset_y = 20)#, high=rgb(0,0,0.75),low="gray")
p1<-gheatmap(p1, rotodfpn[,c(5:6)], offset=0.15,  width=.1, high=rgb(0,.75,0.5),low="gray")#, colnames_angle=90, colnames_offset_y = 20
p2 <- p1 + new_scale_fill()
p2<-gheatmap(p2, rotodfpn[,c(13:14)], offset=0.225, width=.1 )#, colnames_angle=90, colnames_offset_y = 13)
p2<-gheatmap(p2, rotodfpn[,c(15:16)], offset=0.3,   width=.1, high=rgb(0.5,0.25,1),low="gray") #colnames_angle=90, colnames_offset_y = 21, 
p3 <- p2 + new_scale_fill()
p3<-gheatmap(p3, rotodfpn[,c(9:10)], offset=0.375,  width=.1 )#, colnames_angle=90, colnames_offset_y = 25)#,  high=rgb(0,0,0.75),low="gray") 
p3<-gheatmap(p3, rotodfpn[,c(11:12)], offset=0.45,  width=.1,  high=rgb(1.0,0.54,0),low="gray") # legend_title="interval")   , colnames_angle=90, colnames_offset_y = 25,      
p3   #
dev.off()

pdf("ggtreeheat_interval_shortphylo_bio.pdf",width=5,height=5)
gtree <- ggtree(hier.phynodenams,ladderize=F,color="gray") + geom_hilight(data=filldat, mapping=aes(node=id, fill=type,extend=9),show.legend=FALSE)  + geom_nodelab() + scale_fill_manual(values =
			c(rgb(0,0.1,0.2),rgb(0.1,0,0.26),rgb(0,0,0.32),rgb(0,0.1,0.38),rgb(0.1,0,0.44), rgb(0,0,0.5),rgb(0,0.1,0.56),rgb(0.1,0,0.72),rgb(0,0,0.8),rgb(0,0.1,0.86), rgb(0.1,0,0.92)) )
p <- gtree + new_scale_fill()
p<-gheatmap(p, biotodfpn[,c(7:8)], offset=0,  width=.1, )#colnames_angle=90, colnames_offset_y = 20,
p1<-gheatmap(p, biotodfpn[,c(11:12)], offset=0.075,  width=.1,  high=rgb(0,0,0),low="gray")#, colnames_angle=90, colnames_offset_y = 20)#, high=rgb(0,0,0.75),low="gray")
p1
dev.off()


##### Follow ups based on balance analysis.


####Use balances to identify clades of interest
#recover balances most correlated, visually inspect the nodes at which these fall (e.g. nodes underlying highlighted tips in main balance figure)
	reccormat <- sapply(1:nrow(rbind(mixAAsaltm,mixAAnosaltm)), function(z) cor(rbind(mixAAsaltm,mixAAnosaltm)[z,],t(subbalphy)) )
	reccormatO <- reccormat[,c(1,9,2,10,3,11,4,12,5,13,6,14,7,15,8,16)]
tagged <- sapply(1:ncol(reccormatO), function(z) paste("y",which(abs(reccormatO[,z])> quantile(abs(reccormatO[,z]),0.95)) -1 , sep=""))
###REPLACE numbers at HERE with column number for product of interest. does not run inside for loop due to ggtree errors
pdf(paste("check_cor_nodes",1,".pdf", sep=""), width = 9, height=24) ##HERE, currently 1
hier.taggednod <- hier.phy2
hier.taggednod$tip.label <- as.character(treefam)
hier.taggednod$node.label[!(hier.taggednod$node.label%in%tagged[,1])] <- "  " ###HERE
gtree <- ggtree(hier.taggednod,ladderize=F,color="gray") + geom_tiplab()  + geom_nodelab() + scale_fill_manual(values =rep(c(rgb(0,0,0.4),rgb(0,0,0.55),rgb(0,0,0.7),rgb(0,0,0.95)), times=9) )
gtree + xlim(0,25)
dev.off()
#BZT: 1 and 2 share y131 and y147 & nodes below, and for column 1, node y73 is just before these; they also share 105, but is v small clade
#glyBZT: 3 and 4, no interesting tips of significance in balance figure, skip; there are some shared nodes however - y23, y39, y137, y118, y131
#BZT a+ aa: 5 and 6, both share hits in rhizobiaceae, though only y84 itself shared (5 has an additional deeper tree hit, 6 hits tippier, try y22 which splits a big clade, and y10 as origin), and pseudomonadaceae, but mostly not the same nodes try y119 which splits sets of tagged nodes.
#col 7 and 8 are aniline.
#a3ph: col 9 and ten, PRIMARILY Aeromonadaceae hits, either within or leading to. Try y138, which is major clade split, y77, which splits from entero*. y122 shared but one branch a single, fairly rare taxon.  y13 also appears in both (in both early chitin* and cytoph*)
#phz: col 11 shares aeromonadaceae hits and y13 with 9 and ten, whereas 12 has no tips passing thesholds in our figure
#mBZT: for 13,  many of the highest ones are in xanthomonadaceae. might consider pseudomonadaceae; column 14 nearly no thresholds passed according to full figure; y195 in aero shared, y94 in pseudos, col14 has bunch of hits in pseudomonadaceae
#moxyBZT: skip, not much of interest in balance figure. only cluster passing threshold occurs at less than 1% and in only 3 inocula

###
#get sums of taxa subtending families:
phylopseudos <- colSums(mastertab.s[hier.phy$tip.label[unlist(Descendants(hier.phy,which(hier.phy$node.label=="y32")+250,type="tips"))],])
phyloaeros <- colSums(mastertab.s[hier.phy$tip.label[unlist(Descendants(hier.phy,which(hier.phy$node.label=="y97")+250,type="tips"))],])
phylorhizos <- colSums(mastertab.s[hier.phy$tip.label[unlist(Descendants(hier.phy,which(hier.phy$node.label=="y22")+250,type="tips"))],])
phyloxanthos <- colSums(mastertab.s[hier.phy$tip.label[unlist(Descendants(hier.phy,which(hier.phy$node.label=="y15")+250,type="tips"))],])
phyloenteros <- colSums(mastertab.s[hier.phy$tip.label[unlist(Descendants(hier.phy,which(hier.phy$node.label=="y96")+250,type="tips"))],])
tolumonasclade <- colSums(mastertab.s[hier.phy$tip.label[unlist(Descendants(hier.phy,which(hier.phy$node.label=="y154")+250,type="tips"))],])
aeromonasclade <- colSums(mastertab.s[hier.phy$tip.label[unlist(Descendants(hier.phy,which(hier.phy$node.label=="y155")+250,type="tips"))],])
caulobact <- colSums(mastertab.s[hier.phy$tip.label[unlist(Descendants(hier.phy,which(hier.phy$node.label=="y21")+250,type="tips"))],])
rhizoout <- colSums(mastertab.s[hier.phy$tip.label[unlist(Descendants(hier.phy,which(hier.phy$node.label=="y39")+250,type="tips"))],])
rhizoin <- colSums(mastertab.s[hier.phy$tip.label[unlist(Descendants(hier.phy,which(hier.phy$node.label=="y38")+250,type="tips"))],])
pveron <- colSums(mastertab.s[hier.phy$tip.label[unlist(Descendants(hier.phy,which(hier.phy$node.label=="y133")+250,type="tips"))],])
pother <- colSums(mastertab.s[hier.phy$tip.label[unlist(Descendants(hier.phy,which(hier.phy$node.label=="y134")+250,type="tips"))],])

#deterimine balance direction:
cor(subbalphy[139,],aeromonasclade); cor(subbalphy[139,],tolumonasclade)#as balance increases, relative aeromonas increases
cor(subbalphy[11,],phylorhizos); cor(subbalphy[11,],caulobact) #as balance increases, relative caulobacter decreases
cor(subbalphy[78,],phyloenteros); cor(subbalphy[78,],phyloaeros) # as balance increases, relative entero. decrease and relative aero. increase
cor(subbalphy[23,],rhizoout); cor(subbalphy[23,],rhizoin) # as balance increases, outgroup rhizobia relatively increase, and ingroup rhizobiaceae somewhat relatively decrease
cor(subbalphy[120,],pveron); cor(subbalphy[120,],pother) # as balance increases, relative veronii decrease



# rhizobiaceae and pseudomonadaceae with  bztaa, across salt
BZTaaBacdat <- data.frame(BZTaa = c(mixAAsaltm[3,],mixAAnosaltm[3,] ),
						y22=rep(subbalphy[23,],times=2),y10=rep(subbalphy[11,],times=2), y17=rep(subbalphy[18,],times=2), y119=rep(subbalphy[120,],times=2), 
						phylorhizos = rep(phylorhizos,times=2), #rhizos= rep(colSums(mastertab.s[which(mastertax.s[,5]=="f__Rhizobiaceae"),]),times=2),
						phylopseudos = rep(phylopseudos,times=2), #pseudos= rep(colSums(mastertab.s[which(mastertax.s[,5]=="f__Pseudomonadaceae"),]),times=2),
						rddens = rep(tapply(transf$density,transf$Genotype,mean),times=2)/1000) 
mgenos <- transf$Genotype[ transf$Microbes=="Yes"]
mindex <- sapply(mgenos, function(z) which(colnames(subbalphy)==z)) 
#full replicate transformation products data
BZTaaBacdatf <- data.frame(BZTaa = mixAA[transf$Microbes=="Yes",3], 
					phylorhizos = phylorhizos[ mindex  ], #rhizos= colSums(mastertab.s[which(mastertax.s[,5]=="f__Rhizobiaceae"),])[mindex],
					phylopseudos = phylopseudos[ mindex  ], #pseudos= colSums(mastertab.s[which(mastertax.s[,5]=="f__Pseudomonadaceae"),])[mindex],
					y22 = subbalphy[23,mindex], y10= subbalphy[11,mindex], y17 = subbalphy[18,mindex], y119 = subbalphy[120,mindex],
					rddens = transf$density[transf$Microbes=="Yes"]/1000)
summary(MCMCglmm(BZTaa~rddens ,data=BZTaaBacdatf,nitt=100000,burnin=10000,thin=10,verbose=F))
summary(MCMCglmm(BZTaa~log(phylorhizos + min(phylorhizos[phylorhizos>0]) ) ,data=BZTaaBacdatf,nitt=100000,burnin=10000,thin=10,verbose=F))
summary(MCMCglmm(BZTaa~log(phylorhizos + min(phylorhizos[phylorhizos>0]) )+rddens ,data=BZTaaBacdatf,nitt=100000,burnin=10000,thin=10,verbose=F))
summary(MCMCglmm(BZTaa~y22 ,data=BZTaaBacdatf,nitt=100000,burnin=10000,thin=10,verbose=F))
summary(MCMCglmm(BZTaa~y22+rddens ,data=BZTaaBacdatf,nitt=100000,burnin=10000,thin=10,verbose=F))
summary(MCMCglmm(BZTaa~y10 ,data=BZTaaBacdatf,nitt=100000,burnin=10000,thin=10,verbose=F))
summary(MCMCglmm(BZTaa~y10+rddens ,data=BZTaaBacdatf,nitt=100000,burnin=10000,thin=10,verbose=F))
summary(MCMCglmm(BZTaa~y119 ,data=BZTaaBacdatf,nitt=100000,burnin=10000,thin=10,verbose=F))
summary(MCMCglmm(BZTaa~y119+rddens ,data=BZTaaBacdatf,nitt=100000,burnin=10000,thin=10,verbose=F))
summary(MCMCglmm(BZTaa~phylopseudos ,data=BZTaaBacdatf,nitt=100000,burnin=10000,thin=10,verbose=F))
summary(MCMCglmm(BZTaa~phylopseudos+rddens ,data=BZTaaBacdatf,nitt=100000,burnin=10000,thin=10,verbose=F))
#y10 best. y22, y119 and phylorhizos also sig. phylorhizos sig pos.
#y10 is split between Rhizobacteriaceae+close families clade and caulobacteriaceae.
#may want to refit simplest to undestand for example plot:
bztaaBrhizo <- (MCMCglmm(BZTaa~ log(phylorhizos + min(phylorhizos[phylorhizos>0]) ),data=BZTaaBacdatf,nitt=100000,burnin=10000,thin=10,verbose=F))#refit simple to understand
rhizo.s <- seq(from=min(log(BZTaaBacdatf$phylorhizos+ min(phylorhizos[phylorhizos>0]))),to=max(log(BZTaaBacdatf$phylorhizos+ min(phylorhizos[phylorhizos>0]))),length.out=1000)#
#may want to refit best
#bztaaBrhizo <- (MCMCglmm(BZTaa~ y10,data=BZTaaBacdatf,nitt=100000,burnin=10000,thin=10,verbose=F))#refit simple to understand
#rhizo.s <- seq(from=min(BZTaaBacdatf$y10),to=max(BZTaaBacdatf$y10),length.out=1000)#
solbztaaR <- bztaaBrhizo$Sol
prd.bztaaR <- sapply(rhizo.s, function(z) exp(mean(solbztaaR[,1] + solbztaaR[,2]*z )) ) 
hpdi.bztaaR <- sapply(rhizo.s, function(z) exp(HPDi(solbztaaR[,1] + solbztaaR[,2]*z ,.95)) )

#aeromonadaceae with a3ph, and phenazine -- phz only without salt.
data3ph <- data.frame(amino3ph = c(mixAAsaltm[5,],mixAAnosaltm[5,] ), phyloaeros=rep(phyloaeros,times=2),
					 y138 = rep(subbalphy[139,],times=2), y77=rep( subbalphy[78,], times=2), rddens = rep(tapply(transf$density,transf$Genotype,mean),times=2)/1000)
data3phf <- data.frame(amino3ph = mixAA$amino_3_phenol[transf$Microbes=="Yes"], phyloaeros = phyloaeros[ mindex  ], 
					y138 = subbalphy[139,mindex], y77 = subbalphy[78,mindex],
					rddens = transf$density[transf$Salt==0.8 & transf$Microbes=="Yes"]/1000)
summary(MCMCglmm(amino3ph~rddens ,data=data3phf,nitt=100000,burnin=10000,thin=10,verbose=F))
summary(MCMCglmm(amino3ph~phyloaeros ,data=data3phf,nitt=100000,burnin=10000,thin=10,verbose=F))
summary(MCMCglmm(amino3ph~phyloaeros + rddens ,data=data3phf,nitt=100000,burnin=10000,thin=10,verbose=F))
summary(MCMCglmm(amino3ph~y77 ,data=data3phf,nitt=100000,burnin=10000,thin=10,verbose=F))
summary(MCMCglmm(amino3ph~y77 + rddens ,data=data3phf,nitt=100000,burnin=10000,thin=10,verbose=F))
summary(MCMCglmm(amino3ph~y138 ,data=data3phf,nitt=100000,burnin=10000,thin=10,verbose=F))
summary(MCMCglmm(amino3ph~y138 + rddens ,data=data3phf,nitt=100000,burnin=10000,thin=10,verbose=F))
#y138 sig and best, note phyloaeros is n.s.
#refit best model
a3phAeros <- MCMCglmm(amino3ph~y138,data=data3phf,nitt=100000,burnin=10000,thin=10,verbose=F)
aer.s <- seq(from=min(data3ph$y138),to=max(data3ph$y138),length.out=1000)
sola3paos <- a3phAeros$Sol
prd.a3paos <- sapply(aer.s, function(z) exp(mean(sola3paos[,1] + sola3paos[,2]*z)))
hpdi.a3paos <- sapply(aer.s, function(z) exp(HPDi(sola3paos[,1] + sola3paos[,2]*z,.95)))


datphz <- data.frame(phz = c(mixAAsaltm[6,]), phyloaeros=phyloaeros,
					 y138 = subbalphy[139,], y77= subbalphy[78,],rddens = tapply(transf$density,transf$Genotype,mean)/1000)
smgenos <- transf$Genotype[transf$Salt==0.8 & transf$Microbes=="Yes"]
smindex <- sapply(smgenos, function(z) which(colnames(subbalphy)==z)) 
datphzf <- data.frame(phz = mixAA$phenazine[transf$Microbes=="Yes" & transf$Salt==0.8], phyloaeros = phyloaeros[ smindex  ], 
					y138 = subbalphy[139,smindex], y77 = subbalphy[78,smindex],
					rddens = transf$density[transf$Salt==0.8 & transf$Microbes=="Yes"]/1000)
summary(MCMCglmm(phz~rddens ,data=datphzf,nitt=100000,burnin=10000,thin=10,verbose=F))
summary(MCMCglmm(phz~phyloaeros ,data=datphzf,nitt=100000,burnin=10000,thin=10,verbose=F))
summary(MCMCglmm(phz~phyloaeros + rddens ,data=datphzf,nitt=100000,burnin=10000,thin=10,verbose=F))
summary(MCMCglmm(phz~y77 ,data=datphzf,nitt=100000,burnin=10000,thin=10,verbose=F))
summary(MCMCglmm(phz~y77 + rddens ,data=datphzf,nitt=100000,burnin=10000,thin=10,verbose=F))
summary(MCMCglmm(phz~y138 ,data=datphzf,nitt=100000,burnin=10000,thin=10,verbose=F))
summary(MCMCglmm(phz~y138 + rddens ,data=datphzf,nitt=100000,burnin=10000,thin=10,verbose=F))
#phyloaeros and y138 both sig alone, phyloaeros best by a little bit
phzAeros <- MCMCglmm(phz~y138,data=datphzf,nitt=100000,burnin=10000,thin=10,verbose=F)
aer.s <- seq(from=min(datphz$y138),to=max(datphz$y138),length.out=1000)
solphz <- phzAeros$Sol
prd.phz <- sapply(aer.s, function(z) exp(mean(solphz[,1] + solphz[,2]*z)))
hpdi.phz <- sapply(aer.s, function(z) exp(HPDi(solphz[,1] + solphz[,2]*z,.95)))

# xantho/pseduomonadaceae with 1mBZT, only with salt.
datmbzt <- data.frame(mBZT = mixAAsaltm[7,], phylopseudos = phylopseudos, phyloxanthos = phyloxanthos,
						rddens = tapply(transf$density,transf$Genotype,mean)/1000)
datmbztf <- data.frame(mBZT = mixAA$methylBZT[transf$Salt==0.8 & transf$Microbes=="Yes"], 
					phylopseudos = phylopseudos[ smindex  ], phyloxanthos = phyloxanthos[smindex],
					 rddens = transf$density[transf$Salt==0 & transf$Microbes=="Yes"]/1000)
summary(MCMCglmm(scale(mBZT)~ rddens,data=datmbztf,nitt=100000,burnin=10000,thin=10,verbose=F))
summary(MCMCglmm(scale(mBZT)~phylopseudos,data=datmbztf,nitt=100000,burnin=10000,thin=10,verbose=F))#
summary(MCMCglmm(scale(mBZT)~phylopseudos + rddens,data=datmbztf,nitt=100000,burnin=10000,thin=10,verbose=F))
summary(MCMCglmm(scale(mBZT)~phyloxanthos,data=datmbztf,nitt=100000,burnin=10000,thin=10,verbose=F))#best, sig and pos, rddens never sig
summary(MCMCglmm(scale(mBZT)~phyloxanthos + rddens,data=datmbztf,nitt=100000,burnin=10000,thin=10,verbose=F))
mBZTsaltX <- (MCMCglmm(scale(mBZT)~phylopseudos,data=datmbztf,nitt=100000,burnin=10000,thin=10,verbose=F))#refit best
pseud.s <- seq(from=min(datmbztf$phylopseudos),to=max(datmbztf$phylopseudos),length.out=1000)
solmbztX <- mBZTsaltX$Sol
prd.mbztX <- sapply(pseud.s, function(z) mean( (solmbztX[,1] + solmbztX[,2]*z)*sd(datmbztf$mBZT) + mean(datmbztf$mBZT))) 
hpdi.mbztX <- sapply(pseud.s, function(z) HPDi((solmbztX[,1] + solmbztX[,2]*z)*sd(datmbztf$mBZT) + mean(datmbztf$mBZT),.95))


#plot selected results, remove any log, or scale, and set units to nanomolar
pdf("exploratorymicrobegroups_short.pdf",width=6,height=6)
par(mfrow=c(2,2))
par(mar=c(5,5,1,1))
  plot(exp(BZTaa)~log(BZTaaBacdat$phylorhizos+ min(phylorhizos[phylorhizos>0])),data=BZTaaBacdat,pch=NA,xlab="",ylab="",cex=1.5,ylim=bufferX(exp(BZTaaBacdat$BZTaa),.25))
  	polygon(c(rhizo.s,rev(rhizo.s)),y=c(hpdi.bztaaR[1,],rev(hpdi.bztaaR[2,])),border=NA,col=rgb(0,0,0,alpha=0.25))
  	lines(prd.bztaaR~rhizo.s)
  	points(exp(BZTaa)~log(BZTaaBacdat$phylorhizos+ min(phylorhizos[phylorhizos>0])),data=BZTaaBacdat,pch=16,cex=1.5)
  	mtext("log of Rhizobacteriaceae+",side=1,line=2.5)
  	mtext("relative abundance",side=1,line=4)
	mtext("BZTa + BZTaa nM",side=2,line=2)
plot(mBZT~phylopseudos,data=datmbzt,xlab="",ylab="",cex=1.5,pch=NA,ylim=bufferX(datmbzt$mBZT,.05))#,col=rgb(range01(datPseudo$rddens),0,0))
	polygon(c(pseud.s,rev(pseud.s)),y=c(hpdi.mbztX[1,],rev(hpdi.mbztX[2,])),border=NA,col=rgb(0,0,0,alpha=0.25))
	lines(prd.mbztX~pseud.s)
	points(mBZT~phylopseudos,data=datmbzt,pch=16,cex=1.5)
	mtext("1-methylBZT nM",side=2,line=2)
	mtext("Pseudomonadaceae",side=1,line=2.5)
	mtext("relative abundance",side=1,line=4)
plot(exp(amino3ph)~y138,data=data3ph,pch=NA,xlab="",ylab="",cex=1.5,ylim=bufferX(exp(data3ph$amino3ph),.25))
	polygon(c(aer.s,rev(aer.s)),y=c(hpdi.a3paos[1,],rev(hpdi.a3paos[2,])),border=NA,col=rgb(0,0,0,alpha=0.25))
	lines(prd.a3paos~aer.s)
	points(exp(amino3ph)~y138,data=data3ph,pch=16,cex=1.5)
	mtext("amino-3-phenol nM",side=2,line=2)
	mtext("Balance of",side=1,line=2.5)
	mtext("Aeromonas/Tolumonas",side=1,line=4)
plot(exp(phz)~y138,data=datphz,pch=NA,xlab="",ylab="",cex=1.5,ylim=bufferX(exp(datphz$phz),.25))
	polygon(c(aer.s,rev(aer.s)),y=c(hpdi.phz[1,],rev(hpdi.phz[2,])),border=NA,col=rgb(0,0,0,alpha=0.25))
	lines(prd.phz~aer.s)
	points(exp(phz)~y138,data=datphz,pch=16,cex=1.5)
	mtext("phenazine nM",side=2,line=2)
	mtext("Balance of",side=1,line=2.5)
	mtext("Aeromonas/Tolumonas",side=1,line=4)
dev.off()



####################Diversity-Function Analysis

#test relationship between microbiome diversity and remaining benzotriazole
divdf <- data.frame(community = colnames(mastertab.s),richness=transfdivF,shannon = transfdivS,faithphylodiv=transfdivPD)
mrichindex <- sapply(mixAAtrtbio$Genotype, function(z) which(divdf$community == z) )
datrich <- data.frame(cbind(mixAAtrtbio[mixAAtrtbio$Microbes=="Yes",]), 
					divdf[mrichindex,][mixAAtrtbio$Microbes=="Yes",])
#even though there is more growth data observations in the treatments without benzotriazole, 
	#the point is to compare to the relationship between total transformed BZT and microbial community richness
	#therefore, the same slightly reduced dataset should be used. 
	#commented out code below shows similar results including additional treatments without BZT. 
	#results are very similar, with model for optical density becoming significant but slope very similar.
# 	mrichindexbio <- sapply(subbio$genotype, function(z) which(divdf$community == z) )
# 	datrichgrow <- data.frame(cbind(subbio[subbio$microbe=="+",]), 
# 						divdf[mrichindexbio,][subbio$microbe=="+",])
# 	datrichOD <- datrichgrow[!is.na(datrichgrow$lnOD600),]
# 	summary(MCMCglmm(deltasqmm~richness,data=datrichgrow,nitt=100000,burnin=10000,thin=10,verbose=F)) 
# 	summary(MCMCglmm(lnOD600~richness,data=datrichOD,nitt=100000,burnin=10000,thin=10,verbose=F)) 


#primary reported measure of richness
BZTr_richness <- (MCMCglmm(conc_corrected~richness,data=datrich,nitt=100000,burnin=10000,thin=10,verbose=F)) # marginally sig pos
summary(BZTr_richness)
rich.s <- seq(from = min(datrich$richness),to=max(datrich$richness),length.out=1000)
BZTrsol <- BZTr_richness$Sol
BZTr.prd <- sapply(rich.s, function(z) mean(BZTrsol[,1] + BZTrsol[,2]*z   ) )
BZTr.hpdi <- sapply(rich.s, function(z) HPDi(BZTrsol[,1] + BZTrsol[,2]*z, 0.95   ) )

#test relationship between microbiome diversity and growth metrics, 1 mean per treatment
duckgrowthrichness <- (MCMCglmm(deltasqmm~richness,data=datrich,nitt=100000,burnin=10000,thin=10,verbose=F)) 
micrgrowthrichness <- (MCMCglmm(lnOD600~richness,data=datrich,nitt=100000,burnin=10000,thin=10,verbose=F)) 
# sig neg for dsqm
summary(duckgrowthrichness)
summary(micrgrowthrichness)#n.s.
duckgrrsol <- duckgrowthrichness$Sol
duckgr.prd <- sapply(rich.s, function(z) mean(duckgrrsol[,1] + duckgrrsol[,2]*z   ) )
duckgr.hpdi <- sapply(rich.s, function(z) HPDi(duckgrrsol[,1] + duckgrrsol[,2]*z, 0.95   ) )

richBmn <- tapply(datrich$conc_corrected,datrich$community,mean)
richBse <- tapply(datrich$conc_corrected,datrich$community,std.error)
richOmn <- tapply(datrich$lnOD600,datrich$community,mean)
richOse <- tapply(datrich$lnOD600,datrich$community,std.error)
richDmn <- tapply(datrich$deltasqmm,datrich$community,mean)
richDse <- tapply(datrich$deltasqmm,datrich$community,std.error)
richlev <- tapply(datrich$richness,datrich$community,mean)

###diversity plots
pdf("diversity_bzt.pdf",width=3,height=3)
par(mar=c(4,4,1,1))
plot(richBmn~richlev,pch=NA,ylim=bufferX(c(richBmn-richBse,richBmn+richBse),p=0.05), ylab="", xlab="")
	polygon(c(rich.s,rev(rich.s)), c(BZTr.hpdi[1,],rev(BZTr.hpdi[2,])), border=NA,col=rgb(0,0,0,alph=0.25))
	lines(BZTr.prd~rich.s)
	points(richBmn~richlev,pch=16,cex=1.5)
	arrows(richlev, y0=richBmn-richBse, y1=richBmn+richBse, length=0)
	mtext("%BZT removed",side=2,line=2.5)
	text(24,81,"p < 0.01")
	mtext("microbial community richness",side=1,line=2.5)
dev.off()


pdf("diversity_growth.pdf",width=6,height=3)
par(mfrow=c(1,2))
par(mar=c(5,4,2,0))
par(oma=c(0,1,1,1))
plot(richDmn~richlev,pch=NA, ylim=bufferX(c(richDmn-richDse,richDmn+richDse),p=0.05), ylab="",xlab="")
	polygon(c(rich.s,rev(rich.s)), c(duckgr.hpdi[1,],rev(duckgr.hpdi[2,])), border=NA,col=rgb(0,0,0,alph=0.25))
	lines(duckgr.prd~rich.s)
	points(richDmn~richlev,pch=16,cex=1.5)
	arrows(richlev, y0=richDmn-richDse, y1=richDmn+richDse, length=0)
	mtext(bquote("plant growth, mm"^2),side=2,line=2.5)
	mtext("p < 0.05",side=3,line=1)
	mtext("microbial community richness",side=1,line=2.5, adj=-7)
plot(richOmn~richlev,pch=NA, ylim=bufferX(c(richOmn-richOse,richOmn+richOse),p=0.05), ylab="", xlab="")
	points(richOmn~richlev,pch=16,cex=1.5)
	arrows(richlev, y0=richOmn-richOse, y1=richOmn+richOse, length=0)
	mtext("microbe growth ln(OD)",side=2,line=2.5)
	mtext("n.s.",side=3,line=1)
dev.off()

#microbiome richness and urban location
commrddens <- tapply(biodat$rddns,biodat$genotype,mean)
rddensSub <- commrddens[names(commrddens)%in%divdf$community]
urbanrich <- data.frame(rddens = rddensSub, richness = divdf$richness[order(divdf$community)], shannon = divdf$shannon[order(divdf$community)], faithphylodiv = divdf$faithphylodiv[order(divdf$community)])
richrddens <- MCMCglmm(richness~rddens,data=urbanrich,nitt=100000,burnin=10000,thin=10,verbose=F)
summary(richrddens)
#no relationship
