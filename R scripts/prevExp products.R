##analysis of additional transformation products data for a previously published experiment
#@article{o2020resilience,
#  title={Resilience to multiple stressors in an aquatic plant and its microbiome},
#  author={O'Brien, Anna M and Yu, Zhu Hao and Luo, Dian-ya and Laurich, Jason and Passeport, Elodie and Frederickson, Megan E},
#  journal={American journal of botany},
#  volume={107},
#  number={2},
#  pages={273--285},
#  year={2020},
#  publisher={Wiley Online Library}
#}

library(MCMCglmm)
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

products <- read.csv("experiment A products.csv",header=T, stringsAsFactors=F) 
molmassTPs <- c(136 ,207,249 ,282, 166,150,134,94) # hydroxyBZT, bztalanine, bzt acetyl alanine, glycosylated bzt, pthalic acid
mgperL <- sapply(1:8, function(z) molmassTPs[z]*products[,z+5]*(1000/2.5))  #molarmass * mols per 2.5 mL * 1000/2.5 = grams /L
colnames(mgperL) <- colnames(products)[6:13]
prdngL <- data.frame(cbind(products[,c(1:4,14)],mgperL*1e06))
prdngL$bztaplusa <- prdngL$BZTalanine + prdngL$BZTacetylalanine

cbind(sapply(6:14, function(z) shapiro.test(prdngL[,z])$statistic) ,
	sapply(6:14, function(z) shapiro.test(log(prdngL[,z] + min(prdngL[,z][which(prdngL[,z]>0)])))$statistic))
logimproves <- sapply(6:14, function(z) shapiro.test(prdngL[,z])$statistic) < 
	sapply(6:14, function(z) shapiro.test(log(prdngL[,z] + min(prdngL[,z][which(prdngL[,z]>0)])))$statistic)	
prdngLmix <- prdngL
for(i in 1:9){
	if(logimproves[i]){
		prdngLmix[,i+5] <- log(prdngL[,i+5] + min(prdngL[,i+5][which(prdngL[,i+5]>0)]) )
	}
}

growth<- read.csv("BZTinBZS1_trts_within_geno_plusgrowth.csv",header=T, stringsAsFactors=F)#sorted identically to products
prdngLmix$logOD <- growth$od.mn
prdngLmix$ducksize <- growth$px.mn

#MODELS
##only 24 datapoints, so no interaction terms

#BZT decrease analysis exists in O'Brien et al 2020, repeated here only for completeness for reader
BrmBSM <- MCMCglmm(BZT.percent.decrease ~ BZT_init + salt + microbes , random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F)
BrmPXOD <- MCMCglmm(BZT.percent.decrease ~ logOD + ducksize , random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F)
tapply(prdngLmix$BZT.percent.decrease,prdngLmix$BZT_init,mean); tapply(prdngLmix$BZT.percent.decrease,prdngLmix$BZT_init,std.error)
tapply(prdngLmix$BZT.percent.decrease,prdngLmix$salt,mean); tapply(prdngLmix$BZT.percent.decrease,prdngLmix$salt,std.error)

#phytotransformation
summary(MCMCglmm(bztaplusa ~ BZT_init + salt + microbes , random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))# dic back to similar to full model
summary(MCMCglmm(bztaplusa ~ BZT_init + salt  , random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))# dic back to similar to full model
#
summary(MCMCglmm(bztaplusa ~ logOD + ducksize, random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))
summary(MCMCglmm(bztaplusa ~ ducksize, random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))
#
summary(MCMCglmm(glycosylatedBZT ~ BZT_init + salt + microbes, random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))
summary(MCMCglmm(glycosylatedBZT ~ BZT_init + salt, random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))
#
summary(MCMCglmm(glycosylatedBZT ~ logOD + ducksize, random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))
summary(MCMCglmm(glycosylatedBZT ~ ducksize, random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))
#
BZTaaBS <- MCMCglmm(bztaplusa ~ BZT_init + salt  , random = ~ genotype, data=prdngLmix, nitt=100000,verbose=F)
BZTaaPX <- MCMCglmm(bztaplusa ~ ducksize, random = ~ genotype, data=prdngLmix, nitt=100000,verbose=F)
glyBS <- MCMCglmm(glycosylatedBZT ~ BZT_init + salt, random = ~ genotype, data=prdngLmix, nitt=100000,verbose=F)
glyPX <- MCMCglmm(glycosylatedBZT ~ ducksize, random = ~ genotype, data=prdngLmix, nitt=100000,verbose=F)
tapply(prdngL$bztaplusa, prdngL$BZT_init,mean); tapply(prdngL$bztaplusa, prdngL$BZT_init,std.error)
tapply(prdngL$glycosylatedBZT, prdngL$BZT_init,mean); tapply(prdngL$glycosylatedBZT, prdngL$BZT_init,std.error)
tapply(prdngL$bztaplusa, prdngL$salt,mean); tapply(prdngL$bztaplusa, prdngL$salt,std.error)
tapply(prdngL$glycosylatedBZT, prdngL$salt,mean); tapply(prdngL$glycosylatedBZT, prdngL$salt,std.error)

#microbial transformation
summary(MCMCglmm(methoxyBZT ~ BZT_init + salt + microbes, random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))
summary(MCMCglmm(methoxyBZT ~ BZT_init + salt , random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))
summary(MCMCglmm(methoxyBZT ~ BZT_init , random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))
#
summary(MCMCglmm(methoxyBZT ~ logOD + ducksize, random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))
summary(MCMCglmm(methoxyBZT ~ logOD , random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))
summary(MCMCglmm(methoxyBZT ~ 1, random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))
#
summary(MCMCglmm(methylBZT ~ BZT_init + salt + microbes, random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))
summary(MCMCglmm(methylBZT ~ BZT_init + salt , random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))
summary(MCMCglmm(methylBZT ~ BZT_init  , random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))
summary(MCMCglmm(methylBZT ~ 1  , random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))
#
summary(MCMCglmm(methylBZT ~ logOD + ducksize, random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))
summary(MCMCglmm(methylBZT ~ logOD , random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))#ns. but rm fits worse
#
mxyB <- MCMCglmm(methoxyBZT ~ BZT_init , random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F)
tapply(prdngL$methoxyBZT, prdngL$BZT_init,mean); tapply(prdngL$methoxyBZT, prdngL$BZT_init,std.error)


#mixed pathway compounds
summary(MCMCglmm(aniline ~ BZT_init + salt + microbes, random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))
summary(MCMCglmm(aniline ~ BZT_init + salt, random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))
#
summary(MCMCglmm(aniline ~ logOD + ducksize, random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))
summary(MCMCglmm(aniline ~ logOD, random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))#marginal
#
summary(MCMCglmm(hydroxyBZT ~ BZT_init + salt + microbes, random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))#slightly worse
summary(MCMCglmm(hydroxyBZT ~ BZT_init + microbes, random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))
summary(MCMCglmm(hydroxyBZT ~ BZT_init, random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))#slightly worse
#
summary(MCMCglmm(hydroxyBZT ~ logOD + ducksize, random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))
summary(MCMCglmm(hydroxyBZT ~ logOD , random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))
summary(MCMCglmm(hydroxyBZT ~ 1, random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))
#
summary(MCMCglmm(pthalic_acid ~ BZT_init + salt + microbes, random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))
summary(MCMCglmm(pthalic_acid ~ BZT_init + salt , random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))
summary(MCMCglmm(pthalic_acid ~ BZT_init , random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))
summary(MCMCglmm(pthalic_acid ~ 1 , random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))
#
summary(MCMCglmm(pthalic_acid ~ logOD + ducksize, random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))
summary(MCMCglmm(pthalic_acid ~ ducksize, random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))
summary(MCMCglmm(pthalic_acid ~ 1, random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))
#


###BZTa and BZTaa separated
summary(MCMCglmm(BZTalanine ~ BZT_init + salt + microbes , random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))# dic back to similar to full model
summary(MCMCglmm(BZTalanine ~ BZT_init + salt  , random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))# dic back to similar to full model
#
summary(MCMCglmm(BZTalanine ~ logOD + ducksize, random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))
summary(MCMCglmm(BZTalanine ~ ducksize, random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))
#
summary(MCMCglmm(BZTacetylalanine ~ BZT_init + salt + microbes , random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))# dic back to similar to full model
#
summary(MCMCglmm(BZTacetylalanine ~ logOD + ducksize, random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))
summary(MCMCglmm(BZTacetylalanine ~ ducksize, random = ~ genotype, data=prdngLmix, nitt=50000,verbose=F))

			