library(MCMCglmm)


range01=function(x){
newnums=(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))
return(newnums)
}

bw <- colorRampPalette(c(rgb(1,1,1),rgb(0,0,0)))

#a remove missing data function:
getfull <- function(dat){
	whichfull <- which(sapply(1:nrow(dat), function(z) any(is.na(dat[z,]) ) )==FALSE)
	return(whichfull)
}

#shortcut for coda HPDinterval 
HPDi <- function(vect,prob){
	int <- HPDinterval(as.mcmc(vect),prob=prob)
	return(int)
}

bufferX <- function(x,p) { 
	r<- range(x,na.rm=T)
	add <- c(-1,1)*p*(r[2]-r[1])
	return(r+add)
	}	
	

std.error <- function(dat, na.rm=TRUE) {sd(dat,na.rm=na.rm)/sqrt(length(dat))}#defaults to na.rm=T

MapToWellsBZ <- function(dat,map,firstcol, sumcols, meancols){ #columns labeled "plate" "number",
		#roi is in the BZT2 files the coordinates of the shape, and "number" is the simple ROI number
	n <- nrow(map)
	outdata <- matrix(NA,nrow=n, ncol=(length(sumcols)+length(meancols)+1))
	for(i in 1:n){
		rois <- map[i,firstcol:ncol(map)] #what happens to NA values
		p <- map$plate[i]
		welldat <- dat[dat$plate==p & dat$number%in%rois,]
		welldat.sums <- colSums(welldat[,sumcols])
		welldat.means <- colMeans(welldat[,meancols],na.rm=T)
		wellstats <- c(nrow(welldat),welldat.sums,welldat.means)
		outdata[i,] <- wellstats
		}
	mappeddata <- cbind(map[,1:(firstcol-1)],outdata)
	colnames(mappeddata) <- c(colnames(map)[1:(firstcol-1)],"particles",colnames(dat)[sumcols],colnames(dat)[meancols])
	return(mappeddata)
}


#################
##Data input and processing
################

bzsinfo <- read.csv("locations_info_bzs2.csv")

#OD dat
bzs2.raw <- read.csv("AO BZS.2.ODall.csv",header=T,stringsAsFactors=F)
#first, need to add coordinates for experimental plate. 2 48-well plates were run for optical density in 1 96-well plate
#first plate in cols 1-6, next in 7-12 two samples from each well so A went to both A&B in the 96-well ,
bzs2.raw$plate <- as.numeric(sapply(1:nrow(bzs2.raw), function(z) unlist(strsplit(bzs2.raw$OD.plate[z],split="[.]"))[1] ))
bzs2.raw$pair <- sapply(1:nrow(bzs2.raw), function(z)   ifelse(bzs2.raw$OD.col[z] <=6, "I","II"))
bzs2.raw$platepair <- paste(bzs2.raw$plate,bzs2.raw$pair,sep=".")
bzs2.raw$row <- sapply(bzs2.raw$OD.row, function(z) 
	  if(z == "A" | z == "B") "A" else if(z == "C" | z == "D") "B" else if(z == "E" | z== "F") "C" else "D" )
bzs2.raw$col <- sapply(bzs2.raw$OD.col, function(z) ifelse(z<= 6, z, z-6))	  
#take the average/midpoint for each well (wells have two measures each)
fullloc <- paste(bzs2.raw$plate,bzs2.raw$pair,bzs2.raw$row,bzs2.raw$col,sep=".") 
od600 <- tapply(bzs2.raw$OD600.raw, fullloc, mean,na.rm=T)
od750 <- tapply(bzs2.raw$OD750.raw, fullloc, mean,na.rm=T)
plate <- tapply(bzs2.raw$plate, fullloc, mean,na.rm=T)
pair <-  c("I","II")[tapply(as.numeric(as.factor(bzs2.raw$pair)), fullloc, mean,na.rm=T)]
cols <- tapply(bzs2.raw$col, fullloc, mean,na.rm=T)
rows <- c("A","B","C","D")[tapply(as.numeric(as.factor(bzs2.raw$row)), fullloc, mean,na.rm=T)]
#
od_datRawu <- data.frame(plate = plate, pair = pair, row=rows, col = cols,od600 = od600, od750 = od750)
od_datRaw <- od_datRawu[order(od_datRawu$plate, od_datRawu$pair,od_datRawu$row,od_datRawu$col,decreasing=F ),] #note the tapply(), then this re-sorts so to sort by plate, then ROW then column, which matches clara's frond data sorting
#sorting still required, because otherwise it's alphabetical by the tapply levels, and plates 10-19 come before plates 2-9 & etc

#frond data
frondsRaw <- read.csv("BZS2 Frond Clara.csv",stringsAsFactors=F,header=T)

#area data, this has extensive processing necessary, continues until line ~250
endmapu <- read.csv("end_map_errorchecked.csv",header=T)
endmap <- endmapu[order( as.numeric(endmapu$plate), endmapu$row, endmapu$column  ),]
endfulldat <- read.csv("end_data_errorchecked.csv",header=T) # note column X.area is %area from imageJ

day5mapu <- read.csv("Feb18-19_map_errorchecked.csv",header=T)
day5map <- day5mapu[order( as.numeric(day5mapu$plate), day5mapu$row, day5mapu$column  ),]
day5fulldat <- read.csv("Feb18-19_data_errorchecked.csv",header=T) # note column X.area is %area from imageJ

day7mapu <- read.csv("Feb20-21map_errorchecked.csv",header=T)
day7map <- day7mapu[order( as.numeric(day7mapu$plate), day7mapu$row, day7mapu$column  ),]
day7fulldat <- read.csv("Feb20-21dat.csv",header=T) # note column X.area is %area from imageJ

startmapu <- read.csv("start map ALL.csv",header=T)
startmap <- startmapu[order( as.numeric(startmapu$plate), startmapu$row, startmapu$column  ),]
startfulldat <- read.csv("start data ALL.csv",header=T) # note column X.area is %area from imageJ

#there is a plate column, but not image name column in the map file. we therefore need to extract plate from the dat file
colors <- c("(blue)","(green)","(red)")
# END PLATES
colorrows <- lapply(colors,function(z) grep(z,endfulldat$label))
endfulldat2 <- endfulldat[-unlist(colorrows),]
endfulldat2$redraw <- endfulldat[colorrows[[3]],]$mean
endfulldat2$greenraw <- endfulldat[colorrows[[2]],]$mean
endfulldat2$blueraw <- endfulldat[colorrows[[1]],]$mean
##CHECK THAT THEY LINE UP, should be 1/4 nrow(endfulldat):
nrow(endfulldat) == 4*length(colorrows[[3]])
sum(endfulldat2$area == endfulldat[colorrows[[3]],]$area) == 0.25*nrow(endfulldat)
sum(endfulldat2$area == endfulldat[colorrows[[2]],]$area)== 0.25*nrow(endfulldat)
sum(endfulldat2$area == endfulldat[colorrows[[1]],]$area)== 0.25*nrow(endfulldat)
endfulldat2$perred <- endfulldat2$redraw/(3*endfulldat2$mean)
endfulldat2$pergreen <- endfulldat2$greenraw/(3*endfulldat2$mean)
endfulldat2$perblue <- endfulldat2$blueraw/(3*endfulldat2$mean)
#
platepix <- sapply(1:nrow(endfulldat2), function(z) strsplit(
									strsplit(as.character(endfulldat2$label[z]),"plate",fixed=T)[[1]][[2]],
										".J") [[1]][[1]]										)
endfulldat2$plate <- platepix
enddatU <- MapToWellsBZ(endfulldat2,endmap,firstcol=6,sumcols = c(5,7),meancols =c(6,8,10:12,16:18) )
enddatc <- enddatU[order( as.numeric(enddatU$plate), enddatU$row, enddatU$column  ),]

#FEB20-21, day 7
colorrows <- lapply(colors,function(z) grep(z,day7fulldat$label))
day7fulldat2 <- day7fulldat[-unlist(colorrows),]
day7fulldat2$redraw <- day7fulldat[colorrows[[3]],]$mean
day7fulldat2$greenraw <- day7fulldat[colorrows[[2]],]$mean
day7fulldat2$blueraw <- day7fulldat[colorrows[[1]],]$mean
##CHECK THAT THEY LINE UP, should be 1/4 nrow(day7fulldat):
nrow(day7fulldat) == 4*length(colorrows[[3]])
sum(day7fulldat2$area == day7fulldat[colorrows[[3]],]$area) == 0.25*nrow(day7fulldat)
sum(day7fulldat2$area == day7fulldat[colorrows[[2]],]$area)== 0.25*nrow(day7fulldat)
sum(day7fulldat2$area == day7fulldat[colorrows[[1]],]$area)== 0.25*nrow(day7fulldat)
day7fulldat2$perred <- day7fulldat2$redraw/(3*day7fulldat2$mean)
day7fulldat2$pergreen <- day7fulldat2$greenraw/(3*day7fulldat2$mean)
day7fulldat2$perblue <- day7fulldat2$blueraw/(3*day7fulldat2$mean)
#
platepix <- sapply(1:nrow(day7fulldat2), function(z) strsplit(
									strsplit(as.character(day7fulldat2$label[z]),"plate",fixed=T)[[1]][[2]],
										".J") [[1]][[1]]										)
day7fulldat2$plate <- platepix
day7datU <- MapToWellsBZ(day7fulldat2,day7map,firstcol=6,sumcols = c(4,6),meancols =c(5,7,9:11,15:17) )
day7datc <- day7datU[order( as.numeric(day7datU$plate), day7datU$row, day7datU$column  ),]

#FEB18-19, day 5
colorrows <- lapply(colors,function(z) grep(z,day5fulldat$label))
day5fulldat2 <- day5fulldat[-unlist(colorrows),]
day5fulldat2$redraw <- day5fulldat[colorrows[[3]],]$mean
day5fulldat2$greenraw <- day5fulldat[colorrows[[2]],]$mean
day5fulldat2$blueraw <- day5fulldat[colorrows[[1]],]$mean
##CHECK THAT THEY LINE UP, should be 1/4 nrow(day5fulldat):
nrow(day5fulldat) == 4*length(colorrows[[3]])
sum(day5fulldat2$area == day5fulldat[colorrows[[3]],]$area) == 0.25*nrow(day5fulldat)
sum(day5fulldat2$area == day5fulldat[colorrows[[2]],]$area)== 0.25*nrow(day5fulldat)
sum(day5fulldat2$area == day5fulldat[colorrows[[1]],]$area)== 0.25*nrow(day5fulldat)
day5fulldat2$perred <- day5fulldat2$redraw/(3*day5fulldat2$mean)
day5fulldat2$pergreen <- day5fulldat2$greenraw/(3*day5fulldat2$mean)
day5fulldat2$perblue <- day5fulldat2$blueraw/(3*day5fulldat2$mean)
#
platepix <- sapply(1:nrow(day5fulldat2), function(z) strsplit(
									strsplit(as.character(day5fulldat2$label[z]),"plate",fixed=T)[[1]][[2]],
										".J") [[1]][[1]]										)
day5fulldat2$plate <- platepix
day5datU <- MapToWellsBZ(day5fulldat2,day5map,firstcol=5,sumcols = c(4,6),meancols =c(5,7,9:11,15:17) )
day5datc <- day5datU[order( as.numeric(day5datU$plate), day5datU$row, day5datU$column  ),]

# START PLATES
colorrowsS <- lapply(colors,function(z) grep(z,startfulldat$label))
startfulldat2 <- startfulldat[-unlist(colorrowsS),]
startfulldat2$redraw <- startfulldat[colorrowsS[[3]],]$mean
startfulldat2$greenraw <- startfulldat[colorrowsS[[2]],]$mean
startfulldat2$blueraw <- startfulldat[colorrowsS[[1]],]$mean
nrow(startfulldat) == 4*length(colorrowsS[[3]])
sum(startfulldat2$area == startfulldat[colorrowsS[[3]],]$area) == 0.25*nrow(startfulldat)
sum(startfulldat2$area == startfulldat[colorrowsS[[2]],]$area)== 0.25*nrow(startfulldat)
sum(startfulldat2$area == startfulldat[colorrowsS[[1]],]$area)== 0.25*nrow(startfulldat)
startfulldat2$perred <- startfulldat2$redraw/(3*startfulldat2$mean)
startfulldat2$pergreen <- startfulldat2$greenraw/(3*startfulldat2$mean)
startfulldat2$perblue <- startfulldat2$blueraw/(3*startfulldat2$mean)
platepixS <- sapply(1:nrow(startfulldat2), function(z) strsplit(
									strsplit(as.character(startfulldat2$label[z]),"plate",fixed=T)[[1]][[2]],
										".J") [[1]][[1]]										)
startfulldat2$plate <- platepixS
startdatU <- MapToWellsBZ(startfulldat2,startmap,firstcol=5,sumcols = c(4,6),meancols =c(5,7,9:11,15:17) )
startdatc <- startdatU[order( as.numeric(startdatU$plate), startdatU$row, startdatU$column  ),]

#### convert to square mm of tissue at same time
#experiment was set up in 2 times steps, plates 1-20; then 21-50
standarddim <- read.csv("wellsizeBZS2photos.csv")
#conversion factors by day
pxpmm <- standarddim$pixlength / standarddim$mm.length #standard width between centers of wells, and measures in representative images
dayofphoto <- c(rep(1,times=960),rep(2,times=1440))
startdat <- startdatc; day5dat <- day5datc; day7dat<-day7datc; enddat <-enddatc
startdat$area <- startdatc$area/(pxpmm[dayofphoto])^2 
day5dat$area <- day5datc$area/(pxpmm[dayofphoto+2])^2 
day7dat$area <- day7datc$area/(pxpmm[dayofphoto+4])^2 
enddat$area <-  enddatc$area/(pxpmm[dayofphoto+6])^2   

arabicpairod <- ifelse(od_datRaw$pair=="I",1,2)
arabicppod <- as.numeric (paste(od_datRaw$plate,arabicpairod,sep="."))
# same_sampleid <- paste(arabicplateod,paste(od_datRaw$row,od_datRaw$col,sep=""),sep=".")
arabicplatefr <- sapply(1:nrow(frondsRaw), function(z) strsplit(frondsRaw$plate[z],"[.]")[[1]][1])
pairfr <-  sapply(1:nrow(frondsRaw), function(z) strsplit(frondsRaw$plate[z],"[.]")[[1]][2])
arabicpairfr <-  ifelse(pairfr=="I",1,2)
arabicppfr <- paste(arabicplatefr,arabicpairfr,sep=".")

#checks that data align as hoped across datasets
sum(arabicppfr == arabicppod)
sum(arabicppfr == enddat$plate)
sum(arabicppfr == startdat$plate)
sum(frondsRaw$column == od_datRaw$col)
sum(frondsRaw$column == enddat$column)
sum(frondsRaw$row == od_datRaw$row)
sum(frondsRaw$row == enddat$row)

frondsRaw$arabicpp <- arabicppfr
od_datRaw$arabicpp <- arabicppod

#remove error rows
fronds <- frondsRaw[frondsRaw$unfixed.error != "ERROR",]
od_dat <- od_datRaw[frondsRaw$unfixed.error != "ERROR",]
endsqm <- enddat[frondsRaw$unfixed.error != "ERROR",]
startpix <- startdat[frondsRaw$unfixed.error != "ERROR",]
midsqm <- day5dat[frondsRaw$unfixed.error != "ERROR",]
midsqm2 <- day7dat[frondsRaw$unfixed.error != "ERROR",]

#
colnames(startpix) <- paste("start",colnames(startdat),sep="_")
colnames(midsqm) <- paste("mid",colnames(midsqm),sep="_")
colnames(midsqm2) <- paste("mid2",colnames(midsqm2),sep="_")

#fix issue with fronds file:
#the experiment prep notes in photographs show col6 of platepair 23 should be MNT, col4 of plaitepaire 23 should be SHB
#in the frond data this is swapped (an error was made in swapping labels, there had been 2 cols of SHB in original design, one column was swapped to MNT but the wrong column was changed in the file)
fronds$genotype[od_dat$plate==23 & fronds$col==6] <- "Mnt"
fronds$genotype[od_dat$plate==23 & fronds$col==4] <- "SHB"

biodat <- data.frame(cbind(fronds[1:14],od_dat$od600,od_dat$od750,endsqm[,c(2,6:16)],startpix[,c(1,5:15)], midsqm[,c(1,5:15)],  midsqm2[,c(1,5:15)]))
biodat$microbeyn <- ifelse(biodat$microbe=="+","Y","N")
biodatlocs <- sapply(biodat$genotype, function(x) bzsinfo$km.cityC[which(as.character(bzsinfo$bzs2)==as.character(x))])
biodat$loc <- biodatlocs
biodat$rddns <- sapply(biodat$genotype, function(x) bzsinfo$roaddens[which(as.character(bzsinfo$bzs2)==as.character(x))])
biodat$deltasqmm <- biodat$area - biodat$start_area #note these are all in units of square millimeters of area. but were measured from pixels
biodat$dsqmmid <- biodat$mid_area - biodat$start_area
biodat$dsqm511 <- biodat$area - biodat$mid_area
biodat$dsqmmid2 <- biodat$mid2_area - biodat$start_area
biodat$dsqm711 <- biodat$area - biodat$mid2_area
write.csv(biodat,"biodat.csv",row.names=F)

#################
##ANALYSIS
################

##CHANGE IN PIXEL AREA, from start to end., model selection
summary(MCMCglmm(deltasqmm~rddns +  microbe + BZT + Salt + rddns:microbe + rddns:BZT + rddns:Salt + BZT:microbe + BZT:Salt + microbe:Salt + 
						rddns:microbe:BZT + rddns:BZT:Salt + rddns:microbe:Salt + microbe:Salt:BZT + microbe:Salt:BZT:rddns, data=biodat,verbose=F, nitt=50000,thin=10,burnin=1000,pr=T))
summary( MCMCglmm(deltasqmm~rddns +  microbe + BZT + Salt + rddns:microbe + rddns:BZT + rddns:Salt + BZT:microbe + BZT:Salt + microbe:Salt + 
						rddns:microbe:BZT + rddns:BZT:Salt + rddns:microbe:Salt + microbe:Salt:BZT, data=biodat,verbose=F, nitt=50000,thin=10,burnin=1000,pr=T))
summary( MCMCglmm(deltasqmm~rddns +  microbe + BZT + Salt + rddns:microbe + rddns:BZT + rddns:Salt + BZT:microbe + BZT:Salt + microbe:Salt + 
						rddns:microbe:BZT + rddns:microbe:Salt + microbe:Salt:BZT, data=biodat,verbose=F, nitt=50000,thin=10,burnin=1000,pr=T))
summary(MCMCglmm(deltasqmm~rddns +  microbe + BZT + Salt + rddns:microbe + rddns:BZT + rddns:Salt + BZT:microbe + BZT:Salt + microbe:Salt + 
						 rddns:microbe:Salt + microbe:Salt:BZT, data=biodat,verbose=F, nitt=50000,thin=10,burnin=1000,pr=T))
summary(MCMCglmm(deltasqmm~rddns +  microbe + BZT + Salt + rddns:microbe + rddns:BZT + rddns:Salt + BZT:microbe + BZT:Salt + microbe:Salt + 
						 microbe:Salt:BZT, data=biodat,verbose=F, nitt=50000,thin=10,burnin=1000,pr=T))
summary( MCMCglmm(deltasqmm~rddns +  microbe + BZT + Salt + rddns:microbe + rddns:BZT + rddns:Salt + BZT:microbe + BZT:Salt + microbe:Salt, data=biodat,verbose=F, nitt=50000,thin=10,burnin=1000,pr=T))
#no change in DIC from r4 to r5...
summary( MCMCglmm(deltasqmm~rddns +  microbe + BZT + Salt + rddns:microbe + rddns:BZT  + BZT:microbe + BZT:Salt + microbe:Salt, data=biodat,verbose=F, nitt=50000,thin=10,burnin=1000,pr=T))
#this one is better than r4 or r5; but pMCMC close on microbe:BZT and rddens:BZT in 4 runs, 2 each worst, try both
summary( MCMCglmm(deltasqmm~rddns +  microbe + BZT + Salt + rddns:microbe  + BZT:microbe + BZT:Salt + microbe:Salt, data=biodat,verbose=F, nitt=50000,thin=10,burnin=1000,pr=T))
summary( MCMCglmm(deltasqmm~rddns +  microbe + BZT + Salt + rddns:microbe + rddns:BZT  + BZT:Salt + microbe:Salt, data=biodat,verbose=F, nitt=50000,thin=10,burnin=1000,pr=T))
#each better than prev, and both find the other term rm next
summary(  MCMCglmm(deltasqmm~rddns +  microbe + BZT + Salt + rddns:microbe  + BZT:Salt + microbe:Salt, data=biodat,verbose=F, nitt=50000,thin=10,burnin=1000,pr=T))
summary(  MCMCglmm(deltasqmm~rddns +  microbe + BZT + Salt + BZT:Salt + microbe:Salt, data=biodat,verbose=F, nitt=50000,thin=10,burnin=1000,pr=T))
# DIC removing microbe:Salt consistently worse, even though is always n.s.
#STOP; urban duckweed grow more, duckweed grow more in the local treatment, grow less with salt or BZT both both together not worse than salt alone or maybe slightly better
dsqmRbest <- MCMCglmm(deltasqmm~rddns +  microbe + BZT + Salt + BZT:Salt + microbe:Salt, 
						data=biodat,verbose=F, nitt=100000,thin=10,burnin=10000,pr=T)
summary(dsqmRbest)
rd.s <- seq(from=min(biodat$rddns),to=max(biodat$rddns),length.out=1000)
dpxsolrd <- dsqmRbest$Sol
dpxng.mnrd <- sapply(1:length(rd.s), function(z) mean(dpxsolrd[,1]+ dpxsolrd[,2]*rd.s[z] + dpxsolrd[,3]*0.5 + dpxsolrd[,4]*0.5+ dpxsolrd[,5]*0.5 + dpxsolrd[,6]*0.5 + dpxsolrd[,7]*0.5) )
dpxng.hpdird <- sapply(1:length(rd.s), function(z) HPDi(dpxsolrd[,1]+ dpxsolrd[,2]*rd.s[z] + dpxsolrd[,3]*0.5 + dpxsolrd[,4]*0.5+ dpxsolrd[,5]*0.5 + dpxsolrd[,6]*0.5 + dpxsolrd[,7]*0.5, 0.95) )
dsqmrdmn <- tapply(biodat$deltasqmm,biodat$rddns,mean)
rdrd <- tapply(biodat$rddns,biodat$rddns,mean)
dsqm.mns <- tapply(biodat$deltasqmm, paste(biodat$microbe,biodat$Salt,biodat$BZT),mean)
dsqm.ses <- tapply(biodat$deltasqmm, paste(biodat$microbe,biodat$Salt,biodat$BZT),std.error)


#OPTICAL DENSITY
#testing normality assumption, and taking log to address issue
shapiro.test((biodat$od_dat.od600))
shapiro.test((biodat$od_dat.od750))
shapiro.test(log(biodat$od_dat.od600+min(biodat$od_dat.od600[biodat$od_dat.od600>0],na.rm=T)))
shapiro.test(log(biodat$od_dat.od750+min(biodat$od_dat.od750[biodat$od_dat.od750>0],na.rm=T)))
biodat$lnod600 <- log(biodat$od_dat.od600+ min(biodat$od_dat.od600[biodat$od_dat.od600>0],na.rm=T ) )
biodat$lnod750 <- log(biodat$od_dat.od750+ min(biodat$od_dat.od750[biodat$od_dat.od750>0],na.rm=T ) )

bioMdat <- biodat[biodat$microbe=="+",]

inocModng <- MCMCglmm(lnod600 ~microbe,data=biodat,verbose=F, nitt=50000,thin=10,burnin=1000)
summary(inocModng) 
#therefore + and - are more "disrupted" and "local"
inocmns <- tapply(biodat$lnod600,biodat$microbe,mean,na.rm=T)
inocses <- tapply(biodat$lnod600,biodat$microbe,std.error,na.rm=T)
pdf("BZT2_no_inoc_effect.pdf",height=3,width=3)
par(mar=c(4,4,1,1))
plot(inocmns~c(1,2),xlim=c(0.5,2.5),xaxt="n",ylim=bufferX(c(inocmns+inocses,inocmns-inocses),0.1) ,cex=2,
	ylab="",xlab=""	)
	arrows(1:2,inocmns-inocses,y1=inocmns+inocses,length=0,lwd=2)
	axis(side=1,at=c(1,2),labels=c("disrupted","local"))
	mtext("microbe treatment",side=1,line=2.2)
	mtext("microbe growth ln(OD)",side=2,line=2.5)
dev.off()

# road dens, model selection
summary(MCMCglmm(lnod600 ~rddns  + BZT + Salt +  rddns:BZT + rddns:Salt + BZT:Salt +
						 rddns:BZT:Salt , data=bioMdat,verbose=F, nitt=50000,thin=10,burnin=1000))
summary(MCMCglmm(lnod600 ~rddns  + BZT + Salt +  rddns:BZT + rddns:Salt + BZT:Salt , data=bioMdat,verbose=F, nitt=50000,thin=10,burnin=1000))
#two terms close, rm 1 at time
summary(MCMCglmm(lnod600 ~rddns  + BZT + Salt +  rddns:BZT + rddns:Salt  , data=bioMdat,verbose=F, nitt=50000,thin=10,burnin=1000))
summary(MCMCglmm(lnod600 ~rddns  + BZT + Salt +  rddns:BZT + BZT:Salt , data=bioMdat,verbose=F, nitt=50000,thin=10,burnin=1000))
#each suggests removal of other in continuing to simpler model
summary(MCMCglmm(lnod600 ~rddns  + BZT + Salt +  rddns:BZT , data=bioMdat,verbose=F, nitt=50000,thin=10,burnin=1000))
summary(MCMCglmm(lnod600 ~rddns  + BZT + Salt  , data=bioMdat,verbose=F, nitt=50000,thin=10,burnin=1000))
summary(MCMCglmm(lnod600 ~rddns  + BZT   , data=bioMdat,verbose=F, nitt=50000,thin=10,burnin=1000))
summary(MCMCglmm(lnod600 ~rddns    , data=bioMdat,verbose=F, nitt=50000,thin=10,burnin=1000))
#stop, refit
mbod600ngrdr6 <- MCMCglmm(lnod600 ~rddns    , data=bioMdat,verbose=F, nitt=100000,thin=10,burnin=10000)
summary(mbod600ngrdr6)
odngRd.mns <- sapply(rd.s,function(z)	mean(mbod600ngrdr6$Sol[,1] + mbod600ngrdr6$Sol[,2]*z))
odngRd.hpdis <- sapply(rd.s,function(z)	HPDi(mbod600ngrdr6$Sol[,1] + mbod600ngrdr6$Sol[,2]*z,0.95))
odrdmns <- tapply(bioMdat$lnod600,bioMdat$rddns,mean)

pdf("dsqmandOD_sigcats_rddens.pdf",width=7,height=2.5)#9.5,3
par(mfrow=c(1,3))
par(mar=c(3.5,4,1.5,0))
par(oma=c(0,0,0,0))
plot(dsqm.mns~c(2,2.1,2.2,2.3,3,3.1,3.2,3.3),xaxt="n",xlab="",ylab="", cex=2, xlim=c(1.8,3.5),
		col = c(rgb(0,0,0),rgb(0,0,0.75),rgb(0.75,0,0),rgb(0.75,0,0.75)), pch=c(0:2,6),
		ylim=c(10,21))
	arrows(c(2,2.1,2.2,2.3,3,3.1,3.2,3.3),dsqm.mns-dsqm.ses,y1=dsqm.mns+dsqm.ses,length=0,
		col = c(rgb(0,0,0),rgb(0,0,0.75),rgb(0.75,0,0),rgb(0.75,0,0.75)), lwd=2)
	arrows(c(2,2.1,2.2,2.3), dsqm.mns[1:4], x1=c(3,3.1,3.2,3.3),dsqm.mns[5:8], length=0, lwd=2,
		col = c(rgb(0,0,0,alpha=0.5),rgb(0,0,0.75,alpha=0.5),rgb(0.75,0,0,alpha=0.5),rgb(0.75,0,0.75,alpha=0.5)),lty=2)
	axis(side=1, at=c(2.15,3.15),labels=c("disrupted","local"),cex.axis=1.25)
	mtext("microbes",side=1,line=2.2)
	mtext(bquote("plant growth, mm"^2),side=2,line=2)
	legend(1.7,22,c("Neither","1 mg/L BZT","0.8 g/L Salt","Both"), bty="n", pch=c(0:2,6),
		 col = c(rgb(0,0,0),rgb(0,0,0.75),rgb(0.75,0,0),rgb(0.75,0,0.75)),cex=1.25)
	mtext("a.",side=3,adj=-0.25,line=0.3)
plot(dsqmrdmn~rdrd,pch=NA,ylab="",xlab="",ylim=c(-3,40))
	polygon(c(rd.s,rev(rd.s)), y=c(dpxng.hpdird[1,],rev(dpxng.hpdird[2,])), col=rgb(0,0,0,alpha=0.5), border=NA)
	lines(dpxng.mnrd~rd.s,lwd=2)	
	points(dsqmrdmn~rdrd,pch=16,cex=1.5)
	mtext(bquote("plant growth, mm"^2),side=2,line=2)
	mtext("road density",side=1,line=2.2) #length per area. - km / km^2
	mtext("b.",side=3,adj=-0.25,line=0.3)
plot(odrdmns~rdrd,pch=16,cex=1.5,xlab="",ylab="",ylim=c(-2.5,-1.6))
	polygon(c(rd.s,rev(rd.s)), y=c(odngRd.hpdis[1,],rev(odngRd.hpdis[2,])), col=rgb(0,0,0,alpha=0.5), border=NA)
	lines(odngRd.mns~rd.s,lwd=2)
	mtext("microbe growth, ln(OD)",side=2,line=2.5)#
	mtext("road density",side=1,line=2.2)
	mtext("c.",side=3,adj=-0.25,line=0.3)
dev.off()

# including disrupted microbes treatments
summary(MCMCglmm(lnod600 ~rddns +  microbe + BZT + Salt + rddns:microbe + rddns:BZT + 
						rddns:Salt + BZT:microbe + BZT:Salt + microbe:Salt + rddns:microbe:BZT +
						rddns:BZT:Salt + rddns:microbe:Salt + microbe:Salt:BZT + rddns:microbe:BZT:Salt, 
						data=biodat,verbose=F, nitt=50000,thin=10,burnin=1000))
summary(MCMCglmm(lnod600 ~rddns +  microbe + BZT + Salt + rddns:microbe + rddns:BZT + 
						rddns:Salt + BZT:microbe + BZT:Salt + microbe:Salt + rddns:microbe:BZT +
						rddns:BZT:Salt + rddns:microbe:Salt + microbe:Salt:BZT , 
						data=biodat,verbose=F, nitt=50000,thin=10,burnin=1000))
summary(MCMCglmm(lnod600 ~rddns +  microbe + BZT + Salt + rddns:microbe + rddns:BZT + 
						rddns:Salt + BZT:microbe + BZT:Salt + microbe:Salt + rddns:microbe:BZT +
						rddns:BZT:Salt + microbe:Salt:BZT , 
						data=biodat,verbose=F, nitt=50000,thin=10,burnin=1000))
summary(MCMCglmm(lnod600 ~rddns +  microbe + BZT + Salt + rddns:microbe + rddns:BZT + 
						rddns:Salt + BZT:microbe + BZT:Salt + microbe:Salt + rddns:microbe:BZT +
						rddns:BZT:Salt, data=biodat,verbose=F, nitt=50000,thin=10,burnin=1000))
summary(MCMCglmm(lnod600 ~rddns +  microbe + BZT + Salt + rddns:microbe + rddns:BZT + 
						rddns:Salt + BZT:microbe + BZT:Salt + microbe:Salt + rddns:microbe:BZT, data=biodat,verbose=F, nitt=50000,thin=10,burnin=1000))
summary(MCMCglmm(lnod600 ~rddns +  microbe + BZT + Salt + rddns:microbe + rddns:BZT + 
						rddns:Salt + BZT:microbe + BZT:Salt + microbe:Salt, data=biodat,verbose=F, nitt=50000,thin=10,burnin=1000))
summary(MCMCglmm(lnod600 ~rddns +  microbe + BZT + Salt + rddns:microbe +
						rddns:Salt + BZT:microbe + BZT:Salt + microbe:Salt, data=biodat,verbose=F, nitt=50000,thin=10,burnin=1000))
summary(MCMCglmm(lnod600 ~rddns +  microbe + BZT + Salt + rddns:microbe +
						 BZT:microbe + BZT:Salt + microbe:Salt, data=biodat,verbose=F, nitt=50000,thin=10,burnin=1000))
summary(MCMCglmm(lnod600 ~rddns +  microbe + BZT + Salt + rddns:microbe + BZT:microbe + microbe:Salt, data=biodat,verbose=F, nitt=50000,thin=10,burnin=1000))
summary(MCMCglmm(lnod600 ~rddns +  microbe + BZT + Salt + rddns:microbe + microbe:Salt, data=biodat,verbose=F, nitt=50000,thin=10,burnin=1000))
summary(MCMCglmm(lnod600 ~rddns +  microbe + Salt + rddns:microbe + microbe:Salt, data=biodat,verbose=F, nitt=50000,thin=10,burnin=1000))
#stop, refit
mbod600fdr10 <- MCMCglmm(lnod600 ~rddns +  microbe + Salt + rddns:microbe + microbe:Salt, data=biodat,verbose=F, nitt=100000,thin=10,burnin=10000)
summary(mbod600fdr10)

###FROND GROWTH
frondsbytime <- data.frame(fronds = c(biodat$X0.FN,biodat$X2.FN,biodat$X4.FN,biodat$X6.FN,biodat$X8.FN,biodat$X10.FN),
							id = as.character(rep(1:nrow(biodat),times=6)), days= rep(c(0,3,5,7,9,11),each=nrow(biodat)),
							BZT = rep(ifelse(biodat$BZT=="B",1,0),times=6), Salt = rep(ifelse(biodat$Salt=="S",1,0), times=6), 
							microbe = rep(ifelse(biodat$microbe=="+",1,0),times=6), loc = rep(biodat$loc,times=6),
							rddns = rep(biodat$rddns,times=6), genotype = rep(biodat$genotype,times=6) )
 mbgrTIMEidintONLYRD <- MCMCglmm(fronds~ days + days:rddns +  days:microbe + days:BZT + days:Salt +
					 days:rddns:microbe + days:rddns:BZT + days:rddns:Salt + days:BZT:microbe + days:BZT:Salt + days:microbe:Salt + 
						days:rddns:microbe:BZT + days:rddns:BZT:Salt + days:rddns:microbe:Salt + days:microbe:Salt:BZT + days:rddns:microbe:BZT:Salt , 
				random = ~ id, 
				data=frondsbytime,verbose=F, nitt=100000,thin=10,burnin=1000,pr=T)
 mbgrTIMEidintONLYRDr1 <- MCMCglmm(fronds~ days + days:rddns +  days:microbe + days:BZT + days:Salt +
					 days:rddns:microbe + days:rddns:BZT + days:rddns:Salt + days:BZT:microbe + days:BZT:Salt + days:microbe:Salt + 
						days:rddns:microbe:BZT + days:rddns:BZT:Salt + days:rddns:microbe:Salt + days:microbe:Salt:BZT , 
				random = ~ id, 
				data=frondsbytime,verbose=F, nitt=100000,thin=10,burnin=1000,pr=T)
 mbgrTIMEidintONLYRDr2 <- MCMCglmm(fronds~ days + days:rddns +  days:microbe + days:BZT + days:Salt +
					 days:rddns:microbe + days:rddns:BZT + days:rddns:Salt + days:BZT:microbe + days:BZT:Salt + days:microbe:Salt + 
						days:rddns:BZT:Salt + days:rddns:microbe:Salt + days:microbe:Salt:BZT , 
				random = ~ id, 
				data=frondsbytime,verbose=F, nitt=100000,thin=10,burnin=1000,pr=T)
#stop			
#already fit with enough iterations
summary(mbgrTIMEidintONLYRDr2)
FgrowSolNGRd <- mbgrTIMEidintONLYRDr2$Sol
forfrondtimeNGRd <- function(b, S, m) {
				mnslp <- sapply(rd.s , function(z) mean(FgrowSolNGRd[,2] +  FgrowSolNGRd[,3]*z +  FgrowSolNGRd[,4]*m +  FgrowSolNGRd[,5]*b + FgrowSolNGRd[,6]*S +
													    FgrowSolNGRd[,7]*z*m +  FgrowSolNGRd[,8]*z*b +  FgrowSolNGRd[,9]*z*S +  FgrowSolNGRd[,10]*b*m +
													   FgrowSolNGRd[,11]*b*S +  FgrowSolNGRd[,12]*m*S +  FgrowSolNGRd[,13]*z*b*S + FgrowSolNGRd[,14]*z*m*S +  FgrowSolNGRd[,15]*m*b*S) )
				hpdislp <- sapply(rd.s , function(z) HPDi(FgrowSolNGRd[,2] +  FgrowSolNGRd[,3]*z +  FgrowSolNGRd[,4]*m +  FgrowSolNGRd[,5]*b + FgrowSolNGRd[,6]*S +
													    FgrowSolNGRd[,7]*z*m +  FgrowSolNGRd[,8]*z*b +  FgrowSolNGRd[,9]*z*S +  FgrowSolNGRd[,10]*b*m +
													   FgrowSolNGRd[,11]*b*S +  FgrowSolNGRd[,12]*m*S +  FgrowSolNGRd[,13]*z*b*S + FgrowSolNGRd[,14]*z*m*S +  FgrowSolNGRd[,15]*m*b*S, 0.95) )
													   return(list(mnslp,hpdislp))
}
slopengr.mnBSM	<- forfrondtimeNGRd(1,1,1)											   
slopengr.mnBSx	<- forfrondtimeNGRd(1,1,0)											   
slopengr.mnBxM	<- forfrondtimeNGRd(1,0,1)											   
slopengr.mnBxx	<- forfrondtimeNGRd(1,0,0)											   
slopengr.mnxSM	<- forfrondtimeNGRd(0,1,1)											   
slopengr.mnxSx	<- forfrondtimeNGRd(0,1,0)											   
slopengr.mnxxM	<- forfrondtimeNGRd(0,0,1)											   
slopengr.mnxxx	<- forfrondtimeNGRd(0,0,0)											   
slopes <- sapply(1:nrow(biodat), function(z) lm(unlist(biodat[z,9:14])~as.vector(c(0,2,4,6,8,10)))$coef[2] )
genoslptrt <- t(sapply(1:length(unique(biodat$genotype)), function(geno) tapply(
		slopes[ biodat$genotype == sort(unique(biodat$genotype))[geno] ],paste( biodat$BZT,biodat$Salt,biodat$microbe )[biodat$genotype == sort(unique(biodat$genotype))[geno] ],mean )
	))

pdf("frondsperday_nogenoRd.pdf",height=4,width=4)
par(mfrow=c(2,2))
par(mar=c(2,2,0,1))
par(oma=c(3,3,2.5,1))
plot(genoslptrt[,1]~bzsinfo$roaddens,pch=1,ylab="Growth rate",xlab="",ylim=c(-0.25,1.2))
	polygon(c(rd.s,rev(rd.s)), y=c(slopengr.mnxxx[[2]][1,],rev(slopengr.mnxxx[[2]][2,])), col=rgb(0,0,0,alpha=0.5), border=NA)
	polygon(c(rd.s,rev(rd.s)), y=c(slopengr.mnxSx[[2]][1,],rev(slopengr.mnxSx[[2]][2,])), col=rgb(1,0,0,alpha=0.5), border=NA)
	points(genoslptrt[,3]~bzsinfo$roaddens,pch=1,col=rgb(1,0,0))
 lines(slopengr.mnxxx[[1]]~rd.s)
 lines(slopengr.mnxSx[[1]]~rd.s,col=rgb(1,0,0))
 mtext("Fronds/day",side=2,line=3,adj=-1.5)
 mtext("0 mg/L BZT",side=2,line=2)
 mtext("disrupted",side=3,line=1.5)
 mtext("microbes",side=3,line=0.5)
plot(genoslptrt[,2]~bzsinfo$roaddens,pch=1,col=rgb(0,0,0),ylab="",xlab="",ylim=c(-0.25,1.2))
	polygon(c(rd.s,rev(rd.s)), y=c(slopengr.mnxxM[[2]][1,],rev(slopengr.mnxxM[[2]][2,])), col=rgb(0,0,0,alpha=0.5), border=NA)
	polygon(c(rd.s,rev(rd.s)), y=c(slopengr.mnxSM[[2]][1,],rev(slopengr.mnxSM[[2]][2,])), col=rgb(1,0,0,alpha=0.5), border=NA)
	points(genoslptrt[,4]~bzsinfo$roaddens,pch=1,col=rgb(1,0,0))
 lines(slopengr.mnxxM[[1]]~rd.s)
 lines(slopengr.mnxSM[[1]]~rd.s,col=rgb(1,0,0))
 mtext("local",side=3,line=1.5)
 mtext("microbes",side=3,line=0.5)
plot(genoslptrt[,5]~bzsinfo$roaddens,pch=1,col=rgb(0,0,0),ylab="Growth rate",xlab="",ylim=c(-0.25,1.2))
	polygon(c(rd.s,rev(rd.s)), y=c(slopengr.mnBxx[[2]][1,],rev(slopengr.mnBxx[[2]][2,])), col=rgb(0,0,0,alpha=0.5), border=NA)
	polygon(c(rd.s,rev(rd.s)), y=c(slopengr.mnBSx[[2]][1,],rev(slopengr.mnBSx[[2]][2,])), col=rgb(1,0,0,alpha=0.5), border=NA)
	points(genoslptrt[,7]~bzsinfo$roaddens,pch=1,col=rgb(1,0,0))
 mtext("1 mg/L BZT",side=2,line=2)
 lines(slopengr.mnBxx[[1]]~rd.s,col=rgb(0,0,0))
 lines(slopengr.mnBSx[[1]]~rd.s,col=rgb(1,0,0))
plot(genoslptrt[,6]~bzsinfo$roaddens,pch=1,col=rgb(0,0,0),ylab="",xlab="",ylim=c(-0.25,1.2))
	polygon(c(rd.s,rev(rd.s)), y=c(slopengr.mnBxM[[2]][1,],rev(slopengr.mnBxM[[2]][2,])), col=rgb(0,0,0,alpha=0.5), border=NA)
	polygon(c(rd.s,rev(rd.s)), y=c(slopengr.mnBSM[[2]][1,],rev(slopengr.mnBSM[[2]][2,])), col=rgb(1,0,0,alpha=0.5), border=NA)
	points(genoslptrt[,8]~bzsinfo$roaddens,pch=1,col=rgb(1,0,0))
 lines(slopengr.mnBxM[[1]]~rd.s,col=rgb(0,0,0))
 lines(slopengr.mnBSM[[1]]~rd.s,col=rgb(1,0,0))
 mtext("road density",side=1,line=2,adj=-3)
 legend(0,1.35,c("0 g/L NaCl","0.8 g/L NaCl"),fill=c(rgb(0,0,0),rgb(1,0,0) ),bty="n")
dev.off()
slopeng.mns <- tapply(slopes,paste( biodat$BZT,biodat$Salt,biodat$microbe ),mean)
slopeng.ses <- tapply(slopes,paste( biodat$BZT,biodat$Salt,biodat$microbe ),std.error)
pdf("frondsperday_nogenoRd_append.pdf",height=3,width=4)
par(mar=c(4,4,1,1))
plot(slopeng.mns ~ c(1,2,1,2,1,2,1,2), col = rgb( c(0,0,1,1,0,0,1,1) ,0, c(0,0,0,0,1,1,1,1) ), 
	xlim=c(0.5,4.5),xaxt="n",xlab="",ylab="", cex=2,ylim=c(0.18,0.28) )
 mtext("Fronds/day",side=2,line=2.5)
	axis(side=1, at=c(1,2),labels=c("disrupted","local"))
	arrows( c(1,2,1,2,1,2,1,2),slopeng.mns-slopeng.ses, y1= slopeng.mns + slopeng.ses,length=0,
	 col = rgb( c(0,0,1,1,0,0,1,1) ,0, c(0,0,0,0,1,1,1,1) ) )
	arrows( c(1,1,1,1), slopeng.mns[c(1,3,5,7)], x1=c(2,2,2,2), y1= slopeng.mns[c(2,4,6,8)],length=0,
	 col = rgb( c(0,1,0,1) ,0, c(0,0,1,1),alpha=0.5 ),lwd=2,lty=2 )
	legend(2.25,0.28,c("Neither","1 mg/L BZT","0.8 g/L Salt","Both"), bty="n",
		 fill = c(rgb(0,0,0),rgb(0,0,0.75),rgb(0.75,0,0),rgb(0.75,0,0.75)))
dev.off()
