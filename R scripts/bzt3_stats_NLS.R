#import libraries
library(readxl)
library(MCMCglmm) # substituting lm() in the place of the one instance of MCMCglmm() provides identical results
#read data
bzt = read_excel("bzt3_concs_original_full.xlsx")
bzt = as.data.frame(bzt)
bzt$Treatment <- factor(bzt$Treatment, levels = c("Abiotic Control",
                                                  "Microbes",
                                                  "Duckweed + Microbes + Salt",
                                                  "Duckweed + Microbes",
                                                  "Algae",
                                                  "Duckweed",
                                                  "Duckweed + Algae"))

nls_Algae <- nls(Concentration ~ (A0 - A_inf)*exp(-k*Day)+A_inf,
              data = bzt[bzt$Treatment=='Algae',],
              start = list(A0 = 5, A_inf = 0.6, k=0.1),
              lower=c(0,0,0), upper=c(10, 5, 1),
              algorithm="port"
              )

nls_Duckweed <- nls(Concentration ~ (A0 - A_inf)*exp(-k*Day)+A_inf,
                 data = bzt[bzt$Treatment=='Duckweed',],
                 start = list(A0 = 5, A_inf = 0.6, k=0.1),
                 lower=c(0,0,0), upper=c(10, 5, 1),
                 algorithm="port")

nls_Duckweed_Algae <- nls(Concentration ~ (A0 - A_inf)*exp(-k*Day)+A_inf,
                    data = bzt[bzt$Treatment=='Duckweed + Algae',],
                    start = list(A0 = 5, A_inf = 0.6, k=0.1),
                    lower=c(0,0,0), upper=c(10, 5, 1),
                    algorithm="port")

nls_Duckweed_Microbes <- nls(Concentration ~ (A0 - A_inf)*exp(-k*Day)+A_inf,
                          data = bzt[bzt$Treatment=='Duckweed + Microbes',],
                          start = list(A0 = 5, A_inf = 0.6, k=0.1),
                          lower=c(0,0,0), upper=c(10, 5, 1),
                          algorithm="port")

nls_Duckweed_Microbes_Salt <- nls(Concentration ~ (A0 - A_inf)*exp(-k*Day)+A_inf,
                             data = bzt[bzt$Treatment=='Duckweed + Microbes + Salt',],
                             start = list(A0 = 5, A_inf = 0.6, k=0.1),
                             lower=c(0,0,0), upper=c(10, 5, 1),
                             algorithm="port")

nls_Microbes <- nls(Concentration ~ (A0 - A_inf)*exp(-k*Day)+A_inf,
                                  data = bzt[bzt$Treatment=='Microbes',],
                                  start = list(A0 = 5, A_inf = 0.6, k=0.1),
                                  lower=c(0,0,0), upper=c(10, 5, 1),
                                  algorithm="port")

nls_Abiotic_Control <- nls(Concentration ~ (A0 - A_inf)*exp(-k*Day)+A_inf,
                    data = bzt[bzt$Treatment=='Abiotic Control',],
                    start = list(A0 = 5, A_inf = 0.6, k=0.1),
                    lower=c(0,0,0), upper=c(10, 5, 1),
                    algorithm="port")

k_values_nls <- data.frame("Treatment" = c("Abiotic Control","Microbes","Duckweed + Microbes + Salt",
                                       "Duckweed + Microbes","Algae","Duckweed","Duckweed + Algae"), 
                       "k_value" = c(as.vector(coefficients(nls_Abiotic_Control)[3]),
                                     as.vector(coefficients(nls_Microbes)[3]),
                                     as.vector(coefficients(nls_Duckweed_Microbes_Salt)[3]),
                                     as.vector(coefficients(nls_Duckweed_Microbes)[3]),
                                     as.vector(coefficients(nls_Algae)[3]),
                                     as.vector(coefficients(nls_Duckweed)[3]),
                                     as.vector(coefficients(nls_Duckweed_Algae)[3])),
                       "lowerCI" = c(as.vector(confint(nls_Abiotic_Control)[3,1]),
                                     as.vector(confint(nls_Microbes)[3,1]),
                                     as.vector(confint(nls_Duckweed_Microbes_Salt)[3,1]),
                                     as.vector(confint(nls_Duckweed_Microbes)[3,1]),
                                     as.vector(confint(nls_Algae)[3,1]),
                                     as.vector(confint(nls_Duckweed)[3,1]),
                                     as.vector(confint(nls_Duckweed_Algae)[3,1])),
                       "upperCI" = c(as.vector(confint(nls_Abiotic_Control)[3,2]),
                                     as.vector(confint(nls_Microbes)[3,2]),
                                     as.vector(confint(nls_Duckweed_Microbes_Salt)[3,2]),
                                     as.vector(confint(nls_Duckweed_Microbes)[3,2]),
                                     as.vector(confint(nls_Algae)[3,2]),
                                     as.vector(confint(nls_Duckweed)[3,2]),
                                     as.vector(confint(nls_Duckweed_Algae)[3,2])))
k_values_nls$Treatment <- factor(k_values_nls$Treatment, levels = c("Abiotic Control",
                                                            "Microbes",
                                                            "Duckweed + Microbes",
                                                            "Duckweed + Microbes + Salt",
                                                            "Algae",
                                                            "Duckweed",
                                                            "Duckweed + Algae"))
DT50_values_nls <- log(2)/k_values_nls[,2:4]
colnames(DT50_values_nls)[1] <- "DT50"
DT50_values_nls$Treatment <- k_values_nls$Treatment

#setting up variables for predictions and figures
day.s <- seq(from = 0, to =21, length.out=1000)
modls <- list(nls_Abiotic_Control, nls_Microbes, nls_Duckweed_Microbes_Salt, 
				nls_Duckweed_Microbes, 
				nls_Algae, nls_Duckweed, nls_Duckweed_Algae)
profiles <- list(profile(nls_Abiotic_Control), profile(nls_Microbes), profile(nls_Duckweed_Microbes_Salt), 
				profile(nls_Duckweed_Microbes), 
				profile(nls_Algae), profile(nls_Duckweed), profile(nls_Duckweed_Algae))
colors <- c("black","blue","tan4","turquoise3","darkgoldenrod1","forestgreen","yellowgreen")				
dfs <- list(bzt[bzt$Treatment=='Abiotic Control',], bzt[bzt$Treatment=='Microbes',], bzt[bzt$Treatment=='Duckweed + Microbes + Salt',], 
			bzt[bzt$Treatment=='Duckweed + Microbes',], 
			bzt[bzt$Treatment=='Algae',], bzt[bzt$Treatment=='Duckweed',], bzt[bzt$Treatment=='Duckweed + Algae',])
modnames <- c("Abiotic Control","Microbes","Duckweed + Microbes + Salt", "Duckweed + Microbes", "Algae","Duckweed","Duckweed + Algae")

#predictions
predtime<- sapply(modls, function(m) predict(m, newdata=data.frame(Day=day.s)))

#new result
summary(MCMCglmm(Concentration~Treatment, data = bzt[bzt$Day==21,],verbose=F,nitt=501000,burnin=1000,thin=10))
#

std.error <- function(dat, na.rm=TRUE) {sd(dat,na.rm=na.rm)/sqrt(length(dat))}#defaults to na.rm=T
means <- sapply(dfs, function(m) tapply(m$Concentration,m$Day,mean))
ses <- sapply(dfs, function(m) tapply(m$Concentration,m$Day,std.error))
#range in percent removed
range((5-means[6,])/5)
100*(5-means[6,])/5
100*ses[6,]/5

pdf("kinetics_meansSEs_preds.pdf",width=7,height=5)
layout(matrix(c(1,1,1,2,3,4,5,2,6:9),byrow=T,ncol=4),heights=c(1,0.4,0.4))
par(oma=c(4,4,1.5,1))
par(mar=c(4,0,0,1))
plot(bzt$Concentration~bzt$Day,pch=NA,xlab="",ylab="",ylim=c(0,6.60))
	mtext("a.",side=3,line=0,at=-3)
	mtext("Day",side=1,line=2.2)
	mtext(expression("Benzotriazole "*mu*"g/L"),side=2,line=2)
	abline(h=0,lty=3)
	for(i in 1:7){
		xjitter <- jitter(c(0,3,7,11,17,21),factor=0.5)
		points(means[,i]~xjitter, type="b",col=colors[i],pch=1,lwd=1.5)
		arrows(x0=xjitter,y0=means[,i]-ses[,i],y1=means[,i]+ ses[,i],length=0,col=colors[i])
	}
  	legend(2,6.95,modnames, fill=colors,bty="n",ncol=2)
par(mar=c(13,4.5,0,0))
plot(k_values_nls$k_value ~ c(1:7),xaxt="n",pch=16,col=colors,ylim=c(0,0.27),ylab="",xlab="")
	mtext("b.",side=3,line=0,at=-3)
	arrows(c(1:7),y0=k_values_nls$upperCI, y1= k_values_nls$lowerCI,length=0,col=colors)
	mtext(expression(First-order~Kinetic~Constant~(day^{-1})),side=2,line=2,at=-0.01)
	axis(at=c(1:7),labels=(modnames),las=2,side=1)
par(mar=c(0,0,0,0))
for(i in 1:ncol(predtime)){
plot(Concentration~Day,data=dfs[[i]],ylim=c(0,5.75),pch=NA,xlab="",ylab="",xaxt="n",yaxt="n")
	lines(predtime[,i]~day.s,col=colors[i])
	abline(h=0,lty=3)
	points(Concentration~Day,data=dfs[[i]],col=colors[i])
	text(10.5,5.5,modnames[i],adj=c(0.5,0.5))#,side=3,line=0.5)
	if(i%in%c(1,4)){axis(side=2)}
	if(i%in%c(4:8)){axis(side=1)}
	if(i==5){mtext("Day",side=1,line=2.2,at=21)}
	if(i==4){mtext(expression("Benzotriazole "*mu*"g/L"),side=2,at=6.5,line=2)}	
	if(i==1){mtext("c.",side=3,line=0,at=-7)}
}
dev.off()