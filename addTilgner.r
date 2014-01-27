## use data from Tilgner et al. 2012 Genome Research

rm(list=ls())
library(gplots)
library(RColorBrewer)
library(Hmisc)
source("/home/jeremy/ExonPipeline/ExonPipeline_functions.r")
source("/home/jeremy/ExonPipeline/ExonPipeline_analysis_and_plotting_Functions.r")
source("/home/home/jeremy/Code/useful_R/useful_R.r")

setwd("/home/jeremy/ExonPipeline/hg19")
#load("PlottedData6.RData")
load("Tilgner_mergedData.RData")
#load("Tilgner_exons3.RData")
load("Median_intronsRemaining.RData")



exons$length = abs(exons$start-exons$end)
### Keep only the transcripts that house the exon whose terminus is most downstream for that gene
# Also keep only the GENES where there is exon 1 there ...
exons2= mostDownstream(exons[exons$UniqueID %in% names(which(tapply(exons$FeatureCount,exons$UniqueID,min)==1)),])   
## Figure out the end of the gene (poly A).
exons2$terminus = exons2$end * ifelse(exons2$strand=='-',-1,1)
exons2$begin = exons2$start * ifelse(exons2$strand=='-',-1,1)
txEnd=abs(tapply(exons2$terminus, exons2$UniqueID,max)); TXEND = data.frame(UniqueID = names(txEnd), txEnd=txEnd)
txStart=abs(tapply(exons2$begin, exons2$UniqueID,min)); TXSTART = data.frame(UniqueID = names(txStart), txStart=txStart)

## Merge with exons
exons2$energy_ID = paste(">",paste(exons2$UniqueID,exons2$FeatureCount,exons2$isLast,sep='_'),sep='')
exons2$acceptor_ID = paste(exons2$UniqueID,exons2$FeatureCount-1,sep='_')
exons2$donor_ID = paste(exons2$UniqueID,exons2$FeatureCount,sep='_')

energyData2 = energyData; names(energyData2) = c("energy_ID","Stability")
donorData2 = donorData; donorData2$ID = sub("_([0|1])$","",donorData2$V1) ; names(donorData2) = c("","Donor","donor_ID") 
acceptorData2 = acceptorData; acceptorData$ID = sub("_([0|1])$","",acceptorData2$V1) ; names(acceptorData2) = c("","Acceptor","acceptor_ID")

mergedData = merge(merge(merge(merge(merge(exons2,energyData2,by='energy_ID'),donorData2[,-1],by='donor_ID',all=T),acceptorData2[,-1],by='acceptor_ID',all=T),TXEND,by='UniqueID'),TXSTART,by='UniqueID')

## Add HK status
mergedData$isHK = toupper(mergedData$GeneName) %in% HKlist$Gene.name
mergedData$Dist2End = abs(mergedData$txEnd - mergedData$end)
mergedData$Dist2Start = abs(mergedData$txStart - mergedData$start) + 1


## Add the length of the up and downstream introns

mergedData$UpIntronLength = mergedData$DownIntronLength = NA
mergedData = mergedData[order(mergedData$UniqueID,mergedData$FeatureCount),]
diffs = diff(mergedData$FeatureCount)
upstreamCompare = apply(mergedData[,c('strand','start','end')],1, function(x)ifelse(x[1]=='+',x[2],x[3])) # find the location where the upstream intron starts
downstreamCompare = apply(mergedData[,c('strand','start','end')],1, function(x)ifelse(x[1]=='+',x[3],x[2])) # find the location where the downstream intron starts
up_target = 0:(nrow(mergedData)-1)
up_target[mergedData$FeatureCount==1] <- NA
down_target = c(2:(nrow(mergedData)),NA)
down_target[diffs != 1] <- NA

mergedData$UpIntronLength = abs(as.numeric(upstreamCompare) - as.numeric(downstreamCompare[up_target]))
mergedData$DownIntronLength = abs(as.numeric(downstreamCompare) - as.numeric(upstreamCompare[down_target]))

## Merge with coSI_ID
mergedData$coSI_ID = paste(mergedData$chr,rowMins(mergedData[,c('start','end')])+1,rowMaxes(mergedData[,c('start','end')]),mergedData$strand,sep='_')
makeZ = function(x)(x-mean(x,na.rm=T))/sd(x,na.rm=T)
mergedData$UpIntronZ = makeZ(log(mergedData$UpIntronLength))
mergedData$DownIntronZ = makeZ(log(mergedData$DownIntronLength))

U = mergedData$UniqueID*mergedData$isLast
mergedData$MaxFeature = mergedData$FeatureCount[match(mergedData$UniqueID, U)]
g=head
#save(mergedData,file="Tilgner_mergedData.RData")

# Read in exons file
coSI = read.delim("~/RNAseq/Tilgner/v3c-exons-to-use-sixFractions-tab",head=F,stringsAsFactors=F)
coSI = coSI
co.chr = subset(coSI,V2=='NUCpolyAm');dim(co.chr) # subset for chrom-assoc data #NUCpolyAm #CHtot
co.chr2 = subset(coSI,V2=='CHtot');dim(co.chr) # subset for chrom-assoc data #NUCpolyAm #CHtot

length(setdiff(co.chr$V1,mergedData$coSI_ID))
names(co.chr)[3] = 'coSI'
names(co.chr2)[3] = 'coSI'
#(setdiff(co.chr$V1,exons$coSI_ID))
exons3 = merge(mergedData,co.chr,by.x='coSI_ID',by.y=1)
exons3$coSI2 = log_or_NA(1-exons3$coSI)*log(10) # hacky!
exons3$Group = splitByVector(exons3$Dist2End,quantile(exons3$Dist2End,na.rm=T,(1:99)/100))

exons3ch = merge(mergedData,co.chr2,by.x='coSI_ID',by.y=1)
exons3$coSI2 = log_or_NA(1-exons3$coSI)*log(10) # hacky!
exons3$Group = splitByVector(exons3$Dist2End,quantile(exons3$Dist2End,na.rm=T,(1:99)/100))
exons3ch$Group = splitByVector(exons3ch$Dist2End,quantile(exons3ch$Dist2End,na.rm=T,(1:99)/100))


#save(exons3,exons3ch,file="Tilgner_exons3.RData")



#####################################################################
## Model MEDIANS
#####################################################################

#exons3$Group = splitByVector(exons3$Dist2End,quantile(exons3$Dist2End,na.rm=T,(1:99)/100))
getQuantilMedians = function(myData,N=100){
  myData$Group = splitByVector(myData$Dist2End,quantile(myData$Dist2End,na.rm=T,(1:(N-1))/N))
  data.frame(LengthMedian = tapply(myData$Dist2End,myData$Group,median), coSI_median = tapply(myData$coSI,myData$Group,median))
}
QuantileData  = getQuantilMedians(exons3)
QuantileData2  = getQuantilMedians(exons3ch)
QuantileData_HK = getQuantilMedians(subset(exons3,isHK==1))
QuantileData_non = getQuantilMedians(subset(exons3,isHK==0))


# Tilgner-type fit
fit3 = lm( coSI_median ~ log(LengthMedian), data=QuantileData)

#####################################################################
#########  TRY OUT SEVERAL FITS #########
## Fit to gamma
#g=NA
#n = nls(coSI_median ~ A*pgamma(LengthMedian,g,k) + B, data=QuantileData,start=c(A=.8,k=1/10000,B=-.2,g=1) )
#coefs_n = c(coef(n),g)

# Fit to step function
#StepF = function(x,A,B,C) A* x /(x + C) + B
StepF = function(x,A,B,C) A* sapply(x,function(y)max((y-C) /y, 0)) + B
s = nls(coSI_median ~ StepF(LengthMedian,A,B,C), data=QuantileData,start=c(A=.4,B=.6,C=50000))
coefs_s = coef(s)

# Fit to log (Dist)
Log = function(x,A,B) A* log(x) + B
l = nls(coSI_median ~ Log(LengthMedian,A,B), data=QuantileData,start=c(A=.33,B=0.05))
coefs_l = coef(l)

# Fit to sigmoid function
SigmoidF = function(x,A,B,C) A* x /(x + C) + B
f = nls(coSI_median ~ SigmoidF(LengthMedian,A,B,C), data=QuantileData,start=c(A=.4,B=.6,C=50000))
coefs_f = coef(f)

## Fit to AVERAGE of gamma ... (integral/length)
MeanGamma = function(LengthMedian,g,k,A,B)  sapply(LengthMedian, function(x)A * integrate(pgamma,0,x, shape=g, rate=k)$value/x + B)

g = nls(coSI_median ~ MeanGamma(LengthMedian,g,k,A,B), data=QuantileData,start=c(g=3,k=1/1000,A=.8,B=.2),lower=c(0,0,0,0),upper=c(Inf,Inf,Inf,Inf),algo='port');
g2 = nls(coSI_median ~ MeanGamma(LengthMedian,g,k,A,B), data=QuantileData,start=c(g=1,k=1/1000,A=.24,B=.7),lower=c(1,0,0,0),upper=c(1,Inf,Inf,Inf),algo='port');
g3 = nls(coSI_median ~ MeanGamma(LengthMedian,g,k,A,B), data=QuantileData,start=c(g=2,k=1/1000,A=.24,B=.7),lower=c(2,0,0,0),upper=c(2,Inf,Inf,Inf),algo='port');
coefs_g = coef(g)
coefs_g2 = coef(g2)
coefs_g3 = coef(g3)

## Average of phase-type distribution
CDF_phasetype = function(t,a,b){
  (a - b - a*exp(-b*t) + b*exp(-a*t))/(a - b)
}
MeanPhase = function(LengthMedian,a,d,A,B)  sapply(LengthMedian, function(x)A * integrate(CDF_phasetype,0,x, a=a, b=a+d)$value/x + B)
# Here I am making a + d = b, so that I can limit b to strictly positive
p = nls(coSI_median ~ MeanPhase(LengthMedian,a,d,A,B), data=QuantileData,start=c(a=0.0002178,d=1e-6,A=.25,B=.7),lower=c(0,1e-6,0,0),upper=c(Inf,Inf,Inf,Inf),algo='port');
Phase_objective = function(params)sum((MeanPhase(QuantileData$LengthMedian,a=params[1],d=params[2],A=params[3],B=params[4])-QuantileData$coSI_median)^2)
p = nlminb(start=c(a=0.0002178,d=1e-10,A=.25,B=.7),  Phase_objective, lower=c(0,1e-10,0,0))
coefs_p = p$par

## SIGMOID FUNCTION WORKS THE BEST
# But the gamma of .5 is pretty close, followed by the gamma of 1
# The phase-type distribution ends up just converging to the gamma of 1!!! (with the 2nd reaction extremely fast!)

## Can I incorporate the rate of falling off the template?
#
#exp_decay = function(x,rate,decay, A, B){
#
#  A * integrate(function(x,rate,decay)pgamma(
#}
#

#gammaDecay = function(x,shape,rate,decay)pgamma(x,shape=shape, rate=rate)*exp(-decay*x)
##gammaDecay = function(x,shape,rate,decay)pgamma(x,shape=shape, rate=rate)*(1-x*decay)
#DecayGamma = function(LengthMedian,g,k,A,B,d)  sapply(LengthMedian, function(x)A * integrate(gammaDecay,0,x, shape=g, rate=k, decay=d)$value/x + B)
#g1.decay = nls(coSI_median ~ DecayGamma(LengthMedian,g,k,A,B,d), data=QuantileData,start=c(g=1,k=1/50000,A=.24,B=.7,d=.00001),lower=c(0,0,0,0,0),upper=c(Inf,Inf,Inf,Inf,Inf),algo='port');
#coefs_g1.d = coef(g1.decay)
### This just converges to decay = 0 ...

#DecayGamma_objective = function(params)sum((DecayGamma(QuantileData$LengthMedian,g=params[1],k=params[2],A=params[3],B=params[4],d=params[5])-QuantileData$coSI_median)^2)
#p = nlminb(start=c(g=1,k=1/10000,A=.25,B=.7, d=.0001p),  DecayGamma_objective, lower=c(1,0,0,0,0), upper=c(1,Inf,Inf,Inf,Inf))
#coefs_p = p$par

MeanGamma2 = function(LengthMedian,g,k,A,B)  sapply(LengthMedian, function(x)A * integrate(pgamma,0,x+B, shape=g, rate=k)$value/(x+B))
gM2 = nls(coSI_median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileData,start=c(g=.5,k=5e-5,A=.95,B=100),lower=c(0,0,0,0),upper=c(Inf,Inf,Inf,Inf),algo='port');
coef(gM2)->coef_gM2
gM2_1 = nls(coSI_median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileData,start=c(g=1,k=5e-5,A=.95,B=4000),lower=c(1,0,0,0),upper=c(1,Inf,Inf,Inf),algo='port');
coef(gM2_1)->coef_gM2_1
gM2_2 = nls(coSI_median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileData,start=c(g=2,k=5e-5,A=.95,B=4000),lower=c(2,0,0,0),upper=c(2,Inf,Inf,Inf),algo='port');
coef(gM2_2)->coef_gM2_2

#####################################################################
###### Plot linear ######
pdf("Tilgner_median_fits.pdf",width=12,height=6)
par(mfrow=c(1,2))
plot((QuantileData$LengthMedian), QuantileData$coSI_median, ylab='Median coSI', xlab='Distance to End of Gene', main='Median coSI index: Poly(A) (-)')
lines(seq(1,400000,1000),coef(fit3)[1] + coef(fit3)[2] *log(seq(1,400000,1000)),lty=2)
#lines(seq(1,400000,1000),coefs_n[1] * pgamma(seq(1,400000,1000),coefs_n[4], coefs_n[2]) + coefs_n[3],col='blue')
## Try out other shape constants
#lines(seq(1,400000,1000),coefs_n[1] * pgamma(seq(1,400000,1000),1, coefs_n[2]) + coefs_n[3],col='blue', lty=3)
lines(xx<-seq(1,500000,1000),StepF(xx,coefs_s[1], coefs_s[2], coefs_s[3]),col='green')
#lines(xx,SigmoidF(xx,coefs_f[1], coefs_f[2], coefs_f[3]),col='cyan')
lines((xx),MeanGamma(xx,coefs_g[1], coefs_g[2],coefs_g[3], coefs_g[4]),col='red')
lines((xx),MeanGamma(xx,coefs_g2[1], coefs_g2[2],coefs_g2[3], coefs_g2[4]),col='red',lty=2)
lines((xx),MeanGamma(xx,coefs_g3[1], coefs_g3[2],coefs_g3[3], coefs_g3[4]),col='red',lty=3)
legend('bottomright',col=c('black','green','red','red','red'),lty=c(2,1,1,2,3),leg=c('Log(dist) fit','Step function mechanics','Gamma Mechanics','Exponential Mechanics','Gamma(2) Mechanics'))

###### Plot log ######
plot(log10(QuantileData$LengthMedian), QuantileData$coSI_median, ylab='Median coSI', xlab='LOG10 Distance to End of Gene', main='Median coSI index: Poly(A) (-)')
abline(a=coef(fit3)[1], b= coef(fit3)[2], lty=2)
# Gamma
#lines(log(seq(1,400000,1000)),coefs_n[1] * pgamma(seq(1,400000,1000),coefs_n[4], coefs_n[2]) + coefs_n[3],col='blue')
#lines(log(seq(1,400000,1000)),coefs_n[1] * pgamma(seq(1,400000,1000),1, coefs_n[2]) + coefs_n[3],col='red')
lines(log10(xx<-seq(1,500000,1000)),StepF(xx,coefs_s[1], coefs_s[2], coefs_s[3]),col='green')
#lines(log10(xx),SigmoidF(xx,coefs_f[1], coefs_f[2], coefs_f[3]),col='cyan')
lines(log10(xx),MeanGamma(xx,coefs_g[1], coefs_g[2],coefs_g[3], coefs_g[4]),col='red')
#lines(log(xx),MeanPhase(xx,coefs_p[1], coefs_p[2],coefs_p[3], coefs_p[4]),col='blue',lty=1,lwd=2)
lines(log10(xx),MeanGamma(xx,coefs_g2[1], coefs_g2[2],coefs_g2[3], coefs_g2[4]),col='red',lty=2)
lines(log10(xx),MeanGamma(xx,coefs_g3[1], coefs_g3[2],coefs_g3[3], coefs_g3[4]),col='red',lty=3)
lines(log10(xx),MeanGamma2(xx,coef_gM2[1], coef_gM2[2],coef_gM2[3], coef_gM2[4]),col='darkgrey',lty=1)


legend('bottomright',col=c('black','green','red','red','red'),lty=c(2,1,1,2,3),leg=c('Log(dist) fit','Step function mechanics','Gamma Mechanics','Exponential Mechanics','Gamma(2) Mechanics'))
dev.off()


##################################
## What do data look like for lasts vs nonlasts?
exons3$intronsRemaining = exons3$MaxFeature - exons3$FeatureCount
QuantileData_lasts = getQuantilMedians(subset(exons3,intronsRemaining==1),20)
QuantileData_nonLasts = getQuantilMedians(subset(exons3,intronsRemaining == 2),20)
QuantileData_nonLasts3 = getQuantilMedians(subset(exons3,intronsRemaining == 3),20)
QuantileData_nonLasts4 = getQuantilMedians(subset(exons3,intronsRemaining == 4),20)
QuantileData_nonLasts_g8 = getQuantilMedians(subset(exons3,intronsRemaining > 4),100)
g_nonLasts4 = nls(coSI_median ~ MeanGamma(LengthMedian,g,k,A,B), data=QuantileData_nonLasts_g8,start=c(g=2,k=1/1000,A=.3,B=.7),lower=c(0,0,0,0),upper=c(Inf,Inf,Inf,Inf),algo='port');

QuantileEnergyData_nonlasts2 = getQuantilMediansEnergy(subset(exons3,intronsRemaining==2),20)
QuantileEnergyData_nonlasts3 = getQuantilMediansEnergy(subset(exons3,intronsRemaining==3),20)

plot(QuantileData_lasts, xlim=c(500,80000), ylim=c(0.4,1))
points(QuantileData_nonLasts,col=2)
points(QuantileData_nonLasts3,col=3)

HK = c(1,0)

pdf("Median_coSI_by_remaining_introns_withLines.pdf")
plot(log10(xx),MeanGamma(xx,coefs_g2[1], coefs_g2[2],coefs_g2[3], coefs_g2[4]),type='l',col='black', xlim=log10(c(500,420000)), ylim=c(0.6,.95), ylab='Median coSI', xlab='LOG Distance to End of Gene', main='Median coSI index for # of remaining introns: Nuclear Poly(A) (-)')
cols=brewer.pal(10,'Spectral')
allFits = list()
Residuals = list()
MedLengths = list()
MedCoSI = list()
S=sapply(1:10, function(i){
     if(i<10)
      QD = getQuantilMedians(subset(subset(exons3,isHK %in% HK),intronsRemaining==i),20)
     else
      QD = getQuantilMedians(subset(subset(exons3,isHK %in% HK),intronsRemaining >=i),20)      
     MedLengths[[i]] <<- QD[,1]
     MedCoSI[[i]] <<- QD[,2]      
     FIT <- lm(coSI_median ~ log10(LengthMedian),QD)
     allFits[[i]]<<- FIT
     abline(a=coef(FIT)[1], b=coef(FIT)[2],col=cols[i],lty=2)
     points(log10(QD[,1]),QD[,2],col=cols[i])    
     
     Residuals[[i]] <<- QD$coSI - MeanGamma(QD$LengthMedian,coefs_g2[1], coefs_g2[2],coefs_g2[3], coefs_g2[4])
     
})
legend('bottomright',col=cols,pch=20,leg=c(1:9,'>=10'))
#points(log10(QuantileData_nonLasts[,1]),QuantileData_nonLasts[,2],col=2)  
dev.off()

### Next I am trying to show quantitatively that the last/second to last exons have lower coSI
ResidualsAll <<- subset(exons3,isHK %in% HK)$coSI - MeanGamma(subset(exons3,isHK %in% HK)$Dist2End,coefs_g2[1], coefs_g2[2],coefs_g2[3], coefs_g2[4])
group10 = splitByVector(subset(exons3,isHK %in% HK)$Dist2End,quantile(subset(exons3,isHK %in% HK)$Dist2End,na.rm=T,(1:9)/10))
intronsLeft = ifelse(subset(exons3,isHK %in% HK)$intronsRemaining < 10, subset(exons3,isHK %in% HK)$intronsRemaining, 10) 

### Plot to show how residuals are distributed
MedLengths2 = unlist(MedLengths)
MedCoSI2 = unlist(MedCoSI)
gMed = nls(MedCoSI2 ~ MeanGamma(MedLengths2,g,k,A,B),start=c(g=.5,k=1e-5,A=.28,B=.72),lower=c(0,0,0,0),upper=c(Inf,Inf,Inf,Inf),algo='port');

pdf("Median_coSI_by_remaining_introns.pdf",height=10,width=10)
#pdf("Median_coSI_by_remaining_introns_HK.pdf",height=10,width=10)
par(mfrow=c(2,2))
plot(log10(xx),MeanGamma(xx,coefs_g2[1], coefs_g2[2],coefs_g2[3], coefs_g2[4]),type='l',col='black', xlim=log10(c(500,420000)), ylim=c(0.6,.95), ylab='Median coSI', xlab='LOG10 Distance to End of Gene', main='Nuclear Poly(A) (-):\nMedian coSI index for # of remaining introns')
#plot(log10(xx),MeanGamma(xx,coefs_g2[1], coefs_g2[2],coefs_g2[3], coefs_g2[4]),type='l',col='black', xlim=log10(c(500,420000)), ylim=c(0.6,.95), ylab='Median coSI', xlab='LOG10 Distance to End of Gene', main='Nuclear Poly(A) (-): Housekeeping Genes\nMedian coSI index for # of remaining introns')

S=sapply(1:10, function(i)     points(log10(MedLengths[[i]]),MedCoSI[[i]],col=cols[i],pch=20))
lines(log10(xx),MeanGamma(xx,coef(gMed)[1], coef(gMed)[2],coef(gMed)[3], coef(gMed)[4]),col='black',lty=2)

legend('bottomright',col=cols,pch=20,leg=c(1:9,'>=10'))

boxplot(Residuals, main='Residuals of median data for 10 intron-remaining groups',names=c(1:9,'>=10'), xlab='# of introns downstream',notch=T,col=cols)
abline(h=0)

# as a control, show this same data but sort all the medians first ...
cols2 = rev(brewer.pal(10,'PRGn'))

plot(log10(xx),MeanGamma(xx,coefs_g2[1], coefs_g2[2],coefs_g2[3], coefs_g2[4]),type='l',col='black', xlim=log10(c(500,420000)), ylim=c(0.6,.95), ylab='Median coSI', xlab='LOG10 Distance to End of Gene', main='Median coSI index for # of remaining introns')
points(log10(MedLengths2),MedCoSI2,col=cols2[group10.small],pch=20)
legend('bottomright',col=cols2,pch=20,leg=paste("Quantile",1:10))

group10.small = splitByVector(MedLengths2,quantile(MedLengths2,na.rm=T,(1:9)/10))
boxplot(MedCoSI2 - MeanGamma(MedLengths2,coefs_g2[1], coefs_g2[2],coefs_g2[3], coefs_g2[4]) ~ group10.small, main='Residuals of median data for 10 intron-remaining groups',names=c(1:10),xlab='Quantiles of distance to end of gene', notch=T,col=cols2)
abline(h=0)
dev.off()

##################################################
## Redo fit for just those that are not at all the LAST exon
QuantileData_nonLasts_g4 = getQuantilMedians(subset(exons3,intronsRemaining > 4),100)
#g2_nonLasts_g4 = nls(coSI_median ~ MeanGamma(LengthMedian,g,k,A,B), data=QuantileData_nonLasts_g4,start=c(g=1,k=1/1000,A=.24,B=.7),lower=c(1,0,0,0),upper=c(1,Inf,Inf,Inf),algo='port');
#g_nonLasts_g4 = nls(coSI_median ~ MeanGamma(LengthMedian,g,k,A,B), data=QuantileData_nonLasts_g4,start=c(g=.5,k=1e-5,A=.28,B=.72),lower=c(0,0,0,0),upper=c(Inf,Inf,Inf,Inf),algo='port');
gM2_nonLasts_g4 = nls(coSI_median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileData_nonLasts_g4,start=c(g=.5,k=1e-5,A=.96,B=1000),lower=c(0,0,0,0),upper=c(Inf,Inf,Inf,Inf),algo='port');
coef(g2_nonLasts_g4)
#save.image(file="Median_intronsRemaining.RData")
plot(log10(QuantileData_nonLasts_g4$LengthMedian), QuantileData_nonLasts_g4$coSI_median, ylab='Median coSI', xlab='LOG10 Distance to End of Gene', main='Median coSI index: Exons with > 4 introns downstream')
lines(log10(xx),MeanGamma2(xx,coef(gM2_nonLasts_g4)[1], coef(gM2_nonLasts_g4)[2],coef(gM2_nonLasts_g4)[3], coef(gM2_nonLasts_g4)[4]),col='black',lty=2)
#lines(log10(xx),MeanGamma(xx,coef(g_nonLasts_g4)[1], coef(g_nonLasts_g4)[2],coef(g_nonLasts_g4)[3], coef(g_nonLasts_g4)[4]),col='blue',lty=2)
lines(log10(xx),MeanGamma2(xx,coef_gM2[1], coef_gM2[2],coef_gM2[3], coef_gM2[4]),col='black',lty=1)
points(log10(QuantileData$LengthMedian), QuantileData$coSI_median, col='green',pch=20)



##################################################
### Compare fits in HK genes.
g.hk = nls(coSI_median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileData_HK,start=c(g=.3,k=1/10000,A=.95,B=4000),lower=c(0,0,0,0),upper=c(Inf,Inf,Inf,Inf),algo='port');
g2.hk = nls(coSI_median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileData_HK,start=c(g=1,k=1/10000,A=.95,B=1000),lower=c(1,0,0,0),upper=c(1,Inf,Inf,Inf),algo='port');
g3.hk = nls(coSI_median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileData_HK,start=c(g=2,k=1/10000,A=.95,B=1000),lower=c(2,0,0,0),upper=c(2,Inf,Inf,Inf),algo='port');

g.non = nls(coSI_median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileData_non,start=c(g=.22,k=1/10000,A=.95,B=4000),lower=c(0,0,0,0),upper=c(Inf,Inf,Inf,Inf),algo='port');
g2.non = nls(coSI_median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileData_non,start=c(g=1,k=1/10000,A=.95,B=1000),lower=c(1,0,0,0),upper=c(1,Inf,Inf,Inf),algo='port');
g3.non = nls(coSI_median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileData_non,start=c(g=2,k=1/10000,A=.95,B=1000),lower=c(2,0,0,0),upper=c(2,Inf,Inf,Inf),algo='port');

###############
### Make a table of all the rate/elongation parameters in both HK, Non, ALL genes ...
###############

gammaFits = rbind(All=c(Free=coef(gM2)[2]/coef(gM2)[1]*3000, One= coef(gM2)[2]/coef(gM2)[1]*3000, Two=coef(gM2)[2]/coef(gM2)[1]*3000),
  HK=c(Free=coef(g.hk)[2]/coef(g.hk)[1]*3000, One= coef(g2.hk)[2]/coef(g2.hk)[1]*3000, Two=coef(g3.hk)[2]/coef(g3.hk)[1]*3000),
  Non=c(Free=coef(g.non)[2]/coef(g.non)[1]*3000, One= coef(g2.non)[2]/coef(g2.non)[1]*3000, Two=coef(g3.non)[2]/coef(g3.non)[1]*3000))

  #       Free.k     One.k     Two.k
#All 0.2699642 0.2926816 0.3267481
#HK  0.2009518 0.2504865 0.2882817
#Non 0.2809882 0.2960773 0.3293210

# THis is cool! It means that there is a slightly higher rate in non-HK genes ...
write.table(gammaFits, file='GammaFits2.txt',sep="\t",quote=F)


QuantileEnergyData  = getQuantilMediansEnergy(mergedData,100)
QuantileAcceptorData  = getQuantilMediansAcceptor(mergedData,100)

getQuantilMedians2 = function(myData,N=100,col1='Dist2End', col2='coSI'){
  myData$Group = quantileGroups(myData[,col1],N)
  data.frame(GroupMedian = tapply(myData[,col1],myData$Group,median), DataMedian = tapply(myData[,col2],myData$Group,median))
}

QuantileEnergyVsCoSI =   getQuantilMedians2(exons3,100,'Stability','coSI')
QuantileAcceptorVsCoSI = getQuantilMedians2(exons3,100,'Acceptor','coSI')

getQuantilMediansEnergy = function(myData,N=100){
  myData$Group = splitByVector(myData$Dist2End,quantile(myData$Dist2End,na.rm=T,(1:(N-1))/N))
  data.frame(LengthMedian = tapply(myData$Dist2End,myData$Group,median), Energy_median = tapply(myData$Stability,myData$Group,median))
}
getQuantilMediansAcceptor = function(myData,N=100){
  myData$Group = splitByVector(myData$Dist2End,quantile(myData$Dist2End,na.rm=T,(1:(N-1))/N))
  data.frame(LengthMedian = tapply(myData$Dist2End,myData$Group,median,na.rm=T), Acceptor_median = tapply(myData$Acceptor,myData$Group,median,na.rm=T))
}
### Now look at the ENERGY of exons:
pdf("Median_energy_by_remaining_introns.pdf",height=8,width=8)
mergedData$intronsRemaining = mergedData$MaxFeature - mergedData$FeatureCount
par(mfrow=c(2,2))
cols=brewer.pal(10,'Spectral')
plot(log10(QuantileEnergyData[,1]), QuantileEnergyData[,2],,xlim=log10(c(500,420000)),ylim=c(-1,1),main='Median Stability', ylab='Median Exon Stability', xlab="LOG10 Distance to End of Gene") 

MedLengths = list()
MedEnergy = list()

#coefs_f = coef(f)

S=sapply(1:10, function(i){
     if(i<10)
      QD = getQuantilMediansEnergy(subset(mergedData,intronsRemaining==i),20)
     else
      QD = getQuantilMediansEnergy(subset(mergedData,intronsRemaining >=i),20)      
     MedLengths[[i]] <<- QD[,1]
     MedEnergy[[i]] <<- QD[,2]      
     if(i==1)
      plot(log10(QD[,1]),QD[,2],col=cols[i],xlim=log10(c(500,420000)),ylim=c(-1,1),main='Median Stability by # introns remaining', ylab='Median Exon Stability', xlab="LOG10 Distance to End of Gene")
     else
      points(log10(QD[,1]),QD[,2],col=cols[i])                 })
legend('bottomright',col=cols,pch=20,leg=c(1:9,'>=10'))
#points(log10(QuantileData_nonLasts[,1]),QuantileData_nonLasts[,2],col=2)  
#dev.off()


#pdf("Median_acceptor_by_remaining_introns.pdf")
plot(log10(QuantileAcceptorData[,1]), QuantileAcceptorData[,2],,xlim=log10(c(500,420000)),ylim=c(-1,1),main='Median Acceptor', ylab='Median Exon Stability', xlab="LOG10 Distance to End of Gene") 

cols=brewer.pal(10,'Spectral')

MedLengths = list()
MedAcceptor = list()
S=sapply(1:10, function(i){
     if(i<10)
      QD = getQuantilMediansAcceptor(subset(mergedData,intronsRemaining==i),20)
     else
      QD = getQuantilMediansAcceptor(subset(mergedData,intronsRemaining >=i),20)      
     MedLengths[[i]] <<- QD[,1]
     MedAcceptor[[i]] <<- QD[,2]      
     if(i==1)
      plot(log10(QD[,1]),QD[,2],col=cols[i],xlim=log10(c(500,420000)),ylim=c(-1,1),main='Median Splice Site Strength by # introns remaining', ylab='Median Acceptor Score', xlab="LOG10 Distance to End of Gene")
     else
      points(log10(QD[,1]),QD[,2],col=cols[i])                 })
legend('bottomright',col=cols,pch=20,leg=c(1:9,'>=10'))
#points(log10(QuantileData_nonLasts[,1]),QuantileData_nonLasts[,2],col=2)  
dev.off()





############### Making the really cool looking plots for Energy and Acceptor strength as well
# QuantileAcceptorVsCoSI = getQuantilMedians2(exons3,100,'Acceptor','coSI')
pdf("Tilgner_Median_coSI_by_energy-acceptor.pdf")
par(mfrow=c(2,2))

### First plot the coSI vs Energy ###
QuantileEnergyVsCoSI =   getQuantilMedians2(exons3,100,'Stability','coSI')
EnergyFit = lm(DataMedian ~ GroupMedian,QuantileEnergyVsCoSI)
plot(QuantileEnergyVsCoSI, xlab="Stability Score", ylab="Median coSI",ylim=c(.6, 0.96))
lines(EE <- seq(-3,3,.05), predict(EnergyFit,data.frame(GroupMedian=EE))) 

## Next plot for different introns remaining ##
plot(EE, predict(EnergyFit,data.frame(GroupMedian=EE)), type='l',ylim=c(.6, 0.96), xlab="Stability Score", ylab="Median coSI" , main='Median coSI index for # of remaining introns')
cols=brewer.pal(10,'Spectral')
allFits = list()
Residuals = list()
MedLengths = list()
MedCoSI = list()
S=sapply(1:10, function(i){
     if(i<10)
      QD = getQuantilMedians2(subset(exons3,intronsRemaining==i),20,'Stability','coSI')
     else
      QD = getQuantilMedians2(subset(exons3,intronsRemaining>=i),20,'Stability','coSI')      
     MedLengths[[i]] <<- QD[,1]
     MedCoSI[[i]] <<- QD[,2]      
     FIT <- lm(DataMedian ~ GroupMedian,QD)
     abline(a=coef(FIT)[1], b=coef(FIT)[2],col=cols[i],lty=2)
     points(QD[,1],QD[,2],col=cols[i])    
     
     Residuals[[i]] <<- QD$DataMedian - predict(EnergyFit,QD)      
})
legend('bottomright',col=cols,pch=20,leg=c(1:9,'>=10'))



QuantileAcceptorVsCoSI =   getQuantilMedians2(exons3,100,'Acceptor','coSI')
AcceptorFit = lm(DataMedian ~ GroupMedian,QuantileAcceptorVsCoSI)
plot(QuantileAcceptorVsCoSI, xlab="Acceptor Score", ylab="Median coSI",ylim=c(.6, 0.96))
lines(EE <- seq(-3,3,.05), predict(AcceptorFit,data.frame(GroupMedian=EE))) 

## Next plot for different introns remaining ##
plot(EE, predict(AcceptorFit,data.frame(GroupMedian=EE)), type='l',ylim=c(.6, 0.96), xlab="Acceptor Score", ylab="Median coSI" , main='Median coSI index for # of remaining introns')
cols=brewer.pal(10,'Spectral')
allFits = list()
Residuals = list()
MedLengths = list()
MedCoSI = list()
S=sapply(1:10, function(i){
     if(i<10)
      QD = getQuantilMedians2(subset(exons3,intronsRemaining==i),20,'Acceptor','coSI')
     else
      QD = getQuantilMedians2(subset(exons3,intronsRemaining>=i),20,'Acceptor','coSI')      
     MedLengths[[i]] <<- QD[,1]
     MedCoSI[[i]] <<- QD[,2]      
     FIT <- lm(DataMedian ~ GroupMedian,QD)
     abline(a=coef(FIT)[1], b=coef(FIT)[2],col=cols[i],lty=2)
     points(QD[,1],QD[,2],col=cols[i])    
     
     Residuals[[i]] <<- QD$DataMedian - predict(AcceptorFit,QD)      
})
legend('bottomright',col=cols,pch=20,leg=c(1:9,'>=10'))
dev.off()





### Next I am trying to show quantitatively that the residuals are important
ResidualsAll <<- exons3$coSI - MeanGamma(exons3$Dist2End,coefs_g2[1], coefs_g2[2],coefs_g2[3], coefs_g2[4])
group10 = splitByVector(exons3$Dist2End,quantile(exons3$Dist2End,na.rm=T,(1:9)/10))
intronsLeft = ifelse(exons3$intronsRemaining < 10, exons3$intronsRemaining, 10) 

### Plot to show how residuals are distributed
MedLengths2 = unlist(MedLengths)
MedCoSI2 = unlist(MedCoSI)
gMed = nls(MedCoSI2 ~ MeanGamma(MedLengths2,g,k,A,B),start=c(g=.5,k=1e-5,A=.28,B=.72),lower=c(0,0,0,0),upper=c(Inf,Inf,Inf,Inf),algo='port');

pdf("Median_coSI_by_remaining_introns.pdf",height=10,width=10)
par(mfrow=c(2,2))
plot(log10(xx),MeanGamma(xx,coefs_g2[1], coefs_g2[2],coefs_g2[3], coefs_g2[4]),type='l',col='black', xlim=log10(c(500,420000)), ylim=c(0.6,.95), ylab='Median coSI', xlab='LOG10 Distance to End of Gene', main='Nuclear Poly(A) (-)\nMedian coSI index for # of remaining introns')
#points(log10(MedLengths2),MedCoSI2,col=cols[as.vector(matrix(rep(1:10,10),byrow=T,10))])
S=sapply(1:10, function(i)     points(log10(MedLengths[[i]]),MedCoSI[[i]],col=cols[i],pch=20))
lines(log10(xx),MeanGamma(xx,coef(gMed)[1], coef(gMed)[2],coef(gMed)[3], coef(gMed)[4]),col='black',lty=2)

legend('bottomright',col=cols,pch=20,leg=c(1:9,'>=10'))

boxplot(Residuals, main='Residuals of median data for 10 intron-remaining groups',names=c(1:9,'>=10'), xlab='# of introns downstream',notch=T,col=cols)
abline(h=0)

# as a control, show this same data but sort all the medians first ...
cols2 = rev(brewer.pal(10,'PRGn'))

plot(log10(xx),MeanGamma(xx,coefs_g2[1], coefs_g2[2],coefs_g2[3], coefs_g2[4]),type='l',col='black', xlim=log10(c(500,420000)), ylim=c(0.6,.95), ylab='Median coSI', xlab='LOG10 Distance to End of Gene', main='Median coSI index for # of remaining introns')
points(log10(MedLengths2),MedCoSI2,col=cols2[group10.small],pch=20)
legend('bottomright',col=cols2,pch=20,leg=paste("Quantile",1:10))

group10.small = splitByVector(MedLengths2,quantile(MedLengths2,na.rm=T,(1:9)/10))
boxplot(MedCoSI2 - MeanGamma(MedLengths2,coefs_g2[1], coefs_g2[2],coefs_g2[3], coefs_g2[4]) ~ group10.small, main='Residuals of median data for 10 intron-remaining groups',names=c(1:10),xlab='Quantiles of distance to end of gene', notch=T,col=cols2)
abline(h=0)
dev.off()



















































########################################
### Try to relate other measures than distance
########################################
group10 =quantileGroups(subset(exons3,intronsRemaining==1,Stability),10)
# Plot 
plot(log10(subset(exons3,intronsRemaining==1,Dist2End)[,1]), subset(exons3,intronsRemaining==1,coSI)[,1],col=cols[group10],pch=20)

group10.2 =quantileGroups(rowSums(subset(exons3,intronsRemaining>1, c('Acceptor','Donor'))),10)
# Plot 
plot(log10(subset(exons3,intronsRemaining>1,Dist2End)[,1]), subset(exons3,intronsRemaining>1,coSI)[,1],pch=20,col=cols[group10.2],main="color: Donor+Acceptor")

## How about LAST exons Energy?
lastOnly = subset(mergedData,isLast==1)
exons3$lastEnergy = lastOnly$Stability[match(exons3$UniqueID, lastOnly$UniqueID)]

group10.2 =quantileGroups(exons3$lastEnergy, c('Acceptor','Donor'))),10)
# Plot 
plot(log10(exons3$Dist2End), exons3$coSI,pch=20,col=cols[group10.2],main="color: Donor+Acceptor")



#boxplot(ResidualsAll ~ group10,main='Residuals of all fits as a function of 10 distance bins', notch=T)
#boxplot(ResidualsAll ~ intronsLeft, main='Residuals of all fits as a function of intron-remaining groups',names=c(1:9,'>=10'), notch=T)


#################################
## Confidence intervals for fits:
ConfInts = sapply(allFits,function(x)confint(x)[2,])
Means = sapply(allFits,function(x)coef(x)[2])
MeansInt = sapply(allFits,function(x)coef(x)[1])
errbar(1:10,Means,ConfInts[1,], ConfInts[2,])

## 
MeanGamma(xx,coefs_g2[1], coefs_g2[2],coefs_g2[3], coefs_g2[4])
















# Hmm, interesting. Now remake the groups and test in each group whether the LAST exon has different coSI
group20 = splitByVector(exons3$Dist2End,quantile(exons3$Dist2End,na.rm=T,(1:19)/20))
for (i in 1:20){
  print(wilcox.test(exons3$coSI[exons3$MaxFeature-exons3$FeatureCount==1 & group20 == i], exons3$coSI[exons3$MaxFeature-exons3$FeatureCount>1 & group20 == i])$p.)
}
# highly significant!






#####################################################################
# Attempt to do fits also for chromatin RNA 

g.2 = nls(coSI_median ~ MeanGamma(LengthMedian,g,k,A,B), data=QuantileData2,start=c(g=3,k=1/1000,A=.8,B=.2),lower=c(0,0,0,0),upper=c(Inf,Inf,Inf,Inf),algo='port');
g2.2 = nls(coSI_median ~ MeanGamma(LengthMedian,g,k,A,B), data=QuantileData2,start=c(g=1,k=1/1000,A=.24,B=.7),lower=c(1,0,0,0),upper=c(1,Inf,Inf,Inf),algo='port');
g3.2 = nls(coSI_median ~ MeanGamma(LengthMedian,g,k,A,B), data=QuantileData2,start=c(g=2,k=1/1000,A=.24,B=.7),lower=c(2,0,0,0),upper=c(2,Inf,Inf,Inf),algo='port');
coefs_g.2 = coef(g.2)
coefs_g2.2 = coef(g2.2)
coefs_g3.2 = coef(g3.2)

plot((QuantileData2$LengthMedian), QuantileData2$coSI_median, ylab='Median coSI', xlab='Distance to End of Gene', main='Median coSI index: total Chromatin')
lines((xx),MeanGamma(xx,coefs_g.2[1], coefs_g.2[2],coefs_g.2[3], coefs_g.2[4]),col='red')
lines((xx),MeanGamma(xx,coefs_g2.2[1], coefs_g2.2[2],coefs_g2.2[3], coefs_g2.2[4]),col='red',lty=2)
lines((xx),MeanGamma(xx,coefs_g3.2[1], coefs_g3.2[2],coefs_g3.2[3], coefs_g3.2[4]),col='red',lty=3)
legend('bottomright',col=c('red','red','red'),lty=c(1,2,3),leg=c('Gamma Mechanics','Exponential Mechanics','Gamma(shape 2) Mechanics'))


###### Plot log ######
plot(log(QuantileData2$LengthMedian), QuantileData2$coSI_median, ylab='Median coSI', xlab='LOG Distance to End of Gene', main='Median coSI index: total Chromatin')
lines(log(xx),MeanGamma(xx,coefs_g.2[1], coefs_g.2[2],coefs_g.2[3], coefs_g.2[4]),col='red')
lines(log(xx),MeanGamma(xx,coefs_g2.2[1], coefs_g2.2[2],coefs_g2.2[3], coefs_g2.2[4]),col='red',lty=2)
lines(log(xx),MeanGamma(xx,coefs_g3.2[1], coefs_g3.2[2],coefs_g3.2[3], coefs_g3.2[4]),col='red',lty=3)
legend('bottomright',col=c('red','red','red'),lty=c(1,2,3),leg=c('Gamma Mechanics','Exponential Mechanics','Gamma(shape 2) Mechanics'))


dev.off()



#####################################################################
# 0-bottom the quantile data
QuantileData$coSI_0 = (QuantileData$coSI_median-min(QuantileData$coSI_median))*max(QuantileData$coSI_median)/diff(range(QuantileData$coSI_median))
# Gamma fit
n2 = nls(coSI_0 ~ A*pgamma(LengthMedian,g,k) + B, data=QuantileData,start=c(A=.9,k=1/50000,B=0.1,g=.5) )
coefs_n = c(coef(n2))
g_0 = nls(coSI_0 ~ MeanGamma(LengthMedian,g,k,A,B), data=QuantileData,start=c(g=3,k=1/5000,A=.24,B=.7),lower=c(0,0,0,0),upper=c(Inf,Inf,Inf,Inf),algo='port');
plot(log(QuantileData$LengthMedian), QuantileData$coSI_0)
lines(log(xx),MeanGamma(xx,coef(g_0)[1], coef(g_0)[2],coef(g_0)[3], coef(g_0)[4]),col='red')












### Phase-type distribution:  
SimulateMarkov = function(Q,T){
  
  E = eigen(Q);
  D = E$values
  V = E$vectors
  
  D2 = diag(exp(D*T));
  
  V %*% D2 %*% solve(V); # This is equivalent to expm(Q*T)
}

# Phase-type distribution with rates 2, 5
a=2;b=5; S = matrix(c(-a,0,0,a,-b,0,0,b,0),nr=3);
S0=-S%*% rbind(1,1,1)

## Generate the CDF:
fx = sapply(X<- seq(0,20,.1), function(x)tail(SimulateMarkov(S,x)[1,],1)) # gets the top right value of the matrix

## Here it is from Matlab symbolic:
#
CDF_phasetype = function(a,b,t){
  (a - b - a*exp(-b*t) + b*exp(-a*t))/(a - b)
}













### ANOVA (F-test)
### Test whether coSI is related to GeneName ... NOT correcting for Dist2End yet ...
#dat1=subset(exons3,GeneName %in% exons3$GeneName[duplicated(exons3$GeneName)],c('coSI','GeneName'))
dat1=exons3[,c('coSI','GeneName')]

doF = function(dat1){
  measurement = dat1[,1]
  treatment = dat1[,2]
  datMean = mean(measurement)
  df1 = (length(unique(treatment))-1)  
  df2 = (length(treatment)-1)  
  MSE_treatment = sum( tapply(measurement,treatment,function(x)(mean(x)-datMean)^2) * tapply(measurement,treatment,length) )/df1
  MSE_resdidual = sum( tapply(measurement,treatment,function(x)sum((mean(x)-x)^2)))/df2
  Fstat = MSE_treatment / MSE_resdidual
  Fcum <- pf(Fstat,df1,df2)
  cat(sprintf("F statistic: %f     df1: %d      $df2: %d        Expected F mean:  %f    p-value (F > expected): %f\n", Fstat, df1, df2, df2/(df2-2), 1-Fcum))
} 
doF(dat1)
#F statistic: 5.434691     df1: 4175      $df2: 24769        Expected F mean:  1.000081    p-value (F > expected): 0.000000
## This shows that there is a lot of between-gene variance, so now I will incorporate that into the model.



################################################################################
# Regression analysis
#################################################################################
drawFit = function(fit){
  smoothScatter(fit$model[,1], fit$fitted.values,ylim=range(fit$model[,1]),xlim=range(fit$model[,1]),ylab='Predicted Outcome', xlab='Outcome')
  grid()
}
drawFitNLS = function(fit){
  smoothScatter(tmp<-fit$m$lhs(), fitted(fit),ylim=range(tmp),xlim=range(tmp),ylab='Predicted Outcome', xlab='Outcome')
  grid()
}
drawFit2 = function(fit,dataCol, ...){
  smoothScatter(dataCol, fit$model[,1], ...,xlab='Independent Variable', ylab='Outcome')  
  points(dataCol, predict(fit),col='red')
  lines(L<<-lowess(dataCol,predict(fit),iter=6,f=1/3),col='green')
}
drawFitNLS2 = function(fit,dataCol, ...){
  smoothScatter(dataCol, fit$m$lhs(), ...,xlab='Independent Variable', ylab='Outcome')  
  points(dataCol, predict(fit),col='red')
   lines(L<<-loess(dataCol,predict(fit),iter=6,f=1/3),col='green')
}
plotH = function(fit1){
  #fit1 = lm(x$model[,1] ~ x$model[,2])
 # plot(fit1$model[,1],fit1$residuals)
  plot(fitted(fit1), resid(fit1))
}
convert_coSI = function(x) log(1-x) 
unconvert_coSI = function(y) 1-exp(y)


#################################################################################
# What is effect of adding Gene-specific parameters?
keepNames = names(table(exons3$GeneName))[table(exons3$GeneName) > 4]
#keepNames = exons3$GeneName[duplicated(exons3$GeneName)]
exons3.2 = subset(exons3, GeneName %in% keepNames)
exons3.2$coSI2 = convert_coSI(exons3.2$coSI)
exons3.2$coSI2[exons3.2$coSI==1] <- NA# mean(exons3.2$coSI2[exons3.2$coSI != 1]) # Set these to the mean

#drawFit(aa <- lm(logit(coSI) ~ log(Dist2End), data=exons4))
drawFit(a <- lm((coSI) ~ log(Dist2End), data=exons3.2))
drawFit(b <- lm(coSI ~ log(Dist2End) + GeneName, data=exons3.2))
drawFit(c <- lm(coSI ~ log(Dist2End):GeneName, data=exons3.2))

drawFit(d <- lm(coSI2 ~ Dist2End:GeneName, data=exons3.2))
hist(coef(d)[-1],100)
plot(exons3.2$Dist2End, exons3.2$coSI2, xlim=c(0,50000))
plotH(d)

# Do F-test between HK, nonHK:
COEF = coef(b)[-(1:2)] # get rid of intercept, log(Dist) ...
COEF2 = data.frame(GeneName=sub("GeneName","",names(COEF)),Coef = as.numeric(COEF))
exons3.2.withCOEF = merge(exons3.2,COEF2, by = 'GeneName')
doF(exons3.2.withCOEF[,c('Coef','isHK')])

COEF = coef(d)[-(1)]
COEF2 = data.frame(GeneName=sub("GeneName","",names(COEF)),Coef = as.numeric(COEF))
summary(subset(COEF2, GeneName %in% subset(exons3.2,isHK==T,GeneName) ,Coef))
summary(subset(COEF2, GeneName %in% subset(exons3.2,isHK==F,GeneName) ,Coef))


#####################################################################
### For each gene, fit a line separately.
exons3.3 = subset(exons3.2, !is.na(coSI2))
keepNames2 = names(table(exons3.3$GeneName))[table(exons3.3$GeneName) > 4] #1833 genes
geneFits = t(sapply(keepNames2, function(x) coef(lm(exons3.3$coSI2[exons3.3$GeneName==x]~ exons3.3$Dist2End[exons3.3$GeneName==x]))))
geneFits2 = merge(data.frame(GeneName = rownames(geneFits),Intercept = geneFits[,1], Slope=geneFits[,2]), subset(exons3.3,!duplicated(GeneName),select=c('GeneName','isHK')),by='GeneName')
t.test(Slope ~ isHK, data=geneFits2)
#t = -0.5573, df = 572.127, p-value = 0.5776 Not significantly different

## Plot histograms of coefficients deteremined by two different fitting methods:
d <- lm(coSI2 ~ Dist2End:GeneName, data=exons3.3)

par(mfrow=c(3,2))
plot(exons3.2$Dist2End, exons3.2$coSI2, xlim=c(0,50000))

# common intercept
hist(coef(d)[-1],100, main=sprintf('Slope of fits with a common Intercept\n(%.1f%% are < 0)', 100*length(which(coef(d)[-1] < 0))/length(coef(d)[-1])))
abline(v=0, col=2)

hist(geneFits2$Slope,100,  main=sprintf('Slope of fits with a variable Slope, Intercept\n(%.1f%% are < 0)', 100*length(which(geneFits2$Slope < 0))/length(geneFits2$Slope)))
abline(v=0, col=2)

hist(geneFits2$Intercept,100,  main=sprintf('Intercept of fits with a variable Slope, Intercept\n(%.1f%% are < 0)', 100*length(which(geneFits2$Intercept < 0))/length(geneFits2$Intercept)))
abline(v=0, col=2)

# Relationship between parameters
plot(geneFits[,1], geneFits[,2], xlab='Intercept', ylab='Slope')

# What are the good fits?
geneFits2$AcceptableFit = geneFits2$Slope < 0 & geneFits2$Intercept < 0
chisq.test(geneFits2$isHK, geneFits2$AcceptableFit) # Not significant :( 

### BUT, HK genes DO splice better according to this:
t.test(coSI ~ isHK, data=exons3)
#t = -4.346, df = 6709.831, p-value = 1.407e-05
#mean in group FALSE  mean in group TRUE 
#          0.7895546           0.8011040 

## Try modeling with the acceptor score as well ...
hk1 <- lm((coSI2) ~ Dist2End, data=subset(exons3,isHK==T)); summary(hk1)
non1 <- lm((coSI2) ~ Dist2End, data=subset(exons3,isHK==F)); summary(non1)

drawFit(hk1 <- lm((coSI2) ~ Dist2End, data=subset(exons3.2,isHK==T)))
drawFit(non1 <- lm((coSI2) ~ Dist2End, data=subset(exons3.2,isHK==F)))










#################################################################################
######   Compare Tilgner's linear model with an improved one

whichDo = c(T,F); addReadThrough = 0
exons4 = exons3
exons4$coSI[exons4$coSI==1] <- .99 #
#exons4 = exons4[!is.na(exons4$coSI),]
exons4$coSI[exons4$coSI==0] <- .01
exons4$coSI_logit = logit(exons4$coSI)
exons4$coSI_n1 = normalize(exons4$coSI_logit)
unnorm = function(x) inv.logit(  x*sd(exons4$coSI_logit,na.rm=T) + mean(exons4$coSI_logit))



## test lm
par(mfrow=c(2,2))
#drawFit(a<-lm(coSI ~ exp(-(Dist2End+addReadThrough)/70000), exons4))
#drawFit(b<-lm(coSI_logit ~ log(Dist2End), exons4))
#drawFit(c<-lm(coSI_n1 ~ log(Dist2End), exons4))
#smoothScatter(unnorm(c$model[,1]),unnorm(c$fitted.values),ylim=c(0,1))
#
## test heteroskedasticity
# # test heteroskedasticity
#drawFit(d<-lm(coSI_logit ~ log(Dist2End) +exp(Acceptor)+exp(Donor)+Stability+log(length)+log(Dist2Start)+log(UpIntronLength)+log(DownIntronLength), exons4))
#drawFit(e<-lm(coSI_logit ~ exp(-(Dist2End+addReadThrough)/30000) +(Acceptor)+(Donor)+Stability+log(length)+log(Dist2Start)+log(UpIntronLength)+log(DownIntronLength), exons4))
#drawFit(f<-nls(coSI_logit ~ exp(-Dist2End/30000) , exons4,c(Tau=30000,sdf=3)))
#



# Compare linear and gamma CDF fits:
par(mfrow=c(3,2))
drawFit(B <- lm(coSI ~ log(Dist2End), data=exons4))
drawFit2(B, (exons4$Dist2End), xlim=c(0,50000)); 
print(sum(resid(B)^2))     #671.5591

drawFitNLS(A <- nls(coSI~A*pgamma(Dist2End,shape=S,rate=R) + (1-A),data=exons4, start=c(S=1,R=1/70000,A=.2427),lower=c(1,0,0),upper=c(10,1,1),algo='port'));A
drawFitNLS2(A,(exons4$Dist2End),main="NLS", xlim=c(0,50000)); 
print(sum(resid(A)^2)) 


#### Add other variables
#par(mfrow=c(2,2))
drawFit(B2 <- lm(coSI ~ log(Dist2End) + log(Dist2Start), data=exons3))
drawFit2(B2, log(exons3$Dist2End)); 
print(sum(resid(B2)^2))   #659

drawFitNLS(A2 <- nls(coSI ~ A*pgamma(Dist2End,shape=S,rate=R)*pgamma(Dist2Start,shape=S2,rate=R2) + (1-A),data=exons4, start=c(S=2, S2=2, R=1/2000, R2=1/70000,A=0.25),lower=c(1,1,0,0,0),upper=c(10,10,1/1000,1/1000,1),algo='port'))
drawFitNLS2(A2,(exons4$Dist2End),main="NLS", xlim=c(0,50000)) 
print(sum(resid(A2)^2))



# Tilgner
f=function(){
                    
  b=summary( a<-lm((coSI) ~ log(Dist2End+addReadThrough) + log(Dist2Start),subset(exons4,isHK %in% whichDo)));cat(sprintf("Tilgner pos: %f\n",sqrt(b$r.squared)));
  ###   Use model to imrpove
  b=NULL; for(i in 10000 * (1:7)){
    b[i/10000]=summary( a<-lm(logit(coSI) ~ exp(-(Dist2End+addReadThrough)/i) + log(Dist2Start),  subset(exons4,isHK %in% whichDo)))$r.squared}; 
  cat(sprintf("Model pos: %f\n",max(sqrt(b))))
  
  ### Tilgner: 'struct'
  b=summary( a<<-lm(coSI ~ log(Dist2End+addReadThrough)+exp(Acceptor)+exp(Donor)+Stability+log(length)+log(Dist2Start)+log(UpIntronLength)+log(DownIntronLength),
    subset(exons4,isHK %in% whichDo)));cat(sprintf("Tilgner +Struc: %f\n",sqrt(b$r.squared)));
    
  ### model + 'struct'
  b=NULL; for(i in 10000 * (1:3)){b[i/10000]=summary(
    a<<-lm(coSI ~ exp(-(Dist2End+addReadThrough)/i) + exp(Acceptor) + exp(Donor) + Stability + log(length) + log(Dist2Start)  + log(UpIntronLength)+log(DownIntronLength),
    subset(exons4,isHK %in% whichDo)))$r.squared};   cat(sprintf("Model +Struct: %f\n",max(sqrt(b))))
    
}
f()
#### Things aside from length:
#i = 70000; b=summary( a<-lm(coSI ~ + exp(Acceptor)*exp(-(Dist2End+addReadThrough)/i) + exp(-(Dist2End+addReadThrough)/i)*exp(Donor)+exp(-(Dist2End+addReadThrough)/i)*Stability + length + log(Dist2Start)+exp(-(Dist2End+addReadThrough)/i)*log(UpIntronLength)+exp(-(Dist2End+addReadThrough)/i)*log(DownIntronLength),
#  subset(exons4,isHK %in% whichDo)));print(c(sqrt(b$r.squared)));
#  
exons4$minA_up = apply(exons4[,c('Acceptor','UpIntronZ')],1, function(x) min(x[1],-x[2]))
exons4$minD_down = apply(exons4[,c('Donor','DownIntronZ')],1, function(x) min(x[1],-x[2]))

b=NULL; for(i in 10000 * (1:7)){b[i]=summary(
  a<-lm(coSI ~ exp(-(Dist2End+addReadThrough)/i) + (minD_down + minA_up) + Stability + length + log(Dist2Start),
  subset(exons4,isHK %in% whichDo)))};print(max(sqrt(b$r.squared)));

  
  ## Hmm.. play with data
  G=density(a$model[,1],from=0,to=1)
 summary(G$x)

hist(a$model[,1],100)->G
lines(x=G$mids,y<-dbeta(G$mids,5.5->alpha,B(alpha))*nrow(a$model)/100,col=4)


