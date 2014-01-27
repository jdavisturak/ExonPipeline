rm(list=ls())
library(gplots)
library(RColorBrewer)
library(Hmisc)
library(nls2)
source("/home/jeremy/ExonPipeline/ExonPipeline_functions.r")
source("/home/jeremy/ExonPipeline/ExonPipeline_analysis_and_plotting_Functions.r")
source("/home/home/jeremy/Code/useful_R/useful_R.r")
source("/home/jeremy/ExonPipeline/ExonPipeline_simulate_model.r")

setwd("/home/jeremy/ExonPipeline/hg19")
#load("PlottedData6.RData")
#load("Tilgner2_mergedData.RData")  
#load("Tilgner2_mergedData_withSplicing.RData")
#load("Tilgner2_mergedData_noPseudogenes.RData")  
#load("Tilgner2_mergedData_noPseudogenes_withSplicing.RData")

ChromInfo  = cbind(c("PolyAMinus.ssj", "PolyAPlus.ssj"),c('splicingAMinus','splicingAPlus'))
load("Tilgner2_mergedData_noPseudogenes_withSplicing.RData")
load("PlottedData_new1.RData") 
load("PlottedData_new2.RData") # I should fix this so it loads all the stuff from new1 as well


##########################################################
#### GET SSJ data  (Dmitri D. Pervouchine)
##########################################################
##
### My data is 0-indexed, their is not??
### also, i added a '1' b/c the way i ran bam2ssj, the 1 is added like *_1
##
### Read in bam2ssj file of Chromatin data for Smale  http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32916
##ChromInfo  = cbind(c("PolyAMinus.ssj", "PolyAPlus.ssj"),c('splicingAMinus','splicingAPlus'))
##
#splicingData = list()
#splicingData_Theta5 = list()
#splicingData_Theta3 = list()
#for(i in 1:nrow(ChromInfo)){
#  dat=read.delim(paste("/home/jeremy/RNAseq/Tilgner/Nucleus/",ChromInfo[i,1],sep=""), head=F,stringsAsFactors=F)
#  names(dat) = c('coSI_ID','count53', 'count5X', 'countX3', 'count50', 'count03')
#  print(dim(dat))
#  
#  # Compute Theta5, Theta3  (http://bioinformatics.oxfordjournals.org/content/early/2012/11/21/bioinformatics.bts678.full.pdf+html)
#  dat$Denom5 = dat$count53 + dat$count5X + dat$count50
#  dat$Theta5 = (dat$count53 + dat$count5X) / dat$Denom5
#  dat$Denom3 = dat$count53 + dat$countX3 + dat$count03
#  dat$Theta3 = (dat$count53 + dat$countX3) / dat$Denom3
#
#  # My way:
#  dat$DenomJer = dat$count53 + (dat$count50 + dat$count03)/2
#  dat$SpliceJer = dat$count53 / dat$DenomJer
#  
#  splicingData[[i]] = dat[dat$DenomJer  > 15,]
#  splicingData_Theta5[[i]] = dat[dat$Denom5  > 15,]
#  splicingData_Theta3[[i]] = dat[dat$Denom3  > 15,]  
#}                                                                                                                      cd Exo
#
##
##########################################################
#### Merge with SSJ data 
##########################################################
##
##
####
#mergedData$coSI_ID = paste(mergedData$chr,rowMins(mergedData[,c('start','end')])+1,rowMaxes(mergedData[,c('start','end')]),1,sep='_')
#mergedData$splicingAMinus = splicingData[[1]]$SpliceJer[match(mergedData$coSI_ID, splicingData[[1]]$coSI_ID)]
#mergedData$splicingAPlus = splicingData[[2]]$SpliceJer[match(mergedData$coSI_ID, splicingData[[2]]$coSI_ID)]
#
#mergedData$Theta5_minus = splicingData_Theta5[[1]]$Theta5[match(mergedData$coSI_ID, splicingData_Theta5[[1]]$coSI_ID)]
#mergedData$Theta3_minus = splicingData_Theta3[[2]]$Theta3[match(mergedData$coSI_ID, splicingData_Theta3[[2]]$coSI_ID)]
#mergedData$Theta5_minus = splicingData[[1]]$Theta5[match(mergedData$coSI_ID, splicingData[[1]]$coSI_ID)]
#mergedData$Theta3_minus = splicingData[[1]]$Theta3[match(mergedData$coSI_ID, splicingData[[1]]$coSI_ID)]
#introns3 = mergedData[apply(mergedData[,ChromInfo[,2]],1,function(x)any(!is.na(x))),]
#

#table(apply(introns3[,ChromInfo[,2]],1,function(x)length(which(!is.na(x)))))
## 1     2 
## 5893 24120 
#save(introns3,file="Tilgner2_mergedData_noPseudogenes_withSplicing.RData")


#####################################################################
## Model MEDIANS
#####################################################################

getQuantilMedians = function(myData,column='splicing0',N=100,col1='Dist2End'){
  myData$Group = quantileGroups(myData[,col1],N)
  data.frame(LengthMedian = tapply(myData[,col1],myData$Group,median,na.rm=T), median = tapply(myData[,column],myData$Group,median,na.rm=T))
}


getQuantilMedians2 = function(myData,N=100,col1='Dist2End', col2='coSI'){
  myData$Group = quantileGroups(myData[,col1],N)
  data.frame(GroupMedian = tapply(myData[,col1],myData$Group,median,na.rm=T), DataMedian = tapply(myData[,col2],myData$Group,median,na.rm=T))
}


##############################################################################
### Gamma model fits
##############################################################################
QuantileDataMinus  = getQuantilMedians(introns3,column=ChromInfo[1,2],100) 
QuantileDataPlus  = getQuantilMedians(introns3,column=ChromInfo[2,2],100) 
QuantileDataMinus_HK  = getQuantilMedians(subset(introns3,isHK==1),column=ChromInfo[1,2],100) 
QuantileDataMinus_non  = getQuantilMedians(subset(introns3, isHK==0),column=ChromInfo[1,2],100) 


#introns3_5prime_Distances = introns3; introns3_5prime_Distances$Dist2End = introns3$Dist2End + introns3$length # Make this the actual distance to the 5' end ...
QuantileData5  = getQuantilMedians(introns3,column='Theta5_minus',100) 
QuantileData3  = getQuantilMedians(introns3,column='Theta3_minus',100) 
QuantileDataMinus_HK5  = getQuantilMedians(subset(introns3,isHK==1),column='Theta5_minus',100) 
QuantileDataMinus_non5  = getQuantilMedians(subset(introns3, isHK==0),column='Theta5_minus',100) 
QuantileDataMinus_HK3  = getQuantilMedians(subset(introns3,isHK==1),column='Theta3_minus',100) 
QuantileDataMinus_non3  = getQuantilMedians(subset(introns3, isHK==0),column='Theta3_minus',100) 

QuantileDataMinus_start  = getQuantilMedians2(introns3,col1='Dist2Start',col2=ChromInfo[1,2],100) 
QuantileDataMinus_start2  = getQuantilMedians2(introns3,col1='Dist2Start',col2="Dist2End",100) 
## Dists to start and stop are highly correlated!  Let's see what the LM says
summary(lm(splicingAMinus ~ Dist2End * Dist2Start, data=introns3))

#QuantileDataMinus_long  = getQuantilMedians(subset(introns3,GeneLength > sizeCutoff),column=ChromInfo[1,2],100) 
#QuantileDataMinus_short  = getQuantilMedians(subset(introns3,GeneLength < sizeCutoff),column=ChromInfo[1,2],100)
QuantileDataMinus_long  = getQuantilMedians(subset(introns3,GeneSize ==4),column=ChromInfo[1,2],N=20) 
QuantileDataMinus_short  = getQuantilMedians(subset(introns3,GeneSize ==1),column=ChromInfo[1,2],N=20)
QuantileDataMinus_Medshort  = getQuantilMedians(subset(introns3,GeneSize ==2),column=ChromInfo[1,2],N=20)
QuantileDataMinus_Medlong  = getQuantilMedians(subset(introns3,GeneSize ==3),column=ChromInfo[1,2],N=20)

QuantileEnergyVsCoSI =   getQuantilMedians2(introns3,100,'Stability','splicingAMinus')
QuantileAcceptorVsCoSI =   getQuantilMedians2(introns3,100,'Acceptor','splicingAMinus')
QuantileDonorVsCoSI =   getQuantilMedians2(introns3,100,'Donor','splicingAMinus')

MeanGamma2 = function(LengthMedian,g,k,A,B)  sapply(LengthMedian, function(x)A * integrate(pgamma,0,x+B, shape=g, rate=k)$value/(x+B))

gMinus = nls2(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileDataMinus,start=c(g=1,k=5e-5,A=1,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,Inf,Inf),algo='port');
gMinus_noRT2 = nls2(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileDataMinus,start=c(g=1,k=5e-5,A=1,B=0),lower=c(1,0,0,0),upper=c(1,Inf,Inf,0),algo='port');


gPlus = nls(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileDataPlus,start=c(g=1,k=5e-5,A=.85,B=100),lower=c(1,0,0,1),upper=c(1,Inf,Inf,Inf),algo='port');
gMinus_HK = nls(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileDataMinus_HK,start=c(g=1,k=5e-5,A=.85,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,Inf,Inf),algo='port');
gMinus_non = nls(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileDataMinus_non,start=c(g=1,k=5e-5,A=.85,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,Inf,Inf),algo='port');

g5 = nls2(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileData5,start=c(g=1,k=5e-5,A=1,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,Inf,Inf),algo='port');
g3 = nls2(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileData3,start=c(g=1,k=5e-5,A=1,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,Inf,Inf),algo='port');
gMinus_HK5 = nls(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileDataMinus_HK5,start=c(g=1,k=5e-5,A=.85,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,Inf,Inf),algo='port');
gMinus_non5 = nls(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileDataMinus_non5,start=c(g=1,k=5e-5,A=.85,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,Inf,Inf),algo='port');
gMinus_HK3 = nls(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileDataMinus_HK3,start=c(g=1,k=5e-5,A=.85,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,Inf,Inf),algo='port');
gMinus_non3 = nls(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileDataMinus_non3,start=c(g=1,k=5e-5,A=.85,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,Inf,Inf),algo='port');


gMinus_noRT = nls(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileDataMinus[-1,],start=c(g=1,k=5e-5,A=1,B=0),lower=c(1,0,0,0),upper=c(1,Inf,Inf,0),algo='port');
gMinus_noRT = nls(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileDataMinus[-1,],start=c(g=1,k=5e-5,A=1,B=0),lower=c(1,0,0,0),upper=c(1,Inf,Inf,0),algo='port');

# Fit long/short genes, and fix the A value
#A_fixed = 0.9
#g_short = nls2(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileDataMinus_short,start=c(g=1,k=5e-5,A=A_fixed,B=1000),lower=c(1,0,A_fixed,1),upper=c(1,Inf,A_fixed,Inf),algo='port');
#g_Medshort = nls2(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileDataMinus_Medshort,start=c(g=1,k=5e-5,A=A_fixed,B=1000),lower=c(1,0,A_fixed,1),upper=c(1,Inf,A_fixed,Inf),algo='port');
#g_Medlong = nls2(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileDataMinus_Medlong,start=c(g=1,k=5e-5,A=A_fixed,B=1000),lower=c(1,0,A_fixed,1),upper=c(1,Inf,A_fixed,Inf),algo='port');
#g_long = nls2(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileDataMinus_long,start=c(g=1,k=5e-5,A=A_fixed,B=1000),lower=c(1,0,A_fixed,1),upper=c(1,Inf,A_fixed,Inf),algo='port');
g_fixed = 2
g_short = nls2(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileDataMinus_short,start=c(g=g_fixed,k=5e-5,A=1,B=3000),lower=c(g_fixed,0,0,1),upper=c(g_fixed,Inf,Inf,Inf),algo='port');
g_Medshort = nls2(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileDataMinus_Medshort,start=c(g=g_fixed,k=5e-5,A=1,B=3000),lower=c(g_fixed,0,0,1),upper=c(g_fixed,Inf,Inf,Inf),algo='port');
g_Medlong = nls2(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileDataMinus_Medlong,start=c(g=g_fixed,k=5e-5,A=1,B=3000),lower=c(g_fixed,0,0,1),upper=c(g_fixed,Inf,Inf,Inf),algo='port');
g_long = nls2(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileDataMinus_long,start=c(g=g_fixed,k=5e-5,A=1,B=3000),lower=c(g_fixed,0,0,1),upper=c(g_fixed,Inf,Inf,Inf),algo='port');




summary(gMinus)

gammaMinus_k_splice = c(All=coef(gMinus)[2]*3000, HK=coef(gMinus_HK)[2]*3000, non=coef(gMinus_non)[2]*3000)
gammaMinus_k_splice
#    All.k      HK.k     non.k 
#1.0085132 0.7618832 1.1135246 

gammaMinus_Dist_ReadThru = c(All=coef(gMinus)[4], HK=coef(gMinus_HK)[4], non=coef(gMinus_non)[4])
gammaMinus_Dist_ReadThru
#   All.B     HK.B    non.B 
# 4338.947 7856.583 3373.990

## OK, so this is actually pretty cool. 
# When I fix k at 1, it makes exactly 1 min-1 for splicing non-HK,!!!, .78 min-1 for HK!

gammaFits = list(splicingAMinus=coef(gMinus),splicingAPlus=coef(gPlus),splicingAMinus_noRT = coef(gMinus_noRT))
save(gammaFits,getQuantilMedians, getQuantilMedians2,QuantileDataMinus  ,MeanGamma2,gMinus, file='Gamma2_fits.RData')

fitSplit = nls2(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=rbind(QuantileDataMinus_short,QuantileDataMinus_long,QuantileDataMinus_Medshort,QuantileDataMinus_Medlong),start=c(g=1,k=5e-5,A=1,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,Inf,Inf),algo='port');



###################################################################################
###                      Can I model energy and acceptor along with Length?     ###
###################################################################################
column = 'splicingAMinus'

#################### Fit grouped into 1000 bins ###################################
QuantileDataMinus_A1000  = getQuantilMedians2(introns3,col1='Dist2End',col2="Acceptor",1000) 
QuantileDataMinus_E1000  = getQuantilMedians2(introns3,col1='Dist2End',col2="Stability",1000) 
QuantileDataMinus_splice1000  = getQuantilMedians2(introns3,col1='Dist2End',col2="splicingAMinus",1000) 
QuantileDataMinus_splicePlus1000  = getQuantilMedians2(introns3,col1='Dist2End',col2="splicingAPlus",1000) 

QuantileData1000 = data.frame(splicingAMinus=QuantileDataMinus_splice1000$DataMedian, Dist2End=QuantileDataMinus_splice1000$GroupMedian, Acceptor=QuantileDataMinus_A1000$DataMedian, Stability=QuantileDataMinus_E1000$DataMedian)

### Use 1000 bins
###Model Length to end
gMinus_1000L = nls(splicingAMinus ~ MeanGamma2(Dist2End,g,k,A,B), data=QuantileData1000,start=c(g=1,k=5e-5,A=1,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,Inf,Inf),algo='port');
coef(gMinus_1000L) -> c1000L
gMinus_1000Plus = nls(DataMedian ~ MeanGamma2(GroupMedian,g,k,A,B), data=QuantileDataMinus_splicePlus1000,start=c(g=1,k=5e-5,A=1,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,Inf,Inf),algo='port');
coef(gMinus_1000Plus) -> c1000P

## Now also model it WITHOUT the read-through
gMinus_1000L_noRT = nls(splicingAMinus ~ MeanGamma2(Dist2End,g,k,A,B), data=QuantileData1000,start=c(g=1,k=5e-5,A=1,B=0),lower=c(1,0,1,0),upper=c(1,Inf,1,0),algo='port');
coef(gMinus_1000L_noRT) -> c1000L_noRT
gMinus_1000L_noRT2 = nls(splicingAMinus ~ MeanGamma2(Dist2End,g,k,A,B), data=QuantileData1000,start=c(g=1,k=5e-5,A=1,B=0),lower=c(1,0,0,0),upper=c(1,Inf,Inf,0),algo='port');
coef(gMinus_1000L_noRT2) -> c1000L_noRT2
coef(gMinus)->c100L
coef(gMinus_noRT2)->c100L_noRT2

save(QuantileData1000,gMinus_1000L,c1000L,gMinus_1000Plus,c1000P, gMinus_1000L_noRT, c1000L_noRT,gMinus_1000L_noRT2,c1000L_noRT2,gMinus_noRT2, c100L_noRT2,c100L,gMinus, file='Tilgner2_1000bin_fits.RData')



c1000L_noRT2
#          g           k           A           B 
#1.000000000 0.001621073 0.763746059 0.000000000 
c1000L_noRT2[2]*3000
#       k 
#4.86321 for simulating

## Here I am making the rate be a linear function of k + C*acceptor_score + D*donor_score
#MeanGamma4 = function(myData,g,k,A,B,C1,C2)    apply(myData, 1,
#    function(x)A * integrate(pgamma,0,(x[1]+B), shape=g, rate=(k + C1*x[2]+ C2*x[3])/3000)$value/(x[1]+B)   )
#
### Include acceptor
#MeanGamma3 = function(myData,g,k,A,B,sigma)    apply(myData, 1,
#    function(x)A * integrate(pgamma,0,(x[1]+B), shape=g, rate=k*qlnorm(pnorm(x[2]),0,sigma))$value/(x[1]+B)   )
#
# Here I am making the rate be a linear function of k + C*acceptor score
#MeanGamma3 = function(myData,g,k,A,B,sigma)    apply(myData, 1,
#    function(x)A * integrate(pgamma,0,(x[1]+B), shape=g, rate=(k + sigma*x[2])/3000)$value/(x[1]+B)   )
#gMinus_1000A = nls(splicingAMinus ~ MeanGamma3(cbind(Dist2End,Acceptor),g,k,A,B,sigma), data=QuantileData1000,start=c(g=1,k=.8,A=.91,B=1000,sigma=.0001),lower=c(1,0,0,1,0),upper=c(1,Inf,Inf,Inf,5),algo='port');
#coef(gMinus_1000A) -> c1000A 
#cor(predict(gMinus_1000A,newdata=subset(introns3, !is.na(Acceptor))), subset(introns3, !is.na(Acceptor),'splicingAMinus'),use='c')
#  0.4201027



###################################################################################
###                      PLOTS                    ###
###################################################################################

### Plot Fits of 1000 bin data #
#png(sprintf("Tilgner2_Median1000IntronSplice.png"),height=400,width=400)
pdf(sprintf("Tilgner2_Median1000IntronSplice3.pdf"),height=3,width=3)
par(mai=c(.45,.45,0.24,0.1),cex=0.6, cex.main=1)
plot(log10(QuantileData1000$Dist2End), QuantileData1000$splicingAMinus,ylim=c(0,1), xlim=log10(c(100,420000)), ylab='Median psi-intron', xlab='LOG10 Distance to End of Gene' , main=paste('Nulear poly(A) Minus RNA'),pch=20)
#lines(log10(xx<- 10^seq((0),log10(500000),length.out=1000)),MeanGamma2(xx,c1000L_noRT[1],c1000L_noRT[2],c1000L_noRT[3], c1000L_noRT[4]),type='l',col='red',lwd=2,lty=2)
lines(log10(xx<- 10^seq((0),log10(500000),length.out=1000)),MeanGamma2(xx,c1000L[1],c1000L[2],c1000L[3], c1000L[4]),type='l',col='red',lwd=2)
lines(log10(xx<- 10^seq((0),log10(500000),length.out=1000)),MeanGamma2(xx,c1000L_noRT2[1],c1000L_noRT2[2],c1000L_noRT2[3], c1000L_noRT2[4]),type='l',col='red',lwd=2,lty=2)
legend('topleft',col='red',lwd=2,lty=c(2,1), leg=c('No read-Through','With Read-Through'))
dev.off()

pdf(sprintf("Tilgner2_Median1000IntronSplice_noRT.pdf"),height=3,width=3)
par(mai=c(.45,.45,0.24,0.1),cex=0.6, cex.main=1)
plot(log10(QuantileData1000$Dist2End), QuantileData1000$splicingAMinus,ylim=c(0,1), xlim=log10(c(100,420000)), ylab='Median psi-intron', xlab='LOG10 Distance to End of Gene' , main=paste('Nulear poly(A) Minus RNA'),pch=20)
lines(log10(xx<- 10^seq((0),log10(500000),length.out=1000)),MeanGamma2(xx,c1000L_noRT2[1],c1000L_noRT2[2],c1000L_noRT2[3], c1000L_noRT2[4]),type='l',col='red',lwd=2,lty=2)
dev.off()

#pdf(sprintf("Tilgner2_Median100IntronSplice_noRT.pdf"),height=3,width=3)
#par(mai=c(.45,.45,0.24,0.1),cex=0.6, cex.main=1)
#plot(log10(QuantileDataMinus$LengthMedian), QuantileDataMinus$median,ylim=c(0,1), xlim=log10(c(100,420000)), ylab='Median psi-intron', xlab='LOG10 Distance to End of Gene' , main=paste('Nulear poly(A) Minus RNA'),pch=20)
#lines(log10(xx<- 10^seq((0),log10(500000),length.out=1000)),MeanGamma2(xx,c100L_noRT2[1],c100L_noRT2[2],c100L_noRT2[3], c100L_noRT2[4]),type='l',col='red',lwd=2,lty=2)
#lines(log10(xx<- 10^seq((0),log10(500000),length.out=1000)),MeanGamma2(xx,c100L[1],c100L[2],c100L[3], c100L[4]),type='l',col='red',lwd=2)
#dev.off()


#png(sprintf("Tilgner2_Median1000IntronSplice+.png"),height=400,width=400)
pdf(sprintf("Tilgner2_Median1000IntronSplice+.pdf"),height=3,width=3)
par(mai=c(.45,.45,0.24,0.1),cex=0.6, cex.main=1)
plot(log10(QuantileDataMinus_splicePlus1000$GroupMedian), QuantileDataMinus_splicePlus1000$DataMedian,ylim=c(0,1),, xlim=log10(c(100,420000)), ylab='Median psi-intron', xlab='LOG10 Distance to End of Gene' , main=paste('Nulear poly(A) Plus RNA'),pch=20)
lines(log10(xx<- 10^seq((0),log10(500000),length.out=1000)),MeanGamma2(xx,c1000P[1],c1000P[2],c1000P[3], c1000P[4]),type='l',col='red',lwd=3)
dev.off()


pdf(sprintf("Tilgner2_Medians_100_.pdf"),height=3,width=12)
par(mfrow=c(1,4))
plot(log10(QuantileDataPlus$LengthMedian), QuantileDataPlus$median,ylim=c(0.3,1), xlim=log10(c(100,420000)), ylab=expression(paste('Median ', phi)), xlab='LOG10 Distance to End of Gene' , main=paste('Nulear poly(A) Plus RNA'),pch=20)
plot(log10(QuantileDataMinus$LengthMedian), QuantileDataMinus$median,ylim=c(0.3,1), xlim=log10(c(100,420000)), ylab=expression(paste('Median ', phi)), xlab='LOG10 Distance to End of Gene' , main=paste('Nulear poly(A) Minus RNA'),pch=20)
lines(log10(xx<- 10^seq((0),log10(500000),length.out=1000)),MeanGamma2(xx,c100L[1],c100L[2],c100L[3], c100L[4]),type='l',col='red',lwd=2)
lines(log10(xx<- 10^seq((0),log10(500000),length.out=1000)),MeanGamma2(xx,c100L_noRT2[1],c100L_noRT2[2],c100L_noRT2[3], c100L_noRT2[4]),type='l',col='red',lwd=2,lty=2)
legend('topleft',col='red',lwd=2,lty=c(2,1), leg=c('No read-Through','With Read-Through'))
plot((QuantileDonorVsCoSI$GroupMedian), QuantileDonorVsCoSI$DataMedian,ylim=c(0.3,1), ylab=expression(paste('Median ', phi)), xlab="Splicing Donor Sequence Strength" , main=paste('Nulear poly(A) Minus RNA'),pch=20)
plot((QuantileAcceptorVsCoSI$GroupMedian), QuantileAcceptorVsCoSI$DataMedian,ylim=c(0.3,1), ylab=expression(paste('Median ', phi)), xlab="Splicing Acceptor Sequence Strength" , main=paste('Nulear poly(A) Minus RNA'),pch=20)
dev.off()


png('Tilgner2_median_LongShort2.png',height=400,width=400)

plot(log10(QuantileDataMinus_short$LengthMedian), xlim=log10(c(100,420000)), QuantileDataMinus_short$median, ylab='Median psi-intron', xlab='LOG10 Distance to End of Gene' ,pch=20,ylim=c(0.2,1),col=myCols[1])
points(log10(QuantileDataMinus_long$LengthMedian),  QuantileDataMinus_long$median,pch=20, col=myCols[4])
points(log10(QuantileDataMinus_Medlong$LengthMedian),  QuantileDataMinus_Medlong$median,pch=20, col=myCols[3])
points(log10(QuantileDataMinus_Medshort$LengthMedian),  QuantileDataMinus_Medshort$median,pch=20, col=myCols[2])
lines(log10(xx<- 10^seq((0),log10(500000),length.out=1000)),MeanGamma2(xx,coef(gMinus)[1], coef(gMinus)[2],coef(gMinus)[3], coef(gMinus)[4]),type='l',col='black',lwd=1) 
lines(log10(xx<- 10^seq((0),log10(500000),length.out=1000)),MeanGamma2(xx,coef(g_short)[1], coef(g_short)[2],coef(g_short)[3], coef(g_short)[4]),type='l',col=myCols[1],lwd=1) 
lines(log10(xx<- 10^seq((0),log10(500000),length.out=1000)),MeanGamma2(xx,coef(g_Medshort)[1], coef(g_Medshort)[2],coef(g_Medshort)[3], coef(g_Medshort)[4]),type='l',col=myCols[2],lwd=1) 
lines(log10(xx<- 10^seq((0),log10(500000),length.out=1000)),MeanGamma2(xx,coef(g_Medlong)[1], coef(g_Medlong)[2],coef(g_Medlong)[3], coef(g_Medlong)[4]),type='l',col=myCols[3],lwd=1) 
lines(log10(xx<- 10^seq((0),log10(500000),length.out=1000)),MeanGamma2(xx,coef(g_long)[1], coef(g_long)[2],coef(g_long)[3], coef(g_long)[4]),type='l',col=myCols[4],lwd=1) 

legend('bottomright',col=myCols,pch=20,leg=c('short genes','medium short','medium long','long genes'))
dev.off()

pdf('Tilgner2_median_HKfits.pdf')
plot(log10(QuantileDataMinus_HK$LengthMedian), xlim=log10(c(100,420000)), QuantileDataMinus_HK$median, ylab='Median psi-intron', xlab='LOG10 Distance to End of Gene' ,pch=20,ylim=c(0.4,1),col='red')
points(log10(QuantileDataMinus_non$LengthMedian),  QuantileDataMinus_non$median,pch=20, col='black')
lines(log10(xx<- 10^seq((0),log10(500000),length.out=1000)),MeanGamma2(xx,coef(gMinus_HK)[1], coef(gMinus_HK)[2],coef(gMinus_HK)[3], coef(gMinus_HK)[4]),type='l',col='red',lwd=1) 
lines(log10(xx<- 10^seq((0),log10(500000),length.out=1000)),MeanGamma2(xx,coef(gMinus_non)[1], coef(gMinus_non)[2],coef(gMinus_non)[3], coef(gMinus_non)[4]),type='l',col='black',lwd=1) 
legend('topright',col=c(2,1),lwd=1,leg=c('HK','non'))
dev.off()

pdf(sprintf('%s_HK_dataMining.pdf',settings$CommonName),height=4, width=8)
par(cex=0.5, mai=c(0.6,1,.2,.2), mfrow=c(1,2))
boxplot(Stability ~ isHK, mergedData, notch=T,names=c("Regulated","HK"), border=c(1,2),main='',ylab=sprintf("Stability Score",settings$CommonName))
boxplot(exp(Acceptor) ~ isHK, mergedData, notch=T,names=c("Regulated","HK"), border=c(1,2),main='',ylab=sprintf("Acceptor Score",settings$CommonName))

dev.off()


















HK = c(1,0); HKname = 'All'
#HK = c(1); HKname = 'HK'


## Plot crazy exons graph ...
cols=brewer.pal(11,'Spectral')


pdf(sprintf("MedianIntronSplice_by_remaining_exons_%s_withLines.pdf",HKname), height=8, width=12)
par(mfrow=c(2,3))
ylims = data.frame(c(.45,.89),c(.73,.95)); names(ylims) = ChromInfo[,2]
for(column in ChromInfo[,2]){
  
  ### Plot data split into 100 groups, with exponential fit  ###
  QuantileData  = getQuantilMedians(subset(introns3,isHK %in% HK),column=column,100) 
  plot(log10(QuantileData$LengthMedian), xlim=log10(c(500,420000)), QuantileData$median, ylab='Median psi-intron', xlab='LOG10 Distance to End of Gene' , main=paste('Chromatin',column),pch=20)
  lines(log10(xx<- 10^seq((0),log10(500000),length.out=1000)),MeanGamma2(xx,gammaFits[[column]][1], gammaFits[[column]][2],gammaFits[[column]][3], gammaFits[[column]][4]),type='l',col='black')

  
  #QuantileData  = getQuantilMedians(subset(introns3,isHK %in% HK),column=column,20) 
  #pdf(sprintf("Median_%s_by_remaining_introns_%s_withLines.pdf",column,HKname))
  #plot(log10(QuantileData$LengthMedian), xlim=log10(c(500,420000)), ylim=ylims[,column], QuantileData$median, ylab='Median coSI', xlab='LOG10 Distance to End of Gene' , main=paste('Chromatin',column),cex=.8)
  
  ### Plot data split into 5 groups of # remaining exons, plus again 20 groups per.
  # start with fit line
  plot(log10(xx), MeanGamma2(xx,gammaFits[[column]][1], gammaFits[[column]][2],gammaFits[[column]][3], gammaFits[[column]][4]), xlim=log10(c(500,420000)),type='l', ylim=ylims[,column],ylab='Median coSI', xlab='LOG10 Distance to End of Gene' , main=paste('Chromatin',column))

  # Save residuals
  allFits = list()
  MedLengths = list()
  MedCoSI = list()
  #myGroups = list(1:2,NA,3:4,NA,5:6,NA,7:9,NA,10); AA<-seq(1,10,2); leg=c("1-2","3-4","5-6","7-9",">=10")
  myGroups = as.list(1:10); AA <-1:9;leg=c(1:8,">=9")
  Residuals = list()

  S=sapply(AA, function(i){
       if(i<9)
        QD = getQuantilMedians(subset(subset(introns3,isHK %in% HK),intronsRemaining %in% myGroups[[i]]),column,10)
       else
        QD = getQuantilMedians(subset(subset(introns3,isHK %in% HK),intronsRemaining >=myGroups[[i]]),column,10)      
       MedLengths[[i]] <<- QD[,1]
       MedCoSI[[i]] <<- QD[,2]      
       FIT <- lm(median ~ log10(LengthMedian),QD)
       allFits[[i]]<<- FIT
       abline(a=coef(FIT)[1], b=coef(FIT)[2],col=cols[i],lty=2)
       points(log10(QD[,1]),QD[,2],col=cols[i],pch=20,cex=1.2)  
       Residuals[[i]] <<- QD$median - MeanGamma2(QD$LengthMedian,gammaFits[[column]][1], gammaFits[[column]][2],gammaFits[[column]][3], gammaFits[[column]][4])
                
  })
  print(range(sapply(MedCoSI[AA],range)))
  legend('bottomright',col=cols[AA],pch=20,leg=leg)
  
  ### Plot 
  boxplot(Residuals[AA], main='Residuals of median data for intron-remaining groups',names=leg, xlab='# of exons downstream',notch=T,col=cols[AA])
  abline(h=0)


}

dev.off()


################# Splicing compare to distance from STARTL: ##############
plot(log10(QuantileDataMinus_start$GroupMedian), xlim=log10(c(500,420000)), QuantileDataMinus_start$DataMedian, ylab='Median psi-intron', xlab='LOG10 Distance to START of Gene' , main='A-minus: dist to Start',pch=20)
  

###################################################################################
###                     Show fitting                                  ####
###################################################################################
 

## Demonstrate what happends to fit curve when I change parameters

png('Tilgner2_parameter_demonsrtration.png')
plot(log10(QuantileDataMinus$LengthMedian), xlim=log10(c(10,420000)), QuantileDataMinus$median, ylab='Median psi-intron', xlab='LOG10 Distance to End of Gene' , main='Effect of varying the 3 parameters',pch=20,ylim=c(0,1))
lines(log10(xx<- 10^seq((0),log10(500000),length.out=1000)),MeanGamma2(xx,gammaFits[['splicingAMinus']][1], gammaFits[['splicingAMinus']][2],gammaFits[['splicingAMinus']][3], gammaFits[['splicingAMinus']][4]),type='l',col='black',lwd=2)
lines(log10(xx<- 10^seq((0),log10(500000),length.out=1000)),MeanGamma2(xx,gammaFits[['splicingAMinus']][1], gammaFits[['splicingAMinus']][2],gammaFits[['splicingAMinus']][3], .5*gammaFits[['splicingAMinus']][4]) ,type='l',col='red')
lines(log10(xx<- 10^seq((0),log10(500000),length.out=1000)),MeanGamma2(xx,gammaFits[['splicingAMinus']][1], gammaFits[['splicingAMinus']][2],gammaFits[['splicingAMinus']][3], 2*gammaFits[['splicingAMinus']][4]) ,type='l',col='red')

lines(log10(xx<- 10^seq((0),log10(500000),length.out=1000)),MeanGamma2(xx,gammaFits[['splicingAMinus']][1], gammaFits[['splicingAMinus']][2],1, gammaFits[['splicingAMinus']][4]),type='l',col='green')
lines(log10(xx<- 10^seq((0),log10(500000),length.out=1000)),MeanGamma2(xx,gammaFits[['splicingAMinus']][1], .5*gammaFits[['splicingAMinus']][2],gammaFits[['splicingAMinus']][3], gammaFits[['splicingAMinus']][4]),type='l',col='blue',lty=2)
lines(log10(xx<- 10^seq((0),log10(500000),length.out=1000)),MeanGamma2(xx,gammaFits[['splicingAMinus']][1], 2*gammaFits[['splicingAMinus']][2],gammaFits[['splicingAMinus']][3], gammaFits[['splicingAMinus']][4]),type='l',col='blue')

legend('topleft',lwd=1,col=c('green','red','blue'), leg=c('Saturation level','Read-through Distance', 'k-splicing'))
dev.off()
 
 
## Show improvement over linear model ...
fit1 = lm(median ~ log10(LengthMedian),QuantileDataMinus[-1,])
sum(resid(gMinus)^2)
#[1]  0.04026385
sum(fit1$residuals^2)
#[1]  0.05893246

sum(resid(gMinus)^2)/ sum(fit1$residuals^2)
#[1]  0.6832202 My fit has only 68% of the residuals as his does

png('Tilgner2_parameter_demonstration2_%d.png')

for (i in 1:3){

  plot(log10(QuantileDataMinus$LengthMedian), xlim=log10(c(100,420000)), QuantileDataMinus$median, ylab='Median psi-intron', xlab='LOG10 Distance to End of Gene' ,pch=20,ylim=c(0.4,1))
  if(i >= 2)
    lines(log10(xx), predict(fit1,newdata=data.frame(LengthMedian=xx)),lty=2)
    
  if(i >= 3)
    lines(log10(xx<- 10^seq((0),log10(500000),length.out=1000)),MeanGamma2(xx,gammaFits[['splicingAMinus']][1], gammaFits[['splicingAMinus']][2],gammaFits[['splicingAMinus']][3], gammaFits[['splicingAMinus']][4]),type='l',col='black',lwd=2) 
#  if(i >= 4)
#     lines(log10(xx<- 10^seq((1),log10(500000),length.out=1000)),MeanGamma2(xx,gammaFits[['splicingAMinus_noRT']][1], gammaFits[['splicingAMinus_noRT']][2],gammaFits[['splicingAMinus_noRT']][3], gammaFits[['splicingAMinus_noRT']][4]),type='l',col='black',lwd=1)
#
}
dev.off()


png('Tilgner2_median_HKfits.png',height=400,width=400)

plot(log10(QuantileDataMinus_HK$LengthMedian), xlim=log10(c(100,420000)), QuantileDataMinus_HK$median, ylab='Median psi-intron', xlab='LOG10 Distance to End of Gene' ,pch=20,ylim=c(0.4,1),col='red')
points(log10(QuantileDataMinus_non$LengthMedian),  QuantileDataMinus_non$median,pch=20, col='black')
lines(log10(xx<- 10^seq((0),log10(500000),length.out=1000)),MeanGamma2(xx,coef(gMinus_HK)[1], coef(gMinus_HK)[2],coef(gMinus_HK)[3], coef(gMinus_HK)[4]),type='l',col='red',lwd=1) 
lines(log10(xx<- 10^seq((0),log10(500000),length.out=1000)),MeanGamma2(xx,coef(gMinus_non)[1], coef(gMinus_non)[2],coef(gMinus_non)[3], coef(gMinus_non)[4]),type='l',col='black',lwd=1) 
legend('topright',col=c(2,1),lwd=1,leg=c('HK','non'))
dev.off()


###################################### 
## Plot pretty picture of fitting
png('Tilgner2_simpleFitting.png')
plot(log10(QuantileDataMinus$LengthMedian), xlim=log10(c(100,420000)), QuantileDataMinus$median, ylab='Median psi-intron', xlab='LOG10 Distance to End of Gene' , main='Median intron splicing across 100 bins of introns\nGrouped by distance to poly(A) site',pch=20,ylim=c(0,1))
lines(log10(xx<- 10^seq((0),log10(500000),length.out=1000)),MeanGamma2(xx,gammaFits[['splicingAMinus']][1], gammaFits[['splicingAMinus']][2],gammaFits[['splicingAMinus']][3], gammaFits[['splicingAMinus']][4]),type='l',lwd=2,col='blue')
lines(log10(xx), predict(fit1,newdata=data.frame(LengthMedian=xx)),col='black')
abline(h=gammaFits[['splicingAMinus']][3],col='blue',lty=2)
legend('bottomright',col=c('black','blue','blue'),lty=c(1,1,2), lwd=c(1,2,1),leg=c('Linear fit to log(Distance)','Fit with splicing model','Fitted signal saturation level'))
dev.off()

png('Tilgner2_Theta_fitting.png')
par(mfrow=c(2,2))
plot(log10(QuantileData5$LengthMedian), xlim=log10(c(100,420000)), QuantileData5$median, ylab=expression(paste('Median ',theta,'5',sep='')), xlab='LOG10 Distance to End of Gene' , main=expression(paste('Median ',theta,'5 across 100 bins of introns',sep='')),pch=20,ylim=c(0,1))
lines(log10(xx<- 10^seq((0),log10(500000),length.out=1000)),MeanGamma2(xx,coef(g5)[1], coef(g5)[2],coef(g5)[3], coef(g5)[4]),type='l',lwd=2,col='blue')

plot(log10(QuantileData3$LengthMedian), xlim=log10(c(100,420000)), QuantileData3$median, ylab=expression(paste('Median ',theta,'3',sep='')), xlab='LOG10 Distance to End of Gene' , main=expression(paste('Median ',theta,'3 across 100 bins of introns',sep='')),pch=20,ylim=c(0,1))
lines(log10(xx<- 10^seq((0),log10(500000),length.out=1000)),MeanGamma2(xx,coef(g3)[1], coef(g3)[2],coef(g3)[3], coef(g3)[4]),type='l',lwd=2,col='blue')

plot(log10(QuantileDataMinus_HK5$LengthMedian), xlim=log10(c(100,420000)), QuantileDataMinus_HK5$median, ylab=expression(paste('Median ',theta,'5',sep='')), xlab='LOG10 Distance to End of Gene' ,pch=20,ylim=c(0.4,1),col='red')
points(log10(QuantileDataMinus_non5$LengthMedian),  QuantileDataMinus_non5$median,pch=20, col='black')
lines(log10(xx<- 10^seq((0),log10(500000),length.out=1000)),MeanGamma2(xx,coef(gMinus_HK5)[1], coef(gMinus_HK5)[2],coef(gMinus_HK5)[3], coef(gMinus_HK5)[4]),type='l',col='red',lwd=1) 
lines(log10(xx<- 10^seq((0),log10(500000),length.out=1000)),MeanGamma2(xx,coef(gMinus_non5)[1], coef(gMinus_non5)[2],coef(gMinus_non5)[3], coef(gMinus_non5)[4]),type='l',col='black',lwd=1) 
legend('bottomright',col=c(2,1),lwd=1,leg=c('HK','non'))

plot(log10(QuantileDataMinus_HK3$LengthMedian), xlim=log10(c(100,420000)), QuantileDataMinus_HK3$median, ylab=expression(paste('Median ',theta,'3',sep='')), xlab='LOG10 Distance to End of Gene' ,pch=20,ylim=c(0.4,1),col='red')
points(log10(QuantileDataMinus_non3$LengthMedian),  QuantileDataMinus_non3$median,pch=20, col='black')
lines(log10(xx<- 10^seq((0),log10(500000),length.out=1000)),MeanGamma2(xx,coef(gMinus_HK3)[1], coef(gMinus_HK3)[2],coef(gMinus_HK3)[3], coef(gMinus_HK3)[4]),type='l',col='red',lwd=1) 
lines(log10(xx<- 10^seq((0),log10(500000),length.out=1000)),MeanGamma2(xx,coef(gMinus_non3)[1], coef(gMinus_non3)[2],coef(gMinus_non3)[3], coef(gMinus_non3)[4]),type='l',col='black',lwd=1) 
legend('bottomright',col=c(2,1),lwd=1,leg=c('HK','non'))


dev.off()




###################################### 
## Examine numbers for fitting
for(column in ChromInfo[,2]){

	print(cor(introns3[,column], MeanGamma2(introns3$Dist2End,gammaFits[[column]][1], gammaFits[[column]][2],gammaFits[[column]][3], gammaFits[[column]][4]), use='c'))
}
#[1] 0.4203467
#[1] 0.3034644


fit2 = lm(splicingAMinus ~ log10(Dist2End+1),introns3)
 
cor(predict(fit2), fit2$model$splicingAMinus)
# [1]  0.4176859


###################################################################################
###                      PLOTS for ENERGY /ACCEPTOR                             ###
###################################################################################
 
 
###################  1: Compare to splicing measurements    ##################### 
 
pdf("Tilgner2_Median_spliceMinus_by_energy-acceptor.pdf")
par(mfrow=c(2,2))

### First plot the coSI vs Energy ###
QuantileEnergyVsCoSI =   getQuantilMedians2(introns3,100,'Stability','splicingAMinus')
EnergyFit = lm(DataMedian ~ GroupMedian,QuantileEnergyVsCoSI)
plot(QuantileEnergyVsCoSI, xlab="Stability Score", ylab="Median splicingAMinus",ylim=c(.4, 0.96))
lines(EE <- seq(-3,3,.05), predict(EnergyFit,data.frame(GroupMedian=EE))) 

## Next plot for different introns remaining ##
plot(EE, predict(EnergyFit,data.frame(GroupMedian=EE)), type='l',ylim=c(.4, 0.96), xlab="Stability Score", ylab="Median coSI" , main='Median coSI index for # of remaining introns')
cols=brewer.pal(10,'Spectral')
allFits = list()
Residuals = list()
MedLengths = list()
MedCoSI = list()
S=sapply(1:10, function(i){
     if(i<10)
      QD = getQuantilMedians2(subset(introns3,intronsRemaining==i),20,'Stability','splicingAMinus')
     else
      QD = getQuantilMedians2(subset(introns3,intronsRemaining>=i),20,'Stability','splicingAMinus')      
     MedLengths[[i]] <<- QD[,1]
     MedCoSI[[i]] <<- QD[,2]      
     FIT <- lm(DataMedian ~ GroupMedian,QD)
     abline(a=coef(FIT)[1], b=coef(FIT)[2],col=cols[i],lty=2)
     points(QD[,1],QD[,2],col=cols[i])    
     
     Residuals[[i]] <<- QD$DataMedian - predict(EnergyFit,QD)      
})
legend('bottomright',col=cols,pch=20,leg=c(1:9,'>=10'))


QuantileAcceptorVsCoSI =   getQuantilMedians2(introns3,100,'Acceptor','splicingAMinus')
AcceptorFit = lm(DataMedian ~ GroupMedian,QuantileAcceptorVsCoSI)
plot(QuantileAcceptorVsCoSI, xlab="Stability Score", ylab="Median splicingAMinus",ylim=c(.4, 0.96))
lines(EE <- seq(-3,3,.05), predict(AcceptorFit,data.frame(GroupMedian=EE))) 

## Next plot for different introns remaining ##
plot(EE, predict(AcceptorFit,data.frame(GroupMedian=EE)), type='l',ylim=c(.4, 0.96), xlab="Acceptor Score", ylab="Median coSI" , main='Median coSI index for # of remaining introns')
cols=brewer.pal(10,'Spectral')
allFits = list()
Residuals = list()
MedLengths = list()
MedCoSI = list()
S=sapply(1:10, function(i){
     if(i<10)
      QD = getQuantilMedians2(subset(introns3,intronsRemaining==i),20,'Acceptor','splicingAMinus')
     else
      QD = getQuantilMedians2(subset(introns3,intronsRemaining>=i),20,'Acceptor','splicingAMinus')      
     MedLengths[[i]] <<- QD[,1]
     MedCoSI[[i]] <<- QD[,2]      
     FIT <- lm(DataMedian ~ GroupMedian,QD)
     abline(a=coef(FIT)[1], b=coef(FIT)[2],col=cols[i],lty=2)
     points(QD[,1],QD[,2],col=cols[i])    
     
     Residuals[[i]] <<- QD$DataMedian - predict(AcceptorFit,QD)      
})
legend('bottomright',col=cols,pch=20,leg=c(1:9,'>=10'))
dev.off()


################## 2: Compared to Length       ###################################


getQuantilMediansEnergy = function(myData,N=100){
  myData$Group = splitByVector(myData$Dist2End,quantile(myData$Dist2End,na.rm=T,(1:(N-1))/N))
  data.frame(LengthMedian = tapply(myData$Dist2End,myData$Group,median), Energy_median = tapply(myData$Stability,myData$Group,median))
}
getQuantilMediansAcceptor = function(myData,N=100){
  myData$Group = splitByVector(myData$Dist2End,quantile(myData$Dist2End,na.rm=T,(1:(N-1))/N))
  data.frame(LengthMedian = tapply(myData$Dist2End,myData$Group,median,na.rm=T), Acceptor_median = tapply(myData$Acceptor,myData$Group,median,na.rm=T))
}


### Now look at the ENERGY of exons:

refseq$GeneLength = refseq$txEnd - refseq$txStart
qr = quantile(refseq$GeneLength,c(1:3)/4)
sizeCutoff = qr[1]; sizeCutoff2 = qr[2]; sizeCutoff3 = qr[3]
#QuantileEnergyData = getQuantilMediansEnergy(mergedData,20)
QuantileEnergyData = getQuantilMediansEnergy(subset(mergedData,isLast==0,20)
QuantileEnergyData_last = getQuantilMediansEnergy(subset(mergedData,isLast==1 ),20)

nums = 10
QuantileEnergyData2 = getQuantilMediansEnergy(subset(mergedData,GeneLength < sizeCutoff & isLast==0),nums)
QuantileEnergyData3 = getQuantilMediansEnergy(subset(mergedData,GeneLength > sizeCutoff& isLast==0),nums)
QuantileEnergyData4 = getQuantilMediansEnergy(subset(mergedData,GeneLength >= sizeCutoff & GeneLength < sizeCutoff2& isLast==0),nums)
QuantileEnergyData5 = getQuantilMediansEnergy(subset(mergedData,GeneLength >= sizeCutoff2 & GeneLength < sizeCutoff3& isLast==0),nums)
QuantileEnergyData6 = getQuantilMediansEnergy(subset(mergedData,GeneLength >= sizeCutoff3 & isLast==0),nums)

QuantileEnergyData_last2 = getQuantilMediansEnergy(subset(mergedData,GeneLength < sizeCutoff & isLast==1),nums)
QuantileEnergyData_last3 = getQuantilMediansEnergy(subset(mergedData,GeneLength > sizeCutoff& isLast==1),nums)
QuantileEnergyData_last4 = getQuantilMediansEnergy(subset(mergedData,GeneLength >= sizeCutoff & GeneLength < sizeCutoff2& isLast==1),nums)
QuantileEnergyData_last5 = getQuantilMediansEnergy(subset(mergedData,GeneLength >= sizeCutoff2 & GeneLength < sizeCutoff3 & isLast==1),nums)
QuantileEnergyData_last6 = getQuantilMediansEnergy(subset(mergedData,GeneLength >= sizeCutoff3 & isLast==1),nums)


#QuantileAcceptorData2 = getQuantilMediansAcceptor(subset(mergedData,GeneLength < sizeCutoff),100)
#QuantileAcceptorData3 = getQuantilMediansAcceptor(subset(mergedData,GeneLength > sizeCutoff),100)

QuantileEnergyData_HK = getQuantilMediansEnergy(subset(mergedData,isHK==1),20)
QuantileEnergyData2_HK = getQuantilMediansEnergy(subset(mergedData,GeneLength < sizeCutoff & isHK==1),20)
QuantileEnergyData3_HK = getQuantilMediansEnergy(subset(mergedData,GeneLength > sizeCutoff & isHK==1),20)


#QuantileEnergyData2_last = getQuantilMediansEnergy(subset(mergedData,GeneLength < sizeCutoff & isLast==1& isHK==1),100)
#QuantileEnergyData3_last = getQuantilMediansEnergy(subset(mergedData,GeneLength > sizeCutoff & isLast==1& isHK==1),100)

#png("Tilgner2_Median_by_size.png",height=600,width=800)
#par(mfrow=c(2,2))
#plot(log10(QuantileEnergyData[,1]), QuantileEnergyData[,2],xlim=log10(c(80,420000)),ylim=c(-1.24,1.24),main='Median Stability: Non-last exons', ylab='Median Exon Stability', xlab="LOG10 Distance to End of Gene") 
#plot(log10(QuantileEnergyData2[,1]), QuantileEnergyData2[,2],pch=19,xlim=log10(c(80,900000)),ylim=c(-1.24,1.24),main='Median Stability in non-last exons: Long vs short genes', ylab='Median Exon Stability', xlab="LOG10 Distance to End of Gene") 
##points(log10(QuantileEnergyData3[,1]), QuantileEnergyData3[,2])
#points(log10(QuantileEnergyData4[,1]), QuantileEnergyData4[,2],pch=19,col='blue')
#points(log10(QuantileEnergyData5[,1]), QuantileEnergyData5[,2],pch=19,col='red')
#points(log10(QuantileEnergyData6[,1]), QuantileEnergyData6[,2],pch=19,col='green')
#
#plot(log10(QuantileEnergyData_last[,1]), QuantileEnergyData_last[,2],xlim=log10(c(80,420000)),ylim=c(-1.24,1.24),main='Median Stability: Last Exons', ylab='Median Exon Stability', xlab="LOG10 Distance to End of Gene") 
#plot(log10(QuantileEnergyData_last2[,1]), QuantileEnergyData_last2[,2],pch=19,xlim=log10(c(80,900000)),ylim=c(-1.24,1.24),main='Median Stability in Last Exons: Long vs short genes', ylab='Median Exon Stability', xlab="LOG10 Distance to End of Gene") 
##points(log10(QuantileEnergyData3[,1]), QuantileEnergyData3[,2])
#points(log10(QuantileEnergyData_last4[,1]), QuantileEnergyData_last4[,2],pch=19,col='blue')
#points(log10(QuantileEnergyData_last5[,1]), QuantileEnergyData_last5[,2],pch=19,col='red')
#points(log10(QuantileEnergyData_last6[,1]), QuantileEnergyData_last6[,2],pch=19,col='green')
#legend('topright',pch=19,col=c(1,'blue','red','green'),leg=c(sprintf('Gene Length < %d',qr[1]), sprintf('%d <= Gene Length < %d',qr[1], ceiling(qr[2])), sprintf('%d <= Gene Length < %d',ceiling(qr[2]), ceiling(qr[3])), sprintf('%d <= GeneLength',ceiling(qr[3])) ) )
#dev.off()

png("Median_Energy_by_size-last.png",height=400,width=800)
require(RColorBrewer)
myCols = brewer.pal(4,'Dark2')
par(mfrow=c(1,2))
plot(log10(QuantileEnergyData[,1]), QuantileEnergyData[,2],xlim=log10(c(80,420000)),ylim=c(-1.24,1.24),main='Median Stability: Non-last exons', ylab='Median Exon Stability', xlab="LOG10 Distance to End of Gene") 
points(log10(QuantileEnergyData_last[,1]), QuantileEnergyData_last[,2],pch=20)
legend('topleft',pch=c(1,20),leg=c('non-last','last'))

plot(log10(QuantileEnergyData2[,1]), QuantileEnergyData2[,2],col=myCols[1],xlim=log10(c(80,900000)),ylim=c(-1.24,1.24),main='Median Stability in non-last exons: Long vs short genes', ylab='Median Exon Stability', xlab="LOG10 Distance to End of Gene") 
#points(log10(QuantileEnergyData3[,1]), QuantileEnergyData3[,2])
points(log10(QuantileEnergyData4[,1]), QuantileEnergyData4[,2],col=myCols[2])
points(log10(QuantileEnergyData5[,1]), QuantileEnergyData5[,2],col=myCols[3])
points(log10(QuantileEnergyData6[,1]), QuantileEnergyData6[,2],col=myCols[4])

points(log10(QuantileEnergyData_last2[,1]), QuantileEnergyData_last2[,2],pch=20,col=myCols[1])
#points(log10(QuantileEnergyData3[,1]), QuantileEnergyData3[,2])
points(log10(QuantileEnergyData_last4[,1]), QuantileEnergyData_last4[,2],pch=20,col=myCols[2])
points(log10(QuantileEnergyData_last5[,1]), QuantileEnergyData_last5[,2],pch=20,col=myCols[3])
points(log10(QuantileEnergyData_last6[,1]), QuantileEnergyData_last6[,2],pch=20,col=myCols[4])
legend('topright',pch=19,col=myCols,leg=c(sprintf('Gene Length < %d',qr[1]), sprintf('%d <= Gene Length < %d',qr[1], ceiling(qr[2])), sprintf('%d <= Gene Length < %d',ceiling(qr[2]), ceiling(qr[3])), sprintf('%d <= GeneLength',ceiling(qr[3])) ) )
dev.off()



#plot(log10(QuantileEnergyData_HK[,1]), QuantileEnergyData_HK[,2],xlim=log10(c(500,420000)),ylim=c(-1.24,1.24),main='Median Stability', ylab='Median Exon Stability', xlab="LOG10 Distance to End of Gene") 
#plot(log10(QuantileEnergyData2_HK[,1]), QuantileEnergyData2_HK[,2],pch=17,xlim=log10(c(500,420000)),ylim=c(-1.24,1.24),main='Median Stability: Long vs short genes', ylab='Median Exon Stability', xlab="LOG10 Distance to End of Gene") 
#points(log10(QuantileEnergyData3_HK[,1]), QuantileEnergyData3_HK[,2])#,xlim=log10(c(500,420000)),ylim=c(-1.24,1.24),main='Median Stability: LONG GENES', ylab='Median Exon Stability', xlab="LOG10 Distance to End of Gene") 
#
#par(mfrow=c(1,2))
#plot(log10(QuantileEnergyData_last[,1]), QuantileEnergyData_last[,2],xlim=log10(c(500,420000)),ylim=c(-1.24,1.24),main='Median Stability', ylab='Median Exon Stability', xlab="LOG10 Distance to End of Gene") 
#plot(log10(QuantileEnergyData2_last[,1]), QuantileEnergyData2_last[,2],pch=17,xlim=log10(c(500,420000)),ylim=c(-1.24,1.24),main='Median Stability: Long vs short genes', ylab='Median Exon Stability', xlab="LOG10 Distance to End of Gene") 
#points(log10(QuantileEnergyData3_last[,1]), QuantileEnergyData3_last[,2])#,xlim=log10(c(500,420000)),ylim=c(-1.24,1.24),main='Median Stability: LONG GENES', ylab='Median Exon Stability', xlab="LOG10 Distance to End of Gene") 



cols2=cols[c(1,3,7,9,10)]
#mergedData$intronCountGroups = quantileGroups(mergedData$MaxFeature,4)
mergedData$intronCountGroups = splitByVector(mergedData$MaxFeature,c(2.5,4.5,7.5,19.5))

for(HK in list(1,0)){
  myText = ifelse(length(HK)==2,'All',ifelse(HK==1,'HK','non'))
  
  lists=  list(subset(mergedData,GeneLength < qr[1] & isHK %in% HK), subset(mergedData, GeneLength > qr[3]  & isHK %in% HK))
  mains=c("SHORT","LONG")
  
  png(sprintf("Stability_size_numIntrons_%s.png",myText),width=800,height=400)
  par(mfrow=c(1,2))
  for(L in 1:2){
    #groups = list((1:3),4:6,7:9,10:12,13:15,16:18)
    #groups = lapply(tapply(mergedData$MaxFeature,mergedData$intronCountGroups,range),function(x)x[1]:x[2])
    lists[[1]]$intronCountGroups = quantileGroups(lists[[1]]$MaxFeature,4)
    groups = lapply(tapply(lists[[1]]$MaxFeature,lists[[1]]$intronCountGroups,range),function(x)x[1]:x[2])
    if(L==2)
      groups[[length(groups)+1]] = (tail(groups,1)[[1]][2]+1):max(lists[[2]]$MaxFeature)
    for(i in 1:length(groups)){
      QED = getQuantilMediansEnergy(subset(lists[[L]], MaxFeature %in% groups[[i]]),10)
       if(i==1)
        plot(log10(QED[,1]),QED[,2],pch=19,col=cols2[i],xlim=log10(c(500,420000)),ylim=c(-1.25,1.25),main=paste(mains[L], 'Genes: Median Stability by # Total Exons '), ylab='Median Exon Stability', xlab="LOG10 Distance to End of Gene")
       else
        points(log10(QED[,1]),QED[,2],pch=19,col=cols2[i])                
  
    } 
    legend('topright',pch=19,col=cols2, leg=sub(":","-",paste(groups,'introns')))
  }
  dev.off()
  
  png(sprintf("Acceptor_size_numIntrons_%s.png",myText),width=800,height=400)
  par(mfrow=c(1,2))
  for(L in 1:2){
    #groups = list((1:3),4:6,7:9,10:12,13:15,16:18)
    #groups = lapply(tapply(mergedData$MaxFeature,mergedData$intronCountGroups,range),function(x)x[1]:x[2])
    lists[[1]]$intronCountGroups = quantileGroups(lists[[1]]$MaxFeature,4)
    groups = lapply(tapply(lists[[1]]$MaxFeature,lists[[1]]$intronCountGroups,range),function(x)x[1]:x[2])
    if(L==2)
      groups[[length(groups)+1]] = (tail(groups,1)[[1]][2]+1):max(lists[[2]]$MaxFeature)
    for(i in 1:length(groups)){
      QED = getQuantilMediansAcceptor(subset(lists[[L]], MaxFeature %in% groups[[i]]),10)
       if(i==1)
        plot(log10(QED[,1]),QED[,2],pch=19,col=cols2[i],xlim=log10(c(500,420000)),ylim=c(-1.25,1.25),main=paste(mains[L], 'Genes: Median Acceptor Strength by # Total Exons '), ylab='Median Acceptor Strength', xlab="LOG10 Distance to End of Gene")
       else
        points(log10(QED[,1]),QED[,2],pch=19,col=cols2[i])                
  
    } 
    legend('topright',pch=19,col=cols2, leg=sub(":","-",paste(groups,'introns')))
  }
  dev.off()
}  


pdf("Median_energy_by_remaining_introns.pdf",height=6,width=9)
par(mfcol=c(2,3))
cols=brewer.pal(10,'Spectral')

### ALL GENES
plot(log10(QuantileEnergyData[,1]), QuantileEnergyData[,2],xlim=log10(c(500,420000)),ylim=c(-1.24,1.24),main='Median Stability', ylab='Median Exon Stability', xlab="LOG10 Distance to End of Gene") 
MedLengths = list()
MedEnergy = list()

S=sapply(1:10, function(i){
     if(i<10)
      QD = getQuantilMediansEnergy(subset(mergedData,ExonsRemaining==i),20)
     else
      QD = getQuantilMediansEnergy(subset(mergedData,ExonsRemaining >=i),20)      
     MedLengths[[i]] <<- QD[,1]
     MedEnergy[[i]] <<- QD[,2]      
     if(i==1)
      plot(log10(QD[,1]),QD[,2],col=cols[i],xlim=log10(c(500,420000)),ylim=c(-1.24,1.24),main='Median Stability by # exons remaining  \n(including current exon)', ylab='Median Exon Stability', xlab="LOG10 Distance to End of Gene")
     else
      points(log10(QD[,1]),QD[,2],col=cols[i])                 })
legend('bottomright',col=cols,pch=20,leg=c(1:9,'>=10'))

### SHORT GENES

plot(log10(QuantileEnergyData2[,1]), QuantileEnergyData2[,2],,xlim=log10(c(500,420000)),ylim=c(-1.24,1.24),main='Median Stability: SHORT GENES', ylab='Median Exon Stability', xlab="LOG10 Distance to End of Gene") 

MedLengths = list()
MedEnergy = list()

S=sapply(1:10, function(i){
     if(i<10)
      QD = getQuantilMediansEnergy(subset(mergedData,ExonsRemaining==i & GeneLength < 10^4.5),20)
     else
      QD = getQuantilMediansEnergy(subset(mergedData,ExonsRemaining >=i & GeneLength < 10^4.5),20)      
     MedLengths[[i]] <<- QD[,1]
     MedEnergy[[i]] <<- QD[,2]      
     if(i==1)
      plot(log10(QD[,1]),QD[,2],col=cols[i],xlim=log10(c(500,420000)),ylim=c(-1.24,1.24),main='Median Stability by # exons remaining  \n(including current exon)', ylab='Median Exon Stability', xlab="LOG10 Distance to End of Gene")
     else
      points(log10(QD[,1]),QD[,2],col=cols[i])                 })
legend('bottomright',col=cols,pch=20,leg=c(1:9,'>=10'))



### LONG GENES

plot(log10(QuantileEnergyData3[,1]), QuantileEnergyData3[,2],,xlim=log10(c(500,420000)),ylim=c(-1.24,1.24),main='Median Stability: LONG GENES', ylab='Median Exon Stability', xlab="LOG10 Distance to End of Gene") 

MedLengths = list()
MedEnergy = list()

S=sapply(1:10, function(i){
     if(i<10)
      QD = getQuantilMediansEnergy(subset(mergedData,ExonsRemaining==i & GeneLength > 10^4.5),20)
     else
      QD = getQuantilMediansEnergy(subset(mergedData,ExonsRemaining >=i & GeneLength > 10^4.5),20)      
     MedLengths[[i]] <<- QD[,1]
     MedEnergy[[i]] <<- QD[,2]      
     if(i==1)
      plot(log10(QD[,1]),QD[,2],col=cols[i],xlim=log10(c(500,420000)),ylim=c(-1.24,1.24),main='Median Stability by # exons remaining \n (including current exon)', ylab='Median Exon Stability', xlab="LOG10 Distance to End of Gene")
     else
      points(log10(QD[,1]),QD[,2],col=cols[i])                 })
legend('bottomright',col=cols,pch=20,leg=c(1:9,'>=10'))

dev.off()


##### HK Energy ############
pdf("Median_energy_by_remaining_introns_HK.pdf",height=6,width=9)
par(mfcol=c(2,3))
cols=brewer.pal(10,'Spectral')

### ALL GENES
plot(log10(QuantileEnergyData_HK [,1]), QuantileEnergyData_HK [,2],xlim=log10(c(500,420000)),ylim=c(-1.24,1.24),main='HK: Median Stability', ylab='Median Exon Stability', xlab="LOG10 Distance to End of Gene") 
MedLengths = list()
MedEnergy = list()

S=sapply(1:10, function(i){
     if(i<10)
      QD = getQuantilMediansEnergy(subset(mergedData,ExonsRemaining==i & isHK==1),20)
     else
      QD = getQuantilMediansEnergy(subset(mergedData,ExonsRemaining >=i & isHK==1),20)      
     MedLengths[[i]] <<- QD[,1]
     MedEnergy[[i]] <<- QD[,2]      
     if(i==1)
      plot(log10(QD[,1]),QD[,2],col=cols[i],xlim=log10(c(500,420000)),ylim=c(-1.24,1.24),main='Median Stability by # exons remaining \n(including current exon)', ylab='Median Exon Stability', xlab="LOG10 Distance to End of Gene")
     else
      points(log10(QD[,1]),QD[,2],col=cols[i])                 })
legend('bottomright',col=cols,pch=20,leg=c(1:9,'>=10'))

### SHORT GENES

plot(log10(QuantileEnergyData2[,1]), QuantileEnergyData2[,2],,xlim=log10(c(500,420000)),ylim=c(-1.24,1.24),main='Median Stability: SHORT GENES', ylab='Median Exon Stability', xlab="LOG10 Distance to End of Gene") 

MedLengths = list()
MedEnergy = list()

S=sapply(1:10, function(i){
     if(i<10)
      QD = getQuantilMediansEnergy(subset(mergedData,ExonsRemaining==i & GeneLength < 10^4.5  & isHK==1),20)
     else
      QD = getQuantilMediansEnergy(subset(mergedData,ExonsRemaining >=i & GeneLength < 10^4.5  & isHK==1),20)      
     MedLengths[[i]] <<- QD[,1]
     MedEnergy[[i]] <<- QD[,2]      
     if(i==1)
      plot(log10(QD[,1]),QD[,2],col=cols[i],xlim=log10(c(500,420000)),ylim=c(-1.24,1.24),main='Median Stability by # exons remaining  \n(including current exon)', ylab='Median Exon Stability', xlab="LOG10 Distance to End of Gene")
     else
      points(log10(QD[,1]),QD[,2],col=cols[i])                 })
legend('bottomright',col=cols,pch=20,leg=c(1:9,'>=10'))



### LONG GENES

plot(log10(QuantileEnergyData3[,1]), QuantileEnergyData3[,2],,xlim=log10(c(500,420000)),ylim=c(-1.24,1.24),main='Median Stability: LONG GENES', ylab='Median Exon Stability', xlab="LOG10 Distance to End of Gene") 

MedLengths = list()
MedEnergy = list()

S=sapply(1:10, function(i){
     if(i<10)
      QD = getQuantilMediansEnergy(subset(mergedData,ExonsRemaining==i & GeneLength > 10^4.5  & isHK==1),20)
     else
      QD = getQuantilMediansEnergy(subset(mergedData,ExonsRemaining >=i & GeneLength > 10^4.5  & isHK==1),20)      
     MedLengths[[i]] <<- QD[,1]
     MedEnergy[[i]] <<- QD[,2]      
     if(i==1)
      plot(log10(QD[,1]),QD[,2],col=cols[i],xlim=log10(c(500,420000)),ylim=c(-1.24,1.24),main='Median Stability by # exons remaining  \n(including current exon)', ylab='Median Exon Stability', xlab="LOG10 Distance to End of Gene")
     else
      points(log10(QD[,1]),QD[,2],col=cols[i])                 })
legend('bottomright',col=cols,pch=20,leg=c(1:9,'>=10'))

dev.off()




### ACCEPTORS ###########
pdf("Median_acceptor_by_remaining_introns.pdf",height=6,width=9)
par(mfcol=c(2,3))
cols=brewer.pal(10,'Spectral')

### ALL GENES
plot(log10(QuantileAcceptorData[,1]), QuantileAcceptorData[,2],xlim=log10(c(500,420000)),ylim=c(-1,1),main='Median Acceptor Strength', ylab='Median Exon Acceptor Strength', xlab="LOG10 Distance to End of Gene") 
MedLengths = list()
MedAcceptor = list()

S=sapply(1:10, function(i){
     if(i<10)
      QD = getQuantilMediansAcceptor(subset(mergedData,ExonsRemaining==i),20)
     else
      QD = getQuantilMediansAcceptor(subset(mergedData,ExonsRemaining >=i),20)      
     MedLengths[[i]] <<- QD[,1]
     MedAcceptor[[i]] <<- QD[,2]      
     if(i==1)
      plot(log10(QD[,1]),QD[,2],col=cols[i],xlim=log10(c(500,420000)),ylim=c(-1,1),main='Median Acceptor Strength by # exons remaining \n(including current exon)', ylab='Median Exon Acceptor Strength', xlab="LOG10 Distance to End of Gene")
     else
      points(log10(QD[,1]),QD[,2],col=cols[i])                 })
legend('bottomright',col=cols,pch=20,leg=c(1:9,'>=10'))

### SHORT GENES

plot(log10(QuantileAcceptorData2[,1]), QuantileAcceptorData2[,2],,xlim=log10(c(500,420000)),ylim=c(-1,1),main='Median Acceptor Strength: SHORT GENES', ylab='Median Exon Acceptor Strength', xlab="LOG10 Distance to End of Gene") 

MedLengths = list()
MedAcceptor = list()

S=sapply(1:10, function(i){
     if(i<10)
      QD = getQuantilMediansAcceptor(subset(mergedData,ExonsRemaining==i & GeneLength < 10^4.5),20)
     else
      QD = getQuantilMediansAcceptor(subset(mergedData,ExonsRemaining >=i & GeneLength < 10^4.5),20)      
     MedLengths[[i]] <<- QD[,1]
     MedAcceptor[[i]] <<- QD[,2]      
     if(i==1)
      plot(log10(QD[,1]),QD[,2],col=cols[i],xlim=log10(c(500,420000)),ylim=c(-1,1),main='Median Acceptor Strength by # exons remaining \n(including current exon)', ylab='Median Exon Acceptor Strength', xlab="LOG10 Distance to End of Gene")
     else
      points(log10(QD[,1]),QD[,2],col=cols[i])                 })
legend('bottomright',col=cols,pch=20,leg=c(1:9,'>=10'))



### LONG GENES

plot(log10(QuantileAcceptorData3[,1]), QuantileAcceptorData3[,2],,xlim=log10(c(500,420000)),ylim=c(-1,1),main='Median Acceptor Strength: LONG GENES', ylab='Median Exon Acceptor Strength', xlab="LOG10 Distance to End of Gene") 

MedLengths = list()
MedAcceptor = list()

S=sapply(1:10, function(i){
     if(i<10)
      QD = getQuantilMediansAcceptor(subset(mergedData,ExonsRemaining==i & GeneLength > 10^4.5),20)
     else
      QD = getQuantilMediansAcceptor(subset(mergedData,ExonsRemaining >=i & GeneLength > 10^4.5),20)      
     MedLengths[[i]] <<- QD[,1]
     MedAcceptor[[i]] <<- QD[,2]      
     if(i==1)
      plot(log10(QD[,1]),QD[,2],col=cols[i],xlim=log10(c(500,420000)),ylim=c(-1,1),main='Median Acceptor Strength by # exons remaining \n(including current exon)', ylab='Median Exon Acceptor Strength', xlab="LOG10 Distance to End of Gene")
     else
      points(log10(QD[,1]),QD[,2],col=cols[i])                 })
legend('bottomright',col=cols,pch=20,leg=c(1:9,'>=10'))

dev.off()


### ACCEPTORS ###########
#pdf("Median_acceptor_by_remaining_introns.pdf",height=8,width=8)
#par(mfrow=c(2,2))
#plot(log10(QuantileAcceptorData[,1]), QuantileAcceptorData[,2],,xlim=log10(c(500,420000)),ylim=c(-1,1),main='Median Acceptor', ylab='Median Exon Acceptor Strength', xlab="LOG10 Distance to End of Gene") 
#
#cols=brewer.pal(10,'Spectral')
#
#MedLengths = list()
#MedAcceptor = list()
#S=sapply(1:10, function(i){
#     if(i<10)
#      QD = getQuantilMediansAcceptor(subset(mergedData,ExonsRemaining==i),20)
#     else
#      QD = getQuantilMediansAcceptor(subset(mergedData,ExonsRemaining >=i),20)      
#     MedLengths[[i]] <<- QD[,1]
#     MedAcceptor[[i]] <<- QD[,2]      
#     if(i==1)
#      plot(log10(QD[,1]),QD[,2],col=cols[i],xlim=log10(c(500,420000)),ylim=c(-1,1),main='Median Acceptor by # exons remaining (including current exon)', ylab='Median Acceptor Score', xlab="LOG10 Distance to End of Gene")
#     else
#      points(log10(QD[,1]),QD[,2],col=cols[i])                 })
#legend('bottomright',col=cols,pch=20,leg=c(1:9,'>=10'))
#dev.off()

################### 3: Last exons: compared to Splicing but separated for different groups    ##################### 
 
introns3.lasts = subset(introns3, intronsRemaining==0) 
introns3.lasts$lengthGroups10 = quantileGroups(introns3.lasts$Dist2End) 
introns3.lasts$lengthGroups5 = quantileGroups(introns3.lasts$Dist2End,5) 
introns3.lasts$splicingGroups5 = quantileGroups(introns3.lasts$splicingAMinus,5) 
introns3.lasts$splicingGroups10 = quantileGroups(introns3.lasts$splicingAMinus,10) 
QD = getQuantilMedians2(introns3.lasts,10,'Acceptor','splicingAMinus')
plot(QD$GroupMedian, QD$DataMedian, ylab='splicing index', xlab='Acceptor Score', ylim=c(0.3,0.58))
cols2 = cols[seq(1,10,2)]
S=sapply(1:5, function(i){
  QD = getQuantilMedians2(subset(introns3.lasts,lengthGroups5==i),10,'Acceptor','splicingAMinus')
 print(range(QD$DataMedian))
  FIT <- lm(DataMedian ~ GroupMedian,QD)
  abline(a=coef(FIT)[1], b=coef(FIT)[2],col=cols2[i],lty=2)
  points(QD[,1],QD[,2],col=cols2[i],pch=20)    
})
legend('bottomright',col=cols2,pch=20,leg=paste("Group",1:5))



 
 
#pdf("Tilgner2_Median_spliceMinus_by_energy-acceptor.pdf")
par(mfrow=c(2,2))

### First plot the coSI vs Energy ###
QuantileEnergyVsCoSI =   getQuantilMedians2(introns3,100,'Stability','splicingAMinus')
EnergyFit = lm(DataMedian ~ GroupMedian,QuantileEnergyVsCoSI)
plot(QuantileEnergyVsCoSI, xlab="Stability Score", ylab="Median splicingAMinus",ylim=c(.4, 0.96))
lines(EE <- seq(-3,3,.05), predict(EnergyFit,data.frame(GroupMedian=EE))) 

## Next plot for different introns remaining ##
plot(EE, predict(EnergyFit,data.frame(GroupMedian=EE)), type='l',ylim=c(.4, 0.96), xlab="Stability Score", ylab="Median coSI" , main='Median coSI index for # of remaining introns')
cols=brewer.pal(10,'Spectral')
allFits = list()
Residuals = list()
MedLengths = list()
MedCoSI = list()
S=sapply(1:10, function(i){
     if(i<10)
      QD = getQuantilMedians2(subset(introns3,intronsRemaining==i),20,'Stability','splicingAMinus')
     else
      QD = getQuantilMedians2(subset(introns3,intronsRemaining>=i),20,'Stability','splicingAMinus')      
     MedLengths[[i]] <<- QD[,1]
     MedCoSI[[i]] <<- QD[,2]      
     FIT <- lm(DataMedian ~ GroupMedian,QD)
     abline(a=coef(FIT)[1], b=coef(FIT)[2],col=cols[i],lty=2)
     points(QD[,1],QD[,2],col=cols[i])    
     
     Residuals[[i]] <<- QD$DataMedian - predict(EnergyFit,QD)      
})
legend('bottomright',col=cols,pch=20,leg=c(1:9,'>=10'))


QuantileAcceptorVsCoSI =   getQuantilMedians2(introns3,100,'Acceptor','splicingAMinus')
AcceptorFit = lm(DataMedian ~ GroupMedian,QuantileAcceptorVsCoSI)
plot(QuantileAcceptorVsCoSI, xlab="Stability Score", ylab="Median splicingAMinus",ylim=c(.4, 0.96))
lines(EE <- seq(-3,3,.05), predict(AcceptorFit,data.frame(GroupMedian=EE))) 

## Next plot for different introns remaining ##
plot(EE, predict(AcceptorFit,data.frame(GroupMedian=EE)), type='l',ylim=c(.4, 0.96), xlab="Acceptor Score", ylab="Median coSI" , main='Median coSI index for # of remaining introns')
cols=brewer.pal(10,'Spectral')
allFits = list()
Residuals = list()
MedLengths = list()
MedCoSI = list()
S=sapply(1:10, function(i){
     if(i<10)
      QD = getQuantilMedians2(subset(introns3,intronsRemaining==i),20,'Acceptor','splicingAMinus')
     else
      QD = getQuantilMedians2(subset(introns3,intronsRemaining>=i),20,'Acceptor','splicingAMinus')      
     MedLengths[[i]] <<- QD[,1]
     MedCoSI[[i]] <<- QD[,2]      
     FIT <- lm(DataMedian ~ GroupMedian,QD)
     abline(a=coef(FIT)[1], b=coef(FIT)[2],col=cols[i],lty=2)
     points(QD[,1],QD[,2],col=cols[i])    
     
     Residuals[[i]] <<- QD$DataMedian - predict(AcceptorFit,QD)      
})
legend('bottomright',col=cols,pch=20,leg=c(1:9,'>=10'))
dev.off()





############################################################################################################
###### See what happens if I try to model whole data, not just medians:

# Basic fit
global_fit0 = lm(splicingAMinus ~ log10(Dist2End), data=introns3)
global_fit1 = nls(splicingAMinus ~ MeanGamma2(Dist2End,g,k,A,B), data=introns3,start=c(g=1,k=5e-5,A=1,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,Inf,Inf),algo='port');
global_fit1.HK = nls(splicingAMinus ~ MeanGamma2(Dist2End,g,k,A,B), data=subset(introns3,isHK==1),start=c(g=1,k=5e-5,A=1,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,Inf,Inf),algo='port');
global_fit1.non = nls(splicingAMinus ~ MeanGamma2(Dist2End,g,k,A,B), data=subset(introns3,isHK==0),start=c(g=1,k=5e-5,A=1,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,Inf,Inf),algo='port');

# Add acceptor
global_fit0.2 = lm(splicingAMinus ~ log10(Dist2End) + Acceptor, data=introns3)
global_fit0.2.HK = lm(splicingAMinus ~ log10(Dist2End) + Acceptor, data=subset(introns3,isHK==1))
global_fit0.2.non = lm(splicingAMinus ~ log10(Dist2End) + Acceptor, data=subset(introns3,isHK==0))
global_fit2 = nls2(splicingAMinus ~ MeanGamma3(cbind(Dist2End,Acceptor),g,k,A,B,sigma), data=introns3,start=c(g=1,k=.8,A=.91,B=1000,sigma=.0001),lower=c(1,0,0,1,0),upper=c(1,Inf,Inf,Inf,5),algo='port');
global_fit2.HK = nls2(splicingAMinus ~ MeanGamma3(cbind(Dist2End,Acceptor),g,k,A,B,sigma), data=subset(introns3,isHK==1),start=c(g=1,k=.8,A=.91,B=1000,sigma=.0001),lower=c(1,0,0,1,0),upper=c(1,Inf,Inf,Inf,5),algo='port');
global_fit2.non = nls2(splicingAMinus ~ MeanGamma3(cbind(Dist2End,Acceptor),g,k,A,B,sigma), data=subset(introns3,isHK==0),start=c(g=1,k=.8,A=.91,B=1000,sigma=.0001),lower=c(1,0,0,1,0),upper=c(1,Inf,Inf,Inf,5),algo='port');

# Add donor
global_fit0.3.non = lm(splicingAMinus ~ log10(Dist2End) + Acceptor + Donor, data=subset(introns3,isHK==0))
global_fit0.3.HK = lm(splicingAMinus ~ log10(Dist2End) + Acceptor + Donor, data=subset(introns3,isHK==1))
global_fit3.HK = nls2(splicingAMinus ~ MeanGamma4(cbind(Dist2End,Acceptor,Donor),g,k,A,B,C1,C2), data=subset(introns3,isHK==1),start=c(g=1,k=.8,A=.91,B=1000,C1=.06,C2=0.01),lower=c(1,0,0,1,0,0),upper=c(1,Inf,Inf,Inf,5,5),algo='port');
global_fit3.non= nls2(splicingAMinus ~ MeanGamma4(cbind(Dist2End,Acceptor,Donor),g,k,A,B,C1,C2), data=subset(introns3,isHK==0),start=c(g=1,k=.8,A=.91,B=1000,C1=.06,C2=0.01),lower=c(1,0,0,1,0,0),upper=c(1,Inf,Inf,Inf,5,5),algo='port');

###################### 
cor(predict(global_fit0), global_fit0$model$splicingAMinus)
# 0.4173612
cor(predict(global_fit0.2), global_fit0.2$model$splicingAMinus)
# 0.4199824  Linear with Acceptor

cor(introns3$splicingAMinus, MeanGamma2(introns3$Dist2End,coef(global_fit1)[1], coef(global_fit1)[2],coef(global_fit1)[3], coef(global_fit1)[4]), use='c')
# 0.4209286
Model 
cor(subset(introns3,!is.na(Acceptor),c('splicingAMinus')), m3<-MeanGamma3(subset(introns3,!is.na(Acceptor),c('Dist2End','Acceptor')),coef(global_fit2)[1], coef(global_fit2)[2],coef(global_fit2)[3], coef(global_fit2)[4], coef(global_fit2)[5]), use='c')
# 0.4230337 Model with Acceptor

###################### 
### Test all models on HK data only:
c0 = cor(predict(global_fit0,newdata=subset(introns3,isHK==1)), subset(introns3,isHK==1,'splicingAMinus'),use='c')
# 0.4242964 non, .4066297 HK
c0.2 = sapply(c(0,1),function(x)cor(predict(global_fit0.2,newdata=subset(introns3,isHK==x)), subset(introns3,isHK==x,'splicingAMinus'),use='c' ))
#  0.4267795 non, 0.4090903 HK  Linear with Acceptor

cor(predict(global_fit1,newdata=subset(introns3 ,isHK==1 & !is.na(Acceptor))), subset(introns3,isHK==1 & !is.na(Acceptor),'splicingAMinus'),use='c' )
#  0.4267975 Non, 0.410721 HK
cor(predict(global_fit2,newdata=subset(introns3,isHK==1 & !is.na(Acceptor))), subset(introns3,isHK==1 & !is.na(Acceptor),'splicingAMinus'),use='c' )
# Model with Acceptor: 0.428770 non,  4118466 HK

###################### 
## Test Hk-specific models
# Basic linear fit
c0.2cor(predict(global_fit0.2.HK,newdata=subset(introns3,isHK==1)), subset(introns3,isHK==1,'splicingAMinus'),use='c')
#  0.4091151
cor(predict(global_fit0.2.non,newdata=subset(introns3,isHK==0)), subset(introns3,isHK==0,'splicingAMinus'),use='c' )
#     0.4267843

# Fit with acceptor
cor(predict(global_fit2.HK,newdata=subset(introns3 ,isHK==1 & !is.na(Acceptor))), subset(introns3,isHK==1 & !is.na(Acceptor),'splicingAMinus'),use='c' )
#   0.4141876
cor(predict(global_fit2.non,newdata=subset(introns3,isHK==0 & !is.na(Acceptor))), subset(introns3,isHK==0 & !is.na(Acceptor),'splicingAMinus'),use='c' )
#0.4288815


cor(predict(global_fit0.3.HK,newdata=subset(introns3,isHK==1)), subset(introns3,isHK==1,'splicingAMinus'),use='c')
#[1,]      0.4051831
cor(predict(global_fit0.3.non,newdata=subset(introns3,isHK==0)), subset(introns3,isHK==0,'splicingAMinus'),use='c' )
#[1,]      0.4271821


## Include Donor Score
cor(predict(global_fit3.HK,newdata=subset(introns3 ,isHK==1 & !is.na(Acceptor) & !is.na(Donor))), subset(introns3,isHK==1 & !is.na(Acceptor) & !is.na(Donor),'splicingAMinus'),use='c' )
#  0.4114617
cor(predict(global_fit3.non,newdata=subset(introns3,isHK==0 & !is.na(Acceptor) & !is.na(Donor))), subset(introns3,isHK==0 & !is.na(Acceptor) & !is.na(Donor),'splicingAMinus'),use='c' )
#  0.4290831


############
# 
range(cbind(1,subset(introns3,isHK==1)$Acceptor, subset(introns3,isHK==1)$Donor)%*%coef(global_fit3.HK)[c(2,5,6)] ,na.rm=T)
#sd((cbind(1,subset(introns3,isHK==1)$Acceptor, subset(introns3,isHK==1)$Donor)%*%coef(global_fit3.HK)[c(2,5,6)])[,1] ,na.rm=T)
#0.5337450 0.773925
range(cbind(1,subset(introns3,isHK==0)$Acceptor, subset(introns3,isHK==0)$Donor)%*%coef(global_fit3.non)[c(2,5,6)] ,na.rm=T)
# 0.4022572 1.1298271
#sd((cbind(1,subset(introns3,isHK==0)$Acceptor, subset(introns3,isHK==0)$Donor)%*%coef(global_fit3.non)[c(2,5,6)])[,1] ,na.rm=T)




########################################################################################################################
########################################################################################################################
## Can I add in some info from GRO-seq data published by Danko 2013???

MCF7_rates = read.delim("MCF7.10-40m.regressionRate.txt")
AC16_rates = read.delim("AC16.regressionRate.txt")
 
library("GenomicRanges")
library("rtracklayer")
MCF7 <- GRanges(seqnames = MCF7_rates$chrom, ranges = IRanges(MCF7_rates$chromStart, end=MCF7_rates$chromEnd), strand=MCF7_rates$strand, uid=1:nrow(MCF7_rates))
refseq.h19.GR = GRanges(seqnames=refseq$chrom, ranges=IRanges(refseq$txStart, end=refseq$txEnd), strand=refseq$strand, id=refseq$name)

# liftover to hg19:
hg18To19 = import.chain("hg18ToHg19.over.chain")
MCF7_hg19 = liftOver(MCF7, hg18To19)

#A=match(MCF7_hg19,refseq.h19.GR)
#MCF7_rates$UniqueID = refseq$UniqueID[A]

# Match IDs to introns3
#B=match(MCF7_hg19, GRanges(seqnames=introns3$chr, ranges=IRanges(introns3$txStart,introns3$txEnd), strand=introns3$strand))
#MCF7_rates$UniqueID = introns3$UniqueID[B]
B=match(GRanges(seqnames=introns3$chr, ranges=IRanges(introns3$txStart,introns3$txEnd), strand=introns3$strand),MCF7_hg19)

# Get the data back
introns3$MCF7_elong = MCF7_rates$rate[B]



################################################################################################################
################################################################################################################
#### Simulate model!! 
################################################################################################################
################################################################################################################

load('Multiple_CellTypes_splicing1000.RData')  
chrom_coef1000 <- t(sapply(FitData,coef))[2,]
chrom_noRT_coef1000 = t(sapply(FitData2,coef))[2,]#0.001577696

#fitData_to_model = c1000L
fitData_to_model = chrom_coef1000
fitData_to_model2 = chrom_noRT_coef1000

mergedData = mergedData[order(mergedData$UniqueID,mergedData$FeatureCount),]

# Using the fitted parameters

# Assign each gene the median stability from its group:
S=with(subset(mergedData,isLast==0), tapply(Stability,list(GeneSize,intronCountGroups),median, na.rm=T))
S2 = data.frame(S=as.vector(S),GeneSize=rep(1:4,5),intronCountGroups=  rep(1:5, each=4))
SS=with(subset(mergedData,isLast==0), tapply(Stability,list(GeneSize,intronCountGroups),mean, na.rm=T))
#S2 = data.frame(S=as.vector(SS),GeneSize=rep(1:4,5),intronCountGroups=  rep(1:5, each=4))
mergedData$GroupStability = S2$S[match(with(mergedData,paste(GeneSize,intronCountGroups )), with(S2,paste(GeneSize,intronCountGroups )) )]


Elong = 1/fitData_to_model[2]
readThrough = fitData_to_model[4]
Elong2 = 1/fitData_to_model2[1]
#readThrough = 3840

readThrough_for_pauseSplice = readThrough#3840

Stability_exp=1
Stability_div = .85
K_splice = 1
Acceptor_div = 4
extraDist = rep(80,nrow(mergedData))
extraDist =  mergedData$exonLength
extraDist = rep(147,nrow(mergedData))

# Try out new way of computing stability parameters:
Stability_sdlog = 1.2
Stability_MinOffset = 4.5
#log(dlnorm(mergedData$Stability+5,sdlog=sdlog))+MinOffset

# For simplified version: multiplication factor = 1 + 1/Acceptor_div
mult_factor = 2
Acceptor_div = 1/(mult_factor-1) 

# Use K_adjusted to balance out the Acceptor_div 
K_adjusted = K_splice * nrow(mergedData)/(nrow(mergedData) + (mult_factor-1)*length(unique(mergedData$UniqueID)))
K_adjusted=K_splice

# Check: this should equal K_splice
K_splice==K_adjusted*( (nrow(mergedData)-length(unique(mergedData$UniqueID))) + mult_factor*length(unique(mergedData$UniqueID)))/nrow(mergedData)



#myEfficiencies = tapply(1:nrow(mergedData), mergedData$UniqueID, function(x)simulateGene2(mergedData[x,], extraDist=extraDist,Elong=Elong,readThrough=readThrough, K_splice=K_adjusted,Stability_MinOffset=Stability_MinOffset,Stability_sdlog=Stability_sdlog,Acceptor_div=Acceptor_div))
#myEfficiencies_noSplice = tapply(1:nrow(mergedData), mergedData$UniqueID, function(x)simulateGene2(mergedData[x,], extraDist=extraDist,Elong=Elong,readThrough=readThrough, K_splice=K_splice,Acceptor_div=Inf,Stability_MinOffset=Stability_MinOffset,Stability_sdlog=Stability_sdlog))
mergedData_forFit = data.frame(Stability=mergedData$GroupStability, Dist2End=mergedData$Dist2End)
myEfficiencies = tapply(1:nrow(mergedData), mergedData$UniqueID, function(x)simulateGene(mergedData_forFit[x    ,], extraDist=extraDist,Elong=Elong,readThrough=readThrough, K_splice=K_adjusted,Stability_exp=Stability_exp,Stability_div=Stability_div,Acceptor_div=Acceptor_div))
myEfficiencies_noSplice = tapply(1:nrow(mergedData), mergedData$UniqueID, function(x)simulateGene(mergedData_forFit[x,], extraDist=extraDist,Elong=Elong,readThrough=readThrough, K_splice=K_splice,Acceptor_div=Inf,Stability_exp=Stability_exp,Stability_div=Stability_div))

myEfficiencies_noPause = tapply(1:nrow(mergedData), mergedData$UniqueID, function(x)simulateGene(mergedData[x,], extraDist=extraDist,Elong=Elong,readThrough=readThrough, K_splice=K_adjusted,Stability_div=Inf,Stability_exp=Stability_exp,Acceptor_div=Acceptor_div))

myEfficiencies_plain = tapply(1:nrow(mergedData), mergedData$UniqueID, function(x)simulateGene(mergedData[x,], extraDist=extraDist,Elong=Elong,readThrough=readThrough, K_splice=K_splice,Stability_div=Inf, Acceptor_div=Inf,Stability_exp=Stability_exp))
myEfficiencies_noRT = tapply(1:nrow(mergedData), mergedData$UniqueID, function(x)simulateGene(mergedData[x,], extraDist=extraDist,Elong=Elong,readThrough=0, K_splice=K_splice,Stability_div=Inf, Acceptor_div=Inf,Stability_exp=Stability_exp))
# Use median Read-Through of 976 bp (call it 1000) for this simulation,
# doing GRO-median of 3840 instead!!
myEfficiencies_plain_RT1000 = tapply(1:nrow(mergedData), mergedData$UniqueID, function(x)simulateGene(mergedData[x,], extraDist=extraDist,Elong=Elong,readThrough=3840, K_splice=K_splice,Stability_div=Inf, Acceptor_div=Inf,Stability_exp=Stability_exp))

# Get other info
genes.isHK = names(myEfficiencies) %in% subset(mergedData,isHK==1,UniqueID)[,1]
genes.Size = mergedData$GeneSize[match(names(myEfficiencies), mergedData$UniqueID)]
genes.IntronCount = mergedData$MaxFeature[match(names(myEfficiencies), mergedData$UniqueID)]
efficiencyData = data.frame(UniqueID=names( myEfficiencies), isHK=genes.isHK,GeneSize=genes.Size,MaxFeature=genes.IntronCount,intronCountGroups =  splitByVector(genes.IntronCount,c(2.5,4.5,7.5,19.5)),
  myEfficiencies_noRT,myEfficiencies_plain_RT1000,myEfficiencies_plain,myEfficiencies_noSplice,myEfficiencies_noPause,myEfficiencies)

# Chrom data:
#save(efficiencyData, Elong,readThrough,Snamestability_exp,Stability_div,K_splice,  Acceptor_div,mult_factor ,K_adjusted,extraDist, file='Tilgner2_parameterizedModel_chrom.RData')

plain = sapply(1:4,function(x)summary(myEfficiencies_plain ~ intronCountGroups,data=subset(efficiencyData,GeneSize==x)))
pause = sapply(1:4,function(x)summary(myEfficiencies_noSplice ~ intronCountGroups,data=subset(efficiencyData,GeneSize==x)))
splice= sapply(1:4,function(x)summary(myEfficiencies_noPause ~ intronCountGroups,data=subset(efficiencyData,GeneSize==x)))
both= sapply(1:4,function(x)summary(myEfficiencies ~ intronCountGroups,data=subset(efficiencyData,GeneSize==x)))

i=1
cbind(plain[[i]], pause[[i]],splice[[i]], both[[i]])

#save(efficiencyData, Elong,readThrough,Stability_exp,Stability_div,K_splice,  Acceptor_div,mult_factor ,K_adjusted,extraDist, file='Tilgner2_parameterizedModel_chrom1000-stabMeans_div.85_147.RData')

# Nuc data:
#save(efficiencyData, Elong,readThrough,Stability_exp,Stability_div,K_splice,  Acceptor_div,mult_factor ,K_adjusted,extraDist, file='Tilgner2_parameterizedModel2.RData')
#load(file='Tilgner2_parameterizedModel2.RData')

## Do by both, for *all* genes:
 

plot.dev("GenomeSimulations_numIntrons_by_geneSize_chromElong_147_noAdjust_mean.85_Extreme_only.pdf",'pdf',height=3,width=8)  

notch=F; wid=.16;bw=0.11; offset1 = 1.5; offset2 = .8; xlim=c(-.8,11); cex=0.6; pch=19
COL = c(brewer.pal(9,"YlGnBu")[5:9],brewer.pal(9,"OrRd")[5:9])
par(cex=cex,mai=c(.5,.4,.1,.1),pch=pch)
      
      
Xshift=0; add=FALSE; xpos=NULL
for (HK in list(0:1)){  
  for (exper in 2:4){    
    for(sizes in c(1,4)){
      txt = paste(ifelse(length(HK)==1,ifelse(HK,'HK:','Regulated'),'All:'))#, c('Short','Med short','Med long','Long')[sizes])
      x<-subset(efficiencyData,GeneSize==sizes &isHK%in%HK)
      myCats = sort(unique(x$intronCountGroups))
      
      main=''#ifelse(sizes==1,paste(txt,switch(exper,'Basal parameters','Add Pausing','Add Splicing variability','Add both')),'')
      # Decide which variable to test
      Form0 = myEfficiencies_plain ~ intronCountGroups
      Form = switch(exper,Form0,myEfficiencies_noSplice ~ intronCountGroups, myEfficiencies_noPause ~ intronCountGroups, myEfficiencies ~ intronCountGroups) 
    
      at = 2*(myCats)*wid + Xshift
      xpos = c(xpos,mean(at))
      
      b=boxplot(Form0, data=x,main=main,border=COL[myCats],staplewex=0,range=0.00001,outline=F,ylim=c(0,1),lwd=2,boxwex=bw*.5,at=at-wid*.7,xlim=xlim,names=NA,xaxt='n',notch=notch,add=add,cex=cex,pars=list(pch=pch,cex=cex))
      b=boxplot(Form, data=x,main=main,col=COL[myCats],boxwex=bw,at=at,names=NA,notch=notch,add=T,xaxt='n')
      add=TRUE 
      Xshift = offset1 + Xshift      
    }
      Xshift = offset2 + Xshift          
  }
 abline(h=0,col='gray')
  legend('left',fill=COL,leg=c('1-2 introns','3-4 introns','5-7 introns','8-19 introns','20+ introns'))#sprintf(c('1-2 introns (n=%d)','3-4 introns (n=%d)','5-7 introns (n=%d)','8-19 introns (n=%d)','20+ introns (n=%d)')[myCats], table(subset(efficiencyData,GeneSize==sizes&isHK%in%HK)$intronCountGroups)))
  axis(1,at=xpos, lab = rep(c('Short','Long'),3)) 

}
plot.off()
 
 
 


## Make plot for 'Fig 3', looking at effect of adding read-through

plot.dev("GenomeSimulations_readThrough_comparisons_chromElong_147_noAdjust_mean.85_2only.pdf",'pdf',width=7,height=1.5)
par(mfcol=c(1,1),cex=0.5,mai=c(0.3,0.2,0.01,0.01))
notch=F; wid=.21;bw=0.07; offset = 1.4; xlim=c(-1,9.2)
cols = cols=c(brewer.pal(9,"YlGnBu")[5:9],brewer.pal(9,"OrRd")[5:9])
for (HK in list(0:1)){

  for(sizes in 1:4){
    COL = cols[c(1:5) + HK[1]*5]
    txt = paste(ifelse(length(HK)==1,ifelse(HK,'HK:','Regulated'),'All:'))#, c('Short','Med short','Med long','Long')[sizes])
    x<-subset(efficiencyData,GeneSize==sizes &isHK%in%HK)
    myCats = sort(unique(x$intronCountGroups))

    main=''
    # Decide which variable to test
    Form1 = myEfficiencies_noRT ~ intronCountGroups
    Form3 = myEfficiencies_plain ~ intronCountGroups

    Xshift = offset*(sizes-1 )
    at = 2*(myCats*wid+Xshift) - (sizes-1)*wid*2
    lwd=1
    b=boxplot(Form1, data=x,main=main,border=COL[myCats],ylim=c(0,1),lwd=1.5,staplewex=0,range=0.00001,outline=F,boxwex=bw,at=at-wid*.5,xlim=xlim,names=NA,yaxt=ifelse(sizes==1,'s','n'),notch=notch,add=sizes-1,xaxt='n')
    if (sizes==1) abline(h=0,col='gray')
    b=boxplot(Form3, data=x,main=main,col=COL[myCats],border='black',boxwex=bw,at=at+wid*.5,names=NA,notch=notch,add=T,xaxt='n',lwd=lwd,yaxt='n')
  }
  
  legend('topleft',fill=COL,leg=c('1-2 introns','3-4 introns','5-7 introns','8-19 introns','20+ introns'))
  legend('bottomleft',fill=c(NA,COL[1]), border=c(COL[1],gray.colors(2)[1],'black'),leg=c('No Read-Thru',sprintf('%.1f kb Read-Thru',readThrough_for_pauseSplice/1000)))
 axis(1,at=2*(3*wid+(0:3)*offset)-(0:3)*wid*2, lab = c('Short','Med short','Med long','Long')) 
}
plot.off()


 
 
 
 
plot.dev("GenomeSimulations_numIntrons_by_geneSize_chromElong_147_noAdjust_mean.85.pdf",'pdf',height=10,width=10)  
par(mfcol=c(4,1))

notch=F; wid=.21;bw=0.16; offset = 1.25; xlim=c(-.8,8.2)

for (HK in list(0:1)){
  
  for (exper in 2:4){
    for(sizes in c(1:4)){

      COL = brewer.pal(11,'Spectral')[c(1:5) + HK[1]*6]
      COL = cols=c(brewer.pal(9,"YlGnBu")[5:9],brewer.pal(9,"OrRd")[5:9])
      #cols = cols=c(brewer.pal(9,"YlGnBu")[5:9],brewer.pal(9,"OrRd")[5:9])

      txt = paste(ifelse(length(HK)==1,ifelse(HK,'HK:','Regulated'),'All:'))#, c('Short','Med short','Med long','Long')[sizes])
      x<-subset(efficiencyData,GeneSize==sizes &isHK%in%HK)
      myCats = sort(unique(x$intronCountGroups))
            
      Xshift = offset*(sizes-1)
      main=ifelse(sizes==1,paste(txt,switch(exper,'Basal parameters','Add Pausing','Add Splicing variability','Add both')),'')
      # Decide which variable to test
      Form0 = myEfficiencies_plain ~ intronCountGroups
      Form = switch(exper,Form0,myEfficiencies_noSplice ~ intronCountGroups, myEfficiencies_noPause ~ intronCountGroups, myEfficiencies ~ intronCountGroups) 
    
      if (exper ==1){
        b=boxplot(myEfficiencies_plain ~ intronCountGroups, data=x,ylim=c(0,1),main=main,col=COL[myCats],boxwex=bw,at=(myCats)*wid+Xshift,xlim=xlim,names=NA,xaxt='n',notch=notch,add=(sizes-1))
        abline(h=0,col='gray')
      }else{
      at = 2*((myCats)*wid+Xshift) - (sizes-1)*wid*2
       b=boxplot(Form0, data=x,main=main,border=COL[myCats],staplewex=0,range=0.00001,outline=F,ylim=c(0,1),lwd=2,boxwex=bw*.5,at=at-wid*.7,xlim=xlim,names=NA,xaxt='n',notch=notch,add=(sizes-1))
       abline(h=0,col='gray')
       b=boxplot(Form, data=x,main=main,col=COL[myCats],boxwex=bw,at=at,names=NA,notch=notch,add=T,xaxt='n')
      }
    }
  legend('left',fill=COL,leg=c('1-2 introns','3-4 introns','5-7 introns','8-19 introns','20+ introns'))#sprintf(c('1-2 introns (n=%d)','3-4 introns (n=%d)','5-7 introns (n=%d)','8-19 introns (n=%d)','20+ introns (n=%d)')[myCats], table(subset(efficiencyData,GeneSize==sizes&isHK%in%HK)$intronCountGroups)))
  #abline(h=c(0)); grid(nx=0,ny=5)
  axis(1,at=2.5*wid+(0:3)*offset, lab = c('Short','Med short','Med long','Long')) 
  }
}
plot.off()
 
 
 

## Make plot for 'Fig 3', looking at effect of adding read-through

plot.dev("GenomeSimulations_readThrough_comparisons_chromElong_147_noAdjust_mean.85.pdf",'pdf',width=7,height=1.5)
par(mfcol=c(1,1),cex=0.5,mai=c(0.3,0.2,0.01,0.01))
notch=F; wid=.21;bw=0.07; offset = 1.4; xlim=c(-1,9.2)
cols = cols=c(brewer.pal(9,"YlGnBu")[5:9],brewer.pal(9,"OrRd")[5:9])
for (HK in list(0:1)){

  for(sizes in 1:4){
    COL = cols[c(1:5) + HK[1]*5]
    txt = paste(ifelse(length(HK)==1,ifelse(HK,'HK:','Regulated'),'All:'))#, c('Short','Med short','Med long','Long')[sizes])
    x<-subset(efficiencyData,GeneSize==sizes &isHK%in%HK)
    myCats = sort(unique(x$intronCountGroups))

    main=''
    # Decide which variable to test
    Form1 = myEfficiencies_noRT ~ intronCountGroups
    Form2 = myEfficiencies_plain_RT1000 ~ intronCountGroups
    Form3 = myEfficiencies_plain ~ intronCountGroups

    Xshift = offset*(sizes-1 )
    at = 2*(myCats*wid+Xshift) - (sizes-1)*wid*2
    lwd=1
    b=boxplot(Form1, data=x,main=main,border=COL[myCats],ylim=c(0,1),lwd=1.5,staplewex=0,range=0.00001,outline=F,boxwex=bw,at=at-wid*.5,xlim=xlim,names=NA,yaxt=ifelse(sizes==1,'s','n'),notch=notch,add=sizes-1,xaxt='n')
    if (sizes==1) abline(h=0,col='gray')
    b=boxplot(Form2, data=x,main=main,col=COL[myCats],border=gray.colors(2)[1],staplewex=0,range=0.00001,outline=F,boxwex=bw,at=at,names=NA,notch=notch,add=T,xaxt='n',lty=1,lwd=lwd,yaxt='n')
    b=boxplot(Form3, data=x,main=main,col=COL[myCats],border='black',boxwex=bw,at=at+wid*.5,names=NA,notch=notch,add=T,xaxt='n',lwd=lwd,yaxt='n')
  }
  
  legend('topleft',fill=COL,leg=c('1-2 introns','3-4 introns','5-7 introns','8-19 introns','20+ introns'))
  legend('bottomleft',fill=c(NA,COL[1],COL[1]), border=c(COL[1],gray.colors(2)[1],'black'),leg=c('No Read-Thru', '3840 bp Read-Thru',sprintf('%.1f kb Read-Thru',readThrough_for_pauseSplice/1000)))
 axis(1,at=2*(3*wid+(0:3)*offset)-(0:3)*wid*2, lab = c('Short','Med short','Med long','Long')) 
}
plot.off()



plot.dev("GenomeSimulations_numIntrons_by_geneSize.pdf",'pdf',height=12,width=12)         
par(mfcol=c(4,2))
notch=F

for (HK in list(0:1,1)){
 for(sizes in 1:4){
    COL = brewer.pal(11,'Spectral')[c(1:5) + HK[1]*6]
    txt = paste(ifelse(length(HK)==1,ifelse(HK,'HK:','Regulated'),'All:'), c('Short','Med short','Med long','Long')[sizes])
    x<-subset(efficiencyData,GeneSize==sizes &isHK%in%HK)
    myCats = sort(unique(x$intronCountGroups))
     
    wid=.21;bw=0.16; offset = 1.25
    for(i in 0:1){
      b=boxplot(myEfficiencies_plain ~ intronCountGroups, data=x,main=txt,col=COL[myCats],boxwex=bw,at=(myCats)*wid,xlim=c(-1,4.8),names=NA,xaxt='n',notch=notch,add=i)
      if(i==0) abline(h=b$stats[3,], lwd=1,lty=2,col=COL)
    }
    boxplot(myEfficiencies_noSplice ~ intronCountGroups,data=x,col=COL[myCats],boxwex=bw,at=(myCats)*wid+1*offset,names=NA,xaxt='n',add=T,notch=notch)
    boxplot(myEfficiencies_noPause ~ intronCountGroups,data=x,col=COL[myCats],boxwex=bw,at=(myCats)*wid+2*offset,names=NA,xaxt='n',add=T,notch=notch)   
    boxplot(myEfficiencies ~ intronCountGroups,data=x,col=COL[myCats],boxwex=bw,at=(myCats)*wid+3*offset,names=NA,xaxt='n',add=T,notch=notch)
    axis(1,at=2.5*wid+(0:3)*offset,lab=c('Basal parameters','Add Pausing','Add Splicing variability','Add both'))
    legend('topleft',fill=COL,leg=sprintf(c('1-2 introns (n=%d)','3-4 introns (n=%d)','5-7 introns (n=%d)','8-19 introns (n=%d)','20+ introns (n=%d)')[myCats], table(subset(efficiencyData,GeneSize==sizes&isHK%in%HK)$intronCountGroups)))
  }
 }
 plot.off()
 
plot.dev("GenomeSimulations_geneSize_by_numIntrons.pdf",'pdf',height=12,width=12)        
par(mfcol=c(5,2))
notch=F

for (HK in list(0:1,1)){
  for(sizes in 1:5){
    COL = brewer.pal(11,'Spectral')[c(1:4) + HK[1]*6]
    txt = paste(ifelse(length(HK)==1,ifelse(HK,'HK:','Regulated'),'All:'), c('1-2 introns','3-4 introns','5-7 introns','8-19 introns','20+ introns')[sizes])
    x<-subset(efficiencyData,intronCountGroups==sizes &isHK%in%HK)
    myCats = sort(unique(x$GeneSize))
     
    wid=.21;bw=0.16; offset = 1.25
    for(i in 0:1){
      b=boxplot(myEfficiencies_plain ~ GeneSize, data=x,main=txt,col=COL[myCats],boxwex=bw,at=(myCats)*wid,xlim=c(-1,4.8),names=NA,xaxt='n',notch=notch,add=i)
      if(i==0) abline(h=b$stats[3,], lwd=1,lty=2,col=COL)
    }
    boxplot(myEfficiencies_noSplice ~ GeneSize,data=x,col=COL[myCats],boxwex=bw,at=(myCats)*wid+1*offset,names=NA,xaxt='n',add=T,notch=notch)
    boxplot(myEfficiencies_noPause ~ GeneSize,data=x,col=COL[myCats],boxwex=bw,at=(myCats)*wid+2*offset,names=NA,xaxt='n',add=T,notch=notch)   
    boxplot(myEfficiencies ~ GeneSize,data=x,col=COL[myCats],boxwex=bw,at=(myCats)*wid+3*offset,names=NA,xaxt='n',add=T,notch=notch)
    axis(1,at=2.5*wid+(0:3),lab=c('Basal parameters','Add Pausing','Add Splicing variability','Add both'))
    legend('topleft',fill=COL,leg=sprintf(c('Short (n=%d)','Med short (n=%d)','Med long (n=%d)','Long (n=%d)')[myCats], table(subset(efficiencyData,intronCountGroups==sizes&isHK%in%HK)$GeneSize)))
  }
}
plot.off()


  #for(sizes in 1:4){
#    #COL = cols[1:4+colorOffset + (sizes-1)*2]
#    COL = brewer.pal(11,'Spectral')[c(1:4) + HK*7]
#    txt = paste(ifelse(length(HK)==1,ifelse(HK,'HK:','Regulated'),'All:'),c('Short','Med short','Med long','Long')[sizes])
#  
#    wid=.21;bw=0.16
#    for(i in 0:1){
#      b=boxplot(myEfficiencies_plain ~ intronCountGroups, data=subset(efficiencyData,GeneSize==sizes &isHK %in%HK),main=txt,col=COL,boxwex=bw,at = (1:4)*wid,xlim=c(-1,3.8),names=NA,xaxt='n',notch=notch,add=i)
#      if(i==0) abline(h=b$stats[3,], lwd=1,lty=2,col=COL)
#    }
#    boxplot(myEfficiencies_noSplice ~ intronCountGroups,data=subset(efficiencyData,GeneSize==sizes&isHK%in%HK),col=COL,boxwex=bw,at=(1:4)*wid+1,names=NA,xaxt='n',add=T,notch=notch)
#    boxplot(myEfficiencies_noPause ~ intronCountGroups,data=subset(efficiencyData,GeneSize==sizes&isHK%in%HK),col=COL,boxwex=bw,at=(1:4)*wid+2,names=NA,xaxt='n',add=T,notch=notch)   
#    boxplot(myEfficiencies ~ intronCountGroups,data=subset(efficiencyData,GeneSize==sizes&isHK%in%HK),col=COL,boxwex=bw,at=(1:4)*wid+3,names=NA,xaxt='n',add=T,notch=notch)
#    axis(1,at=2.5*wid+(0:3),lab=c('Basal parameters','Add Pausing','Add Splicing variability','Add both'))
#    legend('topleft',fill=COL,leg=sprintf(c('1-3 introns (n=%d)','4-6 introns (n=%d)','7-11 introns (n=%d)','12+ introns (n=%d)'), table(subset(efficiencyData,GeneSize==sizes&isHK==HK)$intronCountGroups)))
#  }
#}
#dev.off()

#
#par(mfcol=c(4,2))
#hist(myEfficiencies_plain,100)
#hist(myEfficiencies_noSplice,100, main='AddPause')
#hist(myEfficiencies_noPause,100, main='AddSplice')
#hist(myEfficiencies,100, main='Add both')
#plot(myEfficiencies_plain,myEfficiencies_noSplice,main='Add Pause')
#plot(myEfficiencies_plain,myEfficiencies_noPause,main='Add Splice')
#plot(myEfficiencies_plain,myEfficiencies, main='Add both')
#plot(myEfficiencies_plain,myEfficiencies,col=cols2[genes.Size])
#legend('topleft',col=cols2[1:4],pch=1,leg=paste('Gene Size',1:4))
#plot(myEfficiencies_plain,myEfficiencies,col=genes.isHK+1)

#par(mfcol=c(4,2)) 
##### Intron count groups
#for(HK in c(0,1)){
#  COL = brewer.pal(11,'Spectral')[c(1,4,8,10)+HK]
#  txt =ifelse(HK,'HK','non')
#
#  wid=.21;bw=0.16
#  for(i in 0:1){
#    b=boxplot(myEfficiencies_plain ~ intronCountGroups, data=subset(efficiencyData,isHK==HK),main=txt,col=COL,boxwex=bw,at = (1:4)*wid,xlim=c(0,4),names=NA,xaxt='n',notch=T,add=i)
#    if(i==0) abline(h=b$stats[3,], lwd=2,lty=2,col=COL)
#  }
#  boxplot(myEfficiencies_noSplice ~ intronCountGroups,data=subset(efficiencyData,isHK==HK),col=COL,boxwex=bw,at=(1:4)*wid+1,names=NA,xaxt='n',add=T,notch=T)
#  boxplot(myEfficiencies_noPause ~ intronCountGroups,data=subset(efficiencyData,isHK==HK),col=COL,boxwex=bw,at=(1:4)*wid+2,names=NA,xaxt='n',add=T,notch=T)
#  boxplot(myEfficiencies ~ intronCountGroups,data=subset(efficiencyData,isHK==HK),col=COL,boxwex=bw,at=(1:4)*wid+3,names=NA,xaxt='n',add=T,notch=T)
#  axis(1,at=2.5*wid+(0:3),lab=c('Basal parameters','Add Pausing','Add Splicing variability','Add both'))
#  legend('topleft',col=COL,pch=1,leg=c('1-3 introns','4-6 introns','7-11 introns','12+ introns'))
#}
#
##### Gene size
#for(HK in c(0,1)){
#
#  COL = brewer.pal(11,'Spectral')[c(1,4,8,10)+HK] 
#  txt =ifelse(HK,'HK','non')
#
#  wid=.21;bw=0.16
#  for(i in 0:1){
#    b=boxplot(myEfficiencies_plain ~ GeneSize, data=subset(efficiencyData,isHK==HK),main=txt,col=COL,boxwex=bw,at = (1:4)*wid,xlim=c(0,4),names=NA,xaxt='n',notch=T,add=i)
#    if(i==0) abline(h=b$stats[3,], lwd=2,lty=2,col=COL)
#  }
#  boxplot(myEfficiencies_noSplice ~ GeneSize,data=subset(efficiencyData,isHK==HK),col=COL,boxwex=bw,at=(1:4)*wid+1,names=NA,xaxt='n',add=T,notch=T)
#  boxplot(myEfficiencies_noPause ~ GeneSize,data=subset(efficiencyData,isHK==HK),col=COL,boxwex=bw,at=(1:4)*wid+2,names=NA,xaxt='n',add=T,notch=T)
#  boxplot(myEfficiencies ~ GeneSize,data=subset(efficiencyData,isHK==HK),col=COL,boxwex=bw,at=(1:4)*wid+3,names=NA,xaxt='n',add=T,notch=T)
#  axis(1,at=2.5*wid+(0:3),lab=c('Basal parameters','Add Pausing','Add Splicing variability','Add both'))
#  legend('topleft',col=COL,pch=1,leg=c('Short','Med short','Med long','Long'))
#}



###############################################################################
###############################################################################
#### 8/26/13
###############################################################################
###############################################################################

# Adding AGE of gene ... 
#ENS2refseq = read.delim("../ensemble_to_refseq_human.txt",stringsAsFactors=F)
#ens2refseq = with(ENS2refseq,data.frame(ensTx = sub("\\.[0-9]*$", "", X.hg19.wgEncodeGencodeBasicV17.name), ensGene = sub("\\.[0-9]*$", "", hg19.wgEncodeGencodeAttrsV17.geneId), refseq = hg19.wgEncodeGencodeRefSeqV17.rnaAcc, stringsAsFactors=F))
#ens2refseq = subset(ens2refseq, refseq!= 'n/a')
#
#ens2refseq2 = sapply(1:nrow(ens2refseq),function(x) cbind(ens2refseq[x,1:2],refseq=sub("\\.[0-9]*$","",unlist(strsplit(ens2refseq[x,3],",")))),simplify = F )
#ens2refseq3 = matrix(unlist(lapply((ens2refseq2), t)),nc=3,byrow=T, dimnames=list(NULL,c('ensTx','ensGene','refseq')))
#ens2refseq3=ens2refseq3[!duplicated(ens2refseq3),]
#write.delim(ens2refseq3, File='../ensemble_to_refseq_simplified_human.txt')

ref2ENS = read.delim('../ensemble_to_refseq_simplified_human.txt', stringsAsFactors=F)
ref2ENS = ref2ENS[!duplicated(ref2ENS[,2:3]),]
ensAge = read.delim('../Gene_age_Zhang2010_human.txt', stringsAsFactors=F)
ref2ENS$AgeCluster = ensAge$branch[match(ref2ENS$ensGene, ensAge$gene)]

## Figure out if there are any refseqs that got assigned more than one age
table(a<-tapply(ref2ENS$AgeCluster, ref2ENS$refseq, function(x)length(unique(x[!is.na(x)]))))
ref2ENS = subset(ref2ENS,! refseq%in% names(a)[a>1] & !is.na(AgeCluster))

##### NOW APPLY THIS TO OTHER DATASETS ############
refseq$AgeCluster = ref2ENS$AgeCluster[match(refseq$name, ref2ENS$refseq)]
mergedData$AgeCluster = refseq$AgeCluster[match(mergedData$UniqueID, refseq$UniqueID)]

sapply(0:12,function(x)t.test(Stability ~ isLast, data=subset(mergedData,AgeCluster==x))$p.va)->stab_p
sapply(0:12,function(x)t.test(Acceptor ~ isLast, data=subset(mergedData,AgeCluster==x))$p.va)->stab_a

###############################################################################
# Compare differences in LAST exon stability/ intron strength as a function of gene age:
###############################################################################
acceptorComparison_by_age = sapply(0:12, function(x)with(subset(mergedData,AgeCluster==x), median(Acceptor[isLast==1], na.rm=T)-median(Acceptor[isLast==0], na.rm=T) ))
acceptorLast_by_age = sapply(0:12, function(x)with(subset(mergedData,AgeCluster==x), median(Acceptor[isLast==1], na.rm=T)))
stabilityComparison_by_age = sapply(0:12, function(x)with(subset(mergedData,AgeCluster==x), median(Stability[isLast==1], na.rm=T)-median(Stability[isLast==0], na.rm=T) ))
stabilityLast_by_age = sapply(0:12, function(x)with(subset(mergedData,AgeCluster==x), median(Stability[isLast==1], na.rm=T)))
stabilityNonLast_by_age = sapply(0:12, function(x)with(subset(mergedData,AgeCluster==x), median(Stability[isLast==0], na.rm=T)))

numIntrons_by_age = tapply(refseq$exonCount,refseq$AgeCluster,median, na.rm=T)

#acceptorComparison_by_age = sapply(0:12, function(x)with(subset(mergedData,AgeCluster==x), median(Acceptor[isLast==1], na.rm=T)-median(Acceptor[isLast==0], na.rm=T) ))

# I got these numbers by using Engauge digitizer on fig 2 of Zhang 2010:
myr=c(-480.979,-392.321,-324.815,-243.678,-178.203,-127.053,-84.8575,-76.5389,-54.5904,-29.3928,-17.749, -10.6576,-2.74751 )

cor.test(myr,stabilityComparison_by_age)
cor.test(myr,acceptorComparison_by_age)

# Other: get data from nuc. occupancy ...

pdf('GeneAge_medians_comparisons3.pdf')
par(mfrow=c(2,2), cex=0.8)

#### First plot both non-last and last stability
#plot(myr,stabilityNonLast_by_age,cex.main=1,ylab= 'Median Nucleosome Stability Score', ylim=c(-.7,.25),pch=19, xlab='Million Years that gene has existed', main='Stability Scores in Exons')
#points(myr,stabilityLast_by_age,col='green',pch=19)
#legend('bottom',pch=19,col=c('black','green'),leg=c('Internal exons','Last exons'))

## Then plot Delta Stability
plot(myr,stabilityComparison_by_age,cex.main=.8, main='Median Nucleosome Stability Scores in last - internal exons',ylab= expression(paste(Delta,'Stability Score')), ylim=c(-.25,.25),pch=19, xlab='Million Years that gene has existed')
abline(h=0)
points(myr,stabilityComparison_by_age+diff(par('usr')[3:4])/30,pch=ifelse(stab_p < 0.01,'*',NA),cex=1.2)

## Then plot Delta Acceptor
plot(myr,acceptorLast_by_age,cex.main=.8, main='Median Acceptor Scores in last - nonlast introns',ylim=c(.15,.55),ylab= expression(paste(Delta,'Acceptor Score')),pch=19, xlab='Million Years that gene has existed')
abline(h=0)
points(myr,acceptorLast_by_age+diff(par('usr')[3:4])/30,pch=ifelse(stab_a < 0.01,'*',NA),cex=1.2)

## Then plot Number of exons per gene ...
boxplot(MaxFeature ~ AgeCluster,data=mergedData, range=0.001, outline=F, notch=T,xaxt='n',cex.main=1, ylab='# Introns per gene')
axis(1,at=1:13,lab=round(myr), las=2)

dev.off()










###############################################################################
###############################################################################
#### Looking at Pile-up of Pol, PolS2.. August 2013
###############################################################################
###############################################################################
library(rtracklayer)
library(RColorBrewer)
source("/home/jeremy/ExonPipeline/analyzeBigWig.r")
load('ENCODE_PolCHIP_250.RData')

### IMPORT .BED FILES
dist=250
bed_nonLastPlus   = import.bed(sprintf("/home/home/jeremy/ExonPipeline/hg19/Human_exons_-%dTo%d_nonLast+.bed",dist,dist),asRangedData=FALSE)
bed_nonLastMinus  = import.bed(sprintf("/home/home/jeremy/ExonPipeline/hg19/Human_exons_-%dTo%d_nonLast-.bed",dist,dist),asRangedData=FALSE)
bed_lastPlus      = import.bed(sprintf("/home/home/jeremy/ExonPipeline/hg19/Human_exons_-%dTo%d_last+.bed",dist,dist),asRangedData=FALSE)
bed_lastMinus     = import.bed(sprintf("/home/home/jeremy/ExonPipeline/hg19/Human_exons_-%dTo%d_last-.bed",dist,dist),asRangedData=FALSE)
bed_last = c(bed_lastPlus,bed_lastMinus)
bed_nonLast = c(bed_nonLastPlus,bed_nonLastMinus)

bed_last$UniqueID = sub(">([0-9]*)_.*","\\1",bed_last$name)
bed_nonLast$UniqueID = sub(">([0-9]*)_.*","\\1",bed_nonLast$name)


allBeds_selection = make_beds_selection(c(bed_lastPlus, bed_lastMinus, bed_nonLastPlus, bed_nonLastMinus))

## Now do this for the polyA sites of genes only
dist2=1000
bed_endPlus   = import.bed(sprintf("/home/home/jeremy/ExonPipeline/hg19/Human_polyA_-%dTo%d+.bed",dist2,dist2),asRangedData=FALSE)
bed_endMinus  = import.bed(sprintf("/home/home/jeremy/ExonPipeline/hg19/Human_polyA_-%dTo%d-.bed",dist2,dist2),asRangedData=FALSE)
allBeds_selection.ends = make_beds_selection(bed_ends<-c(bed_endPlus, bed_endMinus))
bed_ends$UniqueID = sub(">([0-9]*)_.*","\\1",bed_ends$name)





###############################################################################
### IMPORT .BIGWIG FILES
## Pol ChIP
bw_pol      = import.bw("/home/RNAseq/ENCODE/wgEncodeSydhTfbsK562Pol2IggmusSig.bigWig", selection=allBeds_selection,asRangedData=FALSE)
bw_polS2    = import.bw("/home/RNAseq/ENCODE/wgEncodeSydhTfbsK562Pol2s2IggrabSig.bigWig",selection=allBeds_selection,asRangedData=FALSE)

RangeScores_pol = cbind(as.matrix(ranges(bw_pol)),score(bw_pol));
RangeScores_polS2 = cbind(as.matrix(ranges(bw_polS2)),score(bw_polS2));

# Also for different .bed files 
bw_pol.ends      = import.bw("/home/RNAseq/ENCODE/wgEncodeSydhTfbsK562Pol2IggmusSig.bigWig", selection=allBeds_selection.ends,asRangedData=FALSE)
bw_polS2.ends    = import.bw("/home/RNAseq/ENCODE/wgEncodeSydhTfbsK562Pol2s2IggrabSig.bigWig",selection=allBeds_selection.ends,asRangedData=FALSE)

RangeScores_pol.ends = cbind(as.matrix(ranges(bw_pol.ends)),score(bw_pol.ends));
RangeScores_polS2.ends = cbind(as.matrix(ranges(bw_polS2.ends)),score(bw_polS2.ends));

  
##  Nucleosomes
bw_NucK562 = import.bw('/home/RNAseq/ENCODE/wgEncodeSydhNsomeK562Sig.bigWig',selection=allBeds_selection,asRangedData=F)
bw_NucGm12878 = import.bw('/home/RNAseq/ENCODE/wgEncodeSydhNsomeGm12878Sig.bigWig',selection=allBeds_selection,asRangedData=F)

RangeScores_NucK562 = cbind(as.matrix(ranges(bw_NucK562)),score(bw_NucK562));
RangeScores_NucGm12878 = cbind(as.matrix(ranges(bw_NucGm12878)),score(bw_NucGm12878));
###############################################################################

###############################################################################
## Compute overlaps
###############################################################################

################## Pol ChIP : intron/exon junctions

#date()
#scores_lastPol = scoreBW_bed(bw_pol, bed_last, RangeScores_pol)#scores_lastPol_size->scores_lastPol
#scores_lastPolS2 = scoreBW_bed(bw_polS2, bed_last, RangeScores_polS2)
#colnames(scores_lastPol) = colnames(scores_lastPolS2) = bed_last$UniqueID
#
#scores_nonLastPol = scoreBW_bed(bw_pol, bed_nonLast, RangeScores_pol)#scores_lastPol_size->scores_lastPol
#scores_nonLastPolS2 = scoreBW_bed(bw_polS2, bed_nonLast, RangeScores_polS2)
#colnames(scores_nonLastPol) = colnames(scores_nonLastPolS2) = bed_nonLast$UniqueID
#
#date()
#
##save.image(file='CHIP.RData')
##save(scores_lastPol, scores_lastPolS2,scores_nonLastPol, scores_nonLastPolS2,file=sprintf('ENCODE_PolCHIP_%d.RData', dist))
# 
#
################### Pol ChIP: 3' ends         
#
#### Overlap with BW files:
#scores_lastPol.ends = scoreBW_bed(bw_pol.ends, bed_ends, RangeScores_pol.ends)
#scores_lastPolS2.ends = scoreBW_bed(bw_polS2.ends, bed_ends, RangeScores_polS2.ends)
#colnames(scores_lastPol.ends) = colnames(scores_lastPolS2.ends) = bed_ends$UniqueID
#save(scores_lastPol.ends, scores_lastPolS2.ends,file='ENCODE_PolCHIP_ends1000.RData')
# 
################ Nucleosomes: intron/exon junctions
#    
## K562
#scores_last_nucK562 = scoreBW_bed(bw_NucK562, bed_last, RangeScores_NucK562)
#scores_nonLast_nucK562 = scoreBW_bed(bw_NucK562, bed_nonLast, RangeScores_NucK562)
# 
#
#nucAverages1 = cbind(last_K562 = rowMeans(scores_last_nucK562),nonLast_K562 = rowMeans(scores_nonLast_nucK562))
#
## Gm12878
#scores_last_nucGm12878 = scoreBW_bed(bw_NucGm12878, bed_last, RangeScores_NucGm12878)
#scores_nonLast_nucGm12878 = scoreBW_bed(bw_NucGm12878, bed_nonLast, RangeScores_NucGm12878)
#
#nucAverages2 = cbind(last_Gm12878 = rowMeans(scores_last_nucGm12878),nonLast_Gm12878 = rowMeans(scores_nonLast_nucGm12878))
#
# colnames(scores_last_nucGm12878) = colnames(scores_last_nucK562)= bed_last$UniqueID
# colnames(scores_nonLast_nucGm12878) = colnames(scores_nonLast_nucK562)= bed_nonLast$UniqueID
##save(scores_last_nucGm12878,scores_nonLast_nucGm12878,scores_last_nucK562,scores_nonLast_nucK562,nucAverages1,nucAverages2,
##  file=sprintf('ENCODE_Nuc_occupancy_%d.RData', dist))
load("ENCODE_Nuc_occupancy_250.RData")

###### GATHER STATISTICS FOR PLOTTING ###########

#### NUCLEOSOMES
# First need to concatenate AND make vector of UniqueIDs
scores_nucK562 = cbind(scores_last_nucK562, scores_nonLast_nucK562)
scores_nucGm12878 = cbind(scores_last_nucGm12878, scores_nonLast_nucGm12878)
ids.bed = c(bed_last$UniqueID, bed_nonLast$UniqueID)

## Average over sizes
averages_nucK562_size = sapply(1:4, function(x)rowMeans(scores_nucK562[,ids.bed %in% subset(mergedData,GeneSize==x)$UniqueID]))
averages_nucGm12878_size = sapply(1:4, function(x)rowMeans(scores_nucGm12878[,ids.bed %in% subset(mergedData,GeneSize==x)$UniqueID]))

## Average over gene size/# introns
averages_nucK562_size_numIntrons = sapply(1:4,function(size)sapply(1:5, function(num)rowMeans(scores_nucK562[,ids.bed  %in% subset(mergedData,GeneSize==
size & intronCountGroups==num)$UniqueID])),simplify=F)

averages_nucGm12878_size_numIntrons = sapply(1:4,function(size)sapply(1:5, function(num)rowMeans(scores_nucGm12878[,ids.bed  %in% subset(mergedData,GeneSize==
size & intronCountGroups==num)$UniqueID])),simplify=F)

averages_nucK562_HK_last = sapply(0:1, function(x)rowMeans(scores_last_nucK562[,colnames(scores_last_nucK562) %in% subset(mergedData,isHK==x)$UniqueID]))
averages_nucGm12878_HK_last = sapply(0:1, function(x)rowMeans(scores_last_nucGm12878[,colnames(scores_last_nucGm12878) %in% subset(mergedData,isHK==x)$UniqueID]))

averages_nucK562_HK = sapply(0:1, function(x)rowMeans(scores_nucK562[,ids.bed %in% subset(mergedData,isHK==x)$UniqueID]))
averages_nucGm12878_HK = sapply(0:1, function(x)rowMeans(scores_nucGm12878[,ids.bed %in% subset(mergedData,isHK==x)$UniqueID]))


############# POL
#load('CHIP.RData')

## Average over gene sizes:
averages_lastPol_size = sapply(1:4, function(x)rowMeans(scores_lastPol[,colnames(scores_lastPol) %in% subset(mergedData,GeneSize==x)$UniqueID])) 
averages_lastPolS2_size = sapply(1:4, function(x)rowMeans(scores_lastPolS2[,colnames(scores_lastPolS2) %in% subset(mergedData,GeneSize==x)$UniqueID]))
averages_nonLastPol_size = sapply(1:4, function(x)rowMeans(scores_nonLastPol[,colnames(scores_nonLastPol) %in% subset(mergedData,GeneSize==x &)$UniqueID])) 
averages_nonLastPolS2_size = sapply(1:4, function(x)rowMeans(scores_nonLastPolS2[,colnames(scores_nonLastPolS2) %in% subset(mergedData,GeneSize==x)$UniqueID]))


## Average over gene size/# introns
averages_lastPol_size_numIntrons = sapply(1:4,function(size)sapply(1:5, function(num)rowMeans(scores_lastPol[,colnames(scores_lastPol) %in% subset(mergedData,GeneSize==size & intronCountGroups==num)$UniqueID])),simplify=F) 
averages_lastPolS2_size_numIntrons = sapply(1:4,function(size)sapply(1:5, function(num)rowMeans(scores_lastPolS2[,colnames(scores_lastPolS2) %in% subset(mergedData,GeneSize==size & intronCountGroups==num)$UniqueID])),simplify=F) 
averages_nonLastPol_size_numIntrons = sapply(1:4,function(size)sapply(1:5, function(num)rowMeans(scores_nonLastPol[,colnames(scores_nonLastPol) %in% subset(mergedData,GeneSize==size & intronCountGroups==num)$UniqueID])),simplify=F) 
averages_nonLastPolS2_size_numIntrons = sapply(1:4,function(size)sapply(1:5, function(num)rowMeans(scores_nonLastPolS2[,colnames(scores_nonLastPolS2) %in% subset(mergedData,GeneSize==size & intronCountGroups==num)$UniqueID])),simplify=F) 

averages_lastPol_size_numIntrons2 = sapply(1:4,function(size)sapply(1:5, function(num)rowMeans(scores_lastPol[,colnames(scores_lastPol) %in% subset(mergedData,GeneSize==size & intronCountGroups==num)$UniqueID & bed_last$name %in% subset(mergedData,FeatureCount>1)$energy_ID])),simplify=F) 
averages_lastPolS2_size_numIntrons2 = sapply(1:4,function(size)sapply(1:5, function(num)rowMeans(scores_lastPolS2[,colnames(scores_lastPolS2) %in% subset(mergedData,GeneSize==size & intronCountGroups==num)$UniqueID & bed_last$name %in% subset(mergedData,FeatureCount>1)$energy_ID])),simplify=F) 
averages_nonLastPol_size_numIntrons2 = sapply(1:4,function(size)sapply(1:5, function(num)rowMeans(scores_nonLastPol[,colnames(scores_nonLastPol) %in% subset(mergedData,GeneSize==size & intronCountGroups==num)$UniqueID & bed_nonLast$name %in% subset(mergedData,FeatureCount>1)$energy_ID])),simplify=F) 
averages_nonLastPolS2_size_numIntrons2 = sapply(1:4,function(size)sapply(1:5, function(num)rowMeans(scores_nonLastPolS2[,colnames(scores_nonLastPolS2) %in% subset(mergedData,GeneSize==size & intronCountGroups==num)$UniqueID & bed_nonLast$name %in% subset(mergedData,FeatureCount>1)$energy_ID])),simplify=F) 


# Total genes:
scores_pol = cbind(scores_lastPol, scores_nonLastPol)
scores_polS2 = cbind(scores_lastPolS2, scores_nonLastPolS2)
bed_all = c(bed_last, bed_nonLast)

############# POL: 3' END
averages_lastPol_size.ends = sapply(1:4, function(x)rowMeans(scores_lastPol.ends[,colnames(scores_lastPol.ends) %in% subset(mergedData,GeneSize==x)$UniqueID])) 
averages_lastPolS2_size.ends = sapply(1:4, function(x)rowMeans(scores_lastPolS2.ends[,colnames(scores_lastPolS2.ends) %in% subset(mergedData,GeneSize==x)$UniqueID]))

averages_lastPol_size_numIntrons.ends = sapply(1:4,function(size)sapply(1:5, function(num)rowMeans(scores_lastPol.ends[,colnames(scores_lastPol.ends) %in% subset(mergedData,GeneSize==size & intronCountGroups==num)$UniqueID])),simplify=F) 
averages_lastPolS2_size_numIntrons.ends = sapply(1:4,function(size)sapply(1:5, function(num)rowMeans(scores_lastPolS2.ends[,colnames(scores_lastPolS2.ends) %in% subset(mergedData,GeneSize==size & intronCountGroups==num)$UniqueID])),simplify=F) 

averages_lastPolS2_HK.ends = sapply(0:1, function(x)rowMeans(scores_lastPolS2.ends[,colnames(scores_lastPolS2.ends) %in% subset(mergedData,isHK==x)$UniqueID]))



###### PLOTTING ###########
myCols = brewer.pal(10,'Paired')
myCols4 = brewer.pal(4,'Set2')
myX = (1:(2*dist)) - dist
myCols5 = cols=(brewer.pal(9,"YlGnBu"))[5:9]
        
plotSize_numIntrons = function(SIZE, SIZE_NUM, main='',dist=250,xlab="bp from 3' splice site",ylim=NULL,...){
  par(mfrow=c(2,3),cex=0.75)

  myX = (1:(2*dist)) - dist
	  
  matplot(myX, SIZE,pch=19,type='l',lwd=4,col=myCols4, lty=1,xlab=xlab,ylab='average read coverage',ylim=ylim, ...)
  abline(v=0)
  
  plot(1,col='white',xaxt='n',yaxt='n',xlab='',ylab='', bty='n',main=main)
  legend('topleft',col=myCols4,lwd=3,leg=c('Short','Med. short','Med. long','Long'))
  legend('bottomright',col=myCols5,lwd=3,leg=c('1-2 introns','3-4 introns','5-7 introns','8-19 introns','20+ introns'))
  
  invisible(sapply(1:4, function(size){
    matplot(myX, SIZE_NUM[[size]],ylim=ylim,xlab=xlab,ylab='average read coverage',main=c('Short','Med. Short','Med. Long','Long')[size],pch=19,type='l',lwd=4,col=myCols5, lty=1)
    grid()
    abline(v=0)
  }))

}
dist2=1000

pdf("ENCODE_pol_ChIP_last.pdf", height=8,width=12)
plotSize_numIntrons(averages_lastPol_size,averages_lastPol_size_numIntrons,'Pol ChIP: \nCoverage over final intron-exon junction',ylim=c(0,4.1))
dev.off()

pdf("ENCODE_polS2_ChIP_last.pdf", height=8,width=12)
plotSize_numIntrons(averages_lastPolS2_size,averages_lastPolS2_size_numIntrons,'PolS2 ChIP: \nCoverage over final intron-exon junction', ylim=c(0,35))
dev.off()

pdf("ENCODE_pol_ChIP_nonLast.pdf", height=8,width=12)
plotSize_numIntrons(averages_nonLastPol_size,averages_nonLastPol_size_numIntrons,'Pol ChIP: \nCoverage over non-last intron-exon junctions',ylim=c(0,4.1))
dev.off()

pdf("ENCODE_polS2_ChIP_nonLast.pdf", height=8,width=12)
plotSize_numIntrons(averages_nonLastPolS2_size,averages_nonLastPolS2_size_numIntrons,'PolS2 ChIP: \nCoverage over non-last intron-exon junctions', ylim=c(0,35))
dev.off()


plotSize_numIntrons2 = function(SIZE_NUM_last, SIZE_NUM_nonLast, main='',dist=250,xlab="bp from 3' splice site",ylim=NULL,...){
  par(cex=0.75)
  layout(heights=c(1,1,1),widths=c(1,1,1),mat=matrix(c(1,2,5,3,4,5),nr=2,byrow=T))

  myX = (1:(2*dist)) - dist
  myTest=c(1:5)
  
  invisible(sapply(1:4, function(size){
    matplot(myX, SIZE_NUM_last[[size]][,myTest],ylim=ylim,xlab=xlab,ylab='average read coverage',main=c('Short','Med. Short','Med. Long','Long')[size],pch=19,type='l',lwd=3,col=myCols5[myTest], lty=1)
    matplot(myX, SIZE_NUM_nonLast[[size]][,myTest],ylim=ylim,xlab=xlab,ylab=,pch=19,type='l',lwd=2,col=myCols5[myTest], lty=2,add=T)
    grid()
    abline(v=0)
    
  }))

   plot(1,col='white',xaxt='n',yaxt='n',xlab='',ylab='', bty='n',main=main)
   legend('topleft',cex=1.2,col=rep(myCols5,each=2),lwd=3,lty = rep(c(1,2),5),leg=c('1-2 introns',('(non-last)'),'3-4 introns',('(non-last)'),'5-7 introns',('(non-last)'),'8-19 introns',('(non-last)'),'20+ introns',('(non-last)')))      


}


## Plot both last/nonlast PolS2 signals, with for both last/internal exon junctions
pdf("ENCODE_pol_ChIP_all.pdf", height=8,width=12)
plotSize_numIntrons2(averages_lastPol_size_numIntrons2,averages_nonLastPol_size_numIntrons2,   'Pol ChIP: \nCoverage over intron-exon junctions',ylim=c(0,6))
dev.off()

pdf("ENCODE_polS2_ChIP_all.pdf", height=8,width=12)
plotSize_numIntrons2(averages_lastPolS2_size_numIntrons2,averages_nonLastPolS2_size_numIntrons2,'PolS2 ChIP: \nCoverage over intron-exon junction  for internal exons', ylim=c(0,35))
dev.off()

## Here I'm allowing use of exon/intron 1
pdf("ENCODE_pol_ChIP_ends.pdf", height=8,width=12)
plotSize_numIntrons(averages_lastPol_size.ends,averages_lastPol_size_numIntrons.ends,'Pol ChIP: \nCoverage over the poly(A) site',xlab='',ylim=c(0,10),dist=dist2)
dev.off()

pdf("ENCODE_polS2_ChIP_ends.pdf", height=8,width=12)
plotSize_numIntrons(averages_lastPolS2_size.ends,averages_lastPolS2_size_numIntrons.ends,'PolS2 ChIP: \nCoverage over the poly(A) site',xlab='', ylim=c(0,70),dist=dist2)
dev.off()

plot.dev("ENCODE_polS2_ChIP_ends_HK_comparison.pdf",'pdf',height=3,width=3)
par(mai=c(.8,.5,0.1,0.3))
dist=1000;myX = (1:(2*dist)) - dist
matplot(myX, averages_lastPolS2_HK.ends,ylab='PolS2 read coverage', ylim=c(0,80),xlab='bp from start of final exon',pch=19,type='l',col=1:2, lty=1,lwd=2) 
lines(myX, averages_lastPolS2_HK.ends[,1]*30/ averages_lastPolS2_HK.ends[1000,1],lty=2,lwd=2)
legend('topleft',col=c(1,1:2),lwd=2,lty=c(1,2,1),leg=c('non-HK','normalized non-HK','HK'))      
plot.off()




## Do this also for HK vs. non-HK genes:


for(HK in c(0,1)){
  lab=ifelse(HK,'HK','non')
  lastPol = sapply(1:4,function(size)sapply(1:5, function(num)rowMeans(scores_lastPol[,bed_last$name %in% subset(mergedData,FeatureCount>1 & GeneSize==size & intronCountGroups==num & isHK==HK)$energy_ID])),simplify=F) 
  lastPolS2 = sapply(1:4,function(size)sapply(1:5, function(num)rowMeans(scores_lastPolS2[,bed_last$name %in% subset(mergedData,FeatureCount>1 & GeneSize==size & intronCountGroups==num & isHK==HK)$energy_ID])),simplify=F) 
  nonLastPol = sapply(1:4,function(size)sapply(1:5, function(num)rowMeans(scores_nonLastPol[,bed_nonLast$name %in% subset(mergedData,FeatureCount>1 & GeneSize==size & intronCountGroups==num & isHK==HK)$energy_ID])),simplify=F) 
  nonLastPolS2 = sapply(1:4,function(size)sapply(1:5, function(num)rowMeans(scores_nonLastPolS2[,bed_nonLast$name %in% subset(mergedData,FeatureCount>1 & GeneSize==size & intronCountGroups==num & isHK==HK)$energy_ID])),simplify=F) 
  
  pdf(sprintf("ENCODE_pol_ChIP_all_%s.pdf",lab), height=8,width=12)
  plotSize_numIntrons2(lastPol,nonLastPol,   'Pol ChIP: \nCoverage over intron-exon junctions',ylim=c(0,13))
  dev.off()
  
  pdf(sprintf("ENCODE_polS2_ChIP_all_%s.pdf",lab), height=8,width=12)
  plotSize_numIntrons2(lastPolS2,nonLastPolS2,'PolS2 ChIP: \nCoverage over intron-exon junction  for internal exons', ylim=c(0,50))
  dev.off()
}

## Now take the means of the COLUMNS, which gives me average signal for each gene
# also, ONLY take the first 250 bp, NOT the upstream intronic part
getStats = function(mat){
  #cbind(apply(mat,2,mean,na.rm=T), apply(mat,2,sd,na.rm=T)/sqrt(nrow(mat)))
  cbind(sapply(mat,mean,na.rm=T), sapply(mat,sd,na.rm=T)/sqrt(sapply(mat,length)))
  #cbind(sapply(mat,mean,na.rm=T), sapply(mat,sd,na.rm=T),sapply(mat,length))
}
propagateDivision = function(A,B){
  cbind(A[,1]/B[,1], sqrt(A[,2]^2 + (B[,2]*A[,1])^2)/B[,1] )
}

propagateSubtraction = function(A,B){
  cbind(A[,1]-B[,1], sqrt(A[,2]^2 + B[,2]^2))
}


extent = dist+(1:dist)
for(HK in list(0,1,c(0,1))){

  lab = ifelse(length(HK)==2,'All',ifelse(HK==1,'HK','non'))
  featureStart =1
  lastPol_ = sapply(1:4,function(size)sapply(1:5, function(num)colMeans(scores_lastPol[extent,bed_last$name %in% subset(mergedData,FeatureCount>=featureStart & GeneSize==size & intronCountGroups==num & isHK %in% HK)$energy_ID])),simplify=F) 
  lastPolS2_ = sapply(1:4,function(size)sapply(1:5, function(num)colMeans(scores_lastPolS2[extent,bed_last$name %in% subset(mergedData,FeatureCount>=featureStart & GeneSize==size & intronCountGroups==num & isHK %in% HK)$energy_ID])),simplify=F) 
  nonLastPol_ = sapply(1:4,function(size)sapply(1:5, function(num)colMeans(scores_nonLastPol[extent,bed_nonLast$name %in% subset(mergedData,FeatureCount>=featureStart & GeneSize==size & intronCountGroups==num & isHK %in% HK)$energy_ID])),simplify=F) 
  nonLastPolS2_ = sapply(1:4,function(size)sapply(1:5, function(num)colMeans(scores_nonLastPolS2[extent,bed_nonLast$name %in% subset(mergedData,FeatureCount>=featureStart & GeneSize==size & intronCountGroups==num & isHK %in% HK)$energy_ID])),simplify=F) 
 
  p.polS2 = sapply(1:4,function(size)sapply(1:5, function(num){
    as.numeric(try(t.test(lastPolS2_[[size]][[num]],nonLastPolS2_[[size]][[num]])$p.val))
 }))
 p.pol = sapply(1:4,function(size)sapply(1:5, function(num){
    as.numeric(try(t.test(lastPol_[[size]][[num]],nonLastPol_[[size]][[num]])$p.val))
 }))
  
 pdf(sprintf("ENCODE_ChIP_%s_last-v-nonLast_exonsOnly.pdf",lab), height=6, width=9)
 at = (1:23)[-c(6,12,18)]
 boxplot(unlist(lastPol_, recursive=F),range=0.001,outline=F, ylab='Average Read Coverage',col=myCols5,at=at,boxwex=.4, notch=T, ylim=c(0,6), main=sprintf("Pol ChIP: %s genes ", lab),xaxt='n')
 boxplot(unlist(nonLastPol_, recursive=F),range=0.001,outline=F, border=myCols5,col='white',lwd=2,at=at-.4,boxwex=.25,add=T,xaxt='n', notch=T)
 points(at,rep(0,length(at)),pch=ifelse(p.pol<0.01,"*",NA)) 
 
 b = boxplot(unlist(lastPolS2_, recursive=F),range=0.001,outline=F, ylab='Average Read Coverage', col=myCols5, at=at,boxwex=.4, ylim=c(0,50), notch=T, main=sprintf("PolS2 ChIP: %s genes ", lab),xaxt='n')
 boxplot(unlist(nonLastPolS2_, recursive=F),range=0.001,outline=F, border=myCols5,col='white',lwd=2,at=at-.4,boxwex=.25,add=T, notch=T,xaxt='n')
 points(at,rep(0,length(at)),pch=ifelse(p.polS2<0.01,"*",NA)) 
 dev.off()
 
}
  
#  # Get stats
#  pol_stats = sapply(1:4, function(x){
#    stats_lastPol = getStats(lastPol_[[x]])
#    stats_nonLastPol = getStats(nonLastPol_[[x]])    
#    cbind(propagateDivision(stats_lastPol,stats_nonLastPol)    ,stats_lastPol,stats_nonLastPol)
#  })#,simplify=F)
#  
#  polS2_stats = sapply(1:4, function(x){
#    stats_lastPolS2 = getStats(lastPolS2_[[x]])
#    stats_nonLastPolS2 = getStats(nonLastPolS2_[[x]]) 
#    cbind(propagateDivision(stats_lastPolS2,stats_nonLastPolS2) , propagateSubtraction(stats_lastPolS2,stats_nonLastPolS2),stats_lastPolS2,stats_nonLastPolS2)    
#  })#,simplify=F)
# 

####################################
### Nucleosomes:
####################################

pdf("ENCODE_nuc_K562.pdf", height=8,width=12)
plotSize_numIntrons(averages_nucK562_size,averages_nucK562_size_numIntrons,'Nucleosomes in K562: \nCoverage over final intron-exon junction',ylim=c(0.5,1.5))
dev.off()

pdf("ENCODE_nuc_Gm12878.pdf", height=8,width=12)
plotSize_numIntrons(averages_nucGm12878_size,averages_nucGm12878_size_numIntrons,'Nucleosomes in Gm12878: \nCoverage over final intron-exon junction', ylim=c(0.5,1.5))
dev.off()


plot.dev('ENCODE_nucleosomes_lastNonlast_withHK.pdf','pdf',height=12,width=8)
dist=250;myX = (1:(2*dist)) - dist
par(mfrow=c(3,2))

cols2=gray.colors(3)[1:2]#[1:2]
matplot(myX,nucAverages1,main='K562',col=cols2,type='l',lty=1,lwd=3,ylab='average read coverage',xlab="bp from exon start")
abline(v=0)
matplot(myX,nucAverages2,main='Gm12878',col=cols2,type='l',lty=1,lwd=3,ylab='average read coverage',xlab="bp from exon start")
abline(v=0)
legend('topleft',col=rev(cols2),lwd=3,leg=c('internal','last'))

matplot(myX, averages_nucK562_HK_last,ylab='MNase-seq read coverage',  main='K562', xlab='bp from start of final exon',pch=19,type='l',col=1:2, lty=1,lwd=2) 
matplot(myX, averages_nucGm12878_HK_last,ylab='MNase-seq read coverage', main='Gm12878',xlab='bp from start of final exon',pch=19,type='l',col=1:2, lty=1,lwd=2) 

legend('topright',col=1:2,lwd=2,leg=c('non','HK'))      

matplot(myX, averages_nucK562_HK,ylab='MNase-seq read coverage',  main='K562', xlab='bp from start of exon',pch=19,type='l',col=1:2, lty=1,lwd=2) 
matplot(myX, averages_nucGm12878_HK,ylab='MNase-seq read coverage', main='Gm12878',xlab='bp from start of exon',pch=19,type='l',col=1:2, lty=1,lwd=2) 


plot.off()


# Plot Nuc.
matplot(myX, nucAverages1[,1:2],pch=19,type='b',col=myCols[8:7],main='Nuc. Occupancy (K562)')
legend('bottom',col=myCols[6:5],lwd=3,leg=c('non-last','last'))

matplot(myX, nucAverages2[,1:2],pch=19,type='b',col=myCols[8:7],main='Nuc. Occupancy (Gm12878)')
legend('bottom',col=myCols[8:7],lwd=3,leg=c('non-last','last'))


########################################################################
### Plot Pol ChIP signals as a function of intron count
########################################################################

### Add this stuff to mergedData

pol_concat = data.frame(Pol = colMeans(scores_pol) , PolS2 = colMeans(scores_polS2))
rownames(pol_concat)=bed_all$name

chipData = cbind(mergedData, pol_concat[match(mergedData$energy_ID, rownames(pol_concat)),])

boxplot(Pol ~ FeatureCount, data=chipData)

cor(chipData[,c('Pol','PolS2','Dist2End','Dist2Start','isHK','ExonsRemaining','MaxFeature', 'GeneLength','Stability','Acceptor','length','FeatureCount','isLast')])[,1:2]





