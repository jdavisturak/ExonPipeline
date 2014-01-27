rm(list=ls())
library(gplots)
library(RColorBrewer)
library(Hmisc)
#library(nls2)
source("/home/jeremy/ExonPipeline/ExonPipeline_functions.r")
source("/home/jeremy/ExonPipeline/ExonPipeline_analysis_and_plotting_Functions.r")
source("/home/home/jeremy/Code/useful_R/useful_R.r")
source("/home/jeremy/ExonPipeline/ExonPipeline_simulate_model.r")

setwd("/home/jeremy/ExonPipeline/hg19")

#load("Tilgner2_mergedData_noPseudogenes_withSplicing.RData")
#load("PlottedData_new1.RData") 
#load("PlottedData_new2.RData") # I should fix this so it loads all the stuff from new1 as well
load("PlottedData_new3.RData") 
load("PaperOnly_Fall13_mergedData_noPseudogenes_withSplicing.RData")

getQuantilMedians = function(myData,column='splicing0',N=100,col1='Dist2End'){
  myData$Group = quantileGroups(myData[,col1],N)
  data.frame(LengthMedian = tapply(myData[,col1],myData$Group,median,na.rm=T), median = tapply(myData[,column],myData$Group,median,na.rm=T))
}

#
###########################################################
##### GET SSJ data  (Dmitri D. Pervouchine)
###########################################################
###
### My data is 0-indexed, their is not??
### also, i added a '1' b/c the way i ran bam2ssj, the 1 is added like *_1
#
#ChromInfo  = cbind(c("PolyAMinus.ssj", "Chromatin.ssj","PolyAPlus.ssj"),c('splicingAMinus','Chromatin','splicingAPlus'))
#
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
##   My way:
#  dat$DenomJer = dat$count53 + (dat$count50 + dat$count03)/2
#  dat$SpliceJer = dat$count53 / dat$DenomJer
#  
#  splicingData[[i]] = dat[dat$DenomJer  > 15,]
#  splicingData_Theta5[[i]] = dat[dat$Denom5  > 15,]
#  splicingData_Theta3[[i]] = dat[dat$Denom3  > 15,]  
#}                                                                                                                

#########################################################
#### Merge with SSJ data 
##########################################################

#mergedData$coSI_ID = paste(mergedData$chr,rowMins(mergedData[,c('start','end')])+1,rowMaxes(mergedData[,c('start','end')]),1,sep='_')
#mergedData$splicingAMinus = splicingData[[1]]$SpliceJer[match(mergedData$coSI_ID, splicingData[[1]]$coSI_ID)]
#mergedData$Chromatin = splicingData[[2]]$SpliceJer[match(mergedData$coSI_ID, splicingData[[2]]$coSI_ID)]
#mergedData$splicingAPlus = splicingData[[3]]$SpliceJer[match(mergedData$coSI_ID, splicingData[[3]]$coSI_ID)]
#
#introns3 = mergedData[apply(mergedData[,ChromInfo[,2]],1,function(x)any(!is.na(x))),]
#introns3$improvement = with(introns3,splicingAMinus-Chromatin)
#save(introns3,file="PaperOnly_Fall13_mergedData_noPseudogenes_withSplicing.RData")



#########################################################
#### Only use genes that use the canonical polyA site
##########################################################
#load("/home/home/jeremy/RNAseq/Tilgner/K562_polyA_depleted/ENCODE_K562_cleavedFractions_refseq.RData")
#cyt_polyA_coverage = AllCleavedAverages$CytPlus[match(introns3$UniqueID, rownames(AllCleavedAverages))]
#introns4 = introns3[cyt_polyA_coverage > 0.90,] 

load('Human_multiple_cellTypes_splicingInfo.RData')
cyt_polyA_coverage = myCleavedFractions$K562[match(introns3$UniqueID, names(myCleavedFractions$K562))] 
introns4 = introns3[!(introns3$UniqueID %in% names(myCleavedFractions$K562[myCleavedFractions$K562 < 0.95])), ]


#####################################################################
## Model MEDIANS
#####################################################################

QuantileData  = getQuantilMedians(introns4,col1='Dist2End',column="Chromatin",1000) 
QuantileData_nuc  = getQuantilMedians(introns4,col1='Dist2End',column="splicingAMinus",1000) 
QuantileData_nucPlus = getQuantilMedians(introns4,col1='Dist2End',column="splicingAPlus",1000) 

QuantileData100  = getQuantilMedians(introns4,col1='Dist2End',column="Chromatin",100) 
QuantileData100_nuc  = getQuantilMedians(introns4,col1='Dist2End',column="splicingAMinus",100) 
QuantileData100_nucPlus  = getQuantilMedians(introns4,col1='Dist2End',column="splicingAPlus",100) 
QuantileData100_improvement  = getQuantilMedians(introns4,col1='Dist2End',column="improvement",100) 

QuantileData_HK  = getQuantilMedians(subset(introns4,isHK==1),col1='Dist2End',column="Chromatin",100)  
QuantileData_non  = getQuantilMedians(subset(introns4, isHK==0),col1='Dist2End',column="Chromatin",100) 
QuantileData_HK_nuc  = getQuantilMedians(subset(introns4,isHK==1),col1='Dist2End',column="splicingAMinus",100)  
QuantileData_non_nuc  = getQuantilMedians(subset(introns4, isHK==0),col1='Dist2End',column="splicingAMinus",100) 
                                                                                                                 
MeanGamma2 = function(LengthMedian,g,k,A,B)  sapply(LengthMedian, function(x)A * integrate(pgamma,0,x+B, shape=g, rate=k)$value/(x+B))
callMeanGamma2 = function(dat,coefs)MeanGamma2(dat,coefs[1],coefs[2], coefs[3],coefs[4])
MeanExp = function(LengthMedian,k,A)  sapply(LengthMedian, function(x)A * integrate(pgamma,0,x, shape=1, rate=k)$value/x)

model_1000_RT = nls(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileData,start=c(g=1,k=5e-5,A=.86,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,Inf,Inf),algo='port');
model_1000_noRT = nls(median ~ MeanExp(LengthMedian,k,A), data=QuantileData,start=c(k=2e-3,A=.85),lower=c(0,0),upper=c(Inf,1),algo='port');
model_nuc_1000_RT = nls(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileData_nuc,start=c(g=1,k=5e-5,A=.86,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,Inf,Inf),algo='port');
model_nuc_1000_noRT = nls(median ~ MeanExp(LengthMedian,k,A), data=QuantileData_nuc,start=c(k=2e-3,A=.85),lower=c(0,0),upper=c(Inf,1),algo='port');

model_100_HK = nls(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileData_HK,start=c(g=1,k=5e-5,A=.86,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,Inf,Inf),algo='port');
model_100_non = nls(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileData_non,start=c(g=1,k=5e-5,A=.86,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,Inf,Inf),algo='port');
model_100_HK_nuc = nls(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileData_HK_nuc,start=c(g=1,k=5e-5,A=.86,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,Inf,Inf),algo='port');
model_100_non_nuc = nls(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileData_non_nuc,start=c(g=1,k=5e-5,A=.86,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,Inf,Inf),algo='port');

save(model_1000_RT,model_1000_noRT,model_nuc_1000_RT,model_nuc_1000_noRT,model_100_HK, model_100_non , model_100_HK_nuc,  model_100_non_nuc,file='MedianModels_introns4_cyt_not_lt95.RData')

## Show all parameters of interest:
outputParams = function(x)sprintf("%.f bp per splicing event = ; %.f seconds intron half-life if elongation is 3kb/min",1/x,60*log(2)/(3000*x))

outputParams(coef(model_nuc_1000_noRT)[1])
#"597 bp per splicing event = ; 8 seconds intron half-life if elongation is 3kb/min"

outputParams(coef(model_nuc_1000_RT)[2]); sprintf("Read through is %.f bp", coef(model_nuc_1000_RT)[4])
#3443 bp per splicing event = ; 48 seconds intron half-life if elongation is 3kb/min"
#[1] "Read through is 5344 bp"

## HK vs Non
outputParams(coef(model_100_HK_nuc)[2]); sprintf("Read through is %.f bp", coef(model_100_HK_nuc)[4])
#[1] "3700 bp per splicing event = ; 51 seconds intron half-life if elongation is 3kb/min"
#[1] "Read through is 7535 bp"


outputParams(coef(model_100_non_nuc)[2]);sprintf("Read through is %.f bp", coef(model_100_non_nuc)[4])
#[1] "2961 bp per splicing event = ; 41 seconds intron half-life if elongation is 3kb/min"
#[1] "Read through is 3773 bp"
 
### NOTE: these numbers are SLIGHTLY different than those obtained in the multiple-cell-types analysis.  I'm not sure exactly why ....


### Chromatin:
### Older numbers are from using ALL genes even those that aren't poly-Adenylated in the correct spot. (i.e. using 'introns3' instead of 'introns4')
#outputParams(coef(model_1000_RT)[2]);sprintf("Read through is %.f bp", coef(model_1000_RT)[4])
#[1] "9934 bp per splicing event = ; 138 seconds intron half-life if elongation is 3kb/min"
#[1] "Read through is 16311 bp"
## "7984 bp per splicing event = ; 111 seconds intron half-life if elongation is 3kb/min"P
##[1] "Read through is 11878 bp"
#
#outputParams(coef(model_1000_noRT)[1])
#"387 bp per splicing event = ; 5 seconds intron half-life if elongation is 3kb/min"
##[1] "631 bp per splicing event = ; 9 seconds intron half-life if elongation is 3kb/min"
#



################################################################################################
### PLOT: Chromatin (all genes), 1000
################################################################################################

blueCol = rgb(0,174,239,max=255)
plot.dev("Chromatin_splicing1000.pdf",'pdf', height=3,width=3)
plot(log10(QuantileData[,1]), (QuantileData[,2]), main="Chromatin-associated RNA", pch=19, xlab= 'kb from intron to poly(A) site' , ylim=c(0,1),ylab=expression(paste('median ',Sigma)),cex=.5,cex.main=0.8,xaxt='n',cex.axis=.8,yaxt='n')
formatAxes_kb()
lines(log10(xx<- 10^seq((0),log10(2000000),length.out=1000)),callMeanGamma2(xx,coef(model_1000_RT)),type='l',col=blueCol,lwd=2)  
lines(log10(xx<- 10^seq((0),log10(2000000),length.out=1000)),MeanExp(xx,coef(model_1000_noRT)[1],coef(model_1000_noRT)[2]),type='l',col=blueCol,lwd=2,lty=2)     
plot.off()

plot.dev("NUCminus_splicing1000_1B_2A.pdf",'pdf', height=3,width=3)
plot(log10(QuantileData_nuc[,1]), (QuantileData_nuc[,2]), main="Nuclear poly(A) Minus RNA", pch=19, xlab= 'kb from intron to poly(A) site' , ylim=c(0,1),ylab=expression(paste('median ',Sigma)),cex=.5,xaxt='n',yaxt='n',cex.axis=.8,cex.main=.8)
formatAxes_kb()
lines(log10(xx<- 10^seq((0),log10(2000000),length.out=1000)),callMeanGamma2(xx,coef(model_nuc_1000_RT)),type='l',col=blueCol,lwd=2)  
lines(log10(xx<- 10^seq((0),log10(2000000),length.out=1000)),MeanExp(xx,coef(model_nuc_1000_noRT)[1],coef(model_nuc_1000_noRT)[2]),type='l',col=blueCol,lwd=2,lty=2)     
plot.off()


################################################################################################
### PLOT Housekeeping/non
################################################################################################

plot.dev("NUC_splicing100_HKvNon_4A.pdf",'pdf', height=3, width=3)
par(cex=0.5)
plot(log10(QuantileData_HK_nuc$LengthMedian), xlim=log10(c(100,420000)), QuantileData_HK_nuc$median, ylab=expression(paste('median ',Sigma)), xlab='kb from intron to poly(A) site' ,pch=20,ylim=c(0,1),col='red',main='Nuc (-)',xaxt='n',yaxt='n')
points(log10(QuantileData_non_nuc$LengthMedian),  QuantileData_non_nuc$median,pch=20, col='black',yaxt='n')
formatAxes_kb()
lines(log10(xx<- 10^seq((0),log10(500000),length.out=100)),callMeanGamma2(xx,coef(model_100_HK_nuc)),type='l',col='red',lwd=1) 
lines(log10(xx<- 10^seq((0),log10(500000),length.out=100)),callMeanGamma2(xx,coef(model_100_non_nuc)),type='l',col='black',lwd=1) 
legend('topleft',col=c(2,1),lwd=1,leg=c('HK','non'))
plot.off()

plot.dev("CHROM_splicing100_HKvNon.pdf",'pdf', height=3, width=3)
par(cex=0.5)
plot(log10(QuantileData_HK$LengthMedian), xlim=log10(c(100,420000)), QuantileData_HK$median, ylab=expression(paste('median ',Sigma)), xlab='kb from intron to poly(A) site' ,pch=20,ylim=c(0,1),col='red',main='Chromatin',xaxt='n',yaxt='n')
points(log10(QuantileData_non$LengthMedian),  QuantileData_non$median,pch=20, col='black',yaxt='n')
formatAxes_kb()
lines(log10(xx<- 10^seq((0),log10(500000),length.out=100)),callMeanGamma2(xx,coef(model_100_HK)),type='l',col='red',lwd=1) 
lines(log10(xx<- 10^seq((0),log10(500000),length.out=100)),callMeanGamma2(xx,coef(model_100_non)),type='l',col='black',lwd=1) 
legend('topleft',col=c(2,1),lwd=1,leg=c('HK','non'))
plot.off()






################################################################################################################
################################################################################################################
#### Simulate model!! 
################################################################################################################
################################################################################################################

mergedData = mergedData[order(mergedData$UniqueID,mergedData$FeatureCount),]

# Assign each gene the median stability from its group:
S=with(subset(mergedData,isLast==0), tapply(Stability,list(GeneSize,intronCountGroups),median, na.rm=T))
S2 = data.frame(S=as.vector(S),GeneSize=rep(1:4,5),intronCountGroups=  rep(1:5, each=4))
mergedData$GroupStability = S2$S[match(with(mergedData,paste(GeneSize,intronCountGroups )), with(S2,paste(GeneSize,intronCountGroups )) )]
mergedData_forFit = data.frame(Stability=mergedData$GroupStability, Dist2End=mergedData$Dist2End)
fitData_to_model = coef(model_nuc_1000_RT)
#fitData_to_model2 = chrom_noRT_coef1000


# Using the fitted parameters
Elong = 1/fitData_to_model[2]
Elong_noRT = 1/coef(model_nuc_1000_noRT)[1]
readThrough = fitData_to_model[4]
readThrough_for_pauseSplice = readThrough#3840

Stability_exp=1
Stability_div = max(S,na.rm=T)*1.20
K_splice = 1
Acceptor_div = 4
#extraDist = rep(80,nrow(mergedData))
#extraDist =  mergedData$exonLength
extraDist = rep(147,nrow(mergedData))

# For simplified version: multiplication factor = 1 + 1/Acceptor_div
mult_factor = 2
Acceptor_div = 1/(mult_factor-1) 
K_adjusted=K_splice

### Use K_adjusted to balance out the Acceptor_div 
#K_adjusted = K_splice * nrow(mergedData)/(nrow(mergedData) + (mult_factor-1)*length(unique(mergedData$UniqueID)))
### Check: this should equal K_splice
#K_splice==K_adjusted*( (nrow(mergedData)-length(unique(mergedData$UniqueID))) + mult_factor*length(unique(mergedData$UniqueID)))/nrow(mergedData)
#


myEfficiencies = tapply(1:nrow(mergedData), mergedData$UniqueID, function(x)simulateGene(mergedData_forFit[x ,], extraDist=extraDist,Elong=Elong,readThrough=readThrough, K_splice=K_adjusted,Stability_exp=Stability_exp,Stability_div=Stability_div,Acceptor_div=Acceptor_div))
myEfficiencies_noSplice = tapply(1:nrow(mergedData), mergedData$UniqueID, function(x)simulateGene(mergedData_forFit[x,], extraDist=extraDist,Elong=Elong,readThrough=readThrough, K_splice=K_splice,Acceptor_div=Inf,Stability_exp=Stability_exp,Stability_div=Stability_div))

myEfficiencies_noPause = tapply(1:nrow(mergedData), mergedData$UniqueID, function(x)simulateGene(mergedData[x,], extraDist=extraDist,Elong=Elong,readThrough=readThrough, K_splice=K_adjusted,Stability_div=Inf,Stability_exp=Stability_exp,Acceptor_div=Acceptor_div))

myEfficiencies_plain = tapply(1:nrow(mergedData), mergedData$UniqueID, function(x)simulateGene(mergedData[x,], extraDist=extraDist,Elong=Elong,readThrough=readThrough, K_splice=K_splice,Stability_div=Inf, Acceptor_div=Inf,Stability_exp=Stability_exp))
myEfficiencies_noRT = tapply(1:nrow(mergedData), mergedData$UniqueID, function(x)simulateGene(mergedData[x,], extraDist=extraDist,Elong=Elong,readThrough=0, K_splice=K_splice,Stability_div=Inf, Acceptor_div=Inf,Stability_exp=Stability_exp))
#myEfficiencies_noRT = tapply(1:nrow(mergedData), mergedData$UniqueID, function(x)simulateGene(mergedData[x,], extraDist=extraDist,Elong=Elong_noRT,readThrough=0, K_splice=K_splice,Stability_div=Inf, Acceptor_div=Inf,Stability_exp=Stability_exp))
                     s
# Get other info
genes.isHK = names(myEfficiencies) %in% subset(mergedData,isHK==1,UniqueID)[,1]
genes.Size = mergedData$GeneSize[match(names(myEfficiencies), mergedData$UniqueID)]
genes.IntronCount = mergedData$MaxFeature[match(names(myEfficiencies), mergedData$UniqueID)]
efficiencyData = data.frame(UniqueID=names( myEfficiencies), isHK=genes.isHK,GeneSize=genes.Size,MaxFeature=genes.IntronCount,intronCountGroups =  splitByVector(genes.IntronCount,c(2.5,4.5,7.5,19.5)),
  myEfficiencies_noRT,myEfficiencies_plain,myEfficiencies_noSplice,myEfficiencies_noPause,myEfficiencies)


### PLOT better ####
myCols5 = cols=(brewer.pal(9,"YlGnBu"))[5:9]; pch=16
library(reshape2)
eff2 = melt(efficiencyData[,c(!grepl('(RT)|(plain)',dimnames(efficiencyData)[[2]]))],id=c('UniqueID','isHK','GeneSize','MaxFeature','intronCountGroups'))
plainMedians = sapply(c(1,4),function(size)sapply(1:5,function(num)median(subset(efficiencyData,intronCountGroups==num & GeneSize==size)$myEfficiencies_plain, na.rm=T)))


plot.dev(sprintf("GenomeSimulations_numIntrons_by_geneSize_nucElong_147_noAdjust_median%.3f_Extreme_only_3C.pdf",Stability_div),'pdf',height=2,width=5)  
cex=0.6;width=0.4;
par(cex=cex,mai=c(.02,.3,.02,.02),pch=pch)

at = c(c(1:5,7:11), 13 + c(1:5,7:11), 26  + c(1:5,7:11))
x=outer(c(1,1,NA),at)+rep(.4*c(-width,width,NA),30)
boxplot(value ~ intronCountGroups +GeneSize + variable, data= subset(eff2,GeneSize %in% c(1,4)) , at=at,boxwex=width, col=myCols5,pch=20,xaxt='n')
arrows(x[1,],rep(plainMedians,3),x[2,],rep(plainMedians,3),col='maroon',angle=135,code=3,lwd=2,length=width/Inf)
plot.off()




plot.dev("GenomeSimulations_readThrough_comparisons_chromElong_147_noAdjust_2only_2D.pdf",'pdf',width=5,height=1.5)
#plot.dev("GenomeSimulations_readThrough_comparisons_chromElong_147_noAdjust_diffElong_2only_2D.pdf",'pdf',width=5,height=1.5)
par(mfcol=c(1,1),cex=0.5,mai=c(0.3,0.2,0.01,0.01))
notch=F; wid=.21;bw=0.1085; offset = 1.4; xlim=c(0.5,9.2); widthFraction = 0.4
cols =c(brewer.pal(9,"YlGnBu")[5:9],brewer.pal(9,"OrRd")[5:9])
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
    b=boxplot(Form1, data=x,main=main,border=COL[myCats],ylim=c(0,1),lwd=1.5,staplewex=0,range=0.00001,outline=F,boxwex=bw,at=at-wid*widthFraction,xlim=xlim,names=NA,yaxt=ifelse(sizes==1,'s','n'),notch=notch,add=sizes-1,xaxt='n')
    #if (sizes==1) abline(h=0,col='gray')
    b=boxplot(Form3, data=x,main=main,col=COL[myCats],border='black',boxwex=bw,at=at+wid*widthFraction,names=NA,notch=notch,add=T,xaxt='n',lwd=lwd,yaxt='n')
  }
  
  #legend('topleft',fill=COL,leg=c('1-2 introns','3-4 introns','5-7 introns','8-19 introns','20+ introns'))
  #legend('bottomleft',fill=c(NA,COL[1]), border=c(COL[1],gray.colors(2)[1],'black'),leg=c('No Read-Thru',sprintf('%.1f kb Read-Thru',readThrough_for_pauseSplice/1000)))
 axis(1,at=2*(3*wid+(0:3)*offset)-(0:3)*wid*2, lab = c('Short','Med short','Med long','Long')) 
}
plot.off()

#plot.dev(sprintf("GenomeSimulations_numIntrons_by_geneSize_nucElong_147_noAdjust_median%.3f_Extreme_only.pdf",Stability_div),'pdf',height=3,width=8)  
#
#notch=F; wid=.16;bw=0.11; offset1 = 1.5; offset2 = .8; xlim=c(-.8,11); cex=0.6; pch=19
#COL = c(brewer.pal(9,"YlGnBu")[5:9],brewer.pal(9,"OrRd")[5:9])
#par(cex=cex,mai=c(.5,.4,.1,.1),pch=pch)
#      
#      
#Xshift=0; add=FALSE; xpos=NULL
#for (HK in list(0:1)){  
#  for (exper in 2:4){    
#    for(sizes in c(1,4)){
#      txt = paste(ifelse(length(HK)==1,ifelse(HK,'HK:','Regulated'),'All:'))#, c('Short','Med short','Med long','Long')[sizes])
#      x<-subset(efficiencyData,GeneSize==sizes &isHK%in%HK)
#      myCats = sort(unique(x$intronCountGroups))
#      
#      main=''#ifelse(sizes==1,paste(txt,switch(exper,'Basal parameters','Add Pausing','Add Splicing variability','Add both')),'')
#      # Decide which variable to test
#      Form0 = myEfficiencies_plain ~ intronCountGroups
#      Form = switch(exper,Form0,myEfficiencies_noSplice ~ intronCountGroups, myEfficiencies_noPause ~ intronCountGroups, myEfficiencies ~ intronCountGroups) 
#    
#      at = 2*(myCats)*wid + Xshift
#      xpos = c(xpos,mean(at))
#      
#      b=boxplot(Form0, data=x,main=main,border=COL[myCats],staplewex=0,range=0.00001,outline=F,ylim=c(0,1),lwd=2,boxwex=bw*.5,at=at-wid*.7,xlim=xlim,names=NA,xaxt='n',notch=notch,add=add,cex=cex,pars=list(pch=pch,cex=cex))
#      b=boxplot(Form, data=x,main=main,col=COL[myCats],boxwex=bw,at=at,names=NA,notch=notch,add=T,xaxt='n')
#      add=TRUE 
#      Xshift = offset1 + Xshift      
#    }
#      Xshift = offset2 + Xshift          
#  }
# abline(h=0,col='gray')
#  legend('left',fill=COL,leg=c('1-2 introns','3-4 introns','5-7 introns','8-19 introns','20+ introns'))#sprintf(c('1-2 introns (n=%d)','3-4 introns (n=%d)','5-7 introns (n=%d)','8-19 introns (n=%d)','20+ introns (n=%d)')[myCats], table(subset(efficiencyData,GeneSize==sizes&isHK%in%HK)$intronCountGroups)))
#  axis(1,at=xpos, lab = rep(c('Short','Long'),3)) 
#
#}
#plot.off()
 
 




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
## Compute overlaps
###############################################################################


################## Pol ChIP: 3' ends         

### Overlap with BW files:
#scores_lastPol.ends = scoreBW_bed(bw_pol.ends, bed_ends, RangeScores_pol.ends)
#scores_lastPolS2.ends = scoreBW_bed(bw_polS2.ends, bed_ends, RangeScores_polS2.ends)
#colnames(scores_lastPol.ends) = colnames(scores_lastPolS2.ends) = bed_ends$UniqueID
#save(scores_lastPol.ends, scores_lastPolS2.ends,file='ENCODE_PolCHIP_ends1000.RData')
# 
############### Nucleosomes: intron/exon junctions   
# K562
#scores_last_nucK562 = scoreBW_bed(bw_NucK562, bed_last, RangeScores_NucK562)
#scores_nonLast_nucK562 = scoreBW_bed(bw_NucK562, bed_nonLast, RangeScores_NucK562)
#
# Gm12878
#scores_last_nucGm12878 = scoreBW_bed(bw_NucGm12878, bed_last, RangeScores_NucGm12878)
#scores_nonLast_nucGm12878 = scoreBW_bed(bw_NucGm12878, bed_nonLast, RangeScores_NucGm12878)
#
#colnames(scores_last_nucGm12878) = colnames(scores_last_nucK562)= bed_last$UniqueID
#colnames(scores_nonLast_nucGm12878) = colnames(scores_nonLast_nucK562)= bed_nonLast$UniqueID
#save(scores_last_nucGm12878,scores_nonLast_nucGm12878,scores_last_nucK562,scores_nonLast_nucK562,nucAverages1,nucAverages2,
#  file=sprintf('ENCODE_Nuc_occupancy_%d.RData', dist))



###############################################################################
## Redo with different files
distU=1000
distD=5000
bed_endPlus   = import.bed(sprintf("/home/home/jeremy/ExonPipeline/hg19/Human_polyA_-%dTo%d+.bed",distU,distD),asRangedData=FALSE)
bed_endMinus  = import.bed(sprintf("/home/home/jeremy/ExonPipeline/hg19/Human_polyA_-%dTo%d-.bed",distU,distD),asRangedData=FALSE)
allBeds_selection.ends = make_beds_selection(bed_ends<-c(bed_endPlus, bed_endMinus))
bed_ends$UniqueID = sub(">([0-9]*)_.*","\\1",bed_ends$name)

## Pol ChIP
bw_polS2.ends2    = import.bw("/home/RNAseq/ENCODE/wgEncodeSydhTfbsK562Pol2s2IggrabSig.bigWig",selection=allBeds_selection.ends,asRangedData=FALSE)
RangeScores_polS2.ends = cbind(as.matrix(ranges(bw_polS2.ends2)),score(bw_polS2.ends2));
scores_lastPolS2.ends2 = scoreBW_bed(bw_polS2.ends2, bed_ends, RangeScores_polS2.ends)
colnames(scores_lastPolS2.ends2) = bed_ends$UniqueID
save(scores_lastPolS2.ends2,file=sprintf("ENCODE_PolCHIP_ends%d-%d.RData",distU,distD)) 











load("ENCODE_Nuc_occupancy_250.RData")
load('ENCODE_PolCHIP_ends1000.RData')
load('ENCODE_PolCHIP_ends1000-5000.RData')
nucAverages1 = cbind(last_K562 = rowMeans(scores_last_nucK562),nonLast_K562 = rowMeans(scores_nonLast_nucK562))
nucAverages2 = cbind(last_Gm12878 = rowMeans(scores_last_nucGm12878),nonLast_Gm12878 = rowMeans(scores_nonLast_nucGm12878))

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




### Pol S2: get stats
averages_lastPolS2_size.ends = sapply(1:4, function(x)rowMeans(scores_lastPolS2.ends[,colnames(scores_lastPolS2.ends) %in% subset(mergedData,GeneSize==x)$UniqueID]))
averages_lastPolS2_size_numIntrons.ends = sapply(1:4,function(size)sapply(1:5, function(num)rowMeans(scores_lastPolS2.ends[,colnames(scores_lastPolS2.ends) %in% subset(mergedData,GeneSize==size & intronCountGroups==num)$UniqueID])),simplify=F) 
averages_lastPolS2_HK.ends = sapply(0:1, function(x)rowMeans(scores_lastPolS2.ends[,colnames(scores_lastPolS2.ends) %in% subset(mergedData,isHK==x)$UniqueID]))

averages_lastPolS2_HK_size.ends = sapply(1:4, function(x)rowMeans(scores_lastPolS2.ends[,colnames(scores_lastPolS2.ends) %in% subset(mergedData,GeneSize==x &isHK==1)$UniqueID]))
averages_lastPolS2_non_size.ends = sapply(1:4, function(x)rowMeans(scores_lastPolS2.ends[,colnames(scores_lastPolS2.ends) %in% subset(mergedData,GeneSize==x &isHK==0)$UniqueID]))

# Get Average signal over the 2000 Kb window
totalAverages_lastPolS2.ends = lapply(averages_lastPolS2_size_numIntrons.ends,function(x){x=colMeans(x[1:dist,]);x[is.nan(x)]<-NA;x})


############################################################
## Plotting ChIP-seq data
############################################################

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

plot.dev("ENCODE_polS2_ChIP_ends.pdf", ,'pdf',height=8,width=12)
plotSize_numIntrons(averages_lastPolS2_size.ends,averages_lastPolS2_size_numIntrons.ends,'PolS2 ChIP: \nCoverage over the poly(A) site',xlab='bp from poly(A) site', ylim=c(0,70),dist=dist2)
plot.off()

plot.dev("ENCODE_polS2_ChIP_ends_HK_comparison.pdf",'pdf',height=3,width=3)
par(mai=c(.8,.5,0.1,0.3))
dist=1000;myX = (1:(2*dist)) - dist
matplot(myX, averages_lastPolS2_HK.ends,ylab='PolS2 read coverage', ylim=c(0,80),xlab='bp from poly(A) site',pch=19,type='l',col=1:2, lty=1,lwd=2) 
lines(myX, averages_lastPolS2_HK.ends[,1]*averages_lastPolS2_HK.ends[1000,2]/ averages_lastPolS2_HK.ends[1000,1],lty=2,lwd=2)
legend('topleft',col=c(1,1:2),lwd=2,lty=c(1,2,1),leg=c('non-HK','normalized non-HK','HK'))      
plot.off()

plot.dev("ENCODE_polS2_ChIP_ends_HK_size_comparison.pdf",'pdf',height=3,width=6)
par(mfrow=c(1,2),mai=c(.8,.5,0.1,0.3))
dist=1000;myX = (1:(2*dist)) - dist
matplot(myX, (apply(averages_lastPolS2_HK_size.ends,2,function(x)x/x[1])),ylab='PolS2 read coverage',xlab='bp from poly(A) site',pch=19,type='l',col=myCols4, lty=1,lwd=2) 
matplot(myX, (apply(averages_lastPolS2_non_size.ends,2,function(x)x/x[1])),ylab='PolS2 read coverage',xlab='bp from poly(A) site',pch=19,type='l',col=myCols4, lty=1,lwd=2) 
#lines(myX, averages_lastPolS2_HK.ends[,1]*averages_lastPolS2_HK.ends[1000,2]/ averages_lastPolS2_HK.ends[1000,1],lty=2,lwd=2)
legend('topleft',col=myCols4,lwd=2,leg=1:4)      
plot.off()

plot.dev("ENCODE_polS2_ChIP_ends_HK_size_comparison2.pdf",'pdf',height=3,width=3)
par(mai=c(.8,.5,0.1,0.3))
dist=1000;myX = (1:(2*dist)) - dist
cols = c('red','red','black','black'); lty=c(1,2,1,2)
matplot(myX, (apply(cbind(averages_lastPolS2_HK_size.ends[,c(1,4)],averages_lastPolS2_non_size.ends[,c(1,4)]),2,function(x)x/x[1])),ylab='PolS2 read coverage',xlab='bp from poly(A) site',pch=19,type='l',col=cols, lty=lty,lwd=2) 
#lines(myX, averages_lastPolS2_HK.ends[,1]*averages_lastPolS2_HK.ends[1000,2]/ averages_lastPolS2_HK.ends[1000,1],lty=2,lwd=2)
legend('topleft',col=cols,lwd=2,lty=lty,leg=paste(rep(c('HK','non'),each=2),rep(c('short','long'),2)))      
plot.off()


plot.dev("ENCODE_polS2_ChIP_ends_Normalized_0.pdf",'pdf',height=2,width=10)
par(mfrow=c(1,5), mai=c(.6,.5,0.2,0.1),mgp=c(2,1,0))
dist=1000;myX = (1:(2*dist)) - dist
sapply(1:4,function(size){
  matplot(myX, (apply(averages_lastPolS2_size_numIntrons.ends[[size]],2,function(x)x/x[1])), ylim=c(0.75,2),xlab='bp from poly(A) site',ylab='PolS2 normalized read coverage',main=c('Short','Med. Short','Med. Long','Long')[size],pch=19,type='l',lwd=3,col=myCols5, lty=1,yaxt='n')
  abline(h=1,col='grey',lwd=2)
  axis(2,at=seq(0.8,2,.2),lab=c(NA,1,NA,NA,NA,NA,2.0));grid()
  })
barplot(matrix(unlist(totalAverages_lastPolS2.ends),nc=4),bes=T,col=myCols5,names=c('Short','Med. Short','Med. Long','Long'),ylab='average read coverage in upstream of poly(A) site',space=c(0.2,1),las=2)
plot.off()


plot.dev("ENCODE_nuc_K562.pdf", height=8,width=12)
plotSize_numIntrons(averages_nucK562_size,averages_nucK562_size_numIntrons,'Nucleosomes in K562: \nCoverage over final intron-exon junction',ylim=c(0.5,1.5))
plot.off()

plot.dev("ENCODE_nuc_Gm12878.pdf", height=8,width=12)
plotSize_numIntrons(averages_nucGm12878_size,averages_nucGm12878_size_numIntrons,'Nucleosomes in Gm12878: \nCoverage over final intron-exon junction', ylim=c(0.5,1.5))
plot.off()


plot.dev('ENCODE_nucleosomes_lastNonlast.pdf','pdf',height=4,width=8)
dist=250;myX = (1:(2*dist)) - dist
par(mfrow=c(1,2))

cols2=gray.colors(3)[1:2]#[1:2]
matplot(myX,nucAverages1,main='K562',col=cols2,type='l',lty=1,lwd=3,ylab='average read coverage',xlab="bp from exon start")
abline(v=0)
matplot(myX,nucAverages2,main='Gm12878',col=cols2,type='l',lty=1,lwd=3,ylab='average read coverage',xlab="bp from exon start")
abline(v=0)
legend('topleft',col=rev(cols2),lwd=3,leg=c('internal','last'))

plot.off()


























































#####################################################################################################################
## Not in paper: examining in more detail the separation of the splicing data based on whether or not it's poly-adenylated:



 

#########################################################
## Get the actual polyA ends from ENCODE's processed data
##########################################################
require(rtracklayer)
k562_polyA = import.bed('/home/RNAseq/Tilgner/Nucleus/K562_NucPap_genes2.bed')
k562_polyAm = import.bed('/home/RNAseq/Tilgner/Nucleus/K562_NucPam_genes2.bed')

# Overlap with refseq genes
ref.GRange = GRanges(refseq$chrom, IRanges(refseq$txStart,refseq$txEnd), refseq$strand, UniqueID=refseq[,'UniqueID'])
ref.GRange2 = ref.GRange[ref.GRange$UniqueID %in% mergedData$UniqueID,]
ref.GRange3 = ref.GRange[ref.GRange$UniqueID %in% introns3$UniqueID,]


getExtraTxn = function(bed,ref){
  overlap_plus = as.matrix(findOverlaps(refPlus <- ref[strand(ref)=='+',],bedPlus <- bed[bed$strand=='+',]))
  overlap_minus = as.matrix(findOverlaps(refMinus <- ref[strand(ref)=='-',],bedMinus <- bed[bed$strand=='-',]))
  overlap2 = rbind(cbind(which(strand(ref)=='+')[overlap_plus[,1]], which(bed$strand=='+')[overlap_plus[,2]]),
    cbind(which(strand(ref)=='-')[overlap_minus[,1]], which(bed$strand=='-')[overlap_minus[,2]]))
  dups = subset(as.data.frame(overlap2), V1 %in% V1[duplicated(V1)])
  
  myRanges_plus = refPlus[overlap_plus[,1],]
  myRanges_minus = refMinus[overlap_minus[,1],]
  
  myRanges_plus$END2 = end(bedPlus)[overlap_plus[,2]]
  p = tapply((myRanges_plus$END2), as.factor(myRanges_plus$UniqueID),max)
  myRanges_minus$START2 = start(bedMinus)[overlap_minus[,2]]
  m = tapply((myRanges_minus$START2), as.factor(myRanges_minus$UniqueID),min)
  
  extraTxn_plus = data.frame(UniqueID = names(p), extra = p - refseq$txEnd[match(names(p),refseq$UniqueID)])
  extraTxn_minus = data.frame(UniqueID = names(m), extra = refseq$txStart[match(names(m),refseq$UniqueID)] - m)
  extraTxn = rbind(extraTxn_plus,extraTxn_minus)
  extraTxn
}

extra_plus1 = getExtraTxn(k562_polyA,ref.GRange3)
extra_minus1 = getExtraTxn(k562_polyAm,ref.GRange3)

summary(extra_plus1[,2])
#      Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#-2145000.0    -1117.0      -39.0      539.4     1242.0   341300.0 
summary(extra_minus1[,2])
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-2032000    -3472      -37    -1910     2547   446600 

minus_v_plus_txn = merge(extra_minus1, extra_plus1,by=1)
polyA_diff = with(minus_v_plus_txn,(extra.y-extra.x))

extra_plus2 = getExtraTxn(k562_polyA[k562_polyA$score > 0.1,], ref.GRange3)
extra_minus2 = getExtraTxn(k562_polyAm[k562_polyAm$score > 0.1,], ref.GRange3)
minus_v_plus_txn2 = merge(extra_minus2, extra_plus2,by=1)
polyA_diff2 = with(minus_v_plus_txn2,(extra.y-extra.x))

# For ref.GRange3:
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-232600.0    -859.8    1415.0   12930.0    9611.0  377300.0 

# For ref.GRange2:
summary(polyA_diff2)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-266300   -2069     791   12430   11350  377300
summary(polyA_diff2[minus_v_plus_txn$UniqueID %in% subset(refseq,strand=='+')$UniqueID])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#-223800   -2836     688   11850   11560  242700    2826 
summary(polyA_diff2[minus_v_plus_txn$UniqueID %in% subset(refseq,strand=='-')$UniqueID])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#-266300   -1618    1104   13020   10860  377300    2681 
 

## First use predict
#predictions_nuc_RT = predict(model_nuc_1000_RT, newdata = data.frame(LengthMedian=introns3$'Dist2End'))
## Now look at the relationship between prediction accuracy and polyA coverage
#diff_nuc_RT = introns3$splicingAMinus - predictions_nuc_RT
#predictions_chrom_RT = predict(model_1000_RT, newdata = data.frame(LengthMedian=introns3$'Dist2End'))
#diff_chrom_RT = introns3$Chromatin - predictions_chrom_RT

nucAminus_polyA_coverage = AllCleavedAverages$NucMinus [match(introns3$UniqueID, rownames(AllCleavedAverages))]
nucAplus_polyA_coverage = AllCleavedAverages$NucPlus [match(introns3$UniqueID, rownames(AllCleavedAverages))]
chrom_polyA_coverage = AllCleavedAverages$Chromatin[match(introns3$UniqueID, rownames(AllCleavedAverages))]
cyt_polyA_coverage = AllCleavedAverages$CytPlus[match(introns3$UniqueID, rownames(AllCleavedAverages))]


### The real thing to do would be to plot the difference Quantiles for two populations
QuantileData_nuc_100_nonA = getQuantilMedians(introns3[nucAminus_polyA_coverage < 0.50,],col1='Dist2End',column="splicingAMinus",100) 
QuantileData_nuc_100_A = getQuantilMedians(introns3[nucAminus_polyA_coverage > 0.75,],col1='Dist2End',column="splicingAMinus",100) 
#QuantileData_nuc_1000_nonA = getQuantilMedians(introns3[nucAminus_polyA_coverage < 0.464,],col1='Dist2End',column="splicingAMinus",1000) 
#QuantileData_nuc_1000_A = getQuantilMedians(introns3[nucAminus_polyA_coverage > 0.64,],col1='Dist2End',column="splicingAMinus",1000) 
QuantileData_nucP_100_nonA = getQuantilMedians(introns3[nucAplus_polyA_coverage < 0.50,],col1='Dist2End',column="splicingAPlus",100) 
QuantileData_nucP_100_A = getQuantilMedians(introns3[nucAplus_polyA_coverage > 0.75,],col1='Dist2End',column="splicingAPlus",100) 

QuantileData_chrom_100_nonA = getQuantilMedians(introns3[chrom_polyA_coverage < 0.50,],col1='Dist2End',column="Chromatin",100) 
QuantileData_chrom_100_A = getQuantilMedians(introns3[chrom_polyA_coverage > 0.75,],col1='Dist2End',column="Chromatin",100) 

## This is the important one:
QuantileData_nuc_100_non_cytA = getQuantilMedians(introns3[cyt_polyA_coverage < 0.70,],col1='Dist2End',column="splicingAMinus",100) 
QuantileData_nuc_100_cytA = getQuantilMedians(introns3[cyt_polyA_coverage > 0.9,],col1='Dist2End',column="splicingAMinus",100) 

## Now filter in Both!!!
QuantileData_nuc_100_non_cytA2 = getQuantilMedians(introns3[cyt_polyA_coverage > 0.90 & nucAminus_polyA_coverage < 0.50,],col1='Dist2End',column="splicingAMinus",100) 
QuantileData_nuc_100_cytA2 = getQuantilMedians(introns3[cyt_polyA_coverage > 0.9 &  nucAminus_polyA_coverage > 0.75,],col1='Dist2End',column="splicingAMinus",100) 
QuantileData_chrom_100_non_cytA2 = getQuantilMedians(introns3[cyt_polyA_coverage > 0.90 & chrom_polyA_coverage < 0.50,],col1='Dist2End',column="Chromatin",100) 
QuantileData_chrom_100_cytA2 = getQuantilMedians(introns3[cyt_polyA_coverage > 0.9 &  chrom_polyA_coverage > 0.75,],col1='Dist2End',column="Chromatin",100) 



model_nuc_100_nonA = nls(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileData_nuc_100_nonA,start=c(g=1,k=5e-5,A=.86,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,1,Inf),algo='port');
model_nuc_100_A = nls(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileData_nuc_100_A,start=c(g=1,k=5e-5,A=.86,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,1,Inf),algo='port');
model_nuc_1000_nonA = nls(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileData_nuc_1000_nonA,start=c(g=1,k=5e-5,A=.86,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,1,Inf),algo='port');
model_nuc_1000_A = nls(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileData_nuc_1000_A,start=c(g=1,k=5e-5,A=.86,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,1,Inf),algo='port');
model_chrom_100_nonA = nls(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileData_chrom_100_nonA,start=c(g=1,k=5e-5,A=.86,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,1,Inf),algo='port');
model_chrom_100_A = nls(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileData_chrom_100_A,start=c(g=1,k=5e-5,A=.86,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,1,Inf),algo='port');
model_nucP_100_nonA = nls(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileData_nucP_100_nonA,start=c(g=1,k=5e-5,A=.86,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,1,Inf),algo='port');
model_nucP_100_A = nls(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileData_nucP_100_A,start=c(g=1,k=5e-5,A=.86,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,1,Inf),algo='port');

model_nuc_100_noncytA = nls(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileData_nuc_100_non_cytA,start=c(g=1,k=5e-5,A=.86,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,1,Inf),algo='port');
model_nuc_100_cytA = nls(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileData_nuc_100_cytA,start=c(g=1,k=5e-5,A=.86,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,1,Inf),algo='port');
model_nuc_100_noncytA2 = nls(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileData_nuc_100_non_cytA2,start=c(g=1,k=5e-5,A=.86,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,1,Inf),algo='port');
model_nuc_100_cytA2 = nls(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileData_nuc_100_cytA2,start=c(g=1,k=5e-5,A=.86,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,1,Inf),algo='port');
model_chrom_100_noncytA2 = nls(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileData_chrom_100_non_cytA2,start=c(g=1,k=5e-5,A=.86,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,1,Inf),algo='port');
model_chrom_100_cytA2 = nls(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileData_chrom_100_cytA2,start=c(g=1,k=5e-5,A=.86,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,1,Inf),algo='port');

# Average expression levels
UpAverages = data.frame(NucMinus = rowMeans(Up[,1:4], T), Chromatin = rowMeans(Up[,5:8], T),
  NucPlus = rowMeans(Up[,9:12], T), CytPlus = rowMeans(Up[,13:16], T))
nucAminus_polyA_up = UpAverages$NucMinus [match(introns3$UniqueID, rownames(AllCleavedAverages))]
nucAplus_polyA_up = UpAverages$NucPlus [match(introns3$UniqueID, rownames(AllCleavedAverages))]
chrom_polyA_up = UpAverages$Chromatin[match(introns3$UniqueID, rownames(AllCleavedAverages))]
                  
plot.dev("PolyA_splicing2.pdf",'pdf', height=12,width=4) 
par(mfrow=c(3,1),cex=0.7)
plot(log10(QuantileData_nuc_100_non_cytA[,1]), (QuantileData_nuc_100_non_cytA[,2]),xlim=log10(c(100,420000)), main="Nuc - RNA", pch=19, xlab= 'Distance to poly(A) site' , ylim=c(0,1),ylab=expression(paste('median ',Omega)),cex=.5,cex.main=0.8,xaxt='n',cex.axis=.8,yaxt='n')
points(log10(QuantileData_nuc_100_cytA[,1]), (QuantileData_nuc_100_cytA[,2]), cex=.5, pch=19,col='blue')
formatAxes()
lines(log10(xx<- 10^seq((0),log10(2000000),length.out=1000)),callMeanGamma2(xx,coef(model_nuc_100_noncytA)),type='l',col='black',lwd=2)  
lines(log10(xx<- 10^seq((0),log10(2000000),length.out=1000)),callMeanGamma2(xx,coef(model_nuc_100_cytA)),type='l',col='blue',lwd=2)  
legend('bottomright',col=c(1,'blue'),pch=19,leg=c('Least cleaved  (in Cytoplasm)','most cleaved  (in Cytoplasm)'))

par(cex=0.7)
plot(log10(QuantileData_nuc_100_non_cytA2[,1]), (QuantileData_nuc_100_non_cytA2[,2]),xlim=log10(c(100,420000)), main="Nuc - RNA\nGenes where Cyt cleavage >.9", pch=19, xlab= 'Distance to poly(A) site' , ylim=c(0,1),ylab=expression(paste('median ',Omega)),cex=.5,cex.main=0.8,xaxt='n',cex.axis=.8,yaxt='n',col='green')
points(log10(QuantileData_nuc_100_cytA2[,1]), (QuantileData_nuc_100_cytA2[,2]), cex=.5, pch=19,col='purple')
formatAxes()
lines(log10(xx<- 10^seq((0),log10(2000000),length.out=1000)),callMeanGamma2(xx,coef(model_nuc_100_noncytA2)),type='l',col='green',lwd=2)  
lines(log10(xx<- 10^seq((0),log10(2000000),length.out=1000)),callMeanGamma2(xx,coef(model_nuc_100_cytA2)),type='l',col='purple',lwd=2)  
legend('bottomright',col=c('green','purple'),pch=19,leg=c('Least cleaved  (in Nuc-)','most cleaved  (in Nuc-)'))

par(cex=0.7)
plot(log10(QuantileData_chrom_100_non_cytA2[,1]), (QuantileData_chrom_100_non_cytA2[,2]),xlim=log10(c(100,420000)), main="chrom  RNA\nGenes where Cyt cleavage >.9", pch=19, xlab= 'Distance to poly(A) site' , ylim=c(0,1),ylab=expression(paste('median ',Omega)),cex=.5,cex.main=0.8,xaxt='n',cex.axis=.8,yaxt='n',col='green')
points(log10(QuantileData_chrom_100_cytA2[,1]), (QuantileData_chrom_100_cytA2[,2]), cex=.5, pch=19,col='purple')
formatAxes()
lines(log10(xx<- 10^seq((0),log10(2000000),length.out=1000)),callMeanGamma2(xx,coef(model_chrom_100_noncytA2)),type='l',col='green',lwd=2)  
lines(log10(xx<- 10^seq((0),log10(2000000),length.out=1000)),callMeanGamma2(xx,coef(model_chrom_100_cytA2)),type='l',col='purple',lwd=2)  
legend('bottomright',col=c('green','purple'),pch=19,leg=c('Least cleaved  (in chrom)','most cleaved  (in chrom)'))


 plot.off()


#plot.dev("PolyA_splicing.pdf",'pdf', height=9,width=9) 
#par(mfrow=c(3,3),cex=0.7)
#plot(log10(QuantileData_chrom_100_nonA[,1]), (QuantileData_chrom_100_nonA[,2]),xlim=log10(c(100,420000)), main="Chrom-assoc RNA", pch=19, xlab= 'Distance to poly(A) site' , ylim=c(0,1),ylab=expression(paste('median ',Omega)),cex=.5,cex.main=0.8,xaxt='n',cex.axis=.8,yaxt='n')
#points(log10(QuantileData_chrom_100_A[,1]), (QuantileData_chrom_100_A[,2]), cex=.5, pch=19,col='blue')
#formatAxes()
#lines(log10(xx<- 10^seq((0),log10(2000000),length.out=1000)),callMeanGamma2(xx,coef(model_chrom_100_nonA)),type='l',col='black',lwd=2)  
#lines(log10(xx<- 10^seq((0),log10(2000000),length.out=1000)),callMeanGamma2(xx,coef(model_chrom_100_A)),type='l',col='blue',lwd=2)  
#boxplot(chrom_polyA_coverage ~ introns3$isHK, notch=T, ylab='Fraction poly-Adenylated',main='Chromatin-associated',names=c('non-HK','HK'),las=2)
#boxplot(log(chrom_polyA_up) ~ introns3$isHK, notch=T, ylab="3' Gene expression",main='Chromatin-associated',names=c('non-HK','HK'),las=2)
#
#
#plot(log10(QuantileData_nuc_100_nonA[,1]), (QuantileData_nuc_100_nonA[,2]),xlim=log10(c(100,420000)), main="Nuc- RNA", pch=19, xlab= 'Distance to poly(A) site' , ylim=c(0,1),ylab=expression(paste('median ',Omega)),cex=.5,cex.main=0.8,xaxt='n',cex.axis=.8,yaxt='n')
#points(log10(QuantileData_nuc_100_A[,1]), (QuantileData_nuc_100_A[,2]), cex=.5, pch=19,col='blue')
#formatAxes()
#lines(log10(xx<- 10^seq((0),log10(2000000),length.out=1000)),callMeanGamma2(xx,coef(model_nuc_100_nonA)),type='l',col='black',lwd=2)  
#lines(log10(xx<- 10^seq((0),log10(2000000),length.out=1000)),callMeanGamma2(xx,coef(model_nuc_100_A)),type='l',col='blue',lwd=2)  
#legend('bottomright',col=c(1,'blue'),pch=19,leg=c('Least poly-Adenylated','most poly-adenylated'))
#boxplot(nucAminus_polyA_coverage ~ introns3$isHK, notch=T, ylab='Fraction poly-Adenylated',main='Nuc-',names=c('non-HK','HK'),las=2)
#boxplot(log(nucAminus_polyA_up) ~ introns3$isHK, notch=T, ylab="3' Gene expression",main='Nuc-',names=c('non-HK','HK'),las=2)
#
#plot(log10(QuantileData_nucP_100_nonA[,1]), (QuantileData_nucP_100_nonA[,2]),xlim=log10(c(100,420000)), main="Nuc+ RNA", pch=19, xlab= 'Distance to poly(A) site' , ylim=c(0,1),ylab=expression(paste('median ',Omega)),cex=.5,cex.main=0.8,xaxt='n',cex.axis=.8,yaxt='n')
#points(log10(QuantileData_nucP_100_A[,1]), (QuantileData_nucP_100_A[,2]), cex=.5, pch=19,col='blue')
#formatAxes()
#lines(log10(xx<- 10^seq((0),log10(2000000),length.out=1000)),callMeanGamma2(xx,coef(model_nucP_100_nonA)),type='l',col='black',lwd=2)  
#lines(log10(xx<- 10^seq((0),log10(2000000),length.out=1000)),callMeanGamma2(xx,coef(model_nucP_100_A)),type='l',col='blue',lwd=2)  
#legend('bottomleft',col=c(1,'blue'),pch=19,leg=c('Least poly-Adenylated','most poly-adenylated'))
#boxplot(nucAplus_polyA_coverage ~ introns3$isHK, notch=T, ylab='Fraction poly-Adenylated',main='Nuc+',names=c('non-HK','HK'),las=2)
#boxplot(log(nucAplus_polyA_up) ~ introns3$isHK, notch=T, ylab="3' Gene expression",main='Nuc+',names=c('non-HK','HK'),las=2)
#
#plot.off()




##################################################################################################
### Look at splicing of last introns: subset for everything that I can think of!!
##################################################################################################
# First, use a good binning
introns3$lastLengthsGroup = splitByVector(introns3$Dist2End,quantile(subset(introns3,isLast==1)$Dist2End,na.rm=T,(1:(N<-25))/N))
getLastMedians = function(myData,column='splicing0',col1='Dist2End'){
  data.frame(LengthMedian = tapply(myData[,col1],myData$lastLengthsGroup,median,na.rm=T), median = tapply(myData[,column],myData$lastLengthsGroup,median,na.rm=T))
}


introns4 = introns3[cyt_polyA_coverage > 0.90,] 
column = 'splicingAMinus'
plot.dev(sprintf("SplitSize_lastNonLast_%s.pdf",column),'pdf', height=6,width=6)
par(mfrow=c(2,2))



sapply(1:4,function(x){
  data1 = subset(introns3,GeneSize==x)
  QuantileLast_last = getLastMedians(subset(data1,isLast==1),column)
  QuantileLast_nonLast = getQuantilMedians(subset(data1,isLast==0),col1='Dist2End',column=column,25*x) 
  #QuantileLast_nonLast = getLastMedians(subset(data1,isLast==0),column)
  
  #plot.dev("NUC_splicingLast50_lastVnon.pdf",'pdf', height=3, width=3)
  par(cex=0.5)
  plot(log10(QuantileLast_last$LengthMedian), xlim=log10(c(100,420000)), QuantileLast_last$median, ylab='Median Sigma',
     xlab='Distance to End of Gene' ,pch=20,ylim=c(0,1),col='black',main=paste(column, x),xaxt='n',yaxt='n')
  points(log10(QuantileLast_nonLast$LengthMedian),  QuantileLast_nonLast$median,pch=20, col='red',yaxt='n')
  formatAxes()
  legend('topleft',col=c(2,1),lwd=1,leg=c('Last','nonLast')) 
}); plot.off()

QuantileData200_nuc_nonLast = getQuantilMedians(subset(introns3,isLast==0),col1='Dist2End',column="Chromatin",400) 
model_nuc200_nonLas = nls(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileData200_nuc_nonLast,start=c(g=1,k=5e-5,A=.86,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,1,Inf),algo='port');
model_nuc200_nonLas_noRT = nls(median ~ MeanExp(LengthMedian,k,A), data=QuantileData200_nuc_nonLast,start=c(k=2e-3,A=.85),lower=c(0,0),upper=c(Inf,1),algo='port');

QuantileData20_nuc_last = getQuantilMedians(subset(introns3,isLast==1),col1='Dist2End',column="Chromatin",20) 
plot(log10(QuantileData200_nuc_nonLast$LengthMedian), xlim=log10(c(100,420000)), QuantileData200_nuc_nonLast$median, ylab='Median Sigma',
     xlab='Distance to End of Gene' ,pch=20,ylim=c(0,1),col='red',main=paste('Nuc A -: nonLast'),xaxt='n',yaxt='n')
points(log10(QuantileData20_nuc_last$LengthMedian),  QuantileData20_nuc_last$median,pch=20, col='black',yaxt='n')
lines(log10(xx<- 10^seq((0),log10(2000000),length.out=1000)),callMeanGamma2(xx,coef(model_nuc200_nonLas)),type='l',col='blue',lwd=2)  
lines(log10(xx<- 10^seq((0),log10(2000000),length.out=1000)),MeanExp(xx,coef(model_nuc200_nonLas_noRT)[1],coef(model_1000_noRT)[2]),type='l',col='blue',lwd=2,lty=2)     

formatAxes()
legend('topleft',col=c(2,1),lwd=1,leg=c('Last','nonLast')) 

#QuantileData_chrom_100_cytA2 = getQuantilMedians(introns3[cyt_polyA_coverage > 0.9 &  chrom_polyA_coverage > 0.75,],col1='Dist2End',column="Chromatin",100) 
#
#
#
#model_nuc_100_nonA = nls(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileData_nuc_100_nonA,start=c(g=1,k=5e-5,A=.86,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,1,Inf),algo='port');
#
#
























































































































##QuantileData_HK_nucPlus  = getQuantilMedians(subset(introns3,isHK==1),col1='Dist2End',column="splicingAPlus",100)  
#QuantileData_non_nucPlus  = getQuantilMedians(subset(introns3, isHK==0),col1='Dist2End',column="splicingAPlus",100) 
#QuantileData_HK_non3UTR  = getQuantilMedians(subset(introns3,isHK==1 & downFrame > -1),col1='Dist2End',column="Chromatin",100)  
#QuantileData_non_non3UTR  = getQuantilMedians(subset(introns3, isHK==0  & downFrame > -1),col1='Dist2End',column="Chromatin",100) 
#
#QuantileData_HK_nuc_non3UTR  = getQuantilMedians(subset(introns3,isHK==1 & downFrame > -1),col1='Dist2End',column="splicingAMinus",100)  
#QuantileData_non_nuc_non3UTR  = getQuantilMedians(subset(introns3, isHK==0  & downFrame > -1),col1='Dist2End',column="splicingAMinus",100) 
#QuantileData_HK_nucPlus_non3UTR  = getQuantilMedians(subset(introns3,isHK==1 & downFrame > -1),col1='Dist2End',column="splicingAPlus",100)  
#QuantileData_non_nucPlus_non3UTR  = getQuantilMedians(subset(introns3, isHK==0  & downFrame > -1),col1='Dist2End',column="splicingAPlus",100) 


#plot.dev("Chromatin_splicing100_HKvNon.pdf",'pdf', height=8, width=14)
#par(mfrow=c(2,4))
#
#plot(log10(QuantileData_HK$LengthMedian), xlim=log10(c(100,420000)), QuantileData_HK$median, ylab='Median psi-intron', xlab='LOG10 Distance to End of Gene' ,pch=20,ylim=c(0,1),col='red',main='Chromatin',xaxt='n')
#points(log10(QuantileData_non$LengthMedian),  QuantileData_non$median,pch=20, col='black',yaxt='n')
#formatAxes()
#lines(log10(xx<- 10^seq((0),log10(500000),length.out=100)),MeanGamma2(xx,coef(model_100_HK)[1], coef(model_100_HK)[2],coef(model_100_HK)[3], coef(model_100_HK)[4]),type='l',col='red',lwd=1) 
#lines(log10(xx<- 10^seq((0),log10(500000),length.out=100)),MeanGamma2(xx,coef(model_100_non)[1], coef(model_100_non)[2],coef(model_100_non)[3], coef(model_100_non)[4]),type='l',col='black',lwd=1) 
#legend('topleft',col=c(2,1),lwd=1,leg=c('HK','non'))
#
#plot(log10(QuantileData_HK_nuc$LengthMedian), xlim=log10(c(100,420000)), QuantileData_HK_nuc$median, ylab='Median psi-intron', xlab='LOG10 Distance to End of Gene' ,pch=22,ylim=c(0,1),col='red',main='Nuclear poly(A) -')
#points(log10(QuantileData_non_nuc$LengthMedian),  QuantileData_non_nuc$median,pch=22, col='black')
#lines(log10(xx<- 10^seq((0),log10(500000),length.out=100)),MeanGamma2(xx,coef(model_100_HK_nuc)[1], coef(model_100_HK_nuc)[2],coef(model_100_HK_nuc)[3], coef(model_100_HK_nuc)[4]),type='l',col='red',lwd=1) 
#lines(log10(xx<- 10^seq((0),log10(500000),length.out=100)),MeanGamma2(xx,coef(model_100_non_nuc)[1], coef(model_100_non_nuc)[2],coef(model_100_non_nuc)[3], coef(model_100_non_nuc)[4]),type='l',col='black',lwd=1) 
#legend('topleft',col=c(2,1),lwd=1,leg=c('HK','non'))
#
#
#plot(log10(QuantileData_HK$LengthMedian), xlim=log10(c(100,420000)), QuantileData_HK$median, ylab='Median psi-intron', xlab='LOG10 Distance to End of Gene' ,ylim=c(0,1),col='red',main='HK',pch=20)
#points(log10(QuantileData_HK_nuc$LengthMedian), QuantileData_HK_nuc$median, pch=22,col='red')
#lines(log10(xx<- 10^seq((0),log10(500000),length.out=100)),MeanGamma2(xx,coef(model_100_HK)[1], coef(model_100_HK)[2],coef(model_100_HK)[3], coef(model_100_HK)[4]),type='l',col='red',lwd=1) 
#lines(log10(xx<- 10^seq((0),log10(500000),length.out=100)),MeanGamma2(xx,coef(model_100_HK_nuc)[1], coef(model_100_HK_nuc)[2],coef(model_100_HK_nuc)[3], coef(model_100_HK_nuc)[4]),type='l',col='red',lwd=1) 
#
#plot(log10(QuantileData_non$LengthMedian), xlim=log10(c(100,420000)), QuantileData_non$median, ylab='Median psi-intron', xlab='LOG10 Distance to End of Gene' ,ylim=c(0,1),col='black',main='non',pch=20)
#points(log10(QuantileData_non_nuc$LengthMedian), QuantileData_non_nuc$median, pch=22,col='black')
#lines(log10(xx<- 10^seq((0),log10(500000),length.out=100)),MeanGamma2(xx,coef(model_100_non_nuc)[1], coef(model_100_non_nuc)[2],coef(model_100_non_nuc)[3], coef(model_100_non_nuc)[4]),type='l',col='black',lwd=1) 
#lines(log10(xx<- 10^seq((0),log10(500000),length.out=100)),MeanGamma2(xx,coef(model_100_non)[1], coef(model_100_non)[2],coef(model_100_non)[3], coef(model_100_non)[4]),type='l',col='black',lwd=1) 
#
#plot(log10(QuantileData_HK_nucPlus$LengthMedian), xlim=log10(c(100,420000)), QuantileData_HK_nucPlus$median, ylab='Median psi-intron', xlab='LOG10 Distance to End of Gene' ,pch=24,ylim=c(0,1),col='red',main='Nuclear poly(A) +')
#points(log10(QuantileData_non_nucPlus$LengthMedian),  QuantileData_non_nucPlus$median,pch=24, col='black')
#legend('topleft',col=c(2,1),lwd=1,leg=c('HK','non'))
#
#
#plot(log10(QuantileData_non_nuc$LengthMedian), xlim=log10(c(100,420000)), QuantileData_non_nuc$median, ylab='Median psi-intron', xlab='LOG10 Distance to End of Gene' ,ylim=c(0,1),col='black',main='non',pch=22)
#points(log10(QuantileData_non_nucPlus$LengthMedian), QuantileData_non_nucPlus$median, pch=24,col='black')
#
#plot(log10(QuantileData_HK_nuc_non3UTR$LengthMedian), xlim=log10(c(100,420000)), QuantileData_HK_nuc_non3UTR$median, ylab='Median psi-intron', xlab='LOG10 Distance to End of Gene' ,pch=22,ylim=c(0,1),col='red',main='Nuclear Poly(A) -\nno 3\'UTR introns included')
#points(log10(QuantileData_non_nuc_non3UTR$LengthMedian),  QuantileData_non_nuc_non3UTR$median,pch=22, col='black')
#
#plot(log10(QuantileData_HK_nucPlus_non3UTR$LengthMedian), xlim=log10(c(100,420000)), QuantileData_HK_nucPlus_non3UTR$median, ylab='Median psi-intron', xlab='LOG10 Distance to End of Gene' ,pch=24,ylim=c(0,1),col='red',main='Nuclear Poly(A) +\n no 3\'UTR introns included')
#points(log10(QuantileData_non_nucPlus_non3UTR$LengthMedian),  QuantileData_non_nucPlus_non3UTR$median,pch=24, col='black')
#
#plot.off()


## TESTING for now ...

# Try to model Nuclear as a combination of Chromatin and Nucleoplasmic ...
#
#get_nucleoplasmic1 = function(chrom_data,nuc_data, chrom_frac)
#  (nuc_data - chrom_data*chrom_frac)/(1-chrom_frac)
#
#plot.dev("Nucleoplasm_splicing100.pdf",'pdf', height=8, width=8)
#par(mfrow=c(2,2))
#
## First assume a proportion of the signal that is coming from chromatin
#sapply(c(.5,.6,.7,.8),function(chrom_frac){
#  Nucleoplasm_Qdata <<- get_nucleoplasmic1(QuantileData100[,2],QuantileData100_nuc[,2],chrom_frac)
#  Nucleoplasm_Qdata_non <<- get_nucleoplasmic1(QuantileData_non[,2],QuantileData_non_nuc[,2],chrom_frac)
#  Nucleoplasm_Qdata_HK <<- get_nucleoplasmic1(QuantileData_HK[,2],QuantileData_HK_nuc[,2],chrom_frac)
#  
#  plot(log10(QuantileData100$LengthMedian), xlim=log10(c(100,420000)), QuantileData100$median, ylab='Median psi-intron', xlab='LOG10 Distance to End of Gene' ,pch=20,ylim=c(0,1),col='black',main=sprintf('All genes: Nucleoplasm estimate\n(chromatin fraction=%s)',chrom_frac))
#  points(log10(QuantileData100$LengthMedian), Nucleoplasm_Qdata,pch=8,col='blue')
#  points(log10(QuantileData100_nuc$LengthMedian), QuantileData100_nuc$median,pch=22,cex=0.8)
#  legend('topleft',pch=c(20,24,8),leg=c('chromatin','nucleus','nucleoplasm'))
#
#  invisible(NULL)})
#plot.off()
#
#plot.dev("Nucleoplasm_splicing100_HKvNon.pdf",'pdf', height=10, width=10)
#par(mfcol=c(3,3))
#
## First assume a proportion of the signal that is coming from chromatin
#sapply(c(.5,.6,.7),function(chrom_frac){
#  Nucleoplasm_Qdata <<- get_nucleoplasmic1(QuantileData100[,2],QuantileData100_nuc[,2],chrom_frac)
#  Nucleoplasm_Qdata_non <<- get_nucleoplasmic1(QuantileData_non[,2],QuantileData_non_nuc[,2],chrom_frac)
#  Nucleoplasm_Qdata_HK <<- get_nucleoplasmic1(QuantileData_HK[,2],QuantileData_HK_nuc[,2],chrom_frac)
#  
#  plot(log10(QuantileData_HK$LengthMedian), xlim=log10(c(100,420000)), QuantileData_HK$median, ylab='Median psi-intron', xlab='LOG10 Distance to End of Gene' ,pch=20,ylim=c(0,1),col='red',main=sprintf('HK genes: Nucleoplasm estimate\n(chromatin fraction=%s)',chrom_frac))
#  points(log10(QuantileData_HK$LengthMedian), Nucleoplasm_Qdata_HK,pch=8,col='blue')
#  points(log10(QuantileData_HK_nuc$LengthMedian), QuantileData_HK_nuc$median,pch=22,cex=0.8,col='red')
#  legend('topleft',pch=c(20,24,8),leg=c('chromatin','nucleus','nucleoplasm'))
#
#  plot(log10(QuantileData_non$LengthMedian), xlim=log10(c(100,420000)), QuantileData_non$median, ylab='Median psi-intron', xlab='LOG10 Distance to End of Gene' ,pch=20,ylim=c(0,1),col='black',main=sprintf('non-HK genes: Nucleoplasm estimate\n(chromatin fraction=%s)',chrom_frac))
#  points(log10(QuantileData_non$LengthMedian), Nucleoplasm_Qdata_non,pch=8,col='blue')
#  points(log10(QuantileData_non_nuc$LengthMedian), QuantileData_non_nuc$median,pch=22,cex=0.8)
#  legend('topleft',pch=c(20,24,8),leg=c('chromatin','nucleus','nucleoplasm'))
#  
#  plot(log10(QuantileData_non$LengthMedian), xlim=log10(c(100,420000)), Nucleoplasm_Qdata_non, ylab='Median psi-intron', xlab='LOG10 Distance to End of Gene' ,pch=8,ylim=c(0,1),col='black',main=sprintf('Nucleoplasm estimate\n(chromatin fraction=%s)',chrom_frac))
#  points(log10(QuantileData_HK$LengthMedian), Nucleoplasm_Qdata_HK,pch=,col='red',cex=0.8)
#  legend('topleft',col=c(2,1),leg=c('HK','non'),pch=8)
#
#
#  invisible(NULL)})
#  
#plot.off()
#
#