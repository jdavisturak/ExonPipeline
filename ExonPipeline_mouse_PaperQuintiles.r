## Data mining analysis and FIGURES
# For 'Shorter Version' of the paper: Started 3/12/13


source("/home/jeremy/ExonPipeline/ExonPipeline_functions.r")
source("/home/jeremy/ExonPipeline/ExonPipeline_analysis_and_plotting_Functions.r")
source("/home/jeremy/Code/useful_R/useful_R.r")

setwd("/home/home/jeremy/ExonPipeline/mm9")
load("PlottedData4.RData")
library(gplots)

mergedData$nonLastLength = mergedData$GeneLength - mergedData$length
                                                                      
# 3/27/13
# Split it up into 5 groups BASED ON 5 QUINTILES OF SPLICING PROBABILITIES, GIVEN OUR PARAMETERS
#
k_splice = 2 # 2/minute
elong = 2010 # bp/minute  % From matlab:  1/((1+KP/KU)/(KE+KP)) = 2010
lengths =  -log(1-(1:4)/5)*elong/k_splice; lengths              
#[1]  223.1436  510.8256  916.2907 1609.4379

#######################
## Apply some groupings
mergedData$AcceptorGroup = splitByVector(mergedData$Acceptor,quantile(mergedData$Acceptor,na.rm=T,(1:4)/5)) # split into 5 groups
mergedData$EnergyGroup = splitByVector(mergedData$Energy,quantile(mergedData$Energy,na.rm=T,(1:4)/5)) # split into 5 groups
mergedData$LengthGroup = splitByVector(mergedData$length,quantile(mergedData$length,na.rm=T,(1:4)/5)) # split into 5 groups
mergedData$LengthGroup2 = splitByVector(mergedData$length,quantile(mergedData$length,na.rm=T,(1:49)/50)) 
mergedData$LengthGroup3 = splitByVector(mergedData$length,quantile(mergedData$length,na.rm=T,(1:19)/20)) 
mergedData$LengthGroup4 = splitByVector(mergedData$length,quantile(mergedData$length,na.rm=T,(1:9)/10)) # split into 10 groups
mergedData$FeatureGroup = splitByVector(mergedData$FeatureCount,quantile(mergedData$FeatureCount,na.rm=T,(1:4)/5)) # split into 5 groups
mergedData$FeatureGroup3 = splitByVector(mergedData$FeatureCount,quantile(mergedData$FeatureCount,na.rm=T,(1:19)/20)) 
mergedData$FeatureGroup4 = splitByVector(mergedData$FeatureCount,quantile(mergedData$FeatureCount,na.rm=T,(1:9)/10)) 
mergedData$MultipleIntrons = mergedData$FeatureCount > 2 # FeatureCount counts the number of EXONS

mergedData$LengthGroup_S = splitByVector(mergedData$length,lengths) # split into 5 groups

### Using 5 exon groups based on quintiles of splicing probability
## I also made a version '.2' which used  2010 as the elongation parameter
pdf(sprintf('%s_ShorterPaper_Splits_SpliceGroups.pdf',settings$CommonName),height=10,width=8)
par(mfrow=c(2,2),cex.axis=1.2,cex=1.3,cex.main=1.1,las=3,mar=c(6.1 ,4.1 ,4.1, 2.1))
plotSplits(mergedData,'Energy','', Group='LengthGroup_S',  xlab='Last Exon Quintiles',barcol='grey', ylab='Nucleosome Stability Z-score',gap=0,type='l', ylim=c(-.4,.2),pch=NA)
plotSplits(mergedData,'Acceptor','',Group='LengthGroup_S', xlab='Last Exon Quintiles',barcol='grey', ylab="3' Splice Site Strength Score",gap=0,type='l',ylim=c(0.1,0.4),pch=NA)
plotSplitsHK(mergedData,'Energy',"", Group='LengthGroup_S',  xlab='Last Exon Quintiles',ylab='Nucleosome Stability Z-score')
plotSplitsHK(mergedData,'Acceptor',"", Group='LengthGroup_S',  xlab='Last Exon Quintiles', ylab="3' Splice Site Strength Score")
dev.off()

### 
pdf(sprintf('%s_ShorterPaper_barplots_blackGrey2.pdf',settings$CommonName),height=20,width=4)
par(mfrow=c(4,1),cex.axis=1.2,cex=0.7,cex.main=1.1,mar=c(2,3,2,1),lwd=4)
   errorcol=col=c("grey","black")
  makeBarplot(LENGTH,'',doLog=T,col=c("black","grey"),addStats=F,yaxt='n',border=NA,errorcol=errorcol)
  axis(2,at=0:4,lab=10^(0:4))
  makeBarplot(ENERGY,'',doLog=F,col=c("black","grey"), YLIM=c(-.25,0.25),addStats=F,border=NA,errorcol=errorcol)
  makeBarplot(DONOR,'',doLog=F,col=c("black","grey"), YLIM=c(-.25,0.25),addStats=F,border=NA,errorcol=errorcol)
  makeBarplot(ACCEPTOR,'',doLog=F,col=c("black","grey"), YLIM=c(-.25,0.25),addStats=F,border=NA,errorcol=errorcol)
dev.off()


### GRO-seq data 
load("dataWithGRO.RData")
### Runon as a function of last exon length
dataWithGRO$LengthGroup2 = splitByVector(dataWithGRO$length,quantile(dataWithGRO$length,na.rm=T,(1:49)/50)) # split into equal sized groups
dataWithGRO$LengthGroup3 = splitByVector(dataWithGRO$length,quantile(dataWithGRO$length,na.rm=T,(1:19)/20))
dataWithGRO$LengthGroup4 = splitByVector(dataWithGRO$length,quantile(dataWithGRO$length,na.rm=T,(1:9)/10)) # split into 10 groups
dataWithGRO$FeatureGroup3 = splitByVector(dataWithGRO$FeatureCount,quantile(dataWithGRO$FeatureCount,na.rm=T,(1:19)/20))
dataWithGRO$FeatureGroup4 = splitByVector(dataWithGRO$FeatureCount,quantile(dataWithGRO$FeatureCount,na.rm=T,(1:9)/10))
dataWithGRO$MultipleIntrons = dataWithGRO$FeatureCount > 2 # FeatureCount counts the number of EXONS
dataWithGRO$LengthGroup_S = splitByVector(dataWithGRO$length,lengths) # split into 5 groups

## Redo barplots, but with filled bars and PDF

pdf(sprintf('%s_ShorterPaper_Barplots4.pdf',settings$CommonName),height=8,width=2)
par(mfrow=c(4,1),cex.axis=1.2,cex=0.7,cex.main=1.1,mar=c(2,3,2,1),lwd=4)
   errorcol=rgb(1/2,1/4,0.26)
  makeBarplot(LENGTH,'',doLog=T,col=c("black","grey"),addStats=F,yaxt='n',border=c("black","grey"),errorcol=errorcol)
  axis(2,at=0:4,lab=10^(0:4))
  makeBarplot(ENERGY,'',doLog=F,col=c("black","grey"), YLIM=c(-.25,0.25),addStats=F,border=c("black","grey"),errorcol=errorcol)
  makeBarplot(DONOR,'',doLog=F,col=c("black","grey"), YLIM=c(-.25,0.25),addStats=F,border=c("black","grey"),errorcol=errorcol)
  makeBarplot(ACCEPTOR,'',doLog=F,col=c("black","grey"), YLIM=c(-.4,0.4),addStats=F,border=c("black","grey"),errorcol=errorcol)
dev.off()


## Here I'm going to write a new function to replace 'plotCI' in the call of plotSplits

# We are passed in the splits. Start by plotting first points, then RECTANGLES, and then add the points back!

plotCI = function(means, uiw,  main='',pointCol='black',barcol='grey',ylim=NULL,xaxt=NULL,labels=NULL, ...){
    MIN = min(means-uiw)            
    MAX = max(means+uiw)
    if (is.null(ylim)) ylim=c(MIN,MAX)
    
    xx = 1:length(means)
    plot(xx, means,pch=NA, ylim=ylim, xlim = c(-0.5,length(means) + 0.5), ...);
    
    ## Draw actual rectangles
    #rect(xleft = xx-0.5, ybottom = means-uiw, xright = xx+0.5, ytop=means+uiw, col=barcol, border=NA)
    polygon(c(xx,rev(xx)), y= c(means-uiw, rev(means+uiw)), col=barcol, border=NA)
    
    ## Redraw points on TOP of rectangles
    lines(xx, means, col=pointCol)

}

## Split by QUINTILES: 'LengthGroup'

pdf(sprintf('%s_ShorterPaper_Splits_basic.pdf',settings$CommonName),height=10,width=8)
par(mfrow=c(2,2),cex.axis=1.2,cex=1.3,cex.main=1.1,las=3,mar=c(6.1 ,4.1 ,4.1, 2.1))
plotSplits(mergedData,'Energy','', Group='LengthGroup3',  xlab='Last Exon Quintiles',barcol='grey', ylab='Nucleosome Stability Z-score')
plotSplits(mergedData,'Acceptor','',Group='LengthGroup3', xlab='Last Exon Quintiles',barcol='grey', ylab="3' Splice Site Strength Score")
plotSplits(dataWithGRO,'Runon','', Group='LengthGroup3',  xlab='Last Exon Quantiles',barcol='grey', ylab='Runon Length (bp)')
plotSplits(dataWithGRO,'FeatureCount','', Group='LengthGroup3',  xlab='Last Exon Quantiles',barcol='grey', ylab='# Exons')
dev.off()

## Look at differently by HK genes

pdf(sprintf('%s_ShorterPaper_Splits_HK.pdf',settings$CommonName),height=10,width=8)
par(mfrow=c(2,2),cex.axis=1.2,cex=1.3,cex.main=1.1,las=3,mar=c(6.1 ,4.1 ,4.1, 2.1))
plotSplitsHK(mergedData,'Energy',"", Group='LengthGroup3',  xlab='Last Exon Quintiles',ylab='Nucleosome Stability Z-score')
plotSplitsHK(mergedData,'Acceptor',"", Group='LengthGroup3',  xlab='Last Exon Quintiles', ylab="3' Splice Site Strength Score")
plotSplitsHK(dataWithGRO,'Runon',"", Group='LengthGroup3',  xlab='Last Exon Quintiles', ylab="Runon Length (bp)")
plotSplitsHK(dataWithGRO,'FeatureCount',"", Group='LengthGroup3',  xlab='Last Exon Quintiles', ylab="# Exons")
dev.off()
 
## Look at differently by whether gene has more than 1 intron
pdf(sprintf('%s_ShorterPaper_Splits_MultipleIntrons.pdf',settings$CommonName),height=10,width=8)
par(mfrow=c(2,2),cex.axis=1.2,cex=1.3,cex.main=1.1,las=3,mar=c(6.1 ,4.1 ,4.1, 2.1))
plotSplitsHK(mergedData,'Energy',"", Group='LengthGroup3', checkColumn='MultipleIntrons', ylab='Nucleosome Stability Z-score',xlab='Last Exon Quintiles')
plotSplitsHK(mergedData,'Acceptor',"", Group='LengthGroup3',checkColumn='MultipleIntrons', ylab="3' Splice Site Strength Score",xlab='Last Exon Quintiles')
plotSplitsHK(dataWithGRO,'Runon',"", Group='LengthGroup3',checkColumn='MultipleIntrons',  ylab="Runon Length (bp)",xlab='Last Exon Quintiles')
dev.off()
 
pdf(sprintf('%s_ShorterPaper_Splits4_HK.pdf',settings$CommonName),height=10,width=8)
par(mfrow=c(2,2),cex.axis=1.2,cex=1.3,cex.main=1.1,las=3,mar=c(6.1 ,4.1 ,4.1, 2.1))
plotSplitsHK(mergedData,'Energy',"", Group='LengthGroup4',  xlab='Last Exon Quintiles',ylab='Nucleosome Stability Z-score')
plotSplitsHK(mergedData,'Acceptor',"", Group='LengthGroup4',  xlab='Last Exon Quintiles', ylab="3' Splice Site Strength Score")
plotSplitsHK(dataWithGRO,'Runon',"", Group='LengthGroup4',  xlab='Last Exon Quintiles', ylab="Runon Length (bp)")
plotSplitsHK(mergedData,'FeatureCount',"", Group='LengthGroup4',  xlab='Last Exon Quintiles', ylab="# Exons")
dev.off()

pdf(sprintf('%s_ShorterPaper_Splits_FeatureGroups4_HK.pdf',settings$CommonName),height=10,width=8)
par(mfrow=c(2,2),cex.axis=1.2,cex=1.3,cex.main=1.1,las=3,mar=c(6.1 ,4.1 ,4.1, 2.1))
plotSplitsHK(mergedData,'length','', Group='FeatureGroup4',  xlab='Intron number Quantiles', ylab='Last Exon Length (bp)')
plotSplitsHK(mergedData,'Energy','', Group='FeatureGroup4',  xlab='Intron number Quantiles', ylab='Nucleosome Stability Z-score')
plotSplitsHK(mergedData,'Acceptor','',Group='FeatureGroup4', xlab='Intron number Quantiles', ylab="3' Splice Site Strength Score")
plotSplitsHK(dataWithGRO,'Runon','', Group='FeatureGroup4',  xlab='Intron number Quantiles', ylab='Runon Length')
dev.off() 

pdf(sprintf('%s_ShorterPaper_Splits_FeatureGroups1_HK.pdf',settings$CommonName),height=10,width=8)
par(mfrow=c(2,2),cex.axis=1.2,cex=1.3,cex.main=1.1,las=3,mar=c(6.1 ,4.1 ,4.1, 2.1))
plotSplitsHK(mergedData,'length','', Group='FeatureGroup',  xlab='Intron number Quantiles', ylab='Last Exon Length (bp)')
plotSplitsHK(mergedData,'Energy','', Group='FeatureGroup',  xlab='Intron number Quantiles', ylab='Nucleosome Stability Z-score')
plotSplitsHK(mergedData,'Acceptor','',Group='FeatureGroup', xlab='Intron number Quantiles', ylab="3' Splice Site Strength Score")
plotSplitsHK(dataWithGRO,'Runon','', Group='FeatureGroup',  xlab='Intron number Quantiles', ylab='Runon Length')
dev.off() 

library(RColorBrewer)
library(devEMF)
COLS = col2rgb(brewer.pal(5, "Set1")[c(2,1)])
COLS2 =  c(rgb(COLS[1,1]/256,COLS[2,1]/256,COLS[3,1]/256,alpha=1), rgb(COLS[1,2]/256,COLS[2,2]/256,COLS[3,2]/256,alpha=1))

pdf(sprintf('%s_ShorterPaper_GRO_density.pdf',settings$CommonName))
CEX=0.8; CEX.LAB=1.5
par(cex=CEX*1.2, cex.lab=CEX.LAB/.8,cex.axis=CEX.LAB)
a=plot2Densities(list(log10(dataWithGRO$length),log10(dataWithGRO$Runon+dataWithGRO$length)),'Distribution of last exon lengths','Exon Length',cols=COLS2,lwd=3,leg=c('mRNA','pre-mRNA'),xaxt='n',leg.at='topleft')
axis(1,at=1:6,lab=10^(1:6))
grid()
dev.off() 
 

## Now base things off of the Feature group 
pdf(sprintf('%s_ShorterPaper_Splits_FeatureGroups.pdf',settings$CommonName),height=10,width=8)
par(mfrow=c(2,2),cex.axis=1.2,cex=1.3,cex.main=1.1,las=3,mar=c(6.1 ,4.1 ,4.1, 2.1))
plotSplits(mergedData,'length','', Group='FeatureGroup3',  xlab='Intron number Quantiles',barcol='grey', ylab='Last Exon Length (bp)')
plotSplits(mergedData,'Energy','', Group='FeatureGroup3',  xlab='Intron number Quantiles',barcol='grey', ylab='Nucleosome Stability Z-score')
plotSplits(mergedData,'Acceptor','',Group='FeatureGroup3', xlab='Intron number Quantiles',barcol='grey', ylab="3' Splice Site Strength Score")
plotSplits(dataWithGRO,'Runon','', Group='FeatureGroup3',  xlab='Intron number Quantiles',barcol='grey', ylab='Runon Length')
dev.off()

pdf(sprintf('%s_ShorterPaper_Splits_FeatureGroups_HK.pdf',settings$CommonName),height=10,width=8)
par(mfrow=c(2,2),cex.axis=1.2,cex=1.3,cex.main=1.1,las=3,mar=c(6.1 ,4.1 ,4.1, 2.1))
plotSplitsHK(mergedData,'length','', Group='FeatureGroup3',  xlab='Intron number Quantiles', ylab='Last Exon Length (bp)')
plotSplitsHK(mergedData,'Energy','', Group='FeatureGroup3',  xlab='Intron number Quantiles', ylab='Nucleosome Stability Z-score')
plotSplitsHK(mergedData,'Acceptor','',Group='FeatureGroup3', xlab='Intron number Quantiles', ylab="3' Splice Site Strength Score")
plotSplitsHK(dataWithGRO,'Runon','', Group='FeatureGroup3',  xlab='Intron number Quantiles', ylab='Runon Length')
dev.off() 


# make table of HK genes: 2 introns vs more
write.delim(subset(dataWithGRO,isHK==1 & FeatureCount==2,select='GeneName'),File='HK_1-intron.txt')
write.delim(subset(dataWithGRO,isHK==1 & FeatureCount >= 2,select='GeneName'),File='HK_multiple-introns.txt')

write.delim(subset(dataWithGRO,FeatureCount==2,select='GeneName'),File='All_1-intron.txt')


write.delim(subset(dataWithGRO,FeatureCount==9,select='GeneName'),File='All_9-intron.txt')







 
 
###############################################################
#### 12/4/12
# Going to output stuff for matlab so I can redo Fig 4 with the Z scores: also going to make Fig 5C, which uses the Read-through
lengths_last =  tapply(mergedData$length,mergedData$Group,mean)
strength_nonLast = mean(ACCEPTOR[ACCEPTOR[,2]==0,1],na.rm=T)
strength_lasts_avg = tapply(mergedData$Acceptor,mergedData$Group,mean)
energy_lasts_avg = tapply(mergedData$Energy,mergedData$Group,mean)
runon_lasts_avg = tapply(dataWithGRO$Runon,dataWithGRO$Group,mean,na.rm=T)

write.delim(cbind(lengths_last,strength_lasts_avg,strength_nonLast,energy_lasts_avg,runon_lasts_avg), 'Normalized_averages_byLength_withRunon.txt') 


#############################
###############################################################
#### 3/28/13
# Output for matlab: splicing quantiles 
lengths_last =  tapply(mergedData$length,mergedData$LengthGroup_S,mean)
strength_nonLast = mean(ACCEPTOR[ACCEPTOR[,2]==0,1],na.rm=T)
strength_lasts_avg = tapply(mergedData$Acceptor,mergedData$LengthGroup_S,mean)
energy_lasts_avg = tapply(mergedData$Energy,mergedData$LengthGroup_S,mean)
runon_lasts_avg = tapply(dataWithGRO$Runon,dataWithGRO$LengthGroup_S,mean,na.rm=T)


write.delim(cbind(lengths_last,strength_lasts_avg,strength_nonLast,energy_lasts_avg,runon_lasts_avg), 'Normalized_averages_byLength_SpliceGroups.txt') 

## Alex wants to look at "the ratio of measured last exon vs predicted last exon (= correction factor)"
dataWithGRO$correction_factor = (dataWithGRO$Runon +dataWithGRO$length) / dataWithGRO$length
#dataWithGRO$correction_factor2 = (dataWithGRO$Runon+dataWithGRO$length)
#dataWithGRO$correction_factor_random = sample(dataWithGRO$Runon,nrow(dataWithGRO)) / dataWithGRO$length

pdf(sprintf("%s_ShorterPaper_CorrectionFactor_SpliceGroups.2.pdf",settings$CommonName))
plotSplits(dataWithGRO,column='correction_factor','Correction Factor',Group='LengthGroup_S',ylab='pre-mRNA last exon length',xlab='Last Exon Group', lwd=2,col='black',gap=0,type='l',pch=NA);
dev.off()

pdf("Shorter fig 4.pdf", width=11,height=4)
layout(mat=matrix(1:4,nrow=1),widths=c(2,2,3.5,3.5))

makeBoxplot (data.frame(Value=dataWithGRO$length, isLast=dataWithGRO$isHK), mainText='',T, NAMES=c('non-HK','HK'),ylab='Predicted Exon Length (bp)', addStats=F,at=1:4,labels=10^(1:4),notch=T,lwd=1, border=c("black","gray"),cex=1)
makeBoxplot (data.frame(Value=dataWithGRO$Runon,  isLast=dataWithGRO$isHK), mainText='',T, NAMES=c('non-HK','HK'), ylab='Read-through (bp)', addStats=F,at=1:4,labels=10^(1:4),notch=T,lwd=1, border=c("black","gray"),cex=1)

## Make colors' default changeable!
plotSplitsHK(mergedData,column='Energy','',Group='LengthGroup_S',ylab='Nucleosome Stability Score',xlab='', names=NULL,colors=c("gray","black"),lwd=4)
plotSplitsHK(mergedData,'Acceptor','',Group='LengthGroup_S',ylab='Splice Site Strength Score',xlab='', names=NULL,colors=c("gray","black"),lwd=4)
dev.off()
