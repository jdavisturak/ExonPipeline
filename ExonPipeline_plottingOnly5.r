### Script to plot data after pipeline has run


# **********************************************************************************************************************************
#
# 10/23/2012: Modified this to use the NEW way of using BasicSeqToNuc, which does not average over any part of the intron. 
#  Here, the sequences used to calculate Nucleosome occupancy are totally separate from those used to calculate splice site strength
#
#	11/20/2012: Modified this again to use NORMALIZED values of Energy and Log(splicing scores)
#
#  4/4/13: Modified again to be run on data that used Constrained fixed lengths , i.e. doesn't go out into the intron
#
#**********************************************************************************************************************************

rm(list=ls())
library(gplots)
source("/home/jeremy/ExonPipeline/ExonPipeline_functions.r")
source("/home/jeremy/ExonPipeline/ExonPipeline_analysis_and_plotting_Functions.r")

settings = setup(getwd())

# 1) Comparisons of ALL to LAST
cat("	Reading in data to plot:\n")
#  Get lengths
exons = removeDuplicateStartANDEnds(read.delim("exons.txt",stringsAsFactors=F))
LENGTH = cbind(Value = abs(exons$start-exons$end), isLast = exons$isLast)

# Get Nucleosome energy
energyData = read.delim(settings$nucNormScore_constrained,head=F,stringsAsFactors=F)
ENERGY = formatData(energyData)

## Get splice site strength
donorData = read.delim(settings$DonorsNorm,head=F,stringsAsFactors=F)
DONOR = formatData(donorData)

## Get splice site strength
acceptorData = read.delim(settings$AcceptorsNorm,head=F,stringsAsFactors=F)
ACCEPTOR = formatData(acceptorData)
 
bins = c(1,2,5,10) 
HKlist = read.delim("/home/jeremy/ExonPipeline/hg19/HouseKeepingGenes.txt",skip=1)
mergedData = mergeLastData(exons,energyData,donorData,acceptorData,bins*180,HKlist)

save.image("PlottedData5.RData")
 
#load("PlottedData3.RData")
cat("	Plotting data:\n") 

#### 1) Plot Barplots

palette(colorRampPalette(c('white','orange'))(8))
cols=palette()

pdf(sprintf('%s_Barplots5.pdf',settings$CommonName),height=8,width=2)
par(mfrow=c(4,1),cex.axis=1.2,cex=0.7,cex.main=1.1,mar=c(2,3,2,1),lwd=4)
   errorcol=rgb(1/2,1/4,0.26)
  makeBarplot(LENGTH,'',doLog=T,col=c("black","grey"),addStats=F,yaxt='n',border=c("black","grey"),errorcol=errorcol)
  axis(2,at=1:4,lab=10^(1:4))
  makeBarplot(ENERGY,'',doLog=F,col=c("black","grey"), YLIM=c(-.25,0.25),addStats=F,border=c("black","grey"),errorcol=errorcol)
  makeBarplot(DONOR,'',doLog=F,col=c("black","grey"), YLIM=c(-.25,0.25),addStats=F,border=c("black","grey"),errorcol=errorcol)
  makeBarplot(ACCEPTOR,'',doLog=F,col=c("black","grey"), YLIM=c(-.25,0.25),addStats=F,border=c("black","grey"),errorcol=errorcol)
dev.off()


#### 2) Comparisons WITHIN LASTS (including housekeeping genes)


### Do some plots based on splitting up the data into groups

## SPLIT BY LENGTH
pdf(sprintf('%s_simple_splitByLength5.pdf',settings$CommonName),height=12,width=4)

par(mfrow=c(2,1),cex.axis=1.2,cex=1.3,cex.main=1.1,las=3,mar=c(6.1 ,4.1 ,4.1, 2.1))

xlabs = tapply(mergedData$length,mergedData$Group,function(x)paste(range(x),collapse='-'))

plotSplitsHK(mergedData,'Energy','Nucleosome Stability\nAs a function of last exon length',names=xlabs,ylab='Energy Z Score',xlab='')
plotSplits(mergedData,'Energy','',add=T,type='l',lwd=2,gap=0,sfrac=0.02)
# 
plotSplitsHK(mergedData,'Acceptor',"3' splice site strength\nAs a function of last exon length",names=xlabs, ylab='Splicing Z Score',xlab='')
plotSplits(mergedData,'Acceptor','',add=T,type='l',lwd=2,gap=0,sfrac=0.02)
#
dev.off()

pdf(sprintf('%s_simple2_splitByLength5.pdf',settings$CommonName),height=12,width=6)
par(mfrow=c(2,1))
plotSplits(mergedData,'Energy','',add=F,type='l',lwd=2,gap=0,sfrac=0.02)
par(new=F)
plotSplits(mergedData,'Acceptor','',add=F,type='l',lwd=2,gap=0,sfrac=0.02)
dev.off()

         
### Do with medians
pdf(sprintf('%s_simple_splitByLength5_Median.pdf',settings$CommonName),height=12,width=4)
par(mfrow=c(3,1),cex.axis=1.2,cex=1.3,cex.main=1.1,las=3,mar=c(6.1 ,4.1 ,4.1, 2.1))

xlabs = tapply(mergedData$length,mergedData$Group,function(x)paste(range(x),collapse='-'))

plotSplitsHK(mergedData,'Energy','Nucleosome Stability\nAs a function of last exon length',names=xlabs,avg='median',ylab='Energy Z Score',xlab='')
plotSplits(mergedData,'Energy','',avg='median',add=T,type='l',lwd=2,gap=0,sfrac=0.02)
# 
plotSplitsHK(mergedData,'Acceptor',"3' splice site strength\nAs a function of last exon length",avg='median',names=xlabs, ylab='Splicing Z Score',xlab='')
plotSplits(mergedData,'Acceptor','',avg='median',add=T,type='l',lwd=2,gap=0,sfrac=0.02)
#
dev.off()
