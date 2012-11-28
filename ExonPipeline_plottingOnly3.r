### Script to plot data after pipeline has run


# **********************************************************************************************************************************
#
# 10/23/2012: Modified this to use the NEW way of using BasicSeqToNuc, which does not average over any part of the intron. 
#  Here, the sequences used to calculate Nucleosome occupancy are totally separate from those used to calculate splice site strength
#
#	11/20/2012: Modified this again to use NORMALIZED values of Energy and Log(splicing scores)
# **********************************************************************************************************************************


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
energyData = read.delim(settings$nucNormScore,head=F,stringsAsFactors=F)
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

#write.table( cbind(paste(settings$CommonName,c("HK","nonHK")),table(mergedData$isHK)),file="../Species_HK_list.txt",append=T,row.names=F,col.names=F,quote=F,sep="\t")
save.image("PlottedData3.RData")
 
#load("PlottedData3.RData")
cat("	Plotting data:\n") 

#### 1) Plot Barplots

palette(colorRampPalette(c('white','orange'))(8))
cols=palette()

png(sprintf('%s_Barplots3.png',settings$CommonName),height=1200,width=300)
par(mfrow=c(4,1),cex.axis=1.5,cex=1.1,cex.main=1.1)
  makeBarplot(LENGTH,'Exon Length',doLog=T,col=cols[2])
  makeBarplot(DONOR,'Donor Splice Strength Z score',doLog=F,col=cols[4])
  makeBarplot(ACCEPTOR,'Acceptor Splice Strength Z score',doLog=F,col=cols[6])
  makeBarplot(ENERGY,'Nucleosome Energy Z score',doLog=F,col=cols[8])
dev.off()


#### 2) Comparisons WITHIN LASTS (including housekeeping genes)


### Do some plots based on splitting up the data into groups
mergedData$FeatureGroup = splitByVector(mergedData$FeatureCount,quantile(mergedData$FeatureCount,na.rm=T,(1:4)/5))
mergedData$LengthGroup = splitByVector(mergedData$length,quantile(mergedData$length,na.rm=T,(1:4)/5))

## SPLIT BY LENGTH
png(sprintf('%s_simple_splitByLength3.png',settings$CommonName),height=1200,width=400)
par(mfrow=c(3,1),cex.axis=1.2,cex=1.3,cex.main=1.1,las=3,mar=c(6.1 ,4.1 ,4.1, 2.1))

xlabs = tapply(mergedData$length,mergedData$Group,function(x)paste(range(x),collapse='-'))

plotSplitsHK(mergedData,'Energy','Nucleosome Stability\nAs a function of last exon length',names=xlabs,ylab='Energy Z Score',xlab='')
plotSplits(mergedData,'Energy','',add=T,type='b',lwd=2)
# 
plotSplitsHK(mergedData,'Acceptor',"3' splice site strength\nAs a function of last exon length",names=xlabs, ylab='Splicing Z Score',xlab='')
plotSplits(mergedData,'Acceptor','',add=T,type='b',lwd=2)
#

plotSplitsHK(mergedData,'FeatureCount',"Number of Exons in Gene\nAs a function of last exon length",names=xlabs,ylab='# Exons',xlab='')
plotSplits(mergedData,'FeatureCount','',add=T,type='b',lwd=2)
#
dev.off()

## SPLIT BY LENGTH GROUP, which is by length quintiles
png(sprintf('%s_simple_splitByLengthGroup3.png',settings$CommonName),height=1200,width=400)
par(mfrow=c(3,1),cex.axis=1.2,cex=1.3,cex.main=1.1,las=3,mar=c(6.1 ,4.1 ,4.1, 2.1))

xlabs = tapply(mergedData$length,mergedData$LengthGroup,function(x)paste(range(x),collapse='-'))

plotSplitsHK(mergedData,'Energy','Nucleosome Stability\nAs a function of last exon length',Group='LengthGroup',names=xlabs,ylab='Energy Score',xlab='')
plotSplits(mergedData,'Energy','',add=T,type='b',lwd=2,Group='LengthGroup',xaxt='n')
# 
plotSplitsHK(mergedData,'Acceptor',"3' splice site strength\nAs a function of last exon length",Group='LengthGroup',names=xlabs, ylab='PWM Score',xlab='')
plotSplits(mergedData,'Acceptor','',add=T,type='b',lwd=2,Group='LengthGroup')
#

plotSplitsHK(mergedData,'FeatureCount',"Number of Exons in Gene\nAs a function of last exon length",Group='LengthGroup',names=xlabs,ylab='# Exons',xlab='')
plotSplits(mergedData,'FeatureCount','',add=T,type='b',lwd=2,Group='LengthGroup')
#
dev.off()
