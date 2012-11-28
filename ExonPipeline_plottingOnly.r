### Script to plot data after pipeline has run

library(gplots)
source("/home/jeremy/ExonPipeline/ExonPipeline_functions.r")
source("/home/jeremy/ExonPipeline/ExonPipeline_analysis_and_plotting_Functions.r")

settings = setup(getwd())
#source("/home/jeremy/ExonPipeline/ExonPipeline.r")
             

cat("Doing Plots\n")

# 1) Comparisons of ALL to LAST

#  Get lengths
exons = removeDuplicateStartANDEnds(read.delim("exons.txt",stringsAsFactors=F))
LENGTH = cbind(Value = abs(exons$start-exons$end), isLast = exons$isLast)

# Get Nucleosome energy
energyData = read.delim(settings$nucTestFileFixed_mean,head=F,stringsAsFactors=F)
ENERGY = formatData(energyData)

## Get splice site strength
donorData = read.delim(settings$DonorScores,head=F,stringsAsFactors=F)
DONOR = formatData(donorData)

## Get splice site strength
acceptorData = read.delim(settings$AcceptorScores,head=F,stringsAsFactors=F)
ACCEPTOR = formatData(acceptorData)
 
 
mergedData = mergeLastData(exons,energyData,donorData,acceptorData,bins*180,HKlist)


#write.table( cbind(paste(settings$CommonName,c("HK","nonHK")),table(mergedData$isHK)),file="../Species_HK_list.txt",append=T,row.names=F,col.names=F,quote=F,sep="\t")
#save.image("PlottedData.RData")
 
load("PlottedData.RData")
 
### Plot Boxplots
#png(sprintf('%s.png',settings$CommonName))
#par(mfrow=c(2,2))
#  makeBoxplot(LENGTH,'Exon Length',doLog=T)
#  makeBoxplot(ENERGY,'Nucleosome Energy',doLog=F)
#  makeBoxplot(DONOR,'Donor Splice Strength',doLog=T)
#  makeBoxplot(ACCEPTOR,'Acceptor Splice Strength',doLog=T)
#dev.off()

### Plot Barplots

palette(colorRampPalette(c('white','orange'))(8))
cols=palette()

png(sprintf('%s_Barplots.png',settings$CommonName),height=1200,width=300)
par(mfrow=c(4,1),cex.axis=1.5,cex=1.1,cex.main=1.1)
  makeBarplot(LENGTH,'Exon Length',doLog=T,col=cols[2])
  makeBarplot(DONOR,'Donor Splice Strength',doLog=F,col=cols[4])
  makeBarplot(ACCEPTOR,'Acceptor Splice Strength',doLog=F,col=cols[6])
  makeBarplot(ENERGY,'Nucleosome Energy',doLog=F,col=cols[8])
dev.off()







## 2) Comparisons WITHIN LASTS (including housekeeping genes)
bins = c(1,2,5,10)
HKlist = read.delim("/home/jeremy/ExonPipeline/hg19/HouseKeepingGenes.txt",skip=1)


### Do some plots based on splitting up the data into groups
mergedData$FeatureGroup = splitByVector(mergedData$FeatureCount,quantile(mergedData$FeatureCount,na.rm=T,(1:4)/5))
mergedData$LengthGroup = splitByVector(mergedData$length,quantile(mergedData$length,na.rm=T,(1:4)/5))


png(sprintf('%s_splits.png',settings$CommonName),height=600,width=200)
par(mfrow=c(3,1))
plotSplits(mergedData,'Energy','Nucleosome Stability')
plotSplits(mergedData,'Donor',"5' splice site strength")
plotSplits(mergedData,'Acceptor',"3' splice site strength")
dev.off()

png(sprintf('%s_simple_splitByLength.png',settings$CommonName),height=1200,width=400)
par(mfrow=c(3,1),cex.axis=1.2,cex=1.3,cex.main=1.1,las=3,mar=c(6.1 ,4.1 ,4.1, 2.1))

xlabs = tapply(mergedData$length,mergedData$Group,function(x)paste(range(x),collapse='-'))

plotSplitsHK(mergedData,'Energy','Nucleosome Stability\nAs a function of last exon length',names=xlabs,ylab='Energy Score',xlab='')
plotSplits(mergedData,'Energy','',add=T,type='b',lwd=2)
# 
plotSplitsHK(mergedData,'Acceptor',"3' splice site strength\nAs a function of last exon length",names=xlabs, ylab='PWM Score',xlab='')
plotSplits(mergedData,'Acceptor','',add=T,type='b',lwd=2)
#

plotSplitsHK(mergedData,'FeatureCount',"Number of Exons in Gene\nAs a function of last exon length",names=xlabs,ylab='# Exons',xlab='')
plotSplits(mergedData,'FeatureCount','',add=T,type='b',lwd=2)
#
dev.off()

source("/home/jeremy/ExonPipeline/ExonPipeline_analysis_and_plotting_Functions.r")

png(sprintf('%s_simple_splitByLengthGroup.png',settings$CommonName),height=1200,width=400)
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
#print(par('mar'))
dev.off()

### Plot length of gene as a function of number of introns

#pdf(sprintf('%s_featureCount.pdf',settings$CommonName), height=8,width=3)
#par(mfrow=c(3,1),las=2)
#boxplot(log10(length) ~ FeatureGroup + isHK,data= mergedData, col=tenCols, ylab='Last Exon Length (log10 bp)', main='Number of Exons\n',xaxt='n')
#boxplot(Energy ~ FeatureGroup + isHK,data= mergedData, col=tenCols, ylab='Nucleosome Stability', main='Number of Exons\n',xaxt='n')
#boxplot(log10(Acceptor) ~ FeatureGroup + isHK,data= mergedData, col=tenCols, ylab='Acceptor Strength (log10)', main='Number of Exons\n',xaxt='n')
#dev.off()

png(sprintf('%s_splitsHK_featureCount.png',settings$CommonName),height=1200,width=400)
par(mfrow=c(3,1))
plotSplitsHK(mergedData,'length','Last Exon Length',Group='FeatureGroup')
plotSplitsHK(mergedData,'Energy','Nucleosome Stability',Group='FeatureGroup')
plotSplitsHK(mergedData,'Acceptor',"3' splice site strength",Group='FeatureGroup')
dev.off()

## Same but split up into HK, non-HK
png(sprintf('%s_splitsHK.png',settings$CommonName),height=600,width=200)
par(mfrow=c(3,1))
plotSplitsHK(mergedData,'Energy','Nucleosome Stability')
plotSplitsHK(mergedData,'Donor',"5' splice site strength")
plotSplitsHK(mergedData,'Acceptor',"3' splice site strength")
dev.off()

### Get rid of donor, add length boxplot
#png(sprintf('%s_splitsHK2.png',settings$CommonName),height=900,width=300)
#par(mfrow=c(3,1),cex.axis=1.1,cex=1.1,cex.main=1.1)
#boxplot((length) ~ isHK, data=mergedData, col = c(rgb(.5,.5,1),rgb(1,.5,.5)),ylab="(bp)",names=c("Inducible","HK"),main="Last Exon Length")
#plotSplitsHK(mergedData,'Energy','Nucleosome Stability',ylab="Nucleosome Stability Score",xlab='Length Bins')
#plotSplitsHK(mergedData,'Acceptor',"3' splice site strength",ylab="PWM score" ,xlab='Length Bins')
#plotSplitsHK(mergedData,'isHK',"HK gene status")
#dev.off()
