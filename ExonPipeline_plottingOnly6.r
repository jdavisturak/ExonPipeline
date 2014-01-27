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
# 5/2/13: Instead of using arbitrary bins, I am using 5 bins based on 'Splice Groups'
# 
#**********************************************************************************************************************************

rm(list=ls())
library(gplots)
source("/home/jeremy/ExonPipeline/ExonPipeline_functions.r")
source("/home/jeremy/ExonPipeline/ExonPipeline_analysis_and_plotting_Functions.r")
source("/home/home/jeremy/Code/useful_R/useful_R.r")


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

save.image("PlottedData6.RData")
 
#load("PlottedData3.RData")
cat("	Plotting data:\n") 

#### 1) Plot Barplots
source("/home/jeremy/ExonPipeline/ExonPipeline_analysis_and_plotting_Functions.r")
pdf(sprintf('%s_Barplots6.pdf',settings$CommonName),height=9,width=2)
par(mfrow=c(4,1),cex.axis=1.2,cex=0.7,cex.main=1.1,mar=c(2,3,2,1),lwd=4)
  #errorcol=rgb(1/2,1/4,0.26)
  errorcol='black';
  barCols = c('gray','gray')
  makeBarplot(LENGTH,'',doLog=T,col=barCols,addStats=F,yaxt='n',border=NA,errorcol=errorcol,YLIM=c(1,4.3), xpd=F)
  axis(2,at=0:5,lab=10^(c(NA,1:4,NA)))
  abline(h=par('usr')[3], lwd=1)
  makeBarplot(ENERGY,'',doLog=F,col=barCols, YLIM=c(-.25,0.25),xaxt='s',addStats=F,border=NA,errorcol=errorcol)
  makeBarplot(DONOR,'',doLog=F,col=barCols, YLIM=c(-.25,0.25),addStats=F,border=NA,errorcol=errorcol)
  makeBarplot(ACCEPTOR,'',doLog=F,col=barCols, YLIM=c(-.25,0.25),addStats=F,border=NA,errorcol=errorcol)
dev.off()


#### 2) Comparisons WITHIN LASTS (including housekeeping genes)


### Do some plots based on splitting up the data into groups

# Split it up into 5 groups BASED ON 5 QUINTILES OF SPLICING PROBABILITIES, GIVEN OUR PARAMETERS
#
k_splice = 2 # 2/minute
elong = 2010 # bp/minute  % From matlab:  1/((1+KP/KU)/(KE+KP)) = 2010
lengths =  -log(1-(1:4)/5)*elong/k_splice; lengths              
mergedData$LengthGroup_S = splitByVector(mergedData$length,lengths) # split into 5 groups

colors = c('gray', 'black'); # For HK, non
## SPLIT BY LENGTH
pdf(sprintf('%s_simple_splitByLength7.pdf',settings$CommonName),height=10,width=8)
# THis one has 147 bp.  plot 6 uses ALL bp, I think.
par(mfcol=c(2,2),cex.axis=1.2,cex=1.3,cex.main=1.1,las=3,mar=c(6.1 ,4.1 ,4.1, 2.1))
plotSplits(mergedData,'Energy','', Group='LengthGroup_S',  xlab='Last Exon Quintiles',barcol=colors[2], ylab='Nucleosome Stability Z-score',gap=0,type='l',pch=NA)
plotSplits(mergedData,'Acceptor','',Group='LengthGroup_S', xlab='Last Exon Quintiles',barcol=colors[2], ylab="3' Splice Site Strength Score",gap=0,type='l',ylim=c(0.1,0.4),pch=NA)
plotSplitsHK(mergedData,'Energy',"", Group='LengthGroup_S',  xlab='Last Exon Group',ylab='Nucleosome Stability Z-score',colors=colors)
plotSplitsHK(mergedData,'Acceptor',"", Group='LengthGroup_S',  xlab='Last Exon Group', ylab="3' Splice Site Strength Score",colors=colors)

dev.off()

         
### Do with medians
pdf(sprintf('%s_simple_splitByLength6_Median.pdf',settings$CommonName),height=10,width=8)
par(mfcol=c(2,2),cex.axis=1.2,cex=1.3,cex.main=1.1,las=3,mar=c(6.1 ,4.1 ,4.1, 2.1))
plotSplits(mergedData,'Energy','Median',avg='median', Group='LengthGroup_S',  xlab='Last Exon Quintiles',barcol='grey', ylab='Nucleosome Stability Z-score',gap=0,type='l',pch=NA)
plotSplits(mergedData,'Acceptor','Median',avg='median',Group='LengthGroup_S', xlab='Last Exon Quintiles',barcol='grey', ylab="3' Splice Site Strength Score",gap=0,type='l',ylim=c(0.1,0.4),pch=NA)
plotSplitsHK(mergedData,'Energy','Median',avg='median',ylab='Nucleosome Stability Z-score',xlab='Last Exon Group')
plotSplitsHK(mergedData,'Acceptor',"Median",avg='median', ylab="3' Splice Site Strength Score",xlab='Last Exon Group')

dev.off()

### OUTPUT enough data to simulate in MATALB
lengths_last =  tapply(mergedData$length,mergedData$LengthGroup_S,mean)
strength_nonLast = mean(ACCEPTOR[ACCEPTOR[,2]==0,1],na.rm=T)
strength_lasts_avg = tapply(mergedData$Acceptor,mergedData$LengthGroup_S,mean)
energy_lasts_avg = tapply(mergedData$Energy,mergedData$LengthGroup_S,mean)

write.delim(cbind(lengths_last,strength_lasts_avg,strength_nonLast,energy_lasts_avg), 'Normalized_averages_byLength_6.txt') 

