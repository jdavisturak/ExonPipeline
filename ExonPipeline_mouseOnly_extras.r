
########## Expected splicing ~ [Expected elongation length] * Rate, where Rate ~ 1/Splicing_Time
# Use 'expected time', which tracks ALL exons, using a conversion from splice site strength to kinetic rate
#dataWithChromGROTime2$expectedSplicing =   (dataWithChromGROTime2$length +  Energy_BP) #/ dataWithChromGROTime2$MinutesToSplice
#dataWithChromGROTime2$expectedSplicing2 =   (dataWithChromGROTime2$length + dataWithChromGROTime2$Runon + Energy_BP) #/ dataWithChromGROTime2$MinutesToSplice
#dataWithChromGROTime2$expectedSplicingSimple =  (dataWithChromGROTime2$Acceptor) * (dataWithChromGROTime2$length + Energy_BP) 
#dataWithChromGROTime2$expectedSplicingSimple2 =  (dataWithChromGROTime2$Acceptor) * (dataWithChromGROTime2$length + dataWithChromGROTime2$Runon + Energy_BP) 
## Attempts to Normalize, get rid of NAs to make them more comparable.
#dataWithChromGROTime2$expectedSplicingSimple[is.na(dataWithChromGROTime2$expectedSplicingSimple2)] <- NA
#dataWithChromGROTime2$expectedSplicingSimple.Z = normalize(dataWithChromGROTime2$expectedSplicingSimple)
#dataWithChromGROTime2$expectedSplicingSimple2.Z = normalize(dataWithChromGROTime2$expectedSplicingSimple2)
#dataWithChromGROTime2$expectedSplicingSimpleGroup = splitByVector(dataWithChromGROTime2$expectedSplicingSimple, quantile(dataWithChromGROTime2$expectedSplicingSimple,na.rm=T,(1:4)/5)) # split into 5 groups
#dataWithChromGROTime2$expectedSplicingSimpleGroup2 = splitByVector(dataWithChromGROTime2$expectedSplicingSimple2, quantile(dataWithChromGROTime2$expectedSplicingSimple2,na.rm=T,(1:4)/5)) # split into 5 groups

### Use Final length AFTER intron, (plus add a little bit for pausing)
#dataWithChromGROTime2$expectedSplicing =   (dataWithChromGROTime2$length +  Energy_BP)/3000
#dataWithChromGROTime2$expectedSplicing2 =   (dataWithChromGROTime2$LengthAfterIntron +  Energy_BP)/3000 
#NAME1 = "Approximate Time available to Splice (min)\n(Length , 6sec/EnergyZ)"
#NAME2 = "Approximate Time available to Splice (min)\n(Length + Runon, 6sec/EnergyZ)"



## Merge and plot
#png('modelSplicingRunon3.png')
#plot(dat[,2],dat[,3],ylab='WithRunon', xlab='No Runon', main='Predicted CTS efficiency')
#dev.off()
#######################################################################################################################################################
### THis code snippet will only work for the original version of the matlab script, where I just compared 2 things
#######################################################################################################################################################
# dat$GRO_improvement = dat$PredSplice2 - dat$PredSplice1
#png(sprintf('%s_simple7_SplitRunonImprovementbyOthers3.png',settings$CommonName),height=600,width=1200)
#par(mfrow=c(1,2),cex.axis=1.2,cex=1.3,cex.main=1.1,las=3,mar=c(6.1 ,4.1 ,4.1, 2.1)) ; options(warn=-1) -> W
#
#xlabs = tapply(mergedData$length,mergedData$Group,function(x)paste(range(x),collapse='-'))
#plotSplitsHK(dataWithChromGRO,'Runon','GRO-seq Read-Through Length\nAs a function of last exon length',Group='LengthGroup',xlab='',names=xlabs,ylab='')
#
#xlabs = tapply(mergedData$length,mergedData$Group,function(x)paste(range(x),collapse='-'))
#plotSplitsHK(dataWithChromGROpredictions,'PredSplice1','Splicing\nAs a function of last exon length',Group='Group',xlab='',ylab='',names=xlabs,YLIM=c(0.0,0.4))
#plotSplitsHK(dataWithChromGROpredictions,'PredSplice2','',Group='Group',add=T,LTY=2,YLIM=c(0.0,0.4),xlab='Last Exon Length',ylbias=-1)
#
#options(W); dev.off()
##png(sprintf('%s_simple4.2_SplitRunonImprovementbyOthers3.png',settings$CommonName),height=800,width=1200)
#par(mfrow=c(2,3),cex.axis=1.2,cex=1.3,cex.main=1.1,las=3,mar=c(6.1 ,4.1 ,4.1, 2.1)) ; options(warn=-1) -> W
#
#xlabs = tapply(mergedData$length,mergedData$Group,function(x)paste(range(x),collapse='-'))
#
#plotSplitsHK(dataWithChromGROpredictions,'PredSplice1','Splicing \nAs a function of last exon length',Group='Group',xlab='',names=xlabs,ylab='',YLIM=c(0.2,0.9))
#plotSplits(dataWithChromGROpredictions,'PredSplice2','',Group='Group',add=T,type='b',lwd=2)
#plotSplitsHK(dataWithChromGROpredictions,'PredSplice2','Splicing with Read-through \nAs a function of last exon length',Group='Group',xlab='',names=xlabs,ylab='',YLIM=c(0.2,0.9))
#plotSplits(dataWithChromGROpredictions,'PredSplice2','',Group='Group',add=T,type='b',lwd=2)
#
#plotSplitsHK(dataWithChromGROpredictions,'PredSplice1','Splicing\nAs a function of last exon length',Group='Group',ylab='',xlab='Last Exon Length',names=xlabs,YLIM=c(0.2,0.8))
#plotSplitsHK(dataWithChromGROpredictions,'PredSplice2','',Group='Group',add=T,LTY=2,YLIM=c(0.2,0.8))
#
#xlabs = tapply(mergedData$length,mergedData$LengthGroup,function(x)paste(range(x),collapse='-'))
#plotSplitsHK(dataWithChromGROpredictions,'PredSplice1','Splicing\nAs a function of LengthGroup',Group='LengthGroup',ylab='',xlab='',names=xlabs,YLIM=c(0.2,0.9))
#plotSplits(dataWithChromGROpredictions,'PredSplice1','',Group='LengthGroup',add=T,type='b',lwd=2)
#plotSplitsHK(dataWithChromGROpredictions,'PredSplice2','Splicing with Read-through \nAs a function of LengthGroup',Group='LengthGroup',ylab='',xlab='',names=xlabs,YLIM=c(0.2,0.9))
#plotSplits(dataWithChromGROpredictions,'PredSplice2','',Group='LengthGroup',add=T,type='b',lwd=2)
#
#plotSplitsHK(dataWithChromGROpredictions,'PredSplice1','Splicing\nAs a function of LengthGroup',Group='LengthGroup',ylab='',xlab='Last Exon Length',names=xlabs,YLIM=c(0.3,0.8))
#plotSplitsHK(dataWithChromGROpredictions,'PredSplice2','',Group='LengthGroup',add=T,LTY=2,YLIM=c(0.3,0.8))
#
#options(W); dev.off()

### Try doing GRO while ALLOWING 0's:#
#
###################################################################################################################
####################  MOUSE GRO-SEQ ###########################
###################################################################################################################
#
### Read in GRO data
#GRO = read.delim("/home/home/jeremy/RNAseq/Glass/refseq_and_subsequent_transcripts_with_tag_counts.txt",stringsAsFactors=F)
#
### Calculate the Read-Through (hereafter referred to as 'Runon') and ***set genes with no measurable read-through to 0***
#NAs = is.na(GRO$post_transcript_start)
#GRO$post_transcript_start[NAs & GRO$strand==0] <- GRO$gene_end + 1
#GRO$post_transcript_start[NAs & GRO$strand==1] <- GRO$gene_start - 1
#GRO$post_transcript_end[NAs] <- GRO$post_transcript_start[NAs]
#
#GRO$Runon = abs(GRO$post_transcript_end - GRO$post_transcript_start)
#
### Annotate by merging with refGene
#refGene =add_UniqueID(loadRefgene(settings))
#GRO2 = merge(GRO,refGene[,c("name","UniqueID","strand")],by=1)
#
#options(scipen=10) # make it so I don't write scientific numbers here
### Omit any GRO-seq data when the Read-Through ends within 1Kb of a gene.
#write.table(cbind(GRO2$chr, GRO2$post_transcript_start-1, GRO2$post_transcript_end+1, GRO2$UniqueID, 0, GRO2$strand.y),file='GRO2.2.bed',sep="\t",row.names=F,col.names=F, quote=F)
#system("bedtools intersect -a GRO2.2.bed -b refseq.bed -v -s -wa > GRO_no_TSS_interuptions2.bed")
#options(scipen=0)
#
#good_GRO =read.delim("GRO_no_TSS_interuptions2.bed", head=F)
#GRO3 <- GRO2[GRO2$UniqueID %in% good_GRO[,4],]
#save(GRO3,file="GRO_with0_Data.RData")
#
### Merge with 'mergedData'
#dataWithGRO2 = merge(mergedData, GRO3 , by='UniqueID',all.x=T)
#dataWithGRO2$RunonGroup = splitByVector(dataWithGRO2$Runon,quantile(dataWithGRO2$Runon,na.rm=T,(1:4)/5)) # split into 5 groups
#save(dataWithGRO2,file='dataWithGRO_with0.RData')
#
###################################################################################################################
#################### Re-plotting the first figure at least.  I should probably redo the simulation then, if it looks worthwhile!
#png(sprintf('%s_simple6_SplitRunonWith0ImprovementbyOthers3.png',settings$CommonName),height=600,width=1200)
#par(mfrow=c(1,2),cex.axis=1.2,cex=1.3,cex.main=1.1,las=3,mar=c(6.1 ,4.1 ,4.1, 2.1)) ; options(warn=-1) -> W
#
#xlabs = tapply(mergedData$length,mergedData$Group,function(x)paste(range(x),collapse='-'))
#plotSplitsHK(dataWithGRO2,'Runon','Read-Through Length\nAs a function of last exon length',Group='LengthGroup',xlab='',names=xlabs,ylab='')
#
#xlabs = tapply(mergedData$length,mergedData$Group,function(x)paste(range(x),collapse='-'))
#plotSplitsHK(dataWithChromGROpredictions,'PredSplice1','Splicing\nAs a function of last exon length',Group='Group',xlab='',ylab='',names=xlabs,YLIM=c(0.2,0.8))
#plotSplitsHK(dataWithChromGROpredictions,'PredSplice2','',Group='Group',add=T,LTY=2,YLIM=c(0.2,0.8),xlab='Last Exon Length',ylbias=-1)
#
#options(W); dev.off()



#plotSplitsHK(dataWithChromGROTime2,'Runon',paste('Read-through distance (bp)','As a function of Estimated time need to splice',sep="\n"),Group='TimeGroup',xlab='',ylab='',names=xlabs)
#plotSplits(dataWithChromGROTime2,'Runon','',Group='TimeGroup',add=T,type='b',lwd=2)
#predictSimple1 = plotSplitsHK(dataWithChromGROTime2,'expectedSplicing',paste(NAME1,'As a function of Acceptor Strength of Last Intron',sep="\n"),Group='AcceptorGroup',xlab='',ylab='',names=xlabs)
#plotSplits(dataWithChromGROTime2,'expectedSplicing','',Group='AcceptorGroup',add=T,type='b',lwd=2)
#xlabs = tapply(dataWithChromGROTime2$MinutesToSplice,dataWithChromGROTime2$TimeGroup2,function(x)paste(sprintf("%.1f",range(x)),collapse='-'))
#predictSimple1 = plotSplitsHK(dataWithChromGROTime2,'expectedSplicing',paste(NAME1,'As a function of Estimated time needed to splice',sep="\n"),Group='TimeGroup2',xlab='',ylab='',names=xlabs)
#plotSplits(dataWithChromGROTime2,'expectedSplicing','',Group='TimeGroup2',add=T,type='b',lwd=2)

 
#xlabs = tapply(dataWithChromGROTime2$PooledSplicing,dataWithChromGROTime2$SpliceAllGroup,function(x)paste(sprintf("%.2f",range(x)),collapse='-'))
#predictSimple1 = plotSplitsHK(dataWithChromGROTime2,'expectedSplicing',paste(NAME1,'As a function of Measured Splicing of Last Intron',sep="\n"),Group='SpliceAllGroup',xlab='',ylab='',names=xlabs)
#plotSplits(dataWithChromGROTime2,'expectedSplicing','',Group='SpliceAllGroup',add=T,type='b',lwd=2)
#
#predictSimple2 = plotSplitsHK(dataWithChromGROTime2,'expectedSplicing2',paste(NAME2,'As a function of Measured Splicing of Last Intron',sep="\n"),Group='SpliceAllGroup',xlab='',ylab='',names=xlabs)
#plotSplits(dataWithChromGROTime2,'expectedSplicing2','',Group='SpliceAllGroup',add=T,type='b',lwd=2)

#predictSimple1 = plotSplitsHK(dataWithChromGROTime2,'expectedSplicing',paste(NAME1,'As a function of Estimated time need to splice',sep="\n"),Group='TimeGroup',xlab='',ylab='',names=xlabs)
#plotSplits(dataWithChromGROTime2,'expectedSplicing','',Group='TimeGroup',add=T,type='b',lwd=2)
#
#predictSimple2 = plotSplitsHK(dataWithChromGROTime2,'expectedSplicing2',paste(NAME2,'As a function of Estimated time need to splice',sep="\n"),Group='TimeGroup',xlab='',ylab='',names=xlabs)
#plotSplits(dataWithChromGROTime2,'expectedSplicing2','',Group='TimeGroup',add=T,type='b',lwd=2)

#plotSplitsHK(dataWithChromGROTime2,'PooledSplicing','Measured Splicing of Last Intron\nAs a function of expected Splicing group \n(Calculated using length of last exon)',Group='expectedSplicingSimpleGroup',xlab='',ylab='', YLIM=c(0.7,1))
#plotSplits(dataWithChromGROTime2,'PooledSplicing','',Group='expectedSplicingSimpleGroup',add=T,type='b',lwd=2)
#
#plotSplitsHK(dataWithChromGROTime2,'PooledSplicing','Measured Splicing of Last Intron\nAs a function of expected Splicing group \n(Calculated using last exon + Read-through length)',Group='expectedSplicingSimpleGroup2',xlab='',ylab='', YLIM=c(0.7,1))
#plotSplits(dataWithChromGROTime2,'PooledSplicing','',Group='expectedSplicingSimpleGroup2',add=T,type='b',lwd=2)

#predictSimple1 = plotSplitsHK(dataWithChromGROTime2,'expectedSplicingSimple.Z','Expected Splicing Z score \n(Calculated using length of last exon)\nAs a function of Measured Splicing of Last Intron',Group='SpliceAllGroup',xlab='',ylab='')
#plotSplits(dataWithChromGROTime2,'expectedSplicingSimple.Z','',Group='SpliceAllGroup',add=T,type='b',lwd=2)

#predictSimple2 = plotSplitsHK(dataWithChromGROTime2,'expectedSplicingSimple2.Z','Expected Splicing Z score\n(Calculated using last exon + Read-through length)\nAs a function of Measured Splicing of Last Intron',Group='SpliceAllGroup',xlab='',ylab='')
#plotSplits(dataWithChromGROTime2,'expectedSplicingSimple2.Z','',Group='SpliceAllGroup',add=T,type='b',lwd=2)
# Make bar plot 
#dat1 = as.matrix(rbind(HK_Low = predictSimple1[[1]][1,], HK_High=predictSimple1[[1]][5,], nonHK_Low = predictSimple1[[2]][1,], nonHK_High=predictSimple1[[2]][5,]))
#dat2 = as.matrix(rbind(HK_Low = predictSimple2[[1]][1,], HK_High=predictSimple2[[1]][5,], nonHK_Low = predictSimple2[[2]][1,], nonHK_High=predictSimple2[[2]][5,]))
#barplot(dat1[,1], names=rownames(dat1), main='Version 1',ylab=c("Average Predicted Splicing Z Score")) #,ylim=c(-.4,.4)
#barplot(dat2[,1], names=rownames(dat2),  main='Version 2',ylab=c("Average Predicted Splicing Z Score"))
#plotSplitsHK(dataWithChromGROTime2,'PooledSplicing','Measured Splicing of Last Intron\nAs a function of expected Splicing group\n(All introns, splicing Factor=5)',Group='expectedSplicingGroup',xlab='',ylab='')
#plotSplits(dataWithChromGROTime2,'PooledSplicing','',Group='expectedSplicingGroup',add=T,type='b',lwd=2)
#plotSplitsHK(dataWithChromGROTime2,'PooledSplicing','Measured Splicing of Last Intron\nAs a function of expected Splicing group\n(All introns, splicing Factor=2)',Group='expectedSplicingGroup2',xlab='',ylab='')
#plotSplits(dataWithChromGROTime2,'PooledSplicing','',Group='expectedSplicingGroup2',add=T,type='b',lwd=2)
#plotSplitsHK(dataWithChromGROTime2,'Runon','Length of Read-through\nAs a function of expected Splicing group\n(All introns, splicing Factor=5)',Group='expectedSplicingGroup',xlab='',ylab='')
#plotSplits(dataWithChromGROTime2,'Runon','',Group='expectedSplicingGroup',add=T,type='b',lwd=2)
#plotSplitsHK(dataWithChromGROTime2,'Runon','Length of Read-through\nAs a function of expected Splicing group\n(All introns, splicing Factor=2)',Group='expectedSplicingGroup2',xlab='',ylab='')
#plotSplits(dataWithChromGROTime2,'Runon','',Group='expectedSplicingGroup2',add=T,type='b',lwd=2)






### Split up based on splicing rates at 0 min:
#dataWithChromGRO$SpliceGroup0 = splitByVector(dataWithChromGRO$Chromatin_0_min,quantile(dataWithChromGRO$Chromatin_0_min,na.rm=T,(1:4)/5)) # split into 5 groups
#dataWithChromGRO$SpliceGroup60 = splitByVector(dataWithChromGRO$Chromatin_60_min,quantile(dataWithChromGRO$Chromatin_60_min,na.rm=T,(1:4)/5)) # split into 5 groups

### Split up based on Run-on length
#dataWithChromGRO$RunonGroup = splitByVector(dataWithChromGRO$Runon,quantile(dataWithChromGRO$Runon,na.rm=T,(1:4)/5)) # split into 5 groups
#dataWithGRO$FeatureGroup = splitByVector(dataWithGRO$FeatureCount,quantile(dataWithGRO$FeatureCount,na.rm=T,(1:4)/5)) # split into 5 groups
#dataWithChromGRO$SplicingGroup = splitByVector(dataWithChromGRO$Acceptor,quantile(dataWithChromGRO$Acceptor,na.rm=T,(1:4)/5)) # split into 5 groups
#dataWithChromGRO$InducedSplicing_30min = dataWithChromGRO$Chromatin_30_min - dataWithChromGRO$Chromatin_0_min
#dataWithChromGRO$GroupInduced30 = splitByVector(dataWithChromGRO$InducedSplicing_30min,quantile(dataWithChromGRO$InducedSplicing_30min,na.rm=T,(1:4)/5)) # split into 5 groups




#### CHROM_ASSOC:

## First just overall distribution
chromBoxplot = function(outData,main=''){
  par(mar=c(5,10,3,1))
  palette(colorRampPalette(c('white','green'))(8))
  boxplot(lapply(outData,function(x)x$PercentSpliced),col=2:6,notch=T, las=1,names=samples, horizontal=T,xlab="Fraction of last intron spliced", main=paste(main) )
}

png("Chrom1.png",width=500,height=300)
chromBoxplot(outData)
dev.off()
             
png("ChromFiltered.png",width=500,height=300)
chromBoxplot(lapply(outData, function(x) x[x$SplicedReads + x$UnsplicedReads > 5,]))
dev.off()

png("ChromFiltered2nds.png",width=500,height=300)
chromBoxplot(lapply(outData2, function(x) x[x$SplicedReads + x$UnsplicedReads > 5,]))
dev.off()

png("ChromAll.png",width=500,height=300)
par(mfrow=c(2,2))
boxplot(combinedChrom$PercentSpliced[combinedChrom$SplicedReads + combinedChrom$UnsplicedReads > 0],main='All genes',ylab='% spliced') #, function(x) x[x$SplicedReads + x$SplicedReads > 5,]
boxplot(combinedChrom$PercentSpliced[combinedChrom$SplicedReads + combinedChrom$UnsplicedReads > 10],main='Genes at > 10 reads',ylab='% spliced') #, function(x) x[x$SplicedReads + x$SplicedReads > 5,]
dev.off()

png(sprintf('%s_CHROMsplits.png',settings$CommonName))
par(mfrow=c(3,2))
plotSplitsHK(dataWithChrom,'Chromatin_0_min','Chrom-Splicing @ 0 min')
plotSplitsHK(dataWithChrom,'Chromatin_15_min',"Chrom-Splicing @ 15 min")
plotSplitsHK(dataWithChrom,'Chromatin_30_min',"Chrom-Splicing @ 30 min")
plotSplitsHK(dataWithChrom,'Chromatin_60_min',"Chrom-Splicing @ 60 min")
plotSplitsHK(dataWithChrom,'Chromatin_120_min',"Chrom-Splicing @ 120 min")
dev.off()



#### Try to find genes that change their splicing patterns
source("/home/jeremy/useful_R.r")
png("Chrom_Heatmap1.png")
D = data.matrix(dataWithChrom[dataWithChrom$isHK==0,15:18]) # Get rid of housekeeping genes
D = D[apply(D,1,function(x)!any(is.na(x))),]  # Get rid of any rows that have at least one NA
D = D[rowMins(D) < 0.9,]  # Get rid of genes that always have 100% splicing
heatmap1(D,Rowv=TRUE,main="Splicing of last intron in Chromatin-associated RNA")
dev.off()



## Split splicing data into 'Extreme' Categories!!! 

pdf("Mouse_Chrom_Boxplots1.pdf",height=3,width=6)
par(mfrow=c(1,3),las=2)
boxplot(Chromatin_30_min ~ AcceptorGroup + isHK,data= dataWithChrom, col=tenCols, main='Acceptor Strength\n',xaxt='n')
boxplot(Chromatin_30_min ~ EnergyGroup + isHK,data= dataWithChrom, col=tenCols, main='Nucleosome Stability\n',xaxt='n')
boxplot(Chromatin_30_min ~ LengthGroup + isHK,data= dataWithChrom, col=tenCols, main='Exon Length\n',xaxt='n')
dev.off()

### Chrom 30 min versus Feature Counts:
dataWithChrom$FeatureGroup = splitByVector(dataWithChrom$FeatureCount,quantile(dataWithChrom$FeatureCount,na.rm=T,(1:4)/5)) # split into 5 groups

pdf("Mouse_Chrom_Boxplots2.pdf",height=3,width=6)
par(mfrow=c(1,3),las=2)
boxplot(Chromatin_30_min ~ FeatureGroup + isHK,data= dataWithChrom, col=tenCols, main='Number of Introns\n',xaxt='n')
#boxplot(Chromatin_30_min ~ EnergyGroup + isHK,data= dataWithChrom, col=tenCols, main='Nucleosome Stability\n',xaxt='n')
#boxplot(Chromatin_30_min ~ LengthGroup + isHK,data= dataWithChrom, col=tenCols, main='Exon Length\n',xaxt='n')
dev.off()

### Extreme doesn't work that well because low sample size
#dataWithChrom.extreme = dataWithChrom[dataWithChrom$AcceptorGroup %in% c(1,5) & dataWithChrom$EnergyGroup %in% c(1,5) & dataWithChrom$LengthGroup %in% c(1,5)  ,]
#png("Mouse_Chrom_extreme2.png",height=800,width=566)
#par(mfrow=c(2,2),las=2)
#boxplot(Chromatin_30_min ~ AcceptorGroup + isHK,data= dataWithChrom[dataWithChrom$LengthGroup < 2, ], col=tenCols, main='Short Last exons: \nIncreasing Acceptor Strength')
#boxplot(Chromatin_30_min ~ EnergyGroup + isHK,data= dataWithChrom[dataWithChrom$LengthGroup < 2, ], col=tenCols, main='Short Last exons: \nIncreasing Nuclseosome Stability')
#boxplot(Chromatin_30_min ~ AcceptorGroup + isHK,data= dataWithChrom[dataWithChrom$LengthGroup > 4, ], col=tenCols, main='Long Last exons: \nIncreasing Acceptor Strength')
#boxplot(Chromatin_30_min ~ EnergyGroup + isHK,data= dataWithChrom[dataWithChrom$LengthGroup > 4, ], col=tenCols, main='Long Last exons: \nIncreasing Nuclseosome Stability')
#dev.off()




#### Do some genes CHANGE their splicing upon stimulation?  ####

png(sprintf('%s_ChromChanges.png',settings$CommonName))
par(mfrow=c(2,2))
boxplot(Chromatin_15_min - Chromatin_0_min ~ isHK,  data=dataWithChrom,col=tenCols[5:6],main='Change in Splicing Rate: 15 min - basal',ylab='Change in % spliced', notch=T)
boxplot(Chromatin_30_min - Chromatin_0_min ~ isHK,  data=dataWithChrom,col=tenCols[5:6],main='Change in Splicing Rate: 30 min - basal',ylab='Change in % spliced', notch=T)
boxplot(Chromatin_60_min - Chromatin_0_min ~ isHK,  data=dataWithChrom,col=tenCols[5:6],main='Change in Splicing Rate: 60 min - basal',ylab='Change in % spliced', notch=T)
boxplot(Chromatin_120_min - Chromatin_0_min ~ isHK,  data=dataWithChrom,col=tenCols[5:6],main='Change in Splicing Rate: 120 min - basal',ylab='Change in % spliced', notch=T)
dev.off()
 
cutoff = 0.8
frac = function(x)length(which(x))/length(x)
#F=sapply(15:19, function(i) 	c(frac(dataWithChrom[dataWithChrom$isHK==0,i] > cutoff), frac(dataWithChrom[dataWithChrom$isHK==1,i] > cutoff)))



png(sprintf('%s_GROsplits.png',settings$CommonName))
par(mfrow=c(2,2))
plotSplitsHK(mergedData,'Energy','Nucleosome Stability')
plotSplitsHK(mergedData,'Donor',"5' splice site strength")
plotSplitsHK(mergedData,'Acceptor',"3' splice site strength")
plotSplitsHK(dataWithGRO,'Runon',"Run-on Length")
dev.off()



 
png(sprintf('%s_splits2_by_0min.png',settings$CommonName))
par(mfrow=c(2,3))
plotSplitsHK(dataWithChromGRO,'length','Length',Group='SpliceGroup0')
plotSplitsHK(dataWithChromGRO,'Energy','Nucleosome Stability',Group='SpliceGroup0')
plotSplitsHK(dataWithChromGRO,'Donor',"5' splice site strength",Group='SpliceGroup0')
plotSplitsHK(dataWithChromGRO,'Acceptor',"3' splice site strength",Group='SpliceGroup0')
plotSplitsHK(dataWithChromGRO,'Runon',"Run-on Length",Group='SpliceGroup0')
dev.off()
png(sprintf('%s_splits2_by_60min.png',settings$CommonName))
par(mfrow=c(2,3))
plotSplitsHK(dataWithChromGRO,'length','Length',Group='SpliceGroup60')
plotSplitsHK(dataWithChromGRO,'Energy','Nucleosome Stability',Group='SpliceGroup60')
plotSplitsHK(dataWithChromGRO,'Donor',"5' splice site strength",Group='SpliceGroup60')
plotSplitsHK(dataWithChromGRO,'Acceptor',"3' splice site strength",Group='SpliceGroup60')
plotSplitsHK(dataWithChromGRO,'Runon',"Run-on Length",Group='SpliceGroup60')
dev.off()




png(sprintf('%s_splits2_by_Runon.png',settings$CommonName),height=800,width=800)
par(mfrow=c(3,3))
plotSplitsHK(dataWithChromGRO,'length','Length',Group='RunonGroup')
plotSplitsHK(dataWithChromGRO,'Energy','Nucleosome Stability',Group='RunonGroup')
plotSplitsHK(dataWithChromGRO,'Donor',"5' splice site strength",Group='RunonGroup')
plotSplitsHK(dataWithChromGRO,'Acceptor',"3' splice site strength",Group='RunonGroup')
plotSplitsHK(dataWithChromGRO,'Chromatin_0_min',"Splicing 0 min",Group='RunonGroup')
plotSplitsHK(dataWithChromGRO,'Chromatin_15_min',"Splicing 15 min",Group='RunonGroup')
plotSplitsHK(dataWithChromGRO,'Chromatin_30_min',"Splicing 30 min",Group='RunonGroup')
plotSplitsHK(dataWithChromGRO,'Chromatin_60_min',"Splicing 60 min",Group='RunonGroup')
plotSplitsHK(dataWithChromGRO,'Chromatin_120_min',"Splicing 120 min",Group='RunonGroup')
dev.off() 

png(sprintf('%s_splits3_by_Runon.png',settings$CommonName),height=800,width=800)
par(mfrow=c(3,3))
plotSplitsHK(dataWithChromGRO,'Chromatin_0_min',"Splicing 0 min",Group='RunonGroup')
plotSplitsHK(dataWithChromGRO,'Chromatin_15_min',"Splicing 15 min",Group='RunonGroup')
plotSplitsHK(dataWithChromGRO,'Chromatin_30_min',"Splicing 30 min",Group='RunonGroup')
plotSplitsHK(dataWithChromGRO,'Chromatin_60_min',"Splicing 30 min",Group='RunonGroup')
plotSplitsHK(dataWithGRO,'length','Length',Group='RunonGroup')
plotSplitsHK(dataWithGRO,'Energy','Nucleosome Stability',Group='RunonGroup')
plotSplitsHK(dataWithGRO,'Acceptor',"3' splice site strength",Group='RunonGroup')
plotSplitsHK(dataWithGRO,'FeatureCount',"Number of Introns",Group='RunonGroup')
plotSplitsHK(dataWithGRO,'Runon',"Runon-Length vs number of introns",Group='FeatureGroup')
dev.off() 


png(sprintf('%s_splits4_by_Runon.png',settings$CommonName),height=800,width=800)
par(mfrow=c(2,2))
boxplot(Runon ~ isHK, data=dataWithGRO,col=tenCols[5:6],main='Runon Length',ylab='bp', notch=T)
boxplot(Chromatin_0_min ~ isHK, data=dataWithChrom,col=tenCols[5:6],main='Splicing Rate - 0 min',ylab='% spliced', notch=T)
boxplot(Chromatin_15_min ~ isHK, data=dataWithChrom,col=tenCols[5:6],main='Splicing Rate - 15 min',ylab='% spliced', notch=T)
boxplot(Chromatin_30_min ~ isHK, data=dataWithChrom,col=tenCols[5:6],main='Splicing Rate - 30 min',ylab='% spliced', notch=T)
#plotSplitsHK(dataWithGRO,'length','Length',Group='RunonGroup')
dev.off() 


png(sprintf('%s_splits2_by_Length.png',settings$CommonName),height=600,width=600)
par(mfrow=c(3,3))
plotSplitsHK(dataWithChromGRO,'Runon','Run-on Length',Group='Group')
plotSplitsHK(dataWithChromGRO,'Energy','Nucleosome Stability',Group='Group')
plotSplitsHK(dataWithChromGRO,'Donor',"5' splice site strength",Group='Group')
plotSplitsHK(dataWithChromGRO,'Acceptor',"3' splice site strength",Group='Group')
plotSplitsHK(dataWithChromGRO,'Chromatin_0_min',"Splicing 0 min",Group='Group')
plotSplitsHK(dataWithChromGRO,'Chromatin_15_min',"Splicing 15 min",Group='Group')
plotSplitsHK(dataWithChromGRO,'Chromatin_30_min',"Splicing 30 min",Group='Group')
plotSplitsHK(dataWithChromGRO,'Chromatin_60_min',"Splicing 60 min",Group='Group')
plotSplitsHK(dataWithChromGRO,'Chromatin_120_min',"Splicing 120 min",Group='Group')
dev.off() 

### Split up based on 3' splice site strength
png(sprintf('%s_splits2_by_Splicing3.png',settings$CommonName),height=800,width=800)
par(mfrow=c(3,3))
plotSplitsHK(dataWithChromGRO,'length','Length',Group='SplicingGroup')
plotSplitsHK(dataWithChromGRO,'Energy','Nucleosome Stability',Group='SplicingGroup')
plotSplitsHK(dataWithChromGRO,'Donor',"5' splice site strength",Group='SplicingGroup')
plotSplitsHK(dataWithChromGRO,'Runon',"Runon Length",Group='SplicingGroup')
plotSplitsHK(dataWithChromGRO,'Chromatin_0_min',"Splicing 0 min",Group='SplicingGroup')
plotSplitsHK(dataWithChromGRO,'Chromatin_15_min',"Splicing 15 min",Group='SplicingGroup')
plotSplitsHK(dataWithChromGRO,'Chromatin_30_min',"Splicing 30 min",Group='SplicingGroup')
plotSplitsHK(dataWithChromGRO,'Chromatin_60_min',"Splicing 60 min",Group='SplicingGroup')
plotSplitsHK(dataWithChromGRO,'Chromatin_120_min',"Splicing 120 min",Group='SplicingGroup')
dev.off() 


### Split up based on Induced GRO (30 min)
png(sprintf('%s_splits2_by_30minInduction.png',settings$CommonName),height=800,width=800)
par(mfrow=c(2,3))
plotSplitsHK(dataWithChromGRO,'length','Length',Group='GroupInduced30')
plotSplitsHK(dataWithChromGRO,'Energy','Nucleosome Stability',Group='GroupInduced30')
plotSplitsHK(dataWithChromGRO,'Donor',"5' splice site strength",Group='GroupInduced30')
plotSplitsHK(dataWithChromGRO,'Acceptor',"3' splice site strength",Group='GroupInduced30')
plotSplitsHK(dataWithChromGRO,'Runon',"Runon Length",Group='GroupInduced30')
dev.off()



############# ############# ############# ############# ############# ############# 
############# LINEAR MODELS TO PREDICT SPLICING EFFICIENCY #############################
## LMs for splicing
for (HK in 0:1){
  cat(sprintf("\n\n**************** Showing linear models for HK == %d\n\n",HK))
  print(summary(lm(Chromatin_0_min ~ log10(length) + Energy + log10(Acceptor)  ,  data=dataWithChromGRO[dataWithChromGRO$isHK==HK,])))
  print(summary(lm(Chromatin_15_min ~ log10(length) + Energy + log10(Acceptor)  ,  data=dataWithChromGRO[dataWithChromGRO$isHK==HK,])))
  print(summary(lm(Chromatin_30_min ~ log10(length) + Energy + log10(Acceptor) ,  data=dataWithChromGRO[dataWithChromGRO$isHK==HK,])))
  print(summary(lm(Chromatin_60_min ~ log10(length) + Energy + log10(Acceptor) ,  data=dataWithChromGRO[dataWithChromGRO$isHK==HK,])))
  print(summary(lm(Chromatin_120_min ~ log10(length) + Energy + log10(Acceptor)  ,  data=dataWithChromGRO[dataWithChromGRO$isHK==HK,])))
}

## Lms for INDUCED splicing changes
for (HK in list(0,1,0:1)){
  cat(sprintf("\n\n**************** Showing linear models for HK == %d\n\n",HK))
  print(summary(lm(Chromatin_15_min / Chromatin_0_min ~ log10(length) + Energy + log10(Acceptor)  ,  data=dataWithChromGRO[dataWithChromGRO$isHK %in% HK,])))
  print(summary(lm(Chromatin_30_min / Chromatin_0_min ~ log10(length) + Energy + log10(Acceptor)  ,  data=dataWithChromGRO[dataWithChromGRO$isHK %in% HK,])))
  print(summary(lm(Chromatin_60_min / Chromatin_0_min ~ log10(length) + Energy + log10(Acceptor)  ,  data=dataWithChromGRO[dataWithChromGRO$isHK %in% HK,])))
  print(summary(lm(Chromatin_120_min / Chromatin_0_min ~ log10(length) + Energy + log10(Acceptor)  ,  data=dataWithChromGRO[dataWithChromGRO$isHK %in% HK,])))
}
############# ############# ############# ############# ############# 
 


############# ############# ############# ############# ############# ############# ############# 
###  Try some boxplots
############# ############# ############# ############# ############# ############# ############# 

png(sprintf('%s_Chrom_Boxplot_byLength.png',settings$CommonName),height=800,width=800)
par(mfrow=c(2,3),las=3)
boxplot(Chromatin_0_min ~ Group + isHK,data=dataWithChrom, main="0 min",xlab="Length Groups",ylab="% spliced", col = c(rep(rgb(.5,.5,1),5), rep(rgb(1,.5,.5),5)))
boxplot(Chromatin_15_min ~ Group+ isHK,data=dataWithChrom,main="15 min",xlab="Length Groups",ylab="% spliced", col = c(rep(rgb(.5,.5,1),5), rep(rgb(1,.5,.5),5)))
boxplot(Chromatin_30_min ~ Group+ isHK,data=dataWithChrom, main="30 min",xlab="Length Groups",ylab="% spliced", col = c(rep(rgb(.5,.5,1),5), rep(rgb(1,.5,.5),5)))
boxplot(Chromatin_60_min ~ Group+ isHK,data=dataWithChrom, main="60 min",xlab="Length Groups",ylab="% spliced", col = c(rep(rgb(.5,.5,1),5), rep(rgb(1,.5,.5),5)))
boxplot(Chromatin_120_min ~ Group+ isHK,data=dataWithChrom, main="120 min",xlab="Length Groups",ylab="% spliced", col = c(rep(rgb(.5,.5,1),5), rep(rgb(1,.5,.5),5)))
dev.off()

png(sprintf('%s_Chrom_Boxplot_byRunon.png',settings$CommonName),height=800,width=800)
par(mfrow=c(2,3),las=3)
boxplot(Chromatin_0_min ~ RunonGroup + isHK,data=dataWithChromGRO, main="0 min",xlab="Length Groups",ylab="% spliced", col = c(rep(rgb(.5,.5,1),5), rep(rgb(1,.5,.5),5)))
boxplot(Chromatin_15_min ~ RunonGroup+ isHK,data=dataWithChromGRO,main="15 min",xlab="Length Groups",ylab="% spliced", col = c(rep(rgb(.5,.5,1),5), rep(rgb(1,.5,.5),5)))
boxplot(Chromatin_30_min ~ RunonGroup+ isHK,data=dataWithChromGRO, main="30 min",xlab="Length Groups",ylab="% spliced", col = c(rep(rgb(.5,.5,1),5), rep(rgb(1,.5,.5),5)))
boxplot(Chromatin_60_min ~ RunonGroup+ isHK,data=dataWithChromGRO, main="60 min",xlab="Length Groups",ylab="% spliced", col = c(rep(rgb(.5,.5,1),5), rep(rgb(1,.5,.5),5)))
boxplot(Chromatin_120_min ~ RunonGroup+ isHK,data=dataWithChromGRO, main="120 min",xlab="Length Groups",ylab="% spliced", col = c(rep(rgb(.5,.5,1),5), rep(rgb(1,.5,.5),5)))
dev.off()



############# ############# ############# ############# ############# ############# ############# 
## Heatmap of correlations between features
############# ############# ############# ############# ############# ############# ############# 

f= function(Matrix){
	CC <- cor(Matrix,use='c')
	CC[CC==1] <- 0
	CC
}

pdf("HeatmapFactors1.pdf")

heatmap.2(f(dataWithChromGRO[,c(9:12,14,17,20)]),Rowv=T, main="All Genes",Colv=T,trace='n')
heatmap.2(f(dataWithChromGRO[dataWithChromGRO$isHK==0, c(9:12,17,20)]), Rowv=T, main="Inducible Genes",Colv=T,trace='n')
heatmap.2(f(dataWithChromGRO[dataWithChromGRO$isHK==1, c(9:12,17,20)]), Rowv=T, main="HK genes",Colv=T,trace='n')
heatmap.2(f(dataWithChromGRO[,c(9:12,14,17,20)]),Rowv=F, main="All Genes",Colv=F,trace='n')
heatmap.2(f(dataWithChromGRO[dataWithChromGRO$isHK==0, c(9:12,17,20)]), Rowv=F, main="Inducible Genes",Colv=F,trace='n')
heatmap.2(f(dataWithChromGRO[dataWithChromGRO$isHK==1, c(9:12,17,20)]), Rowv=F, main="HK genes",Colv=F,trace='n')
dev.off()

############# ############# ############# ############# ############# ############# ############# 
############# Principle components and analysis                                     ############# 
############# ############# ############# ############# ############# ############# ############# 

prcomp(~ length + Energy + Acceptor + FeatureCount + PooledSplicing + Runon,data=dataWithChromGRO, scale=TRUE)->PC 
prcomp(~ length + Energy + Acceptor  + PooledSplicing,data=dataWithChromGRO, scale=TRUE)->PC2 

matrixNON = PC$x[rownames(PC$x) %in% rownames(dataWithChromGRO)[dataWithChromGRO$isHK==0],]
matrixHK = PC$x[rownames(PC$x) %in% rownames(dataWithChromGRO)[dataWithChromGRO$isHK==1],]
        
distNON = dist(matrixNON)
distHK  = dist(matrixHK)

clustNON = hclust(distNON)
clustHK = hclust(distHK)

dendroNON = as.dendrogram(clustNON)
dendroHK = as.dendrogram(clustHK)

### Order by tree then sort by one of the columns
#orderClusters = function(myTree, myData, column){
#  outOrder = which(my
  
#  for (i in 1:max(myTree)){
#    temp.data = 
  
#  }
#}

library(dynamicTreeCut)
treeNON = cutreeDynamic(clustNON,distM=as.matrix(distNON),deepSplit=0)
treeHK= cutreeDynamic(clustHK,distM=as.matrix(distHK),deepSplit=0)

pdf("PC_heatmaps.pdf")
heatmap1(matrixNON,main='Non HK genes', maxVal=5, Rowv=dendroNON, RowSideColors=as.character(treeNON),dendrogram='row')
heatmap1(matrixHK,main='HK genes', maxVal=5 ,  Rowv=dendroHK,RowSideColors=as.character(treeHK),dendrogram='row')
dev.off()
pdf("PC_heatmaps2.pdf")
heatmap1(matrixNON,main='Non HK genes', maxVal=5, Rowv=dendroNON, RowSideColors=as.character(treeNON),dendrogram='none')
heatmap1(matrixHK,main='HK genes', maxVal=5 ,  Rowv=dendroHK,RowSideColors=as.character(treeHK),dendrogram='none')
dev.off()






