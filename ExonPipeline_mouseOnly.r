### This script is for analyzing mouse GRO-seq and chromatin-associated-Seq data
source("/home/jeremy/ExonPipeline/ExonPipeline_functions.r")
source("/home/jeremy/ExonPipeline/ExonPipeline_analysis_and_plotting_Functions.r")
source("/home/jeremy/GitHub/useful_R/useful_R.r")

settings = setup(getwd())

tenCols = c(rep(rgb(.5,.5,1),5), rep(rgb(1,.5,.5),5))


mergedData$nonLastLength = mergedData$GeneLength - mergedData$length

#######################
## Apply some groupings
mergedData$AcceptorGroup = splitByVector(mergedData$Acceptor,quantile(mergedData$Acceptor,na.rm=T,(1:4)/5)) # split into 5 groups
mergedData$EnergyGroup = splitByVector(mergedData$Energy,quantile(mergedData$Energy,na.rm=T,(1:4)/5)) # split into 5 groups
mergedData$LengthGroup = splitByVector(mergedData$length,quantile(mergedData$length,na.rm=T,(1:4)/5)) # split into 5 groups
mergedData$FeatureGroup = splitByVector(mergedData$FeatureCount,quantile(mergedData$FeatureCount,na.rm=T,(1:4)/5)) # split into 5 groups


##################################################################################################################
### Do chrom-assoc seq analysis:
##################################################################################################################

features = loadFeatures(settings,refseq,genome)  # 	This excludes intron-less genes
settings = makeLastExonBed(settings,features)
settings = make2ndExonBed(settings,features)
settings = makeAllExonBed(settings,features)
                            
## map reads to last exon
dataDir = "/home/RNAseq/Smale/" 
folders = 352298:352302 
samples = paste("Chromatin_",c(0,15,30,60,120),"_min",sep="")                         

for(i in 1:length(samples)){
  #mapChromAssocBam(settings, sprintf("%schrom_%d/accepted_hits.bam", dataDir,folders[i]) ,paste(samples[i],".txt",sep=""),lastsFile=settings$lastExonBed)
  #mapChromAssocBam(settings, sprintf("%schrom_%d/accepted_hits.bam", dataDir,folders[i]) ,paste(samples[i],"_secondExon.txt",sep=""), lastsFile=settings$secondExonBed)  
  mapChromAssocBam(settings, sprintf("%schrom_%d/accepted_hits.bam", dataDir,folders[i]) ,paste(samples[i],"_allExons.txt",sep=""), lastsFile=settings$allExonBed)  
}

### Once that is finished, parse the output (this gives CTS numbers for each sample)   
outData = get_chromatin_splicing (samples,"")
chromMatrix = mergeChromSplicing(outData, 10,samples)

#outData2 = get_chromatin_splicing (samples,"_secondExon") # also for the 2nd exon

## Also do for ALL introns!!!
outData3 = get_chromatin_splicing (samples,"_allExons") # also for the 2nd exon
chromMatrix3 = mergeChromSplicing(outData3, 10,samples)

## Ok, this may be a stupid test.  But I want to check in each gene whether pooled splicing of the last intron is on average lower than the median.
# This function should return a tuple of median, last intron splicing.
a = chromMatrix3
M=aggregate(PooledSplicing ~ UniqueID,data=a,median, na.rm=T, na.action=na.pass) # obtain the median splicing of each gene's introns
num=aggregate(PooledSplicing ~ UniqueID,data=a,function(x)length(which(!is.na(x))), na.action=na.pass) # obtain the median splicing of each gene's introns
names(num)[2] ='NotNA'
M = merge(M,num,by=1)
last=aggregate(FeatureCount ~ UniqueID,data=a,function(x)m=max(as.numeric(as.character(x)))) # Obtain the list of which is the last intron for each gene
last2 = merge(a,w,by=c("UniqueID","FeatureCount")) # merge with data to find the spplicing of the last intron only
M2 = merge(M[M$NotNA > 3,],data.frame(UniqueID=last2$UniqueID, LastSplicing=last2$PooledSplicing), by="UniqueID")  # add to this matrix: gonna plot!

pdf("Chrom_assoc_allExons_1.pdf")
plot(M2$PooledSplicing, M2$LastSplicing)
abline(b=1,a=0,col='red')
hist(M2$PooledSplicing-M2$LastSplicing,100,xlab="Median splicing - 3' intron splicing")
qqnorm(M2$PooledSplicing-M2$LastSplicing)
qqline(M2$PooledSplicing-M2$LastSplicing)
dev.off()

##################################################################################################################################
# Permutation test to see whether the introns in any given gene have relatively similar splicing to one another
##################################################################################################################################

# Function to calculate the standard deviation
getSplicingSD = function(chrom) tapply(chrom$PooledSplicing,chrom$UniqueID,sd) 

myData = subset(chromMatrix3, !is.na(PooledSplicing),select=c("UniqueID","PooledSplicing"))
# Keep only those that have > 1 non-NA intron
myData = subset(myData, UniqueID %in% subset(as.data.frame(table(myData$UniqueID)),Freq>1,select="Var1")[,1],select=c("UniqueID","PooledSplicing"))
myData.aov = myData
myData.aov[,1] = log(as.numeric(myData.aov[,1]))
myData.aov[,1] = factor(myData.aov[,1])
a=myData.aov[1:20,]
 summary(aov(PooledSplicing ~ UniqueID, data=myData.aov))
 #              Df Sum Sq Mean Sq F value Pr(>F)    
#UniqueID     4499  416.6 0.09260   3.393 <2e-16 ***
#Residuals   25429  693.9 0.02729                   
#---
#Signif. codes:  0 â***â 0.001 â**â 0.01 â*â 0.05 â.â 0.1 â â 1 


## Get std dev. of splicing:
#S = getSplicingSD(myData)
## Do permutations:
#PermuteMeans = rep(0,1000)
#for (i in 1:length(PermuteMeans)){
#  if(i %% 10==0)print(i)
#  myData2 = myData
#  myData2$UniqueID = sample(myData2$UniqueID,nrow(myData2)) 
#  PermuteMeans[i] = mean(getSplicingSD(myData2))
#}
## Get P-value from permutation:
## NULL hypothesis: mean SD of data is >= permuted means.
#P = length(which(PermuteMeans < mean(S)))/length(PermuteMeans) # 0; means P < 0.0001
### YES, on average there is significant similarity between splicing rates within genes

## Intraclass correlation

##############################################################################################################################
#### Specify what I am guessing the Post-splice site length to be.
load("dataWithGRO.RData")
lastLengthsWithGRO = subset(dataWithGRO,!is.na(Runon), select=c("UniqueID","length","Runon")) ## 7713 bp total.

## Add Runon to length, for length2
lastLengthsWithGRO$length2 = lastLengthsWithGRO$length + lastLengthsWithGRO$Runon

#### Now merge either length or length 2 with intron info.
# My time to splice = coSI*Dist/E, but E is constant.

## This function removes NAs from pooled splicing of introns matrix, then merges it with final lengths of all genes, as calculated by GRO-seq.
# Then for each intron it estimates how long it takes (in DISTANCE) to splice out
estimateSplicingTimes = function(intronData, lengthData, lengthColumn){
  # Merge
  introns = merge(subset(intronData, !is.na(PooledSplicing)) , cbind(UniqueID=lastLengthsWithGRO$UniqueID, FinalLength = lastLengthsWithGRO[,lengthColumn]), by = "UniqueID")
  # Get Total length of each intron

  introns$TotalDist = as.numeric(as.character(introns$Dist)) + as.numeric(as.character(introns$FinalLength))
  # Estimate distance to splice
  introns$ExpectedBP = introns$TotalDist * as.numeric(as.character(introns$PooledSplicing))
  
  # Do non-parametric correlation to the 
  introns$Acceptor = as.numeric(as.character(introns$Acceptor))
  
  print(cor(introns$Acceptor, introns$ExpectedBP, method="spearman"))
  
  introns  
}

splicingTimes_lastExonLength = estimateSplicingTimes(allIntronsWithcoSI,lastLengthsWithGRO,'length') 
splicingTimes_GRO = estimateSplicingTimes(allIntronsWithcoSI,lastLengthsWithGRO,'length2') 
 print(cor(log(splicingTimes_GRO$Acceptor), splicingTimes_GRO$ExpectedBP, method="spearman"))
summary(lm(ExpectedBP ~ Acceptor, data=splicingTimes_GRO)) 
summary(lm(ExpectedBP ~ Acceptor, data=splicingTimes_lastExonLength))

pdf("Expected BP.pdf")
plot(log10(splicingTimes_GRO$Acceptor), log10(splicingTimes_GRO$ExpectedBP))
plot(log10(splicingTimes_lastExonLength$Acceptor), log10(splicingTimes_lastExonLength$ExpectedBP))
dev.off()



######## MERGE WITH GENOMIC DATA ###########

dataWithChrom = merge(mergedData, chromMatrix,by='UniqueID',all.x=T)

length(which(!apply(dataWithChrom,1,function(x)any(is.na(x)))))
# Only 1435 genes in all 5 conditions

length(which(!is.na(dataWithChrom$PooledSplicing)))
# 6987 genes have a value for PooledSplicing at > 10 reads total



## Next, I NEED TO ASK THE INTERESTING QUESTIONS: DOES THIS STUFF CORRELATE TO OTHER FEATURES?

#### Apply some groupings
dataWithChrom$SpliceAllGroup = splitByVector(dataWithChrom$PooledSplicing,quantile(dataWithChrom$PooledSplicing,na.rm=T,(1:4)/5)) # split into 5 groups
save(dataWithChrom,file= 'dataWithChrom.RData')

##################################################################################################################
###################  MOUSE GRO-SEQ ###########################
##################################################################################################################

## Read in GRO data
#GRO = read.delim("/home/home/jeremy/RNAseq/Glass/refseq_and_subsequent_transcripts.txt",stringsAsFactors=F)
GRO = read.delim("/home/home/jeremy/RNAseq/Glass/refseq_and_subsequent_transcripts_with_tag_counts.txt",stringsAsFactors=F)

## Calculate the Read-Through (hereafter referred to as 'Runon') and get rid of genes with no measurable read-through.
GRO$Runon = abs(GRO$post_transcript_end - GRO$post_transcript_start)
GRO = GRO[!is.na(GRO$Runon),]

## Annotate by merging with refGene
refGene =add_UniqueID(loadRefgene(settings))
GRO2 = merge(GRO,refGene[,c("name","UniqueID","strand")],by=1)

## Create .bed file, overlap with refseq.bed
refGene$TSS = refGene$txStart
refGene$TSS[refGene$strand=='-'] <- refGene$txEnd[refGene$strand=='-']

options(scipen=10) # make it so I don't write scientific numbers here

## Omit any GRO-seq data when the Read-Through ends within 1Kb of a gene.
write.table(cbind(refGene$chr, refGene$TSS-1000, refGene$TSS +1000 , refGene$name, 0, refGene$strand),file='refseq.bed',sep="\t",row.names=F,col.names=F, quote=F)
write.table(cbind(GRO2$chr, GRO2$post_transcript_start-1, GRO2$post_transcript_end+1, GRO2$UniqueID, 0, GRO2$strand.y),file='GRO2.bed',sep="\t",row.names=F,col.names=F, quote=F)
system("bedtools intersect -a GRO2.bed -b refseq.bed -v -s -wa > GRO_no_TSS_interuptions.bed")
options(scipen=0)

good_GRO =read.delim("GRO_no_TSS_interuptions.bed", head=F)
GRO3 <- GRO2[GRO2$UniqueID %in% good_GRO[,4],]
save(GRO3,file="GRO_Data.RData")

## Merge with 'mergedData'
dataWithGRO = merge(mergedData, GRO3 , by='UniqueID',all.x=T)
dataWithGRO$RunonGroup = splitByVector(dataWithGRO$Runon,quantile(dataWithGRO$Runon,na.rm=T,(1:4)/5)) # split into 5 groups
save(dataWithGRO,file='dataWithGRO.RData')



##################################################################################################################
### Combine GRO-seq and Chrom-assoc seq
##################################################################################################################
dataWithChromGRO = merge(dataWithChrom,dataWithGRO[,c("UniqueID","Runon","RunonGroup")],by="UniqueID",all=T)

length(which(!apply(dataWithChromGRO,1,function(x)any(is.na(x)))))
# 1160 Genes have ALL data

cor(dataWithChromGRO[,c(10:13,15:16,26,28)],use='c') #-> Run-on is NEGATIVELY correlated to SPLICING

#save(mergedData, dataWithChrom, dataWithGRO, dataWithChromGRO, file="MouseChromGRO2.RData"); 
load("MouseChromGRO2.RData")
####################  PLOTS ########################


### In terms of Story, I'm looking at splicing first.  So, if I split up based on splicing groups, what do I see?
# Alternatively, show splicing, split up by other groups?

#png(sprintf('%s_simple_SplitBySpliceAll.png',settings$CommonName))
#par(mfrow=c(3,2))
#plotSplitsHK(dataWithChromGRO,'length','Last Exon Length',Group='SpliceAllGroup')
#plotSplitsHK(dataWithChromGRO,'Energy','Nucleosome Stability',Group='SpliceAllGroup')
#plotSplitsHK(dataWithChromGRO,'Acceptor',"3' splice site strength",Group='SpliceAllGroup')
#plotSplitsHK(dataWithChromGRO,'FeatureCount','Number of Exons',Group='SpliceAllGroup')
#dev.off()


png(sprintf('%s_simple2_SplitSpliceAllbyOthers.png',settings$CommonName), height=1200,width=800)
par(mfrow=c(3,2),cex.axis=1.2,cex=1.3,cex.main=1.1,las=3,mar=c(6.1 ,4.1 ,4.1, 2.1))

xlabs = tapply(mergedData$length,mergedData$Group,function(x)paste(range(x),collapse='-'))
plotSplitsHK(dataWithChromGRO,'PooledSplicing','% Splicing of last intron\nAs a function of last exon length',ylab='',names=xlabs,xlab='')
plotSplits(dataWithChromGRO,'PooledSplicing','% Splicing of last intron\nAs a function of last exon length',add=T,type='b',lwd=2)
# HK : -log p > 5
plotSplitsHK(dataWithChromGRO,'PooledSplicing','% Splicing of last intron\nAs a function of nucleosome stability',Group='EnergyGroup',ylab='',xlab='Stability Score quintiles')
plotSplits(dataWithChromGRO,'PooledSplicing','% Splicing of last intron\nAs a function of nucleosome stability',Group='EnergyGroup',add=T,type='b',lwd=2)
# ~
plotSplitsHK(dataWithChromGRO,'PooledSplicing','% Splicing of last intron\nAs a function of splice site strength',Group='AcceptorGroup',ylab='',xlab='Splice Site Strength quintiles')
plotSplits(dataWithChromGRO,'PooledSplicing','% Splicing of last intron\nAs a function of splice site strength',Group='AcceptorGroup',add=T,type='b',lwd=2)
# NonHK > 3; All > 4
xlabs = tapply(mergedData$FeatureCount,mergedData$FeatureGroup,function(x)paste(range(x),collapse='-'))
plotSplitsHK(dataWithChromGRO,'PooledSplicing','% Splicing of last intron\nAs a function of # exons',Group='FeatureGroup',ylab='',xlab='', names=xlabs)
plotSplits(dataWithChromGRO,'PooledSplicing','% Splicing of last intron\nAs a function of # exons',Group='FeatureGroup',add=T,type='b',lwd=2)
# NonHK, ALL < 16!
xlabs = tapply(mergedData$Runon,mergedData$RunonGroup,function(x)paste(range(x),collapse='-'))
plotSplitsHK(dataWithChromGRO,'PooledSplicing','% Splicing of last intron\nAs a function of Extra Runon Length',Group='RunonGroup',ylab='',xlab='', names=xlabs)
plotSplits(dataWithChromGRO,'PooledSplicing','% Splicing of last intron\nAs a function of Extra Runon Length',Group='RunonGroup',add=T,type='b',lwd=2)
# ~
xlabs = tapply(mergedData$length,mergedData$LengthGroup,function(x)paste(range(x),collapse='-'))
plotSplitsHK(dataWithChromGRO,'PooledSplicing','% Splicing of last intron\nAs a function of last exon length',Group='LengthGroup',ylab='', names=xlabs)
plotSplits(dataWithChromGRO,'PooledSplicing','% Splicing of last intron\nAs a function of last exon length',Group='LengthGroup',add=T,type='b',lwd=2)
# HK : -log p > 4

dev.off()
#
#png(sprintf('%s_simple_SplitSpliceAllbyOthers_noHK.png',settings$CommonName))
#par(mfrow=c(2,2))
#plotSplits(dataWithChromGRO,'PooledSplicing','% Splicing of last intron\nAs a function of last exon length',Group='LengthGroup')
#plotSplits(dataWithChromGRO,'PooledSplicing','% Splicing of last intron\nAs a function of nucleosome stability',Group='EnergyGroup')
#plotSplits(dataWithChromGRO,'PooledSplicing','% Splicing of last intron\nAs a function of splice site strength',Group='AcceptorGroup')
#plotSplits(dataWithChromGRO,'PooledSplicing','% Splicing of last intron\nAs a function of # exons',Group='FeatureGroup')
#dev.off()



png(sprintf('%s_simple2_SplitRunonbyOthers.png',settings$CommonName),height=800,width=1200)
par(mfrow=c(2,3),cex.axis=1.2,cex=1.3,cex.main=1.1,las=3,mar=c(6.1 ,4.1 ,4.1, 2.1))

xlabs = tapply(mergedData$length,mergedData$Group,function(x)paste(range(x),collapse='-'))
plotSplitsHK(dataWithChromGRO,'Runon','Extra Runon Length\nAs a function of last exon length',Group='LengthGroup',xlab='',names=xlabs,ylab='')
plotSplits(dataWithChromGRO,'Runon','Extra Runon Length\nAs a function of last exon length',Group='LengthGroup',add=T,type='b',lwd=2)
# Each is > 4
plotSplitsHK(dataWithChromGRO,'Runon','Extra Runon Length\nAs a function of nucleosome stability',Group='EnergyGroup',ylab='',xlab='Stability Score quintiles')
plotSplits(dataWithChromGRO,'Runon','Extra Runon Length\nAs a function of nucleosome stability',Group='EnergyGroup',add=T,type='b',lwd=2)
# HK, all > 2

xlabs = tapply(dataWithChromGRO$PooledSplicing,dataWithChromGRO$SpliceAllGroup,function(x)paste(sprintf("%.2f",range(x)),collapse='-'))
print(xlabs)
plotSplitsHK(dataWithChromGRO,'Runon','Extra Runon Length\nAs a function of Last Intron Splicing %',Group='SpliceAllGroup',xlab='',ylab='',names=xlabs)
plotSplits(dataWithChromGRO,'Runon','',Group='SpliceAllGroup',add=T,type='b',lwd=2)
# All > 4
plotSplitsHK(dataWithChromGRO,'Runon','Extra Runon Length\nAs a function of splice site strength',Group='AcceptorGroup',ylab='',xlab='Splice Site Strength quintiles')
plotSplits(dataWithChromGRO,'Runon','Extra Runon Length\nAs a function of splice site strength',Group='AcceptorGroup',add=T,type='b',lwd=2)
# ~
xlabs = tapply(mergedData$FeatureCount,mergedData$FeatureGroup,function(x)paste(range(x),collapse='-'))
plotSplitsHK(dataWithChromGRO,'Runon','Extra Runon Length\nAs a function of # exons',Group='FeatureGroup',ylab='',xlab='',names=xlabs)
plotSplits(dataWithChromGRO,'Runon','Extra Runon Length\nAs a function of # exons',Group='FeatureGroup',add=T,type='b',lwd=2)
# ~All >3, non > 4
dev.off()




######################################################################################################################################################
######################################################################################################################################################
## 10/25/2012

## Bayesian network inference:

#bayesChrom = dataWithChrom[,c('length',  'Energy',  'Acceptor',  'FeatureCount',  'nonLastLength')]
#res = gs(bayesChrom) I don't think this is working
#res=gs(round(dataWithChrom[,c('length',  'Energy',  'Acceptor',  'FeatureCount',  'nonLastLength')]))

bayesChrom = dataWithChrom[,c('length',  'Energy',  'Acceptor',  'FeatureCount',  'nonLastLength')]
bayesChrom[,1] = log(bayesChrom[,1]); bayesChrom[,3] = log(bayesChrom[,3]); bayesChrom[,5] = log(bayesChrom[,5])

res = gs(bayesChrom+0.0, test='mc-cor')
res2 = iamb(bayesChrom+0.0, test='mc-cor')
res3 = gs(bayesChrom+0.0, test='mc-zf')
res4 = fast.iamb(bayesChrom+0.0, test='mc-zf')

pdf('Bayesian Graphs.pdf');par(mfrow=c(2,2))
plot(res,main='gs, mc-cor')
plot(res2,main='iamb, mc-cor')
plot(res3,main='gs, mc-zf')
plot(res4,main='fast iamb, mc-zf')
dev.off()

### Principle Components, to explain the splicing data: 
prcomp(~length + Energy + Acceptor + FeatureCount + nonLastLength,data=dataWithChrom,scale=T)->PC
dataWithChrom2 = cbind(dataWithChrom, PC$x)
summary(lm(PooledSplicing ~ PC1 + PC2 + PC3 + PC4 + PC5,data=dataWithChrom2))

### Try it for GRO data:
prcomp(~length + Energy + Acceptor + FeatureCount + nonLastLength,data=dataWithChromGRO,scale=T)->PC2
dataWithChromGRO2 = cbind(dataWithChromGRO, PC2$x)
summary(lm(Runon ~ PC1 + PC2 + PC3 + PC4 + PC5 + PooledSplicing,data=dataWithChromGRO2))


######################################################################################################################################################
######################################################################################################################################################



######################################################################################################################################################
######################################################################################################################################################
## 10/29/2012

################# EXPECTED TIME TO FINISH SPLICING #################
#### Run Analytical model of CTS for all genes, to calculate expected time of COMPLETE splicing, AFTER the last intron as been synthesized

# Prepare introns for calculating each genes' time
formatIntronsForExpectedTimeModel(intronsFile='introns.txt', scoreFile='intron_acceptors_for_Burset.txt.scores', 'introns_for_expectedTimeModel.txt');

# Call MATLAB model:
calculateExpectedTime_MATLAB ('introns_for_expectedTimeModel.txt',elongation_rate=3000,splicingFactor=5);
calculateExpectedTime_MATLAB ('introns_for_expectedTimeModel.txt',elongation_rate=3000,splicingFactor=2);

### Now compare to other things ...
expectedTimeData = read.delim('introns_for_expectedTimeModel_expectedSplicingTimes_3000.00_5.00.txt',head=F,col.names=c("UniqueID","MinutesToSplice"))
dataWithChromGROTime = merge(dataWithChromGRO,expectedTimeData,all.x=T, by = 'UniqueID')
dataWithChromGROTime$TimeGroup = splitByVector(dataWithChromGROTime$MinutesToSplice,quantile(dataWithChromGROTime$MinutesToSplice,na.rm=T,(1:4)/5)) # split into 5 groups

dataWithChromGROTime$LengthAfterIntron = rowMeans(dataWithChromGROTime[,c("length","Runon")], na.rm=T)
dataWithChromGROTime$LengthAfterIntronGroup = splitByVector(dataWithChromGROTime$LengthAfterIntron,quantile(dataWithChromGROTime$LengthAfterIntron,na.rm=T,(1:4)/5)) # split into 5 groups

### Try some multiple regressions based on talking to Kat.
# Regress Read-through time ~ expected time to splice + Nucleosome occupancy + last exon length
# Turn Energy into unit normal
Emean = mean(dataWithChromGROTime$Energy,na.rm=T)
Esd = sd(dataWithChromGROTime$Energy,na.rm=T)
dataWithChromGROTime$UnitEnergy = (dataWithChromGROTime$Energy-Emean)/Esd

summary(lm(Runon ~ MinutesToSplice + UnitEnergy + length, data=dataWithChromGROTime))



png(sprintf('%s_simple2_SplitTimebyOthers.png',settings$CommonName),height=800,width=1200)
par(mfrow=c(2,3),cex.axis=1.2,cex=1.3,cex.main=1.1,las=3,mar=c(6.1 ,4.1 ,4.1, 2.1))

xlabs = tapply(mergedData$length,mergedData$Group,function(x)paste(range(x),collapse='-'))
plotSplitsHK(dataWithChromGROTime,'MinutesToSplice','Minutes until splicing is complete\nAs a function of last exon length',Group='LengthGroup',xlab='',names=xlabs,ylab='')
plotSplits(dataWithChromGROTime,'MinutesToSplice','',Group='LengthGroup',add=T,type='b',lwd=2)
#
plotSplitsHK(dataWithChromGROTime,'MinutesToSplice','Minutes until splicing is complete\nAs a function of nucleosome stability',Group='EnergyGroup',ylab='',xlab='Stability Score quintiles')
plotSplits(dataWithChromGROTime,'MinutesToSplice','',Group='EnergyGroup',add=T,type='b',lwd=2)
#

xlabs = tapply(dataWithChromGRO$PooledSplicing,dataWithChromGRO$SpliceAllGroup,function(x)paste(sprintf("%.2f",range(x)),collapse='-'))
print(xlabs)
plotSplitsHK(dataWithChromGROTime,'MinutesToSplice','Minutes until splicing is complete\nAs a function of Last Intron Splicing %',Group='SpliceAllGroup',xlab='',ylab='',names=xlabs)
plotSplits(dataWithChromGROTime,'MinutesToSplice','',Group='SpliceAllGroup',add=T,type='b',lwd=2)
#
xlabs = tapply(mergedData$FeatureCount,mergedData$FeatureGroup,function(x)paste(range(x),collapse='-'))
plotSplitsHK(dataWithChromGROTime,'MinutesToSplice','Minutes until splicing is complete\nAs a function of # exons',Group='FeatureGroup',ylab='',xlab='',names=xlabs)
plotSplits(dataWithChromGROTime,'MinutesToSplice','',Group='FeatureGroup',add=T,type='b',lwd=2)
#
xlabs = tapply(dataWithGRO$Runon,dataWithGRO$RunonGroup,function(x)paste(range(x),collapse='-'))
plotSplitsHK(dataWithChromGROTime,'MinutesToSplice','Minutes until splicing is complete\nAs a function of Read Through length',Group='RunonGroup',ylab='',xlab='',names=xlabs)
plotSplits(dataWithChromGROTime,'MinutesToSplice','',Group='RunonGroup',add=T,type='b',lwd=2)
# All > 2, non > 2

xlabs = tapply(dataWithChromGROTime$LengthAfterIntron,dataWithChromGROTime$LengthAfterIntronGroup,function(x)paste(range(x),collapse='-'))
plotSplitsHK(dataWithChromGROTime,'MinutesToSplice','Minutes until splicing is complete\nAs a function of Read Through length + last Exon Length',Group='LengthAfterIntronGroup',ylab='',xlab='',names=xlabs)
plotSplits(dataWithChromGROTime,'MinutesToSplice','',Group='LengthAfterIntronGroup',add=T,type='b',lwd=2)
# All > 23, non > 24!!!
cor.test(dataWithChromGROTime$LengthAfterIntron,dataWithChromGROTime$MinutesToSplice) 

dev.off()


 
##################################################################################################################################
### 11/6/2012 ###
# How do my splice site and actual slicing times align if I make some assumptions about CTS?
##################################################################################################################################

allIntrons_withDists = read.delim("introns_for_expectedTimeModel.txt",head=F, col.names=c("UniqueID","DistToLastSpliceSite","AcceptorScores"),stringsAsFactors=F)
# Make a new data frame that I can merge with

allIntrons2 = apply(allIntrons_withDists,1,function(x){
  dists = unlist(strsplit(x[2],","))
  scores = unlist(strsplit(x[3],","))
  IDs = paste( as.character(sub("^\\s*","",x[1])), 1:length(dists), sep="_" )

  cbind(ID = as.character(IDs), Dist = as.character(dists), Acceptor = as.character(scores))
  } )

allIntronsMatrix = data.frame(matrix(unlist(lapply(allIntrons2,t)), nc=3, byrow=T,dimnames=dimnames(allIntrons2[[1]])))
allIntronsMatrix$UniqueID   = sub("^([0-9]*)_([0-9]*)$","\\1",allIntronsMatrix$ID)
allIntronsMatrix$FeatureCount   = sub("^([0-9]*)_([0-9]*)$","\\2",allIntronsMatrix$ID)

allIntronsWithcoSI = merge(allIntronsMatrix,subset(chromMatrix3,select=c("PooledSplicing","UniqueID","FeatureCount")), by = c("UniqueID","FeatureCount"))



######################################################################################################################################################
######################################################################################################################################################
### 11/7/2012
# For each pair of important variables, make a scatter plot with HK vs NON-hk genes.
dataForScatter = subset(dataWithChromGROTime,select=c("isHK","PooledSplicing","MinutesToSplice", "length", "Runon","FeatureCount","Energy", "Acceptor","Donor", "LengthAfterIntron", "GeneLength","nonLastLength"))
dataForScatter$SplicingProb =  (dataForScatter$Runon * dataForScatter$PooledSplicing) #
for (doLog in c(4,5,8,9,10,11,12,13)){
  dataForScatter[,doLog] = log10(dataForScatter[,doLog])
}
dataForScatter = dataForScatter[order(dataForScatter$isHK),]

cor(dataForScatter)

plotMultiscatter = function(dataIn,filename,nr=3,nc=4,W=10,H=8, ...){
  if(length(grep("png",filename))>0) png(filename,width=W,height=H)
  if(length(grep("pdf",filename))>0) pdf(filename,width=W,height=H)
   
  par(mfrow=c(nr,nc)) 
  cols = sub("0","blue",sub("1","red",as.numeric(dataIn[,1])))              # Make HK its own category, as colors (red, blue for 1,0)  
  NAMES = dimnames(dataIn)[[2]]
  
  for(i in 2:(ncol(dataIn)-1))
    for(j in (i+1):ncol(dataIn))
      plot(dataIn[,i],dataIn[,j],xlab=NAMES[i], ylab=NAMES[j],col=cols, ...)
      
  dev.off()    
}

plotMultiscatter(dataForScatter,'ScatterData1_%d.png',W=1600,H=1200,cex.lab=2,pch=16)


######################################################################################################################################################
######################################################################################################################################################
### 11/8/2012
# Linear models: Try look at extreme categories.

doLM = function(Formula,inData,name){
	cat(sprintf("********%s********:\n",name))
	
	cat(sprintf("*** Housekeeping genes ***:\n"))
	print(summary(lm(Formula, data=inData, subset=(isHK==1)))[c("coefficients","adj.r.squared","r.squared")])

	cat(sprintf("**** Inducible genes ****:\n"))
	print(summary(lm(Formula, data=inData, subset=(isHK==0)))[c("coefficients","adj.r.squared","r.squared")])
	NULL
}

doLM(PooledSplicing ~ Acceptor * (length + Energy), subset(dataWithChromGROTime,SpliceAllGroup==5), 'High Splicers')
doLM(PooledSplicing ~ Acceptor * (length + Energy), subset(dataWithChromGROTime,SpliceAllGroup==1), 'Low Splicers')

doLM(PooledSplicing ~ Acceptor + (LengthAfterIntron + Energy), subset(dataWithChromGROTime,SpliceAllGroup==5), 'High Splicers')
doLM(PooledSplicing ~ Acceptor + (LengthAfterIntron + Energy), subset(dataWithChromGROTime,SpliceAllGroup==1), 'Low Splicers')


doLM(PooledSplicing ~ MinutesToSplice + (LengthAfterIntron + Energy), subset(dataWithChromGROTime,SpliceAllGroup==5), 'High Splicers')
doLM(PooledSplicing ~ MinutesToSplice + (LengthAfterIntron + Energy), subset(dataWithChromGROTime,SpliceAllGroup==1), 'Low Splicers')


doLM(PooledSplicing ~ MinutesToSplice * (LengthAfterIntron ), subset(dataWithChromGROTime,SpliceAllGroup==5), 'High Splicers')
doLM(PooledSplicing ~ MinutesToSplice * (LengthAfterIntron ), subset(dataWithChromGROTime,SpliceAllGroup==1), 'Low Splicers')


doLM(PooledSplicing ~ Acceptor + (LengthAfterIntron + Energy), subset(dataWithChromGROTime,TimeGroup==5), 'High time to splice needed')
doLM(PooledSplicing ~ Acceptor + (LengthAfterIntron + Energy), subset(dataWithChromGROTime,TimeGroup==1), 'Low time to splice needed')


doLM(Runon ~ length + MinutesToSplice + ( Energy), subset(dataWithChromGROTime,SpliceAllGroup==5), 'High Splicers')
doLM(Runon ~ length + MinutesToSplice + ( Energy), subset(dataWithChromGROTime,SpliceAllGroup==1), 'Low Splicers')


######################################################################################################################################################
######################################################################################################################################################
###### 11/8/2012
# Based on Alex's email, look at a cohort of genes that I expect to be not well spliced.  For example, take length/time.  What do I expect the runon to look like?

## Read in NEW matlab test
expectedTimeData2 = read.delim('introns_for_expectedTimeModel_expectedSplicingTimes_3000.00_2.00.txt',head=F,col.names=c("UniqueID","MinutesToSplice2"))
dataWithChromGROTime2 = merge(dataWithChromGROTime,expectedTimeData2,all.x=T, by = 'UniqueID')
dataWithChromGROTime2$TimeGroup2 = splitByVector(dataWithChromGROTime2$MinutesToSplice2,quantile(dataWithChromGROTime2$MinutesToSplice2,na.rm=T,(1:4)/5)) # split into 5 groups

####### Make an arbitrary mapping between energy and elongation time.
# Convert energy to a Z score
dataWithChromGROTime2$EnergyZ = normalize(dataWithChromGROTime2$Energy)
# Let's say that each Z value counts as what, 3 seconds?
BP_PerEnergyZ = 1/10*3000 # This means that 3000 bp/minute *1/10 minute = 300 bp. 
OffsetBP = 3.5 # O means i will penalize some genes, 3.5 means all genes at at least a tiny boost.
Energy_BP = BP_PerEnergyZ*(dataWithChromGROTime2$EnergyZ + OffsetBP)

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

### Estimate time it took to splice ... 
dataWithChromGROTime2$expectedSplicing =   ((dataWithChromGROTime2$LengthAfterIntron +  Energy_BP)/2000)
dataWithChromGROTime2$expectedSplicing2 =   (1-dataWithChromGROTime2$PooledSplicing)*((dataWithChromGROTime2$LengthAfterIntron +  Energy_BP)/2000) 
NAME1 = "Time available to splice (min)\n(Length + Runon, 6sec/EnergyZ)"
NAME2 = "Time took to Splice (min)\n(Length + Runon, 6sec/EnergyZ)*SplicingMeasured"

dataWithChromGROTime2$expectedSplicingGroup = splitByVector(dataWithChromGROTime2$expectedSplicing, quantile(dataWithChromGROTime2$expectedSplicing,na.rm=T,(1:4)/5)) # split into 5 groups
dataWithChromGROTime2$expectedSplicingGroup2 = splitByVector(dataWithChromGROTime2$expectedSplicing2, quantile(dataWithChromGROTime2$expectedSplicing2,na.rm=T,(1:4)/5)) # split into 5 groups

dataWithChromGROTime2$RankEnergyAcceptor = rank(dataWithChromGROTime2$Energy) + rank(dataWithChromGROTime2$Acceptor) + rank(dataWithChromGROTime2$length)
dataWithChromGROTime2$RankEnergyAcceptorGroup = splitByVector(dataWithChromGROTime2$RankEnergyAcceptor, quantile(dataWithChromGROTime2$RankEnergyAcceptor,na.rm=T,(1:4)/5)) # split into 5 groups

dataWithChromGROTime2$RankEnergyAcceptorRunon = rank(dataWithChromGROTime2$Energy) + rank(dataWithChromGROTime2$Acceptor) + rank(dataWithChromGROTime2$length)+ rank(dataWithChromGROTime2$Runon)
dataWithChromGROTime2$RankEnergyAcceptorRunonGroup = splitByVector(dataWithChromGROTime2$RankEnergyAcceptorRunon, quantile(dataWithChromGROTime2$RankEnergyAcceptorRunon,na.rm=T,(1:4)/5)) # split into 5 groups


png(sprintf('%s_simple3_Split_expectedSplicing8.png',settings$CommonName),height=800,width=1800)
par(mfcol=c(2,4),cex.axis=1.2,cex=1.3,cex.main=1,las=3,mar=c(6.1 ,4.1 ,4.1, 2.1)); options(warn=-1) -> W
xlabs = tapply(dataWithChromGROTime2$MinutesToSplice,dataWithChromGROTime2$TimeGroup,function(x)paste(sprintf("%.1f",range(x)),collapse='-'))
plotSplitsHK(dataWithChromGROTime2,'Runon',paste('A: Read-through distance (bp)','As a function of Estimated time needed to splice',sep="\n"),Group='TimeGroup2',xlab='',ylab='',names=xlabs)
plotSplits(dataWithChromGROTime2,'Runon','',Group='TimeGroup2',add=T,type='b',lwd=2)

plotSplitsHK(dataWithChromGROTime2,'PooledSplicing',paste('B: Measured Splicing of Last intron','As a function of Estimated time needed to splice',sep="\n"),Group='TimeGroup2',xlab='',ylab='',names=xlabs)
plotSplits(dataWithChromGROTime2,  'PooledSplicing','',Group='TimeGroup2',add=T,type='b',lwd=2)

predictSimple2 = plotSplitsHK(dataWithChromGROTime2,'expectedSplicing2',paste("C: ",NAME2,'\nAs a function of Estimated time needed to splice',sep=""),Group='TimeGroup2',xlab='',ylab='',names=xlabs)
plotSplits(dataWithChromGROTime2,'expectedSplicing2','',Group='TimeGroup2',add=T,type='b',lwd=2)

xlabs = tapply(dataWithChromGROTime2$Acceptor,dataWithChromGROTime2$AcceptorGroup,function(x)paste(sprintf("%.e",range(x)),collapse='-'))
predictSimple2 = plotSplitsHK(dataWithChromGROTime2,'expectedSplicing2',paste("D: ", NAME2,'\nAs a function of Acceptor Strength of Last Intron'),Group='AcceptorGroup',xlab='',ylab='',names=xlabs)
plotSplits(dataWithChromGROTime2,'expectedSplicing2','',Group='AcceptorGroup',add=T,type='b',lwd=2)

xlabs = tapply(dataWithChromGROTime2$PooledSplicing,dataWithChromGROTime2$SpliceAllGroup,function(x)paste(sprintf("%.2f",range(x)),collapse='-'))
predictSimple2 = plotSplitsHK(dataWithChromGROTime2,'expectedSplicing',paste("E: ", NAME1,'\nAs a function of Measured Splicing of Last Intron'),Group='SpliceAllGroup',xlab='',ylab='',names=xlabs)
plotSplits(dataWithChromGROTime2,'expectedSplicing','',Group='SpliceAllGroup',add=T,type='b',lwd=2)

plotSplitsHK(dataWithChromGROTime2,'Runon',paste('F: Read-through distance (bp)','\nAs a function of Measured Splicing of Last Intron'),Group='SpliceAllGroup',xlab='',ylab='',names=xlabs)
plotSplits(dataWithChromGROTime2,'Runon','',Group='SpliceAllGroup',add=T,type='b',lwd=2)

### Look at interaction of Energy and Acceptor
xlabs = tapply(dataWithChromGROTime2$RankEnergyAcceptor,dataWithChromGROTime2$RankEnergyAcceptorGroup,function(x)paste(sprintf("%.f",range(x)),collapse='-'))
plotSplitsHK(dataWithChromGROTime2,'PooledSplicing',paste('G: Last Intron Splicing','\nAs a function of Predicted Splicing Group'),Group='RankEnergyAcceptorGroup',xlab='',ylab='',names=xlabs)
plotSplits(dataWithChromGROTime2,'PooledSplicing','',Group='RankEnergyAcceptorGroup',add=T,type='b',lwd=2)

xlabs = tapply(dataWithChromGROTime2$RankEnergyAcceptorRunon,dataWithChromGROTime2$RankEnergyAcceptorRunonGroup,function(x)paste(sprintf("%.f",range(x)),collapse='-'))
plotSplitsHK(dataWithChromGROTime2,'PooledSplicing',paste('H: Last Intron Splicing','\nAs a function of Predicted Splicing Group w/ Read-through'),Group='RankEnergyAcceptorRunonGroup',xlab='',ylab='',names=xlabs)
plotSplits(dataWithChromGROTime2,'PooledSplicing','',Group='RankEnergyAcceptorGroup',add=T,type='b',lwd=2)

options(W)
dev.off()


png(sprintf('%s_simple3_Split_expectedSplicing9.png',settings$CommonName),height=800,width=1800)
par(mfcol=c(2,4),cex.axis=1.2,cex=1.3,cex.main=1,las=3,mar=c(6.1 ,4.1 ,4.1, 2.1)); options(warn=-1) -> W

xlabs = tapply(dataWithChromGROTime2$FeatureCount,dataWithChromGROTime2$FeatureGroup,function(x)paste(range(x),collapse='-'))
plotSplitsHK(dataWithChromGROTime2[dataWithChromGROTime2$RunonGroup==1,],'PooledSplicing','% Splicing of last intron\nAs a function of # exons',Group='FeatureGroup',ylab='',xlab='', names=xlabs)
plotSplits(dataWithChromGROTime2[dataWithChromGROTime2$RunonGroup==1,],'PooledSplicing','% Splicing of last intron\nAs a function of # exons',Group='FeatureGroup',add=T,type='b',lwd=2)

plotSplitsHK(dataWithChromGROTime2[dataWithChromGROTime2$RunonGroup==5,],'PooledSplicing','% Splicing of last intron\nAs a function of # exons',Group='FeatureGroup',ylab='',xlab='', names=xlabs)
plotSplits(dataWithChromGROTime2[dataWithChromGROTime2$RunonGroup==5,],'PooledSplicing','% Splicing of last intron\nAs a function of # exons',Group='FeatureGroup',add=T,type='b',lwd=2)

xlabs = tapply(dataWithChromGROTime2$Runon,dataWithChromGROTime2$RunonGroup,function(x)paste(range(x),collapse='-'))
plotSplitsHK(dataWithChromGROTime2[dataWithChromGROTime2$FeatureGroup==1,],'PooledSplicing','Low-intron genes: % Splicing of last intron\nAs a function of Read-Through',Group='RunonGroup',ylab='',xlab='', names=xlabs)
plotSplits(dataWithChromGROTime2[dataWithChromGROTime2$FeatureGroup==1,],'PooledSplicing','% Splicing of last intron\nAs a function of # exons',Group='RunonGroup',add=T,type='b',lwd=2)

plotSplitsHK(dataWithChromGROTime2[dataWithChromGROTime2$FeatureGroup==5,],'PooledSplicing','Many-intron genes: % Splicing of last intron\nAs a function of Read-Through',Group='RunonGroup',ylab='',xlab='', names=xlabs)
plotSplits(dataWithChromGROTime2[dataWithChromGROTime2$FeatureGroup==5,],'PooledSplicing','% Splicing of last intron\nAs a function of # exons',Group='RunonGroup',add=T,type='b',lwd=2)

options(W); dev.off()


######################################################################################################################################################
######################################################################################################################################################
###### 11/26/2012
#### Prepare to simulate splicing of ALL genes, in MATLAB, with the Read-Through length as well.

## First collect the introns' splice site scores: this collects the normalized scores
formatIntronsForExpectedTimeModel(intronsFile='introns.txt', scoreFile='intron_acceptors_for_Burset.txt.NormScores', 'introns_for_analyticalModel.txt');

## Now I need to merge those score with RefSeq, Read-through length, and EnergyZ.
intronsForProcessing = read.delim('introns_for_analyticalModel.txt', head=F,col.names=c("UniqueID","DistsToLastExon","SpliceStrengths"))

load("PlottedData3.RData")  # This updated mergedData object has Z-scores for Energy
load("GRO_Data.RData")
data3WithGRO = merge(mergedData, GRO3 , by='UniqueID',all.x=T) 
temp1 = merge(intronsForProcessing[,c(1,3)], subset(data3WithGRO, !is.na(Runon), select=c("UniqueID", "Energy", "Runon")), by="UniqueID") # 7525 genes     

refGene =add_UniqueID(loadRefgene(settings))
temp2 = merge(temp1, subset(refGene, select= c("UniqueID","txStart","txEnd","exonStarts","exonEnds","strand")), by = "UniqueID")

## Save the data to file
write.delim(temp2, File='genes_for_analyticalModel_withRunon.txt')

####### Call the matlab code that simulates all genes:
dat = calculateSplicingRunon_MATLAB('genes_for_analyticalModel_withRunon.txt',3000,30,'0.5,2,15,60')

## Merge and plot
png('modelSplicingRunon3.png')
plot(dat[,2],dat[,3],ylab='WithRunon', xlab='No Runon', main='Predicted CTS efficiency')
dev.off()

dat$GRO_improvement = dat$PredSplice2 - dat$PredSplice1
summary( dat$PredSplice2)

dataWithChromGROpredictions = merge(dataWithChromGRO,dat, by='UniqueID',all.x=T)

######################################################################################################################################################
source("/home/jeremy/ExonPipeline/ExonPipeline_analysis_and_plotting_Functions.r")
png(sprintf('%s_simple4.2_SplitRunonImprovementbyOthers3.png',settings$CommonName),height=800,width=1200)
par(mfrow=c(2,3),cex.axis=1.2,cex=1.3,cex.main=1.1,las=3,mar=c(6.1 ,4.1 ,4.1, 2.1)) ; options(warn=-1) -> W

xlabs = tapply(mergedData$length,mergedData$Group,function(x)paste(range(x),collapse='-'))

plotSplitsHK(dataWithChromGROpredictions,'PredSplice1','Splicing \nAs a function of last exon length',Group='Group',xlab='',names=xlabs,ylab='',YLIM=c(0.2,0.9))
plotSplits(dataWithChromGROpredictions,'PredSplice2','',Group='Group',add=T,type='b',lwd=2)
plotSplitsHK(dataWithChromGROpredictions,'PredSplice2','Splicing with Read-through \nAs a function of last exon length',Group='Group',xlab='',names=xlabs,ylab='',YLIM=c(0.2,0.9))
plotSplits(dataWithChromGROpredictions,'PredSplice2','',Group='Group',add=T,type='b',lwd=2)

plotSplitsHK(dataWithChromGROpredictions,'PredSplice1','Splicing\nAs a function of last exon length',Group='Group',ylab='',xlab='Last Exon Length',names=xlabs,YLIM=c(0.2,0.8))
plotSplitsHK(dataWithChromGROpredictions,'PredSplice2','',Group='Group',add=T,LTY=2,YLIM=c(0.2,0.8))

xlabs = tapply(mergedData$length,mergedData$LengthGroup,function(x)paste(range(x),collapse='-'))
plotSplitsHK(dataWithChromGROpredictions,'PredSplice1','Splicing\nAs a function of LengthGroup',Group='LengthGroup',ylab='',xlab='',names=xlabs,YLIM=c(0.2,0.9))
plotSplits(dataWithChromGROpredictions,'PredSplice1','',Group='LengthGroup',add=T,type='b',lwd=2)
plotSplitsHK(dataWithChromGROpredictions,'PredSplice2','Splicing with Read-through \nAs a function of LengthGroup',Group='LengthGroup',ylab='',xlab='',names=xlabs,YLIM=c(0.2,0.9))
plotSplits(dataWithChromGROpredictions,'PredSplice2','',Group='LengthGroup',add=T,type='b',lwd=2)

plotSplitsHK(dataWithChromGROpredictions,'PredSplice1','Splicing\nAs a function of LengthGroup',Group='LengthGroup',ylab='',xlab='Last Exon Length',names=xlabs,YLIM=c(0.3,0.8))
plotSplitsHK(dataWithChromGROpredictions,'PredSplice2','',Group='LengthGroup',add=T,LTY=2,YLIM=c(0.3,0.8))

options(W); dev.off()

summary(lm(GRO_improvement ~ log10(length) + FeatureCount, data=dataWithChromGROpredictions))
######################################################################################################################################################

png(sprintf('%s_simple7_SplitRunonImprovementbyOthers3.png',settings$CommonName),height=600,width=1200)
par(mfrow=c(1,2),cex.axis=1.2,cex=1.3,cex.main=1.1,las=3,mar=c(6.1 ,4.1 ,4.1, 2.1)) ; options(warn=-1) -> W

xlabs = tapply(mergedData$length,mergedData$Group,function(x)paste(range(x),collapse='-'))
plotSplitsHK(dataWithChromGRO,'Runon','GRO-seq Read-Through Length\nAs a function of last exon length',Group='LengthGroup',xlab='',names=xlabs,ylab='')

xlabs = tapply(mergedData$length,mergedData$Group,function(x)paste(range(x),collapse='-'))
plotSplitsHK(dataWithChromGROpredictions,'PredSplice1','Splicing\nAs a function of last exon length',Group='Group',xlab='',ylab='',names=xlabs,YLIM=c(0.0,0.4))
plotSplitsHK(dataWithChromGROpredictions,'PredSplice2','',Group='Group',add=T,LTY=2,YLIM=c(0.0,0.4),xlab='Last Exon Length',ylbias=-1)

options(W); dev.off()

### Try doing GRO while ALLOWING 0's:

##################################################################################################################
###################  MOUSE GRO-SEQ ###########################
##################################################################################################################

## Read in GRO data
GRO = read.delim("/home/home/jeremy/RNAseq/Glass/refseq_and_subsequent_transcripts_with_tag_counts.txt",stringsAsFactors=F)

## Calculate the Read-Through (hereafter referred to as 'Runon') and ***set genes with no measurable read-through to 0***
NAs = is.na(GRO$post_transcript_start)
GRO$post_transcript_start[NAs & GRO$strand==0] <- GRO$gene_end + 1
GRO$post_transcript_start[NAs & GRO$strand==1] <- GRO$gene_start - 1
GRO$post_transcript_end[NAs] <- GRO$post_transcript_start[NAs]

GRO$Runon = abs(GRO$post_transcript_end - GRO$post_transcript_start)

## Annotate by merging with refGene
refGene =add_UniqueID(loadRefgene(settings))
GRO2 = merge(GRO,refGene[,c("name","UniqueID","strand")],by=1)

options(scipen=10) # make it so I don't write scientific numbers here
## Omit any GRO-seq data when the Read-Through ends within 1Kb of a gene.
write.table(cbind(GRO2$chr, GRO2$post_transcript_start-1, GRO2$post_transcript_end+1, GRO2$UniqueID, 0, GRO2$strand.y),file='GRO2.2.bed',sep="\t",row.names=F,col.names=F, quote=F)
system("bedtools intersect -a GRO2.2.bed -b refseq.bed -v -s -wa > GRO_no_TSS_interuptions2.bed")
options(scipen=0)

good_GRO =read.delim("GRO_no_TSS_interuptions2.bed", head=F)
GRO3 <- GRO2[GRO2$UniqueID %in% good_GRO[,4],]
save(GRO3,file="GRO_with0_Data.RData")

## Merge with 'mergedData'
dataWithGRO2 = merge(mergedData, GRO3 , by='UniqueID',all.x=T)
dataWithGRO2$RunonGroup = splitByVector(dataWithGRO2$Runon,quantile(dataWithGRO2$Runon,na.rm=T,(1:4)/5)) # split into 5 groups
save(dataWithGRO2,file='dataWithGRO_with0.RData')

##################################################################################################################
################### Re-plotting the first figure at least.  I should probably redo the simulation then, if it looks worthwhile!
png(sprintf('%s_simple6_SplitRunonWith0ImprovementbyOthers3.png',settings$CommonName),height=600,width=1200)
par(mfrow=c(1,2),cex.axis=1.2,cex=1.3,cex.main=1.1,las=3,mar=c(6.1 ,4.1 ,4.1, 2.1)) ; options(warn=-1) -> W

xlabs = tapply(mergedData$length,mergedData$Group,function(x)paste(range(x),collapse='-'))
plotSplitsHK(dataWithGRO2,'Runon','Read-Through Length\nAs a function of last exon length',Group='LengthGroup',xlab='',names=xlabs,ylab='')

xlabs = tapply(mergedData$length,mergedData$Group,function(x)paste(range(x),collapse='-'))
plotSplitsHK(dataWithChromGROpredictions,'PredSplice1','Splicing\nAs a function of last exon length',Group='Group',xlab='',ylab='',names=xlabs,YLIM=c(0.2,0.8))
plotSplitsHK(dataWithChromGROpredictions,'PredSplice2','',Group='Group',add=T,LTY=2,YLIM=c(0.2,0.8),xlab='Last Exon Length',ylbias=-1)

options(W); dev.off()




















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






