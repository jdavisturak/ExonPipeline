### This script is for analyzing mouse GRO-seq and chromatin-associated-Seq data
source("/home/jeremy/ExonPipeline/ExonPipeline_functions.r")
source("/home/jeremy/ExonPipeline/ExonPipeline_analysis_and_plotting_Functions.r")
source("/home/jeremy/Code/useful_R/useful_R.r")

setwd("/home/home/jeremy/ExonPipeline/mm9")
load("PlottedData4.RData")

settings = setup(getwd())

tenCols = getwc(rep(rgb(.5,.5,1),5), rep(rgb(1,.5,.5),5))


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
## Load GRO-seq data
GRO3 = makeGROdata("/home/home/jeremy/RNAseq/Glass/refseq_and_subsequent_transcripts_with_tag_counts.txt", refGene,"GRO_Data_test.RData",1000)


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
## 10/29/2012
######################################################################################################################################################
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
dat = calculateSplicingRunon_MATLAB('genes_for_analyticalModel_withRunon.txt',3000,30,conversionFactor<-'2,20,5,100,1,2')
#conversionFactor =  '0.5,5,5,300';dat =   read.table(sprintf("genes_for_analyticalModel_withRunon_splicingRatesRunon_3000.00_30.00_%s.txt",conversionFactor),col.names=c("UniqueID","PredSplice","PredSplice_strength","PredSplice_pause","PredSplice_runon","PredSplice_ALL")) 
dataWithChromGROpredictions = merge(dataWithChromGRO,dat, by='UniqueID',all.x=T)

######################################################################################################################################################
#### 11/29/2012

# Here I'm simulating 5 times, first with only the structures, then adding one by one the splicing, pausing, runon, then ALL together
plotFig4Like(dataWithChromGROpredictions, '_simple_SplitImprovementbyLength9','Group','last exon length','median')
plotFig4Like(dataWithChromGROpredictions, '_simple_SplitImprovementbyLengthGroup9','LengthGroup','last exon length','median')


############################################################################################################################################
### (12/3/2012): Putting the AVERAGE splicing of that category into the model, not the individual splicing (it's so noisy!)
####### Call the matlab code that simulates all genes:

strength_nonLast = mean(ACCEPTOR[ACCEPTOR[,2]==0,1],na.rm=T)
strength_lasts_avg = tapply(mergedData$Acceptor,mergedData$Group,mean)

##### For all genes with > 1 intron, model with uniform splicing throughout EXCEPT for the last intron.
# Series of merges to create this file
mergedTemp = merge(data.frame(strength_lasts_avg),mergedData[mergedData$FeatureCount > 2,], by.x=0, by.y = 'Group')
mergedTemp2 = data.frame(UniqueID = mergedTemp$UniqueID,SplicingStrength2= paste(sapply(mergedTemp$FeatureCount-2,function(x)paste(rep(strength_nonLast,x),collapse=",")),mergedTemp$strength_lasts_avg,sep=","))
temp3 = merge(temp2, mergedTemp2, by = 'UniqueID')
temp3$SpliceStrengths = temp3$SplicingStrength2
## Save the data to file
write.delim(temp3[,1:9], File='genesUniform_for_analyticalModel_withRunon.txt')

## Run simulation
dat2 = calculateSplicingRunon_MATLAB('genesUniform_for_analyticalModel_withRunon.txt',3000,30,conversionFactor<-'2,20,5,100,1,1')
#conversionFactor =  '2,20,5,100,1,2';dat2 =   read.table(sprintf("genes_for_analyticalModel_withRunon_splicingRatesRunon_3000.00_30.00_%s.txt",conversionFactor),col.names=c("UniqueID","PredSplice","PredSplice_strength","PredSplice_pause","PredSplice_runon","PredSplice_ALL")) 

##  PLOT FIGURE
dataWithChromGROpredictions = merge(dataWithChromGRO,dat2, by='UniqueID',all.x=T)
plotFig4Like(dataWithChromGROpredictions, '_simple_SplitImprovement_UniformByLength3newBoost','Group','last exon length','median')


################################
## Here I'm simulating by starting first with ALL uniform splicing (but for each gene, make the splicing Z be the same WITHIN the gene first, THEN do it separately)
dat3 = calculateSplicingRunon_MATLAB('genes_for_analyticalModel_withRunon.txt',3000,30,conversionFactor<-'2,20,5,100,1,2',VERSION=2)
conversionFactor =  '2,20,5,100,1,2';dat3 =   read.table(sprintf("genes_for_analyticalModel_withRunon_splicingRatesRunon2_3000.00_30.00_%s.txt",conversionFactor),col.names=c("UniqueID","PredSplice","PredSplice_strength","PredSplice_pause","PredSplice_runon","PredSplice_ALL")) 

##  PLOT FIGURE
dataWithChromGROpredictions3 = merge(dataWithChromGRO,dat3, by='UniqueID',all.x=T)
plotFig4Like(dataWithChromGROpredictions3, '_simple_SplitImprovementV2byLength1','Group','last exon length','median')


###############################################################
#### 12/4/12
# Going to output stuff for matlab so I can redo Fig 4 with the Z scores: also going to make Fig 5C, which uses the Read-through
lengths_last =  tapply(mergedData$length,mergedData$Group,mean)
strength_nonLast = mean(ACCEPTOR[ACCEPTOR[,2]==0,1],na.rm=T)
strength_lasts_avg = tapply(mergedData$Acceptor,mergedData$Group,mean)
energy_lasts_avg = tapply(mergedData$Energy,mergedData$Group,mean)
runon_lasts_avg = tapply(dataWithGRO$Runon,dataWithGRO$Group,mean,na.rm=T)

write.delim(cbind(lengths_last,strength_lasts_avg,strength_nonLast,energy_lasts_avg,runon_lasts_avg), 'Normalized_averages_byLength_withRunon.txt') 

#### 12/5/12
## Alex wants to look at "the ratio of measured last exon vs predicted last exon (= correction factor)"
dataWithGRO$correction_factor = (dataWithGRO$Runon +dataWithGRO$length) / dataWithGRO$length
dataWithGRO$correction_factor2 = (dataWithGRO$Runon+dataWithGRO$length)
dataWithGRO$correction_factor_random = sample(dataWithGRO$Runon,nrow(dataWithGRO)) / dataWithGRO$length

## First set of pictures Alex wants
png("GRO_figures1_%d.png",height=400,width=400)
CEX=1.1
par(cex=CEX)
#par(mfrow=c(1,2));hist(log10(dataWithGRO$length),main='last exons length ',xaxt='n', xlab='Distance (bp)'); axis(1,at=0:5,label=10^(0:5)); 
hist(log10(dataWithGRO$Runon), main= "Read-through length", xaxt='n', xlab='Distance (bp)'); axis(1,at=0:5,label=10^(0:5))

par(cex=CEX)
xlabs = tapply(dataWithGRO$length,dataWithGRO$Group,function(x)paste(range(x),collapse='-'))
plotSplits(dataWithGRO,column='correction_factor_random','Ratio of measured over predicted Last Exon Lengths\nfor each Last Exon bin',Group='Group',ylab='Correction Factor',xlab='', xaxt='n',type='b',lwd=2,col='red');
plotSplits(dataWithGRO,'correction_factor','',Group='Group',add=T,ylab='',xlab='', xaxt='n',type='b',lwd=2); axis(1,at=1:5,lab=xlabs,las=3)

emf("Fig5C.emf",width=5,height=5)
par(cex=CEX)
plotSplits(dataWithGRO,column='correction_factor','Correction Factor for each Last Exon bin',Group='Group',ylab='Ratio of measured over predicted Last Exon Lengths',xlab='', xaxt='n',type='b',lwd=2,col='red');
axis(1,at=1:5,lab=xlabs,las=3)
dev.off()

par(cex=CEX)
plotSplits(dataWithGRO,'correction_factor2','Total last exon length\nAs a function of predicted last exon length',Group='Group',ylab='Length(bp)',xlab='', xaxt='n',type='b',lwd=2,ylim=c(0,5500)); axis(1,at=1:5,lab=xlabs,las=3)
par(lty=2,cex=CEX)
plotSplits(dataWithGRO,column='length','',Group='Group',ylab='Length(bp)',xlab='', xaxt='n',type='b',add=T,lwd=2);
#plotSplits(dataWithGRO,'correction_factor2','Read-through+length/exon length\nAs a function of # exons',Group='Group',ylab='Correction Factor',xlab='', xaxt='n',type='b',lwd=2); axis(1,at=1:5,lab=xlabs,las=3)

dev.off()

## Some more figures, this time different size
png("GRO_figures2_%d.png",height=400,width=400)
CEX=0.8; CEX.LAB=1.5
## Scatter plot
par(cex=CEX, cex.lab=CEX.LAB,cex.axis=CEX.LAB)
plot(log10(dataWithGRO$length), log10(dataWithGRO$Runon + dataWithGRO$length),xlab='Predicted Last Exon Length',ylab='Measured Last Exon Length')
abline(a=0,b=1)
par(cex=CEX, cex.lab=CEX.LAB,cex.axis=CEX.LAB)
plot(log10(dataWithGRO$length), log10(dataWithGRO$Runon),xlab='Predicted Last Exon Length',ylab='Read-through Length')
abline(a=0,b=1)

par(cex=CEX*1.2, cex.lab=CEX.LAB/.8,cex.axis=CEX.LAB)
a=plot2Hists(list(log10(dataWithGRO$length),log10(dataWithGRO$Runon+dataWithGRO$length)),'Re-calculation of last exon lengths','Exon Length',leg=c('Predicted','Measured'),N=30,xaxt='n')
dev.off()

### Redo Figures: Black and grey Figure '6' (HK/nonHK)
source("/home/jeremy/ExonPipeline/ExonPipeline_analysis_and_plotting_Functions.r")

png("HK_fig_withGRO.png",height=300,width=1200)
layout(mat=matrix(1:6,nrow=1),widths=c(1,1,1,2,2,2))
par(mar=c(6,3,1,0),cex=1)
makeBoxplot (data.frame(Value=dataWithGRO$length, isLast=dataWithGRO$isHK), mainText ='Predicted Last Exons',T, NAMES=c('non-HK','HK'), addStats=F,at=1:4,labels=10^(1:4),notch=T,lwd=3, border=c("black","gray"))
makeBoxplot (data.frame(Value=dataWithGRO$Runon, isLast=dataWithGRO$isHK), mainText ='Read-through',T, NAMES=c('non-HK','HK'), addStats=F,at=1:4,labels=10^(1:4),notch=T,lwd=3, border=c("black","gray"))
makeBoxplot (data.frame(Value=dataWithGRO$length+dataWithGRO$Runon, isLast=dataWithGRO$isHK), mainText ='Measured Last Exons',T, NAMES=c('non-HK','HK'), addStats=F,at=1:4,labels=10^(1:4),notch=T,lwd=3, border=c("black","gray"))
#dev.off()
#
#png("HK_fig_withGRO2.png",height=300,width=1200)
#par(mfrow=c(1,4))
par(cex=0.8)
a=plot2Hists(list(log10(dataWithGRO$length[dataWithGRO$isHK==1]),log10(dataWithGRO$length[dataWithGRO$isHK==0])),'Predicted Last Exon Length','Exon Length',leg=c('HK','non-HK'),N=50)
a=plot2Hists(list(log10((dataWithGRO$Runon)[dataWithGRO$isHK==1]),log10((dataWithGRO$Runon)[dataWithGRO$isHK==0])),'Read-through','Exon Length',leg=c('HK','non-HK'),N=50)
a=plot2Hists(list(log10((dataWithGRO$length+dataWithGRO$Runon)[dataWithGRO$isHK==1]),log10((dataWithGRO$length+dataWithGRO$Runon)[dataWithGRO$isHK==0])),'Mesaured Last Exon Length','Exon Length',leg=c('HK','non-HK'),N=50)
dev.off()


#### 12/11/12
source("/home/jeremy/GitHub/useful_R/useful_R.r")
library(RColorBrewer)
library(devEMF)
COLS = col2rgb(brewer.pal(5, "Set1")[c(2,1)])

COLS2 =  c(rgb(COLS[1,1]/256,COLS[2,1]/256,COLS[3,1]/256,alpha=1), rgb(COLS[1,2]/256,COLS[2,2]/256,COLS[3,2]/256,alpha=1))

png("GRO_figure5.2.png",height=400,width=400)
emf("GRO_figure5.2.emf",height=6,width=6)
CEX=0.8; CEX.LAB=1.5
par(cex=CEX*1.2, cex.lab=CEX.LAB/.8,cex.axis=CEX.LAB)
a=plot2Densities(list(log10(dataWithGRO$length),log10(dataWithGRO$Runon+dataWithGRO$length)),'Re-calculation of last exon lengths','Exon Length',cols=COLS2,lwd=3,leg=c('Predicted','Measured'),xaxt='n',leg.at='topleft')
axis(1,at=1:6,lab=10^(1:6))
grid()
dev.off()




############################
## 12/12
# Redo Figure 6
source("/home/jeremy/ExonPipeline/ExonPipeline_analysis_and_plotting_Functions.r")

emf("Fig 6.emf",width=8,height=4)
#pdf("Fig 6.pdf",width=8,height=4)
#par(mfrow=c(1,4))
layout(mat=matrix(1:4,nrow=1),widths=c(2,2,3,3))

makeBoxplot (data.frame(Value=dataWithGRO$length, isLast=dataWithGRO$isHK), mainText='',T, NAMES=c('non-HK','HK'),ylab='Predicted Exon Length (bp)', addStats=F,at=1:4,labels=10^(1:4),notch=T,lwd=4, border=c("black","gray"))
makeBoxplot (data.frame(Value=dataWithGRO$Runon,  isLast=dataWithGRO$isHK), mainText='',T, NAMES=c('non-HK','HK'), ylab='Read-through (bp)', addStats=F,at=1:4,labels=10^(1:4),notch=T,lwd=4, border=c("black","gray"))

## Make colors' default changeable!
plotSplitsHK(mergedData,column='Energy','',Group='Group',ylab='Nucleosome Stability Score',xlab='', names=NULL,colors=c("gray","black"),lwd=4)
plotSplitsHK(mergedData,'Acceptor','',Group='Group',ylab='Splice Site Strength Score',xlab='', names=NULL,colors=c("gray","black"),lwd=4)
dev.off()




#### 12/17/2012
exonCategories = getFirstMiddleLast(exons)
exonLengths = lapply(exonCategories,function(x)log10(abs(x$end-x$start)))
emf("Fig 2D.emf")
boxplot(exonLengths,names=c('First','Middle','Last'),yaxt='n', ylab='Exon Lengths', xlab='', main='',las=3)
axis(2,at=0:6,lab=10^(0:6))
dev.off()



##### 12/20/2012
##################################################################################################################
################### REDO  MOUSE GRO-SEQ ###########################
##################################################################################################################

#
refGene = add_UniqueID(loadRefgene(settings))

# Ok, finally going to make this into a function
GRO3.test=makeGROdata("/home/home/jeremy/RNAseq/Glass/refseq_and_subsequent_transcripts_with_tag_counts.txt", refGene,"GRO_Data_test.RData",1000)

## Merge with 'mergedData'
dataWithGRO2 = merge(mergedData, GRO3 , by='UniqueID',all.x=T)
dataWithGRO$RunonGroup = splitByVector(dataWithGRO$Runon,quantile(dataWithGRO$Runon,na.rm=T,(1:4)/5)) # split into 5 groups
save(dataWithGRO,file='dataWithGRO.RData')


### Possible Redo of figures/data for paper:

# Summary of Runon:
summary(dataWithGRO$Runon[dataWithGRO$isHK==0])
summary(dataWithGRO$Runon[dataWithGRO$isHK==1])

# Summary of Final total length:
summary((dataWithGRO$length+dataWithGRO$Runon)[dataWithGRO$isHK==0])
summary((dataWithGRO$length+dataWithGRO$Runon)[dataWithGRO$isHK==1])

# Summary of lengths that were actually added to the Runon:
summary(dataWithGRO$length[!is.na(dataWithGRO$Runon)])



### 12/21/2012
# Writing a supplemental file of GRO-seq data 
GRO_out = data.frame(RefSeq=dataWithGRO$sequence_identifier, Gene=dataWithGRO$GeneName,Chr = dataWithGRO$chr, Strand=dataWithGRO$strand.x, PredictedExonLength=dataWithGRO$length,MeasuredExonLength=dataWithGRO$length + dataWithGRO$Runon, 
 PredictedTerminationLocation=dataWithGRO$end,MeasuredTerminationLocation=dataWithGRO$end + dataWithGRO$Runon*ifelse(dataWithGRO$strand.x=='+',1,-1))
GRO_out=GRO_out[!is.na(GRO_out[,1]),]

write.delim(GRO_out,'Supplemental_Table_2.txt')

#In the abstract we say that we have a database of termination sites.
#So maybe we should also show the location in the genome where we measure the termination site to be.
