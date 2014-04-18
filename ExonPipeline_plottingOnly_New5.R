### Script to plot data after pipeline has run
# **********************************************************************************************************************************
# Major update 7/31/2013
# 10/9/13: Updated to normalize data again AFTER filtering
#**********************************************************************************************************************************

rm(list=ls())
library(gplots)           
library(RColorBrewer)
library(Hmisc)
source("/home/jeremy/ExonPipeline/ExonPipeline_functions.r")
source("/home/jeremy/ExonPipeline/ExonPipeline_analysis_and_plotting_Functions.r")
source("/home/home/jeremy/Code/useful_R/useful_R.r")

settings = setup(getwd())
refseq<-add_UniqueID(loadRefgene(settings))
features = loadFeatures(genome=loadGenome(settings<-setup( getwd())),settings=settings,refseq)

# 1) Comparisons of ALL to LAlsDOST
cat("	Reading in data to plot:\n")
#  Get lengths
exons = removeDuplicateStartANDEnds(read.delim("exons.txt",stringsAsFactors=F))
LENGTH = cbind(Value = abs(exons$start-exons$end), isLast = exons$isLast)


# Get Nucleosome energy
energyData = read.delim("exons_forSeqToNuc_constrained_all_mean.txt",head=F,stringsAsFactors=F)
ENERGY = formatData(energyData)

## Get splice site strength
donorData = read.delim(settings$DonorsNorm,head=F,stringsAsFactors=F)
DONOR = formatData(donorData)

## Get splice site strength
acceptorData = read.delim(settings$AcceptorsNorm,head=F,stringsAsFactors=F)
ACCEPTOR = formatData(acceptorData)
 
bins = c(1,2,5,10) 
HKlist = read.delim("/home/jeremy/ExonPipeline/hg19/HouseKeepingGenes.txt",skip=1)
# mergedData = mergeLastData(exons,energyData,donorData,acceptorData,bins*180,HKlist)

## Get rid of pseudo-genes 
if (!file.exists('refGene_notPseudo.txt')){
  # Make the new file
  if (!file.exists('pseudogenes.txt')){
    # need to download ... try downloading from pseudogene.org
    stop("Try downloading the pseudogenes file from pseudogene.org, and rename it to pseudogenes.txt")
  }
  
  print("making new refGene_notPseudo.txt file ... ")
  system("awk 'NR > 1{printf \"chr%s\\t%s\\t%s\\t%s\\t0\\t%s\\n\",$2,$3,$4,$1,$5 }' pseudogenes.txt > pseudogenes.bed")
  system("awk 'NR >1{printf \"%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n\",$3,$5,$6,$2,$12,$4,$13,$9,$10,$11 }' refGene.txt > refGene.bed")
  system("awk 'NR == 1{printf \"%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n\",$3,$5,$6,$2,$12,$4,$13,$9,$10,$11 }' refGene.txt > refGene_notPseudo.txt")
  system("bedtools intersect -a refGene.bed -b pseudogenes.bed -s -v >> refGene_notPseudo.txt")
  
  system("awk -F \"\\t\" 'NR >1 {printf \"%s\\t%s\\t%s\\t%d\\t%d\\t%d\\t%d\\t%d\\t%s\\t%s\\n\",$4,$1,$6,$2,$3,$2,$3,$8,$9,$10 }' refGene_notPseudo.txt > refGene_notPseudo.genePred")
  system("genePredToGtf 'file' refGene_notPseudo.genePred refGene_notPseudo.gtf")
  system("sh ~/Code/bam2ssj/gtf2cps.sh refGene_notPseudo.gtf > refGene_notPseudo.cps")
    
}

refseq.noPseudogenes = read.delim("refGene_notPseudo.txt")
nonPseudoIDs = subset(refseq, paste(chrom,txStart,txEnd) %in% paste(refseq.noPseudogenes$chr, refseq.noPseudogenes$txStart,refseq.noPseudogenes$txEnd),'UniqueID')[,1]


keepNames = intersect(nonPseudoIDs, names(which(tapply(features$introns$FeatureCount,features$introns$UniqueID,min)==1)))
introns2= mostDownstream(subset(features$introns, UniqueID %in% keepNames))
exons2= mostDownstream(subset(features$exons,UniqueID %in% keepNames))   
introns2$FeatureCount = as.numeric(introns2$FeatureCount)

introns2$length = abs(introns2$start-introns2$end)

## Figure out the end of the gene (poly A).
introns2$terminus = introns2$end * ifelse(introns2$strand=='-',-1,1)
introns2$begin = introns2$start * ifelse(introns2$strand=='-',-1,1)
myMatch = match(sort(unique(introns2$UniqueID)), refseq$UniqueID)

TXINFO = refseq[myMatch,c("UniqueID","txStart","txEnd","exonCount")]; TXINFO$MaxFeature = TXINFO$exonCount - 1

## Merge with exons
introns2$energy_ID = paste(">",paste(introns2$UniqueID,introns2$FeatureCount+1,introns2$isLast,sep='_'),sep='')
introns2$acceptor_ID = paste(introns2$UniqueID,introns2$FeatureCount,sep='_')
introns2$donor_ID = paste(introns2$UniqueID,introns2$FeatureCount,sep='_')
exons2$energy_ID = paste(">",paste(exons2$UniqueID,exons2$FeatureCount,exons2$isLast,sep='_'),sep='')
exons2$exonLength = abs(exons2$start-exons2$end)
introns2$exonLength = exons2$exonLength[match(introns2$energy_ID, exons2$energy_ID)]

energyData2 = energyData; names(energyData2) = c("energy_ID","Stability")
donorData2 = donorData; donorData2$ID = sub("_([0|1])$","",donorData2$V1) ; names(donorData2) = c("","Donor","donor_ID") 
acceptorData2 = acceptorData; acceptorData2$ID = sub("_([0|1])$","",acceptorData2$V1) ; names(acceptorData2) = c("","Acceptor","acceptor_ID")

# Merge all data together (Merge INTRONS this time!)
mergedData = merge(merge(merge(merge(introns2,energyData2,by='energy_ID'),donorData2[,-1],by='donor_ID',all=T),acceptorData2[,-1],by='acceptor_ID',all=T),TXINFO,by='UniqueID')

## New normalizations
#mergedData$Stability = normalize(mergedData$Stability) # MY function, not built-in function
mergedData$Acceptor = normalize(mergedData$Acceptor) # MY function, not built-in function

mergedData$Dist2End = abs(mergedData$txEnd - mergedData$end)
mergedData$Dist2End[mergedData$strand=='-'] = abs(mergedData$txStart[mergedData$strand=='-'] - mergedData$end[mergedData$strand=='-'])
mergedData$Dist2Start = abs(mergedData$txStart - mergedData$start) + 1
mergedData$Dist2Start[mergedData$strand=='-'] = abs(mergedData$txEnd[mergedData$strand=='-'] - mergedData$start[mergedData$strand=='-'])
mergedData$Group = quantileGroups(mergedData$Dist2End,100)
mergedData$exonLength[mergedData$isLast==1] <- mergedData$Dist2End[mergedData$isLast==1]



## Add HK status
HKlist = read.delim("/home/jeremy/ExonPipeline/hg19/HouseKeepingGenes.txt",skip=1)
mergedData$isHK = toupper(mergedData$GeneName) %in% HKlist$Gene.name
mergedData$GeneLength = mergedData$txEnd - mergedData$txStart
mergedData$ExonsRemaining = mergedData$MaxFeature - mergedData$FeatureCount +1
mergedData$intronsRemaining = mergedData$MaxFeature - mergedData$FeatureCount

## For splitting up by length of gene
refseq$GeneLength = refseq$txEnd - refseq$txStart
qr = round(quantile(refseq$GeneLength,c(1:3)/4))
sizeCutoff = qr[1]; sizeCutoff2 = qr[2]; sizeCutoff3 = qr[3]
refseq$GeneSize = splitByVector(refseq$GeneLength,qr)
refseq$intronCountGroups = splitByVector(refseq$exonCount-1,c(2.5,4.5,7.5,19.5))
mergedData$GeneSize = splitByVector(mergedData$GeneLength,qr)
mergedData$intronCountGroups = splitByVector(mergedData$MaxFeature,c(2.5,4.5,7.5,19.5))
mergedData$intronLengthGroups = quantileGroups(mergedData$length,4)

save(mergedData,file="PlottedData_new5.RData")

