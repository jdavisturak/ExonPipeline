rm(list=ls())
library(gplots)
library(RColorBrewer)
library(Hmisc)

source("/home/jeremy/ExonPipeline/ExonPipeline_functions.r")
source("/home/jeremy/ExonPipeline/ExonPipeline_analysis_and_plotting_Functions.r")
source("/home/home/jeremy/Code/useful_R/useful_R.r")

setwd("/home/jeremy/ExonPipeline/mm9")

load("PlottedData_new3.RData")
load("PlottedData_withGRO-redone_new3.RData")     

###load("GRO_Data.RData")
### 9/11/13
### Use Karmel's redone GRO-seq algorithm 
# GRO3 =  makeGROdata2(inputFile="/home/RNAseq/Glass/post_gene_transcripts.txt", subset(refseq, name %in% refseq.noPseudogenes$name),outputFile = "GRO_Data_redone.RData", geneDistance=1000)
# dataWithGRO = merge(mergedData, GRO3 , by='UniqueID',all.x=T)
# dataWithGRO$RunonGroup = quantileGroups(dataWithGRO$Runon,5)
# save(dataWithGRO, file="PlottedData_withGRO-redone_new3.RData")


############################################################################
#dataWithGRO$GeneLengthGroup = quantileGroups(dataWithGRO$GeneLength,4)

dataWithGRO2 = subset(dataWithGRO, !duplicated(UniqueID))
dataWithGRO2$GeneLengthGroup = quantileGroups(dataWithGRO2$GeneLength,5)
dataWithGRO2$ExonGroup = quantileGroups(dataWithGRO2$exonCount,5)
dataWithGRO2$logRunon = log10(dataWithGRO2$Runon)             

dataWithGRO2$GeneSize_hg19 = splitByVector(dataWithGRO2$GeneLength,c(6442, 20252, 57233 ))


summary(dataWithGRO2$Runon) # this is the new data 9/11/13                      
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    231    2441    3925    4911    6275   35250   13036 
summary(subset(dataWithGRO2,isHK==0)$Runon)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    231    2308    3744    4769    6052   35250   12115 
summary(subset(dataWithGRO2,isHK==1)$Runon)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    242    3052    4512    5453    6845   29870     921

t.test(Runon~isHK,dataWithGRO2)
wilcox.test(Runon~isHK,dataWithGRO2)


############################################################################
############################################################################
#### 10/28/13:   
#### Filtering genes by those that have canonical poly(A) usage in Cytoplasm
############################################################################
############################################################################
# 4568 Genes to start
load("/home/RNAseq/Smale/SmaleCyt_cleavedFractions_refseq.RData")

# I got rid of genes that I could determine to be non-canonically poly-adenylated: had at least 5 reads and less than 95% cleavage ratio
dataWithGRO3 = subset(  dataWithGRO2,!(UniqueID  %in% rownames(subset(AllCleavedAverages,Cyt < 0.95))))
length(unique(subset(dataWithGRO3,!is.na(Runon))$UniqueID)) # 3594 genes
save(dataWithGRO3, file='Mouse_dataWithGRO3_polyAfilt.RData')
wilcox.test(Runon~isHK,dataWithGRO3) # very signif

summary(subset(dataWithGRO3,isHK==0)$Runon) 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    285    2223    3561    4448    5597   35040   11107 
summary(subset(dataWithGRO3,isHK==1)$Runon)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    464    3050    4422    5203    6435   29870     780 
summary(dataWithGRO3$Runon)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    285    2353    3766    4609    5853   35040   11887 

########################## PLOTS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

plot.dev("Mouse_GRO-redone_polyAfilt.pdf",height=2.5, width=3.1*(4.1/3.1))
layout(mat=matrix(1:3,nc=3),widths=c(1.8,1.3,1.1))
par(mai=c(0.7,0.6,.2,0.1),cex=0.7)

# Instead of using the gene size categories based on mouse, use those from human

boxplot(log10(Runon) ~ GeneSize_hg19, data=dataWithGRO3, ylim=c(2,5),pch=19,names=c('Short','Med Short','Med Long','Long'), ylab='Read-through (bp)', main='', notch=T,yaxt='n',las=2)
axis(2,at=1:5,lab=sprintf("%d",(10^(1:5))))
#axis(2,at=1:4,lab=c(expression(10^1),expression(10^2),expression(10^3),expression(10^4)),las=2)

## Housekeeping
#boxplot(log10(Runon) ~ isHK, dataWithGRO3, ylim=c(2,5),border=c('black','gray'), ylab='Read-through (bp)',notch=T,pch=19,lwd=1,names=c("Regulated","HK"),las=2,yaxt='n')
boxplot(log10(Runon) ~ isHK, dataWithGRO3,outline=F, staplewex=0,range=0.00000001,border=1:2, ylab='Read-through (kb)',notch=T,pch=19,lwd=1,names=c("Regulated","HK"),las=2,yaxt='n')
#ticks = c(6000,4000,3000,2500)
ticks=1000*6^c(1/3,1/2,2/3,5/6,1)
axis(2,at=log10(ticks),lab=sprintf("%.1f",ticks/1000))

boxplot(log10(dataWithGRO3$Runon) , ylim=c(2,5), ylab='Read-through (bp)',notch=T,pch=19,lwd=1,yaxt='n')
axis(2,at=1:5,lab=sprintf("%d",(10^(1:5))))

plot.off()

plot.dev("Mouse_GRO-redone_density_polyAfilt.pdf",'pdf',height=2.5, width=2.5)
par(mai=c(0.8,.7,.5,.1),mgp=c(3,.5,0))
plot(a<-density(log10(dataWithGRO3$Runon) ,na.rm=T),xlim=c(2,5), xlab='Read-through (bp)', xaxt='n',yaxt='n',main='',cex.main=0.7,lwd=2,ylab='',plot.first=grid())
axis(1,at=1:5,lab=sprintf("%d",(10^(1:5))),las=2)
title(ylab='density',mgp=c(0.1,0,0))
abline(v=V<-median(log10(dataWithGRO3$Runon),na.rm=T),lty=2)
axis(1,at=V,lab=round(10^V),las=2)
plot.off()


### plot Read-through based on my 19 categories ....

par(mfrow=c(2,2))
invisible(sapply(1:4, function(size){
  
  boxplot(log10(Runon) ~ intronCountGroups, data=subset(dataWithGRO3, GeneSize_hg19==size), ylim=c(2,5),pch=19, ylab='Read-through (bp)', notch=T,yaxt='n',las=2, main=c('Short','Med Short','Med Long','Long')[size])
  axis(2,at=1:5,lab=sprintf("%d",(10^(1:5))))

}))



















## Make nice plotting function for splitting things up into 2 categories and plotting the medians of the response variable
split2medians = function(myData,column=column, col1='Dist2End', col2='splicingAMinus', quants2 = 4,numMedians=6, split2=quantileGroups(myData[,col2], quants2), ...){
  numSplits2 <- length(unique(split2));
  print(tapply(myData[,col2],split2,range))
  
  medianData = sapply(unique(sort(split2)), function(x)getQuantilMedians(myData[split2==x,], col1=col1,column=column,numMedians),simplify=F)
  for (x in 1:length(medianData)){    
    names(medianData[[x]]) <- paste(c(col2,col1),x,sep='')
    if(x==1)
      out = medianData[[x]]
    else out = cbind(out,medianData[[x]])
  } 
  out 
}


######
#plot2medians = function(out,cols=brewer.pal(ncol(out)/2,'Set1'),pchs=rep(19,ncol(out)/2), xlab='', ylab='', xaxt='s', main='', ylim= range(out[,seq(2,ncol(out),2)]), xlim= range(out[,seq(1,ncol(out),2)]), ...){
#  for( i in 1:ncol(out)/2){
#    if(i ==1) 
#      plot(out[,1], out[,2],xlab=xlab,ylab=ylab,xaxt=xaxt,main=main,xlim=xlim,ylim=ylim,col=cols[i],pch=pchs[i])
#    else
#      points(out[,2*i-1], out[,2*i], col=cols[i])  
#  }
#}

plot2medians = function(out,pchs=rep(19,ncol(out)/2), cols=brewer.pal(ncol(out)/2,'Set1'), xlab='', ylab='', xaxt='s', main='', ...){
   matplot(out[,seq(1,ncol(out),2)], out[,seq(2,ncol(out),2)], pch=pchs, col=cols,xlab=xlab,ylab=ylab,xaxt=xaxt,main=main, ...)
}         

groNon = split2medians(subset(dataWithGRO2, isHK==0),'Runon',col1='exonCount',col2='GeneLength', quants2=4, numMedians=5)
groHK = split2medians(subset(dataWithGRO2, isHK==1),'Runon',col1='exonCount',col2='GeneLength', quants2=4, numMedians=5)
groNon2 = split2medians(subset(dataWithGRO2, isHK==0),'Runon',col2='exonCount',col1='GeneLength', quants2=4, numMedians=10)
groHK2 = split2medians(subset(dataWithGRO2, isHK==1),'Runon',col2='exonCount',col1='GeneLength', quants2=4, numMedians=10)

par(mfrow=c(2,2))
plot2medians(groHK, main='HK',xlab='Number of introns', ylab='Runon Distance')
plot2medians(groNon, main='Regulated',xlab='Number of introns', ylab='Runon Distance')
legend('right',pch=19,col=brewer.pal(4,'Set1'),leg=c('short genes','medium short','medium long','long genes'))
plot2medians(groHK2, main='HK', xlab='Gene Length', ylab='Runon Distance')
plot2medians(groNon2, main='Regulated',xlab='Gene Length', ylab='Runon Distance')
legend('right',pch=19,col=brewer.pal(4,'Set1'),leg=c('few introns','medium few introns','medium many introns','many introns'))




###########################################################################################################
######### Do GRO-seq traces over poly-A regions 
### 9/30/2013
############################################################################################################
#library(rtracklayer)
#library(RColorBrewer)
#source("/home/jeremy/ExonPipeline/analyzeBigWig.r")
#load('ENCODE_PolCHIP_250.RData')
#


######################################################################################################################################################
### 10/18/13: Use Kristyn's RNAseq data to look at where actual poly-A locations are
######################################################################################################################################################


require(rtracklayer)

bmdm_genes = import.bed('/home/home/jeremy/RNAseq/Kristyn/AllCufflinks/Kristyn_transcripts.bed')


# Overlap with refseq genes
ref.GRange = GRanges(refseq$chrom, IRanges(refseq$txStart,refseq$txEnd), refseq$strand, UniqueID=refseq[,'UniqueID'])
ref.GRange = ref.GRange[ref.GRange$UniqueID %in% mergedData$UniqueID,]

getExtraTxn = function(bed,ref.GRange){
  # Find overlaps: do strands separately               
  overlap_plus = as.matrix(findOverlaps(ref.GRange[strand(ref.GRange)=='+',],bed[bed$strand=='+',]))
  overlap_minus = as.matrix(findOverlaps(ref.GRange[strand(ref.GRange)=='-',],bed[bed$strand=='-',]))
  # Combine the overlaps (ugly subsetting operation, b/c these are indices!)
  overlap2 = rbind(cbind(which(strand(ref.GRange)=='+')[overlap_plus[,1]], which(bed$strand=='+')[overlap_plus[,2]]),
    cbind(which(strand(ref.GRange)=='-')[overlap_minus[,1]], which(bed$strand=='-')[overlap_minus[,2]]))
  #dups = subset(as.data.frame(overlap2), V1 %in% V1[duplicated(V1)]) # genes that have > 1 hit
  
  # Use the indices to subset for the actual genes
  myRanges_plus = ref.GRange[strand(ref.GRange)=='+',][overlap_plus[,1],]
  myRanges_minus = ref.GRange[strand(ref.GRange)=='-',][overlap_minus[,1],]
  
  # Positive strand: find the END of the overlapping beds, then take the MAX bed for each gene, and save the name of the .bed 
  myRanges_plus$END2 = end(bed[bed$strand=='+',])[overlap_plus[,2]]
  myRanges_plus$NAME = bed[bed$strand=='+',]$name[overlap_plus[,2]]
  p = tapply((myRanges_plus$END2), as.factor(myRanges_plus$UniqueID),max)
  p.list = (myRanges_plus$END2)
  names(p.list) = myRanges_plus$NAME
  ff2=function(x){names(x)[which.max(x)]}
  p.name = tapply(p.list,as.factor(myRanges_plus$UniqueID),ff2)
  
  # negative strand: find the START of the overlapping beds, then take the MIN bed for each gene, and save the name of the .bed 
  myRanges_minus$START2 = start(bed[bed$strand=='-',])[overlap_minus[,2]]
  myRanges_minus$NAME = bed[bed$strand=='-',]$name[overlap_minus[,2]]
  m = tapply((myRanges_minus$START2), as.factor(myRanges_minus$UniqueID),min)
  m.list = (myRanges_minus$START2)
  names(m.list) = myRanges_minus$NAME
  ff2=function(x){names(x)[which.min(x)]}
  m.name = tapply(m.list,as.factor(myRanges_minus$UniqueID),ff2)
  
  
  extraTxn_plus = data.frame(UniqueID = names(p), extra = p - refseq$txEnd[match(names(p),refseq$UniqueID)], bedName=p.name)
  extraTxn_minus = data.frame(UniqueID = names(m), extra = refseq$txStart[match(names(m),refseq$UniqueID)] - m, bedName=m.name)
  extraTxn = rbind(extraTxn_plus,extraTxn_minus)
  extraTxn
}

extra1 = getExtraTxn(bmdm_genes,ref.GRange)
export.bed(bmdm_genes[bmdm_genes$name %in% extra1$bedName,], '/home/home/jeremy/RNAseq/Kristyn/AllCufflinks/Kristyn_transcripts_overGenes.bed')


# extra1 now defines presumtive poly-A regions for all expressed genes, according to Cufflinks annotations

######################################################################################################################################################
###  Examine polyA locations and GRO-seq locations for ATAAA signals (use refseq as a positive control?)
######################################################################################################################################################

## Functions to get polyA signals:

# First obtain regions of interest and strands
refseq_ends = merge(dataWithGRO2,extra1,by='UniqueID',all.x=T)
refseq_ends$END = refseq_ends$txEnd
refseq_ends$END[refseq_ends$strand.x=='-'] <- refseq_ends$txStart[refseq_ends$strand.x=='-']
refseq_ends$END_GRO = refseq_ends$post_gene_end
refseq_ends$END_GRO[refseq_ends$strand.x=='-'] <- refseq_ends$post_gene_start[refseq_ends$strand.x=='-']
refseq_ends$END_CUFF = refseq_ends$END + refseq_ends$extra
refseq_ends$END_CUFF[refseq_ends$strand.x=='-'] <- refseq_ends$END[refseq_ends$strand.x=='-'] - refseq_ends$extra[refseq_ends$strand.x=='-'] 

# Control: use TSS
refseq_ends$START = refseq_ends$txStart
refseq_ends$START[refseq_ends$strand.x=='-'] <- refseq_ends$txEnd[refseq_ends$strand.x=='-']


# inputs: [genome], chromosomes, range, strand
# Then take those regions and get sequence up/down stream: input=widths

getRegions = function(chrs,locations,strand,UP=0,DOWN=0){
  sapply(1:length(chrs), function(i){
   print(i)
#   browser()
    D = subseq(genome[[chrs[i]]], locations[i] - ifelse(strand[i]=='+',UP,DOWN), locations[i] + ifelse(strand[i]=='+',DOWN,UP))
    if(strand[i]=='+') D<-reverseComplement(D)
    D
  })
}


### Subset only for genes that have both measurements (Runon and expression)
refseqs_ends_complete = subset(refseq_ends,!is.na(Runon) & !is.na(extra))

# Retreive sequences upstream 
seqs_ref <- getRegions(refseqs_ends_complete$chr, refseqs_ends_complete$END, refseqs_ends_complete$strand.x, 200 ,0)
seqs_rnaseq <- getRegions(refseqs_ends_complete$chr, refseqs_ends_complete$END_CUFF, refseqs_ends_complete$strand.x, 200 ,0)
seqs_gro <- getRegions(refseqs_ends_complete$chr, refseqs_ends_complete$END_GRO, refseqs_ends_complete$strand.x, 200 ,0)
seqs_tss <- getRegions(refseqs_ends_complete$chr, refseqs_ends_complete$START, refseqs_ends_complete$strand.x, 0,200)

seqs_ref_down <- getRegions(refseqs_ends_complete$chr, refseqs_ends_complete$END, refseqs_ends_complete$strand.x, 0,200)
seqs_rnaseq_down <- getRegions(refseqs_ends_complete$chr, refseqs_ends_complete$END_CUFF, refseqs_ends_complete$strand.x, 0 ,200)
seqs_gro_down <- getRegions(refseqs_ends_complete$chr, refseqs_ends_complete$END_GRO, refseqs_ends_complete$strand.x, 0,200)



# Random control
chrs = names(genome)[1:21]
lens = sapply(chrs,function(x)length(genome[[x]]))

seqs_random <- getRegions(rnd.chrs<-sample(chrs,nrow(refseqs_ends_complete),replace=TRUE), sapply(rnd.chrs,function(x)(sample(lens[[x]]-10000,1)+5000)), refseqs_ends_complete$strand.x, 0,200)
seqs_random2 <- getRegions(rnd.chrs<-sample(chrs,nrow(refseqs_ends_complete),replace=TRUE), sapply(rnd.chrs,function(x)(sample(lens[[x]]-10000,1)+5000)), refseqs_ends_complete$strand.x, 0,200)


# Then look for ATAAA
ATAAA = DNAString("ATAAA")
matches_ref = sapply(seqs_ref,function(x) as.matrix(ranges(matchPattern(pattern=ATAAA, x)))[,1])
matches_gro = sapply(seqs_gro,function(x) as.matrix(ranges(matchPattern(pattern=ATAAA, x)))[,1])
matches_rnaseq = sapply(seqs_rnaseq,function(x) as.matrix(ranges(matchPattern(pattern=ATAAA, x)))[,1])
matches_tss = sapply(seqs_tss,function(x) as.matrix(ranges(matchPattern(pattern=ATAAA, x)))[,1])
matches_random = sapply(seqs_random,function(x) as.matrix(ranges(matchPattern(pattern=ATAAA, x)))[,1])
matches_random2 = sapply(seqs_random2,function(x) as.matrix(ranges(matchPattern(pattern=ATAAA, x)))[,1])
matches_ref_down = sapply(seqs_ref_down,function(x) as.matrix(ranges(matchPattern(pattern=ATAAA, x)))[,1])
matches_gro_down = sapply(seqs_gro_down,function(x) as.matrix(ranges(matchPattern(pattern=ATAAA, x)))[,1])
matches_rnaseq_down = sapply(seqs_rnaseq_down,function(x) as.matrix(ranges(matchPattern(pattern=ATAAA, x)))[,1])

# Retreive sequences upstream 
matches_ref100 = sapply(seqs_ref,function(x) as.matrix(ranges(matchPattern(pattern=ATAAA, subseq(x,100,200))))[,1])
matches_gro100 = sapply(seqs_gro,function(x) as.matrix(ranges(matchPattern(pattern=ATAAA, subseq(x,100,200))))[,1])
matches_rnaseq100 = sapply(seqs_rnaseq,function(x) as.matrix(ranges(matchPattern(pattern=ATAAA, subseq(x,100,200))))[,1])
matches_tss100 = sapply(seqs_tss,function(x) as.matrix(ranges(matchPattern(pattern=ATAAA, subseq(x,100,200))))[,1])
matches_random100 = sapply(seqs_random,function(x) as.matrix(ranges(matchPattern(pattern=ATAAA, subseq(x,100,200))))[,1])
matches_random2100 = sapply(seqs_random2,function(x) as.matrix(ranges(matchPattern(pattern=ATAAA, subseq(x,100,200))))[,1])


# Stats
ATAAA_stats = rbind( (stats200<-cbind(length(which(lapply(matches_ref,length) > 0)),
length(which(lapply(matches_rnaseq,length) > 0)),
length(which(lapply(matches_gro,length) > 0)),
length(which(lapply(matches_tss,length) > 0)),
length(which(lapply(matches_random,length) > 0)),
length(which(lapply(matches_random2,length) > 0)),
length(which(lapply(matches_ref_down,length) > 0)),
length(which(lapply(matches_rnaseq_down,length) > 0)),
length(which(lapply(matches_gro_down,length) > 0))))[,1:6],

cbind(length(which(lapply(matches_ref100,length) > 0)),
length(which(lapply(matches_rnaseq100,length) > 0)),
length(which(lapply(matches_gro100,length) > 0)),
length(which(lapply(matches_tss100,length) > 0)),
length(which(lapply(matches_random100,length) > 0)),
length(which(lapply(matches_random2100,length) > 0))))

dimnames(ATAAA_stats) = list(c("200 bp","100 bp"),c("Refseq","rnaseq: cufflinks","GRO","refseq tss","random","random2"))
colnames(stats200) = c("Refseq","rnaseq: cufflinks","GRO","refseq tss","random","random2",'Refseq downstream','cufflinks downstream','GRO downstream')


######################################################################################################################################################
###  Summary of this investigation:

# The median read-through is 3.8 kb
summary(dataWithGRO2$Runon)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    209    2334    3840    4825    6144   75960   10769 


# But the median end as determined by Cufflinks is about 0, meaning the poly(A) site COULD be on average right around the TSS
summary(extra1[,2])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-215100    -318     -72    1938     112  394400 

# AND if we look at how many of these regions have a canonical ATAAA site within 200 bp upstream, the rnaseq blocks have more 
ATAAA_stats   
 #     Refseq rnaseq: cufflinks GRO refseq tss random random2
#200 bp    673               503 365         56    464     456
#100 bp    320               255 190         23    258     243

stats200[1,c(1,7,2,8,3,9,4:6)]
#              Refseq    Refseq downstream    rnaseq: cufflinks cufflinks downstream                  GRO       GRO downstream 
#                 673                  542                  503                  567                  365                  425 
#          refseq tss               random              random2 
#                  56                  464                  456 

 

              
######################################################################################################################################################              






