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
mergedData$Stability = normalize(mergedData$Stability) # MY function, not built-in function
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

save.image("PlottedData_new3.RData")

cat("	Plotting data:\n") 


 
## Plot basic length of exon thing
exons2$isLast=as.numeric(exons2$isLast)
exons2$exonLocationGroups = exons2$isLast + 1 - 1*(exons2$FeatureCount==1)
pdf(sprintf('%s_ExonLengths.pdf',settings$CommonName),height=3, width=3)
par(cex=0.5, mai=c(0.3,.5,.1,.1))
boxplot(log10(abs(start-end)) ~ exonLocationGroups, exons2, outline=F,notch=T, yaxt='n',names=c('First','Internal','Last'), main='',ylab=sprintf("Exon lengths in\n%s genome",settings$CommonName))
axis(2,at=1:4,lab=c(expression(10^1),expression(10^2),expression(10^3),expression(10^4)))
dev.off()

#### 1) Plot Barplots
source("/home/jeremy/ExonPipeline/ExonPipeline_analysis_and_plotting_Functions.r")
pdf(sprintf('%s_Barplots7.pdf',settings$CommonName),height=9,width=2)
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

#### OUTPUT enough data to simulate in MATALB
#lengths_last =  tapply(mergedData$length,mergedData$LengthGroup_S,mean)
#strength_nonLast = mean(ACCEPTOR[ACCEPTOR[,2]==0,1],na.rm=T)
#strength_lasts_avg = tapply(mergedData$Acceptor,mergedData$LengthGroup_S,mean)
#energy_lasts_avg = tapply(mergedData$Energy,mergedData$LengthGroup_S,mean)
#
#write.delim(cbind(lengths_last,strength_lasts_avg,strength_nonLast,energy_lasts_avg), 'Normalized_averages_byLength_6.txt') 
#
#

cat("Writing .BED files:\n")            

refseq$Ends = ifelse(refseq$strand=='+',refseq$txEnd,refseq$txStart)

UP=250
DOWN=250

mergedData$exonStart = mergedData$start - UP
mergedData$exonEnd = mergedData$start + DOWN
dat.Pos = subset(mergedData,strand=='+')
dat.Neg = subset(mergedData,strand=='-')

write.delim(subset(dat.Pos,isLast==0 ,c('chr','exonStart','exonEnd','energy_ID','FeatureCount','strand')),File=sprintf('%s_exons_-%dTo%d_nonLast+.bed',settings$CommonName,UP,DOWN),col.names=F)
write.delim(subset(dat.Neg,isLast==0 ,c('chr','exonStart','exonEnd','energy_ID','FeatureCount','strand')),File=sprintf('%s_exons_-%dTo%d_nonLast-.bed',settings$CommonName,UP,DOWN),col.names=F)
write.delim(subset(dat.Pos,isLast==1 ,c('chr','exonStart','exonEnd','energy_ID','FeatureCount','strand')),File=sprintf('%s_exons_-%dTo%d_last+.bed',settings$CommonName,UP,DOWN),col.names=F)
write.delim(subset(dat.Neg,isLast==1 ,c('chr','exonStart','exonEnd','energy_ID','FeatureCount','strand')),File=sprintf('%s_exons_-%dTo%d_last-.bed',settings$CommonName,UP,DOWN),col.names=F)


## Make .bed files for poly A region

for (intervals in list(
  data.frame(UP=1000,DOWN=1000), data.frame(UP=1000,DOWN=5000) 
  )){
  UP=  intervals$UP
  DOWN=intervals$DOWN
  refseq$polyAstart = refseq$Ends - UP   * ifelse(refseq$strand=='+',  1, -1)
  refseq$polyAend =   refseq$Ends + DOWN * ifelse(refseq$strand=='+', 1,  -1)
  
  write.delim(subset(refseq,strand=='+' ,c('chrom','polyAstart','polyAend','UniqueID','exonCount','strand')),File=sprintf('%s_polyA_-%dTo%d+.bed',settings$CommonName,UP,DOWN),col.names=F)
  write.delim(subset(refseq,strand=='-' ,c('chrom','polyAend','polyAstart','UniqueID','exonCount','strand')),File=sprintf('%s_polyA_-%dTo%d-.bed',settings$CommonName,UP,DOWN),col.names=F)
}  


### This way is different:
### Do this for /Volumes/bigmax.ucsd.edu/Code/ExonPipeline/ExonPipeline_functions.rpoly(A) ends that I mean to use for bigWigAverageOverBed
UP=DOWN=100

refseq$polyAstart = refseq$Ends - UP   * ifelse(refseq$strand=='+', 1, -1)
refseq$polyAend =   refseq$Ends + DOWN * ifelse(refseq$strand=='+', 1, -1)

## write them with names 'up' & 'down'
refseq$PolyA_upName = paste(refseq$UniqueID,'_Up',sep='')
refseq$PolyA_downName = paste(refseq$UniqueID,'_Down',sep='')

refseq2 = subset(refseq, (strand=='+' & polyAstart > 1) | (strand=='-' & polyAend > 1))

# First write 'up' then down
polyAFile = sprintf('%s_polyA_UpDown_%dTo%d.bed',settings$CommonName,UP,DOWN)
write.delim(subset(refseq2,strand=='+' ,c('chrom','polyAstart','Ends','PolyA_upName','exonCount','strand')),File=polyAFile,col.names=F)
write.delim(subset(refseq2,strand=='-' ,c('chrom','Ends','polyAstart','PolyA_upName','exonCount','strand')),File=polyAFile,col.names=F,append=T)
write.delim(subset(refseq2,strand=='+' ,c('chrom','Ends','polyAend','PolyA_downName','exonCount','strand')),File=polyAFile,col.names=F,append=T)
write.delim(subset(refseq2,strand=='-' ,c('chrom','polyAend','Ends','PolyA_downName','exonCount','strand')),File=polyAFile,col.names=F,append=T)


######################################################
HK = list(0,1)
S=list()
a=list()
for (h in HK){
  a[[length(a)+1]] <- with(subset(mergedData,isHK%in%h), tapply(Acceptor,list(GeneSize,intronCountGroups),median, na.rm=T))
  S[[length(S)+1]] <- with(subset(mergedData,isHK%in%h), tapply(Stability,list(GeneSize,intronCountGroups),median, na.rm=T))
}


cat("	Plotting data 2: \n") 




######################################################
pdf(sprintf('%s_isLast_boxplots3.pdf',settings$CommonName),height=6,width=2)
 par(mai=c(.85,.85,0.24,0.0),cex=0.5, cex.main=1, mfrow=c(2,1))
 boxplot(Stability ~ isLast,data=mergedData,outline=F,notch=T, range=0.00001,ylab='Stability score',xaxt='n',ylim=c(-1.25,1.25))
 boxplot(Acceptor ~ isLast,data=mergedData,outline=F,notch=T,range=0.00001,ylab='Acceptor score',xaxt='n',ylim=c(-1.25,1.25))
 dev.off()





#cols=brewer.pal(11,'Spectral')
cols=c(brewer.pal(9,"YlGnBu")[5:9],brewer.pal(9,"OrRd")[5:9])
HK = list(0:1,0,1)
  for (h in HK){
    txt = ifelse(length(h)==1,ifelse(h,'_HK','_Regulated'),'')
    
    plot.dev(sprintf('%s_Stability_by_numIntrons_geneSize%s3.pdf',settings$CommonName,txt),'pdf',height=3,width=4)
    par(mai=c(.45,.45,0.24,0.1),cex=0.6, cex.main=1) 
    
    for(i in 0:1){
      try(boxplot(Stability ~ intronCountGroups+GeneSize,data=subset(mergedData,isHK%in%h),notch=T,staplewex=0,ylim=2*c(-1.25,1.25),col=cols[1:5],at=c(1:5,7:11,13:17,19:23),outline=F,range=0.00001,add=i,ylab='Nucleosome Stability score',xaxt='n'))
      axis(1,at=c(3,9,15,21),lab=c('Short','Med short','Med long','Long'))
      if(!i) 
        abline(h=seq(-4,4,by=0.5),col='gray',lty=2)
    }
    plot.off()
    
    plot.dev(sprintf('%s_Stability_by_lengthIntrons_geneSize%s3.pdf',settings$CommonName,txt),'pdf',height=3,width=4)
    par(mai=c(.45,.45,0.24,0.1),cex=0.6, cex.main=1)
    for(i in 0:1){
      #boxplot(Stability ~ intronLengthGroups+GeneSize,data=subset(mergedData,isHK%in%h),notch=T,staplewex=0,ylim=c(-1.25,1.25),col=cols[1:4],at=c(1:4,6:9,11:14,16:19),outline=F,range=0.00001,add=i,ylab='Acceptor score',xaxt='n')
      #axis(1,at=c(2.5,7.5,12.5,17.5),lab=c('Short','Med short','Med long','Long'))
      boxplot(Stability ~ GeneSize+intronLengthGroups,data=subset(mergedData,isHK%in%h),notch=T,staplewex=0,ylim=c(-1.25,1.25),col=cols[6:9],at=c(1:4,6:9,11:14,16:19),outline=F,range=0.00001,add=i,ylab='Acceptor score',xaxt='n')
      axis(1,at=c(2.5,7.5,12.5,17.5),lab=c('Short introns','Med short introns','Med long introns','Long introns'))      
      if(!i) 
        abline(h=seq(-4,4,by=0.5),col='gray',lty=2)
    }
    plot.off()
    
    plot.dev(sprintf('%s_Acceptor_by_numIntrons_geneSize%s3.pdf',settings$CommonName,txt),'pdf',height=3,width=4)
    par(mai=c(.45,.45,0.24,0.1),cex=0.6, cex.main=1)
    for(i in 0:1){
      boxplot(Acceptor ~ intronCountGroups+GeneSize,data=subset(mergedData,isHK%in%h),notch=T,staplewex=0,ylim=c(-1.25,1.25),col=cols[1:5],at=c(1:5,7:11,13:17,19:23),outline=F,range=0.00001,add=i,ylab='Acceptor score',xaxt='n')
      axis(1,at=c(3,9,15,21),lab=c('Short','Med short','Med long','Long'))
      if(!i) 
        abline(h=seq(-4,4,by=0.5),col='gray',lty=2)
    }
    plot.off()
   
   
   
    #### Acceptor stuff, for supplement  
    plot.dev(sprintf('%s_Acceptor_by_lengthIntrons_geneSize%s3.pdf',settings$CommonName,txt),'pdf',height=3,width=4)
    par(mai=c(.45,.45,0.24,0.1),cex=0.6, cex.main=1)
    for(i in 0:1){
      boxplot(Acceptor ~ GeneSize+intronLengthGroups,data=subset(mergedData,isHK%in%h),notch=T,staplewex=0,ylim=c(-1.25,1.25),col=cols[6:9],at=c(1:4,6:9,11:14,16:19),outline=F,range=0.00001,add=i,ylab='Acceptor score',xaxt='n')
      axis(1,at=c(2.5,7.5,12.5,17.5),lab=c('Short introns','Med short introns','Med long introns','Long introns'))      
      if(!i) 
        abline(h=seq(-4,4,by=0.5),col='gray',lty=2)
    }
    plot.off()
    
    plot.dev(sprintf('%s_Acceptor_by_geneSize_lengthIntrons%s3.pdf',settings$CommonName,txt),'pdf',height=3,width=4)
    par(mai=c(.45,.45,0.24,0.1),cex=0.6, cex.main=1)
    for(i in 0:1){
      boxplot(Acceptor ~ intronLengthGroups+GeneSize,data=subset(mergedData,isHK%in%h),notch=T,staplewex=0,ylim=c(-1.25,1.25),col=cols[6:9],at=c(1:4,6:9,11:14,16:19),outline=F,range=0.00001,add=i,ylab='Acceptor score',xaxt='n')
      axis(1,at=c(2.5,7.5,12.5,17.5),lab=c('Short','Med short','Med long','Long'))
      if(!i) 
        abline(h=seq(-4,4,by=0.5),col='gray',lty=2)
    }
    plot.off()
    
    
    plot.dev(sprintf('%s_Acceptor_by_lenIntr_gSize_isLast%s3.pdf',settings$CommonName,txt),'pdf',height=3,width=5)
    par(mai=c(.45,.45,0.24,0.1),cex=0.6, cex.main=1)
    AT=c(1:8,10:17,19:26,28:35)
    for(i in 0:1){
      boxplot(Acceptor ~ isLast+GeneSize+intronLengthGroups,data=subset(mergedData,isHK%in%h),border=c(rgb(.75,.75,.75),'black'),notch=T,staplewex=0,ylim=c(-1.25,1.25),col=rep(cols[6:9],each=2),at=AT,outline=F,range=0.00001,add=i,ylab='Acceptor score',xaxt='n')
      if(i) axis(1,at=c(4.5,13.5,22.5,31.5),lab=c('Short introns','Med short introns','Med long introns','Long introns'))
      if(!i) 
        abline(h=seq(-4,4,by=0.5),col='gray',lty=2)
        
    }
    plot.off()
    
    plot.dev(sprintf('%s_Stability_by_lenIntr_gSize_isLast%s3.pdf',settings$CommonName,txt),'pdf',height=3,width=4)
    par(mai=c(.45,.45,0.24,0.1),cex=0.6, cex.main=1)
    AT=c(1:8,10:17,19:26,28:35)
    for(i in 0:1){
      boxplot(Stability ~ isLast+GeneSize+intronLengthGroups,data=subset(mergedData,isHK%in%h),notch=T,staplewex=0,ylim=c(-1.25,1.25),col=rep(cols[6:9],each=2),at=AT,outline=F,range=0.00001,add=i,ylab='Acceptor score',xaxt='n')
       if(i) axis(1,at=c(4.5,13.5,22.5,31.5),lab=c('Short introns','Med short introns','Med long introns','Long introns'))
       if(!i) 
        abline(h=seq(-4,4,by=0.5),col='gray',lty=2)
        
    }
    plot.off()
    
    
 
  }
#
### Also split this up between LAST and INTERNAL exons
#HK = list(0,1)
#  for (h in HK){
#    txt = ifelse(length(h)==1,ifelse(h,'_last','_internal'),'')
#  #png(sprintf('%s_Stability_by_numIntrons_geneSize.png',settings$CommonName),height=300,width=400)
#  pdf(sprintf('%s_Stability_by_numIntrons_geneSize%s3.pdf',settings$CommonName,txt),height=3,width=4)
#  par(mai=c(.45,.45,0.24,0.1),cex=0.6, cex.main=1) 
#  
#  for(i in 0:1){
#    boxplot(Stability ~ intronCountGroups+GeneSize,data=subset(mergedData,isLast%in%h),notch=T,staplewex=0,ylim=c(-1.25,1.25),col=cols[1:5],at=c(1:5,7:11,13:17,19:23),outline=F,range=0.00001,add=i,ylab='Nucleosome Stability score',xaxt='n')
#    axis(1,at=c(3,9,15,21),lab=c('Short','Med short','Med long','Long'))
#    if(!i) 
#      abline(h=seq(-4,4,by=0.5),col='gray',lty=2)
#  }
#  dev.off()
#  
#  #png(sprintf('%s_Acceptor_by_numIntrons_geneSize.png',settings$CommonName),height=300,width=400)
#  pdf(sprintf('%s_Acceptor_by_numIntrons_geneSize%s3.pdf',settings$CommonName,txt),height=3,width=4)
#  par(mai=c(.45,.45,0.24,0.1),cex=0.6, cex.main=1)
#  for(i in 0:1){
#    boxplot(Acceptor ~ intronCountGroups+GeneSize,data=subset(mergedData,isLast%in%h),notch=T,staplewex=0,ylim=c(-1.25,1.25),col=cols[1:5],at=c(1:5,7:11,13:17,19:23),outline=F,range=0.00001,add=i,ylab='Acceptor score',xaxt='n')
#    axis(1,at=c(3,9,15,21),lab=c('Short','Med short','Med long','Long'))
#    if(!i) 
#      abline(h=seq(-4,4,by=0.5),col='gray',lty=2)
#  }
#   dev.off()
#  }
#                
#                
#                
#                
##### 10/30/2013: Try to split up by last intron length somehow
#
#HK = as.list(1:4)
#  for (h in HK){
#    txt = switch(h,`1`='Short Introns',`2`='Med short Introns',`3`='med long introns',`4`='long introns')
#  #png(sprintf('%s_Acceptor_by_numIntrons_geneSize.png',settings$CommonName),height=300,width=400)
#  plot.dev(sprintf('%s_A_numInt_gSize_%s3.pdf',settings$CommonName,txt),'pdf',height=3,width=4)
#  par(mai=c(.45,.45,0.24,0.1),cex=0.6, cex.main=1)
#  for(i in 0:1){
#    boxplot(Acceptor ~ intronCountGroups+GeneSize,data=subset(mergedData,intronLengthGroups%in%h),notch=T,staplewex=0,ylim=c(-1.25,1.25),col=cols[1:5],at=c(1:5,7:11,13:17,19:23),outline=F,range=0.00001,add=i,ylab='Acceptor score',xaxt='n')
#    axis(1,at=c(3,9,15,21),lab=c('Short','Med short','Med long','Long'))
#    if(!i) 
#      abline(h=seq(-4,4,by=0.5),col='gray',lty=2)
#  }
#   plot.off()
#  }
#
##  refseq$intronCountGroups = splitByVector(refseq$exonCount-1,c(2.5,4.5,7.5,19.5))
#HK = as.list(1:5)
#  for (h in HK){
#    txt = switch(h,`1`='1-2 introns',`2`='3-4 introns',`3`='5-7 introns',`4`='8-19 introns',`5`='20+introns')
#  #png(sprintf('%s_Acceptor_by_numIntrons_geneSize.png',settings$CommonName),height=300,width=400)
#  plot.dev(sprintf('%s_A_intSize_gSize%s3.pdf',settings$CommonName,txt),'pdf',height=3,width=4)
#  par(mai=c(.45,.45,0.24,0.1),cex=0.6, cex.main=1)
#  for(i in 0:1){
#    boxplot(Acceptor ~ intronLengthGroups+GeneSize,data=subset(mergedData,intronCountGroups%in%h),notch=T,staplewex=0,ylim=c(-1.25,1.25),col=cols[1:5],at=c(1:4,6:9,11:14,16:19),outline=F,range=0.00001,add=i,ylab='Acceptor score',xaxt='n')
#    devaxis(1,at=c(3,9,15,21),lab=c('Short','Med short','Med long','Long'))
#    if(!i) 
#      abline(h=seq(-4,4,by=0.5),col='gray',lty=2)
#  }
#   plot.off()
#  }
#
#HK = as.list(1:4)
#AT = c(1:5,7:11,13:17,19:23)
#  for (h in HK){
#    txt = switch(h,`1`='Short Genes',`2`='Med short Genes',`3`='med long Genes',`4`='long Genes')
#    at = ifelse(h==1,list(AT[-5*(1:4)]),list(AT))[[1]]
#  plot.dev(sprintf('%s_A_numInt_intSize_%s3.pdf',settings$CommonName,txt),'pdf',height=3,width=4)
#  par(mai=c(.45,.45,0.24,0.1),cex=0.6, cex.main=1)
#  for(i in 0:1){
#    boxplot(Acceptor ~ intronCountGroups+intronLengthGroups,data=subset(mergedData,GeneSize%in%h),notch=T,staplewex=0,ylim=c(-1.25,1.25),col=cols[1:5],at=at,outline=F,range=0.00001,add=i,ylab='Acceptor score',xaxt='n')
#    axis(1,at=c(3,9,15,21),lab=c('Short Introns','Med short Introns','Med long Introns','Long Introns'))
#    if(!i) 
#      abline(h=seq(-4,4,by=0.5),col='gray',lty=2)
#  }
#   plot.off()
#  }
#
#

cat("	Plotting data 3:\n") 






pdf(sprintf('%s_HK_dataMining.pdf',settings$CommonName),height=6, width=6)

  par(cex=0.5, mai=c(0.6,1,.2,.2), mfrow=c(2,2))
  boxplot(Stability ~ isHK, mergedData, notch=T,names=c("Regulated","HK"), border=c(1,2),main='',ylab=sprintf("Stability Score",settings$CommonName),pch=19)
  boxplot((Acceptor) ~ isHK, mergedData, notch=T,names=c("Regulated","HK"), border=c(1,2),main='',ylab=sprintf("Acceptor Score",settings$CommonName),pch=19)
  
  boxplot(Stability ~ isHK, mergedData, notch=T,staplewex=0,range=0.00001,outline=F,names=c("Regulated","HK"), border=c(1,2),main='',ylab=sprintf("Stability Score",settings$CommonName),pch=19)
  boxplot((Acceptor) ~ isHK, mergedData, notch=T,staplewex=0,range=0.00001,outline=F, names=c("Regulated","HK"), border=c(1,2),main='',ylab=sprintf("Acceptor Score",settings$CommonName),pch=19)

dev.off()


##### Split up and plot a bunch of stuff!!
getQuantilMedians = function(myData,col1='Dist2End',column='splicing0',N=100){
  myData$Group = quantileGroups(myData[,col1],N)
  data.frame(LengthMedian = tapply(myData[,col1],myData$Group,median,na.rm=T), median = tapply(myData[,column],myData$Group,median,na.rm=T))
}


num1 = 20
num2 = 10

for(column in c('Stability','Acceptor')){
  colText = ifelse(column=='Acceptor','Acceptor Strength','Nucleosome Stability')
  
  #png(sprintf('%s_%s_new1_%%d.png',settings$CommonName,column),height=450,width=450)
  #png(sprintf('%s_%s_new1.png',settings$CommonName,column),height=350,width=350*5)
  pdf(sprintf('%s_%s_new3.pdf',settings$CommonName,column),height=4,width=4*5)  
  par(mfrow=c(1,4),cex=1.1)


  ## First plot feature versus length
  pchLast = c(20,15)
  Qdata = getQuantilMedians(subset(mergedData,isLast==0),N=num1, column=column)
  Qdata_last = getQuantilMedians(subset(mergedData,isLast==1),N=num1, column=column)
  
  plot(log10(Qdata[,1]), Qdata[,2],xaxt='n',pch=pchLast[1],xlim=log10(c(80,420000)),ylim=c(-2.5,2.5),main=colText, ylab=paste('Median',colText), xlab="Distance to End of Gene") 
  points(log10(Qdata_last[,1]), Qdata_last[,2],pch=pchLast[2])
  axis(1,at=1:6,lab=format(10^(c(NA,2:5,NA)),scientific=F))
 if(column=='Stability')
    legend('topright',pch=pchLast,leg=c('non-last','last'))


  ## Split up by gene size:
  myCols = brewer.pal(4,'Dark2')
  for(i in 1:4){
      Qdata = getQuantilMedians(subset(mergedData,isLast==0 & GeneSize==i),N=num2, column=column)
      Qdata_last = getQuantilMedians(subset(mergedData,isLast==1 & GeneSize==i),N=num2, column=column)

      if(i ==1)
        plot(log10(Qdata[,1]), Qdata[,2],xaxt='n',pch=pchLast[1],xlim=log10(c(80,420000)),col=myCols[i],ylim=c(-1.24,1.24),main=colText, ylab=paste('Median',colText), xlab="Distance to End of Gene") 
      else        
        points(log10(Qdata[,1]), Qdata[,2],col=myCols[i],pch=pchLast[1])
        
        points(log10(Qdata_last[,1]), Qdata_last[,2],col=myCols[i],pch=pchLast[2])   
  
 }

 axis(1,at=1:6,lab=format(10^(c(NA,2:5,NA)),scientific=F))

 if(column=='Stability')
    legend('topright',pch=19,col=myCols,leg=c(sprintf('len < %d',qr[1]), sprintf('len < %d', ceiling(qr[2])), sprintf('len < %d', ceiling(qr[3])), sprintf('len >= %d ',ceiling(qr[3])) ) )

  
  ## Split up long and short genes based on # of exons ...
  cols=brewer.pal(10,'Spectral')
  cols2=cols[c(1,3,7,9,10)]
  cols3 = myCols[c(1,4)]
  

  HK = c(1,0)
  # for(HK in list(1,0)){
  myText = ifelse(length(HK)==2,'All',ifelse(HK==1,'HK','non'))
  lists=  list(subset(mergedData,GeneSize==1 & isHK %in% HK), subset(mergedData,GeneSize==4 &  isHK %in% HK))
  pchs = c(16,17)
  pchs2 = c(1,17,8,6,19)
  lists[[1]]$intronCountGroups = quantileGroups(lists[[1]]$MaxFeature,4)
  groups = lapply(tapply(lists[[1]]$MaxFeature,lists[[1]]$intronCountGroups,range),function(x)x[1]:x[2])
  
  # Do for both short and long
  for(L in 1:2){
    
    if(L==2)
      groups[[length(groups)+1]] = (tail(tail(groups,1)[[1]],1)+1):max(lists[[2]]$MaxFeature)
    for(i in 1:length(groups)){
      QED = getQuantilMedians(subset(lists[[L]], MaxFeature %in% groups[[i]]),N=num2,column=column)
      #COLOR = cols2[i]
      #PCH =pchs[L]
      COLOR = cols3[L]
      PCH = pchs2[i]
      
       if(i==1 & L==1)
        plot(log10(QED[,1]),QED[,2],xaxt='n', pch=PCH,col=COLOR,xlim=log10(c(500,420000)),ylim=c(-1.25,1.25),main=colText, ylab=paste('Median',colText), xlab="Distance to End of Gene")
       else
        points(log10(QED[,1]),QED[,2],pch=PCH,col=COLOR)                  
    } 
  }
  axis(1,at=1:6,lab=format(10^(c(NA,2:5,NA)),scientific=F))
  #legend('topright',pch=19,col=cols2, leg=sub(":","-",paste(groups,'introns')))
  #legend('bottomleft',pch=pchs,leg=c(sprintf('len < %d',qr[1]), sprintf('GeneLength >=',ceiling(qr[3])) ))
  if(column=='Stability'){
    legend('topright',pch=pchs2, leg=sub(":","-",paste(groups,'introns')))
    legend('bottomleft',pch=1,col=cols3,leg=c(sprintf('len < %d',qr[1]), sprintf('len >= %d',ceiling(qr[3])) ))
  }
  
  ## Split up based on HK, instead of first/last-ness ...
  hkCols = c('red','black')
  for(i in 1:2){
      Qdata = getQuantilMedians(subset(mergedData,isHK== 2-i & isLast==0),N=num2, column=column)
      Qdata_last = getQuantilMedians(subset(mergedData,isHK==2-i & isLast==1),N=num2, column=column)

      if(i ==1)
        plot(log10(Qdata[,1]), Qdata[,2],xlim=log10(c(80,420000)),xaxt='n',col=hkCols[i],pch=pchLast[1],ylim=c(-1.24,1.24),main=colText, ylab=paste('Median',colText), xlab="Distance to End of Gene") 
      else        
        points(log10(Qdata[,1]), Qdata[,2],col=hkCols[i],pch=pchLast[1])
        
      points(log10(Qdata_last[,1]), Qdata_last[,2],col=hkCols[i],pch=pchLast[2])   
  }
  axis(1,at=1:6,lab=format(10^(c(NA,2:5,NA)),scientific=F))
  if(column=='Stability')
    legend('topright',pch=c(pchLast,pchLast),col=hkCols[c(1,1,2,2)],leg=c("HK non-last","HK last","Regulated non-last","Regulated last"))


 ## plot Acceptor strength vs.  length of intron
 pchIntronLen=c(18,10)
  for(i in 1:4){
      Qdata = getQuantilMedians(subset(mergedData,intronLengthGroups ==1 & GeneSize==i),N=num2, column=column)
      Qdata_last = getQuantilMedians(subset(mergedData,intronLengthGroups ==4 & GeneSize==i),N=num2, column=column)

      if(i ==1)
        plot(log10(Qdata[,1]), Qdata[,2],xaxt='n',pch=pchIntronLen[1],xlim=log10(c(80,420000)),col=myCols[i],ylim=c(-1.24,1.24),main=colText, ylab=paste('Median',colText), xlab="Distance to End of Gene") 
      else        
        points(log10(Qdata[,1]), Qdata[,2],col=myCols[i],pch=pchIntronLen[1])
        
      if(i> 1) points(log10(Qdata_last[,1]), Qdata_last[,2],col=myCols[i],pch=pchIntronLen[2])   
  
 }

 axis(1,at=1:6,lab=format(10^(c(NA,2:5,NA)),scientific=F))
if(column!='Stability'){
    legend('bottomright',pch=19,col=myCols,leg=c(sprintf('len < %d',qr[1]), sprintf('len < %d', ceiling(qr[2])), sprintf('len < %d', ceiling(qr[3])), sprintf('len >= %d ',ceiling(qr[3])) ) )
    legend('topright',pch=pchIntronLen,leg=c('short introns','long introns'))

}


  dev.off()           
}








 







































### OLD

#
## This is set up for 0-indexed data
#
## Determine positions of the upstream:
#Up.left = refseq$Ends + ifelse(input$strand == '+', -(Up + 1), 0)
#Up.right = refseq$Ends + ifelse(input$strand == '+', -1, Up)
#
## Determine positions of the downstream:
#Down.left = refseq$Ends + ifelse(input$strand == '+', 0, -Down-1)
#Down.right = refseq$Ends + ifelse(input$strand == '+', Down,-1)
#
#nonNeg= (Down.left) > 0 & Down.right > 0 & Up.left > 0 & Up.right > 0
#
#Output = rbind(
#cbind(input$chr,Up.left,Up.right,paste(input$name,'Up',sep='_'),0,input$strand)[nonNeg,],
#cbind(input$chr,Down.left,Down.right,paste(input$name,'Down',sep='_'),0,input$strand)[nonNeg,])
#
#write.delim(Output,outputFile,col.names=F)
#write.delim(Output[!duplicated(Output[,4]) ,],paste(outputFile,'noDups',sep='.'),col.names=F)
#
