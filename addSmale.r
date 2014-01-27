## use data from Tilgner et al. 2012 Genome Research

rm(list=ls())
library(gplots)
library(RColorBrewer)
library(Hmisc)
source("/home/jeremy/ExonPipeline/ExonPipeline_functions.r")
source("/home/jeremy/ExonPipeline/ExonPipeline_analysis_and_plotting_Functions.r")
source("/home/home/jeremy/Code/useful_R/useful_R.r")

setwd("/home/jeremy/ExonPipeline/mm10")
#load("PlottedData6.RData")
load("Smale_mergedData.RData")
load("Smale_mergedData_withSplicing.RData")
ChromInfo = cbind(num=352298:352302, times=c(0,15,30,60,120))


# Load the features (exon, intron info)
features = loadFeatures(genome=loadGenome(settings<-setup( getwd())),settings=settings,refseq = add_UniqueID(loadRefgene(settings)))
keepNames = names(which(tapply(features$introns$FeatureCount,features$introns$UniqueID,min)==1))
introns2= mostDownstream(subset(features$introns, UniqueID %in% keepNames))
exons2= mostDownstream(subset(features$exons,UniqueID %in% keepNames))   

introns2$length = abs(introns2$start-introns2$end)
## Figure out the end of the gene (poly A).
introns2$terminus = introns2$end * ifelse(introns2$strand=='-',-1,1)
introns2$begin = introns2$start * ifelse(introns2$strand=='-',-1,1)
txEnd=abs(tapply(introns2$terminus, introns2$UniqueID,max)); 
txStart=abs(tapply(introns2$begin, introns2$UniqueID,min)); 
numExons = table(introns2$UniqueID)
TXINFO = data.frame(UniqueID = names(txStart), txStart=txStart, txEnd=txEnd, MaxFeature=as.vector(numExons))


## Merge with exons
introns2$energy_ID = paste(">",paste(introns2$UniqueID,introns2$FeatureCount+1,introns2$isLast,sep='_'),sep='')
introns2$acceptor_ID = paste(introns2$UniqueID,introns2$FeatureCount,sep='_')
introns2$donor_ID = paste(introns2$UniqueID,introns2$FeatureCount,sep='_')

energyData2 = energyData; names(energyData2) = c("energy_ID","Stability")
donorData2 = donorData; donorData2$ID = sub("_([0|1])$","",donorData2$V1) ; names(donorData2) = c("","Donor","donor_ID") 
acceptorData2 = acceptorData; acceptorData2$ID = sub("_([0|1])$","",acceptorData2$V1) ; names(acceptorData2) = c("","Acceptor","acceptor_ID")

## Merge all data together (Merge INTRONS this time!)
mergedData = merge(merge(merge(merge(introns2,energyData2,by='energy_ID'),donorData2[,-1],by='donor_ID',all=T),acceptorData2[,-1],by='acceptor_ID',all=T),TXINFO,by='UniqueID')

mergedData$Dist2End = abs(mergedData$txEnd - mergedData$end)
mergedData$Dist2Start = abs(mergedData$txStart - mergedData$start) + 1
mergedData$Group = quantileGroups(mergedData$Dist2End,100)

## Add HK status
#HKlist = read.delim("/home/jeremy/ExonPipeline/hg19/HouseKeepingGenes.txt",skip=1)
mergedData$isHK = toupper(mergedData$GeneName) %in% HKlist$Gene.name
save(mergedData,file="Smale_mergedData.RData")

########################################################
## GET SSJ data  (Dmitri D. Pervouchine)
########################################################

# My data is 0-indexed, their is not??
# also, i added a '1' b/c the way i ran bam2ssj, the 1 is added like *_1

mergedData$coSI_ID = paste(mergedData$chr,rowMins(mergedData[,c('start','end')])+1,rowMaxes(mergedData[,c('start','end')]),1,mergedData$strand,sep='_')


# Read in bam2ssj file of Chromatin data for Smale  http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32916
ChromInfo = cbind(num=352298:352302, times=c(0,15,30,60,120))

splicingData = list()
for(i in 1:nrow(ChromInfo)){
  dat=read.delim(paste("/home/RNAseq/Smale/chrom_",ChromInfo[i,1],"/sorted.ssj",sep=""), head=F,stringsAsFactors=F)
  names(dat) = c('coSI_ID','count53', 'count5X', 'countX3', 'count50', 'count03')
  print(dim(dat))
  
  # Compute Theta5, Theta3  (http://bioinformatics.oxfordjournals.org/content/early/2012/11/21/bioinformatics.bts678.full.pdf+html)
  dat$Denom5 = dat$count53 + dat$count5X + dat$count50
  dat$Theta5 = (dat$count53 + dat$count5X) / dat$Denom5
  dat$Denom3 = dat$count53 + dat$countX3 + dat$count03
  dat$Theta3 = (dat$count53 + dat$countX3) / dat$Denom3

  # My way:
  dat$DenomJer = dat$count53 + (dat$count50 + dat$count03)/2
  dat$SpliceJer = dat$count53 / dat$DenomJer
  
  splicingData[[i]] = dat[dat$DenomJer  > 15,]
}


########################################################
## Merge with SSJ data 
########################################################



#hintrons3 = merge(merge(merge(merge(merge(merge(mergedData,splicingData,by = 'coSI_ID'),by = 'coSI_ID'),by = 'coSI_ID'),by = 'coSI_ID'),by = 'coSI_ID')


##
mergedData$coSI_ID = paste(mergedData$chr,rowMins(mergedData[,c('start','end')])+1,rowMaxes(mergedData[,c('start','end')]),1,sep='_')
mergedData$splicing0 = splicingData[[1]]$SpliceJer[match(mergedData$coSI_ID, splicingData[[1]]$coSI_ID)]
mergedData$splicing15 = splicingData[[2]]$SpliceJer[match(mergedData$coSI_ID, splicingData[[2]]$coSI_ID)]
mergedData$splicing30 = splicingData[[3]]$SpliceJer[match(mergedData$coSI_ID, splicingData[[3]]$coSI_ID)]
mergedData$splicing60 = splicingData[[4]]$SpliceJer[match(mergedData$coSI_ID, splicingData[[4]]$coSI_ID)]
mergedData$splicing120= splicingData[[5]]$SpliceJer[match(mergedData$coSI_ID, splicingData[[5]]$coSI_ID)]

introns3 = mergedData[apply(mergedData[,paste('splicing',ChromInfo[,2],sep='')],1,function(x)any(!is.na(x))),]

table(apply(introns3[,paste('splicing',ChromInfo[,2],sep='')],1,function(x)length(which(!is.na(x)))))

#   1    2    3    4    5 
#3531 2107 1702 1683 3598

### Table on all Splicing data before merging:
# No filtering:
#    1     2     3     4     5 
# 5292  4335  4922  6537 20927 
# 20K introns have 5!

# Filtered at 15 reads:
#   1    2    3    4    5 
#3668 2185 1748 1751 3750 
introns3$intronsRemaining = introns3$MaxFeature - introns3$FeatureCount
save(introns3,file="Smale_mergedData_withSplicing.RData")



#####################################################################
## Model MEDIANS
#####################################################################

#exons3$Group = splitByVector(exons3$Dist2End,quantile(exons3$Dist2End,na.rm=T,(1:99)/100))
getQuantilMedians = function(myData,column='splicing0',N=100){
  myData$Group = splitByVector(myData$Dist2End,quantile(myData$Dist2End,na.rm=T,(1:(N-1))/N))
  data.frame(LengthMedian = tapply(myData$Dist2End,myData$Group,median,na.rm=T), median = tapply(myData[,column],myData$Group,median,na.rm=T))
}

QuantileDataMinus  = getQuantilMedians(introns3,column='splicing0',100) 
QuantileDataMinus_HK  = getQuantilMedians(subset(introns3,isHK==1),column='splicing0',100) 
QuantileDataMinus_non  = getQuantilMedians(subset(introns3, isHK==0),column='splicing0',100) 



MeanGamma2 = function(LengthMedian,g,k,A,B)  sapply(LengthMedian, function(x)A * integrate(pgamma,0,x+B, shape=g, rate=k)$value/(x+B))
gMinus = nls(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileDataMinus,start=c(g=1,k=5e-5,A=1,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,Inf,Inf),algo='port');
gMinus_HK = nls(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileDataMinus_HK,start=c(g=1,k=5e-5,A=.85,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,Inf,Inf),algo='port');
gMinus_non = nls(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileDataMinus_non,start=c(g=1,k=5e-5,A=.85,B=1000),lower=c(1,0,0,1),upper=c(1,Inf,Inf,Inf),algo='port');

gammaMinus_k_splice
    All.k      HK.k     non.k 
0.3356356 0.3444153 0.5649482 
> gammaMinus_Dist_ReadThru
    All.B      HK.B     non.B 
 77232.78 130739.15  38222.58 


## Plot crazy exons graph ...
cols=brewer.pal(10,'Spectral')

for(column in c("splicing0","splicing15","splicing30", "splicing60", "splicing120")){
  QuantileData  = getQuantilMedians(introns3,column=column,20) 

  pdf(sprintf("Median_%s_by_remaining_introns_withLines.pdf",column))
  plot(log10(QuantileData$LengthMedian), xlim=log10(c(500,420000)), ylim=c(0.75,.96), QuantileData$median, ylab='Median coSI', xlab='LOG10 Distance to End of Gene' , main=paste('Chromatin',column),pch=20)
  
  allFits = list()
  MedLengths = list()
  MedCoSI = list()
  myGroups = list(1:2,NA,3:4,NA,5:6,NA,7:9,NA,10)
  S=sapply(AA<-seq(1,10,2), function(i){
       if(i<9)
        QD = getQuantilMedians(subset(introns3,intronsRemaining %in% myGroups[[i]]),column,10)
       else
        QD = getQuantilMedians(subset(introns3,intronsRemaining >=myGroups[[i]]),column,10)      
       MedLengths[[i]] <<- QD[,1]
       MedCoSI[[i]] <<- QD[,2]      
       FIT <- lm(median ~ log10(LengthMedian),QD)
       allFits[[i]]<<- FIT
       abline(a=coef(FIT)[1], b=coef(FIT)[2],col=cols[i],lty=2)
       points(log10(QD[,1]),QD[,2],col=cols[i],pch=20,cex=1.2)                  
  })
  print(range(sapply(MedCoSI[AA],range)))
  legend('bottomright',col=cols[AA],pch=20,leg=c("1-2","3-4","5-6","7-9",">=10"))
  dev.off()

}

## Now make one plot for each exon-category...
myGroups = list(1:2,NA,3:4,NA,5:6,NA,7:9,NA,10)
cols2 = rev(brewer.pal(5,'PRGn'))
S=sapply(AA<-seq(1,10,2), function(i){
  pdf(sprintf("Median_%s_by_timepoint.pdf",sub(":","-",myGroups[[i]])))
  
  columns = c("splicing0","splicing15","splicing30", "splicing60", "splicing120")
  for(Col in 1:5){
    column = columns[Col]  

    if(i<10)
      QD = getQuantilMedians(subset(introns3,intronsRemaining %in% myGroups[[i]]),column,10)
    else
       QD = getQuantilMedians(subset(introns3,intronsRemaining >=myGroups[[i]]),column,10)      

    if (Col==1)
      plot(log10(QD[,1]),QD[,2],col=cols2[Col],pch=20,cex=1.2, xlim=log10(c(500,420000)), ylim=c(0.75,.96),ylab='Median coSI', xlab='LOG10 Distance to End of Gene' )                          
    else
      points(log10(QD[,1]),QD[,2],col=cols2[Col],pch=20,cex=1.2)                  
      
    FIT <- lm(median ~ log10(LengthMedian),QD)
    abline(a=coef(FIT)[1], b=coef(FIT)[2],col=cols2[Col],lty=2)   
    
  }     
  legend('bottomright',col=cols2,pch=20,leg=columns)
  dev.off()
})





#########################################################################################################
## GRO-seq analysis
#########################################################################################################
refseq <-add_UniqueID(loadRefgene(settings<-setup( getwd())))
GRO = makeGROdata("/home/home/jeremy/RNAseq/Glass/refseq_and_subsequent_transcripts_with_tag_counts.txt", refseq,"GRO_Data_test.RData",1000)

introns3$GRO_data = GRO$runon[match(paste(introns3$chr,introns3$txStart+1,introns3$txEnd), paste(GORO$chr,GRO$txStart,GRO$txEnd))













