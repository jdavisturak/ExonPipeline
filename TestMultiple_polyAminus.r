rm(list=ls())
library(gplots)
library(RColorBrewer)
library(Hmisc)

source("/home/jeremy/ExonPipeline/ExonPipeline_functions.r")
source("/home/jeremy/ExonPipeline/ExonPipeline_analysis_and_plotting_Functions.r")
source("/home/home/jeremy/Code/useful_R/useful_R.r")
source("/home/jeremy/ExonPipeline/ExonPipeline_simulate_model.r")

setwd("/home/jeremy/ExonPipeline/hg19")
load("PlottedData_new3.RData")

dataDir = '/home/RNAseq/Tilgner/Nucleus/'
ChromInfo  = cbind(c( "Chromatin.ssj","PolyAMinus.ssj","PolyAPlus.ssj","Gm12878.ssj","H1hesc.ssj","Helas3.ssj","Hepg2.ssj","Huvec.ssj","Nhek.ssj"), # ssj file
c('K562_Chromatin','K562_Nucleus A(-)','K562_Nucleus A(+)',"Gm12878","H1hesc","Helas3","Hepg2","Huvec","Nhek"),
c('K562','K562','K562',"Gm12878","H1hesc","Helas3","Hepg2","Huvec","Nhek"))

 
##########################################################
#### GET SSJ data  (Dmitri D. Pervouchine)
##########################################################

splicingData = list()
splicingData_Theta5 = list()
splicingData_Theta3 = list()
SplicingMerged = matrix(NA,nr=nrow(mergedData), nc = nrow(ChromInfo))
mergedData$coSI_ID = paste(mergedData$chr,rowMins(mergedData[,c('start','end')])+1,rowMaxes(mergedData[,c('start','end')]),1,sep='_')

for(i in 1:nrow(ChromInfo)){
  dat=read.delim(paste(dataDir,ChromInfo[i,1],sep=""), head=F,stringsAsFactors=F)
  names(dat) = c('coSI_ID','count53', 'count5X', 'countX3', 'count50', 'count03')
  print(dim(dat))
  
  # Compute Theta5, Theta3  (http://bioinformatics.oxfordjournals.org/content/early/2012/11/21/bioinformatics.bts678.full.pdf+html)
  dat$Denom5 = dat$count53 + dat$count5X + dat$count50
  dat$Theta5 = (dat$count53 + dat$count5X) / dat$Denom5
  dat$Denom3 = dat$count53 + dat$countX3 + dat$count03
  dat$Theta3 = (dat$count53 + dat$countX3) / dat$Denom3

  # My way:
  #dat$DenomJer = dat$count53 + (dat$count50 + dat$count03)/2
  dat$DenomJer = dat$count53 + (dat$count50 + dat$count03)/2 + dat$count5X+dat$countX3
  dat$SpliceJer = dat$count53 / dat$DenomJer
  
  cutoff=15
  #splicingData[[i]] = dat[dat$DenomJer  > 15,]
  splicingData[[i]] = dat[dat$DenomJer  > cutoff & dat$count5X < 1 & dat$countX3 < 1,]
  splicingData_Theta5[[i]] = dat[dat$Denom5  > cutoff,]
  splicingData_Theta3[[i]] = dat[dat$Denom3  > cutoff,]  
  
  SplicingMerged[,i] = splicingData[[i]]$SpliceJer[match(mergedData$coSI_ID, splicingData[[i]]$coSI_ID)]

}                                                                                                                      

mergedData2 = cbind(mergedData, SplicingMerged)
names(mergedData2)[ncol(mergedData) + (1:ncol(SplicingMerged))] = ChromInfo[,2]



#############################################################################
### New 10/29/13: Analyze Cytoplasmic poly(A) RNAseq of these guys first, THEN  do the rest
#############################################################################

dataDir = '/home/RNAseq/Tilgner/Cytosol/'
getNames = function(file)read.delim(file, stringsAsFactors=F,head=F)[,1]
getFileData = function(file)read.delim(file, stringsAsFactors=F,head=F)[,4]

f=function(cell,rep)cbind(sprintf("wgEncodeCshlLongRnaSeq%sCytosolPapPlusRawSigRep%d.bigWig",cell,rep),sprintf("wgEncodeCshlLongRnaSeq%sCytosolPapMinusRawSigRep%d.bigWig",cell,rep))

myLists = list(K562=f("K562",c(1,2)),Gm12878=f("Gm12878",c(1,2)), H1hesc=f("H1hesc",2),Helas3=f("Helas3",c(1,2)),
   Hepg2=f("Hepg2",c(1,2)), Huvec=f("Huvec",c(3,4)), Nhek=f("Nhek",c(3,4)) )
IDs=sub("(.*)_.*","\\1",getNames(paste(dataDir,myLists[[1]][1,1],".refseq.Up",sep="")))
strand = refseq$strand[match(IDs,refseq$UniqueID)]  

myCleavedFractions = list()   

## For each of the lists, read in 'Up' and 'Down' reads to determine overlaps:     
for (i in 1:length(myLists)){
  Down.Plus = sapply(paste(dataDir,myLists[[i]][,1],".refseq.Down",sep=""),getFileData)
  Down.Minus = sapply(paste(dataDir,myLists[[i]][,2],".refseq.Down",sep=""),getFileData)
  Up.Plus = sapply(paste(dataDir,myLists[[i]][,1],".refseq.Up",sep=""),getFileData)
  Up.Minus = sapply(paste(dataDir,myLists[[i]][,2],".refseq.Up",sep=""),getFileData)
  
  rownames(Down.Plus) = rownames(Down.Minus) = rownames(Up.Plus) = rownames(Up.Minus) = IDs
  Up = Up.Plus
  Down = Down.Plus
  Up[strand=='-',] = Up.Minus[strand=='-',]
  Down[strand=='-',] = Down.Minus[strand=='-',]
  
  ## Now I need to filter by some cutoff and then take the average
  AllCleavedFraction = Up/(Up + Down)#1-Down/Up
  AllCleavedFraction2 = AllCleavedFraction
  AllCleavedFraction2[Up + Down <= 500] <- NA # This is equivalent to ~5 reads if length is 100 bp
  myCleavedFractions[[i]] = rowMeans(AllCleavedFraction2)      
}

names(myCleavedFractions) = names(myLists)


## Now for all of the cell types, omit the splicing data from genes that cleaved less than 95%
mergedData3 = mergedData2
# Now loop through ChromInfo
for ( i in 1:nrow(ChromInfo)){
  myColumn = ChromInfo[i,2]
  myCellType = ChromInfo[i,3]
  
  mergedData3[,myColumn][mergedData3$UniqueID %in% names(myCleavedFractions[[myCellType]][myCleavedFractions[[myCellType]] < 0.95])] <-NA  
}

# This summarizes how many introns passed the cutoffs
sapply(ChromInfo[,2],function(x)c(length(which(!is.na(mergedData2[,x]))),length(which(!is.na(mergedData3[,x])))))
#     K562_Chromatin K562_Nucleus A(-) K562_Nucleus A(+) Gm12878 H1hesc Helas3
#[1,]          27332             23790             27918   23080  14521  22624
#[2,]          21513             18505             22116   16179   9099  15870
#     Hepg2 Huvec  Nhek
#[1,] 25661 24013 20791
#[2,] 19117 16602 15189


#############################################################################
### New 10/31/13: Analyze NUCLEAR poly(A) RNAseq, 
# and just report the % of genes saved above that have high cleavage levels in this sample
############################################################################

dataDir = '/home/RNAseq/Tilgner/Nucleus/'
myNuclear_CleavedFractions = list()   

for (i in 1:nrow(ChromInfo)){
  Name = sub(".ssj","",ChromInfo[i,1])
                  
  Down.Plus = sapply(paste(dataDir,Name,".refseq.A.plus.Down",sep=""),getFileData)
  Down.Minus = sapply(paste(dataDir,Name,".refseq.A.minus.Down",sep=""),getFileData)
  Up.Plus = sapply(paste(dataDir,Name,".refseq.A.plus.Up",sep=""),getFileData)
  Up.Minus = sapply(paste(dataDir,Name,".refseq.A.minus.Up",sep=""),getFileData)
   
  rownames(Down.Plus) = rownames(Down.Minus) = rownames(Up.Plus) = rownames(Up.Minus) = IDs
  Up = Up.Plus
  Down = Down.Plus
  Up[strand=='-',] = Up.Minus[strand=='-',]
  Down[strand=='-',] = Down.Minus[strand=='-',]
  
  ## Now I need to filter by some cutoff and then take the average
  AllCleavedFraction = Up/(Up + Down)#1-Down/Up
  AllCleavedFraction2 = AllCleavedFraction
  AllCleavedFraction2[Up + Down <= 500] <- NA # This is equivalent to ~5 reads if length is 100 bp
  myNuclear_CleavedFractions[[i]] = rowMeans(AllCleavedFraction2)      
}

names(myNuclear_CleavedFractions) = names(myLists)

save(myCleavedFractions,myNuclear_CleavedFractions,mergedData2,mergedData3, file='Spring2014_Human_multiple_cellTypes_splicingInfo.RData')
print("Finished Saving the first part")


## For all genes in mergedData3 (those that have evidence of correct splicing in Cytoplasmic)
## measure the distribution of the cleavage Ratio of those genes in THIS sample
myPolyA = list()
for ( i in 1:nrow(ChromInfo)){
   myColumn = ChromInfo[i,2] 
   Name = sub(".ssj","",ChromInfo[i,1])    
   theseGenes = mergedData3$UniqueID[!is.na(mergedData3[, myColumn])]
   myPolyA[[i]] = myNuclear_CleavedFractions[[i]][names(myNuclear_CleavedFractions[[i]]) %in% theseGenes]      
}

# Obtain a list of genes that have cleavage Ratio of 0.8 or higher
polyA80 = sapply(1:9,function(i)try(length(which(myPolyA[[i]] > 0.8))/length(myPolyA[[i]])))





#############################################################################
### Now make median quantile plots for all of these
#############################################################################

getQuantilMedians = function(myData,col1='Dist2End',column='splicing0',N=100,func=median){
  myData$Group = quantileGroups(myData[,col1],N)
  data.frame(LengthMedian = tapply(myData[,col1],myData$Group,median,na.rm=T), median = tapply(myData[,column],myData$Group,func,na.rm=T))
}
MeanGamma2 = function(LengthMedian,g,k,A,B)  sapply(LengthMedian, function(x)A * integrate(pgamma,0,x+B, shape=g, rate=k)$value/(x+B))
MeanExp = function(LengthMedian,k,A)  sapply(LengthMedian, function(x)A * integrate(pgamma,0,x, shape=1, rate=k)$value/x)

# gammaDecay = function(x,shape,rate,decay) pgamma(x,shape=shape, rate=rate)*dexp(x,decay) # needs to be normalized by area under the decay curve
# MeanGamma3 = function(LengthMedian,g,k,A,B)  sapply(LengthMedian, function(x)integrate(gammaDecay,0,x+B, shape=g, rate=k,decay=A)$value/pexp(x+B,A))
# 


#########################################################################
## Fitting                                                               
#########################################################################
sem=function(x,na.rm=T){
  LEN = ifelse(na.rm,length(x[!is.na(x)]), length(x))  
  sd(x,na.rm=na.rm)/sqrt(LEN)
}

for(N in c(100,1000)) {
  A_min = 0;A_max=1
  #A_min=0.9; A_max=0.9
  FitData = list()
  FitData2 = list()
  FitData_noRT2 = list()
  FitData_RT2 = list()
  QuantileData = list()
  SDData = list()
  
  for (i in 1:nrow(ChromInfo)){                                                        
  
    print(i)
    starts = c(g=1,k=5e-3,A=A_max,B=1000)
    if(ChromInfo[i,2] !='Nhek')
      starts = c(g=1,k=5e-5,A=A_max,B=1000)                                       

    QuantileData[[i]] = getQuantilMedians(mergedData3[!is.na(mergedData3[,ChromInfo[i,2]]),],column=ChromInfo[i,2],N=N)
    SDData[[i]] = getQuantilMedians(mergedData3[!is.na(mergedData3[,ChromInfo[i,2]]),],column=ChromInfo[i,2],N=N, func=sem)
    
    # Fit to exponential + read-through (3 free parameters)
    FitData[[i]] =  nls(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileData[[i]],start=starts,lower=c(1,0,A_min,1),upper=c(1,Inf,A_max,Inf),algo='port');
    
    # Fit to exponential with no read-through (2 free parameters)
    #FitData2[[i]] =  nls(median ~ MeanExp(LengthMedian,k,A), data=QuantileData[[i]],start=starts[2:3],lower=c(0,A_min),upper=c(Inf,A_max),algo='port');
    FitData2[[i]] = nls(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileData[[i]],start=c(starts[1],k=0.01,starts[3],B=0), lower=c(1,0,A_min,0),upper=c(1,Inf,A_max,0),algo='port');
    
    # Fit to Gamma with g=2, RT (3 params)
    FitData_RT2[[i]] =  nls(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileData[[i]],start=c(g=2,starts[-1]), lower=c(2,0,A_min,1),upper=c(2,Inf,A_max,Inf),algo='port');
    
    if(i < 9){
      # Fit to Gamma with g=2, no RT (2 params)
      FitData_noRT2[[i]] =  nls(median ~ MeanGamma2(LengthMedian,g,k,A,B), data=QuantileData[[i]],start=c(g=2,starts[2:3],B=0), lower=c(2,0,A_min,0),upper=c(2,Inf,A_max,0),algo='port');
    }
    
    #FitData3[[i]] =  nls(median ~ MeanGamma3(LengthMedian,g,k,A,B), data=QuantileData[[i]],start=c(g=1,k=5e-5,A=1/10000,B=5000),lower=c(1,0,1/100000000,1),upper=c(1,Inf,1/1000,Inf),algo='port');
  }
                                                                         
  
  #########################################################################
  ## Plotting                                                               
  #########################################################################
  
  plot.dev(sprintf("Spring2014_Human_Multiple_CellTypes_splicing%d.pdf",N),'pdf', height=9, width=9)  
  blueCol = rgb(0,174,239,max=255)
  purpleCol = 'purple'
  par(mfrow=c(3,3))
  for( x in 1:nrow(ChromInfo)){
    plotCI(log10(QuantileData[[x]][,1]), (QuantileData[[x]][,2]),SDData[[x]][,2],sfrac=0.01, gap=0, main=ChromInfo[x,2], pch=19, xlab= 'kb from intron to poly(A) site' , ylim=c(0,1),ylab=expression(paste('median ',Sigma)),cex=.5,cex.main=1.2,xaxt='n',cex.axis=.8,yaxt='n', barcol=grey.colors(1,start=0.8))
    formatAxes_kb()
    
    # Plot fits
    lines(log10(xx<- 10^seq((0),log10(2000000),length.out=1000)),MeanGamma2(xx,coef(FitData[[x]])[1],coef(FitData[[x]])[2],coef(FitData[[x]])[3], coef(FitData[[x]])[4]),type='l',col=blueCol,lwd=2)  
    
    if (x< 9){
        lines(log10(xx<- 10^seq((0),log10(2000000),length.out=1000)),MeanGamma2(xx,coef(FitData2[[x]])[1],coef(FitData2[[x]])[2],coef(FitData[[x]])[3], coef(FitData2[[x]])[4]),type='l',col=blueCol,lwd=2,lty=2)    
      
      lines(log10(xx),MeanGamma2(xx,coef(FitData_RT2[[x]])[1],coef(FitData_RT2[[x]])[2],coef(FitData_RT2[[x]])[3], coef(FitData_RT2[[x]])[4]),type='l',col=purpleCol,lwd=1, lty=3)  
      
      lines(log10(xx),MeanGamma2(xx,coef(FitData_noRT2[[x]])[1],coef(FitData_noRT2[[x]])[2],coef(FitData_noRT2[[x]])[3], coef(FitData_noRT2[[x]])[4]),type='l',col=purpleCol,lwd=1,lty=2)  
    }
    
    # Put some text up:
    text(4,0.3,sprintf("kb/splice: %.0f\n RT: %.0f kb\n A: %.2f\n\n%.1f%% of genes highly\n cleaved",
      1/coef(FitData[[x]])[2],coef(FitData[[x]])[4],coef(FitData[[x]])[3], polyA80[x]*100), adj=0)
    
   }
  plot.off()
     
  t(sapply(FitData,coef))->fitCoefs
  t(sapply(FitData2,coef))->fitCoefs2
  rownames(fitCoefs2)= rownames(fitCoefs) <- ChromInfo[,2]
  save(QuantileData, mergedData3, FitData,FitData2,  fitCoefs,fitCoefs2, file=sprintf('Multiple_CellTypes2_splicing%d.RData',N))
  fitCoefs
  fitCoefs2
}
#1.000e+00 2.904e-04 8.628e-01 5.344e+03 
  
### 1000

   #               g            k         A         B
#K562_Chromatin    1 0.0001147089 0.8153261 13425.457
#K562_Nucleus A(-) 1 0.0002967239 0.8613706  5119.237
#K562_Nucleus A(+) 1 0.0008087497 0.9521556  5483.777
#Gm12878           1 0.0001600168 0.8408256  9508.821
#H1hesc            1 0.0003198619 0.8794459 15690.865
#Helas3            1 0.0004634191 0.8651435  5004.494
#Hepg2             1 0.0003420755 0.8490805  5516.401
#Huvec             1 0.0001375649 0.8305603 21923.495
#Nhek              1 0.0052768552 0.9198924  5375.325
#>   fitCoefs2
#                            k         A
#K562_Chromatin    0.001482321 0.6458452
#K562_Nucleus A(-) 0.001413117 0.7745113
#K562_Nucleus A(+) 0.016036956 0.8869726
#Gm12878           0.001165324 0.7145271
#H1hesc            0.018981906 0.8037911
#Helas3            0.003670126 0.7889294
#Hepg2             0.002130068 0.7659167
#Huvec             0.006018350 0.6985735
#Nhek              0.109483474 0.9109959
#
#
  
### 100:
#K562_Chromatin    1 0.0001119933 0.8235266 13462.346
#K562_Nucleus A(-) 1 0.0003029538 0.8640188  5036.963
#K562_Nucleus A(+) 1 0.0008254840 0.9542015  5708.252
#Gm12878           1 0.0001655685 0.8445809  8915.585
#H1hesc            1 0.0002718087 0.8940079 19376.124
#Helas3            1 0.0004709385 0.8713034  5006.756
#Hepg2             1 0.0003563707 0.8494094  5377.424
#Huvec             1 0.0001376448 0.8341237 22323.350
#Nhek              1 0.0052806117 0.9231211  6665.803
#>   fitCoefs2
#                            k         A
#K562_Chromatin    0.001412625 0.6505984
#K562_Nucleus A(-) 0.001442470 0.7776923
#K562_Nucleus A(+) 0.016607332 0.8916269
#Gm12878           0.001104280 0.7219635
#H1hesc            0.015812655 0.8163814
#Helas3            0.003854257 0.7946783
#Hepg2             0.002245259 0.7677393
#Huvec             0.006580220 0.7015959
#Nhek              0.109941763 0.9152730
#           
           
           
           
  
  
##t(sapply(FitData3,coef))->fitCoefs3
#
#
#
#Chrom_A = fitCoefs[2,3]
#
#Fit_ChromA_noRT = nls(median ~ MeanExp(LengthMedian,k,A), data=QuantileData[[1]],start=c(k=2e-3,A=0.8134891),lower=c(0,0.8134891),upper=c(Inf,0.8134891),algo='port');
#
#i=1
#plot(log10(QuantileData[[x]][,1]), (QuantileData[[x]][,2]), main=ChromInfo[x,2], pch=19, xlab= 'LOG10 Distance to poly(A) site' , ylim=c(0,1),ylab=expression(paste('median',phi)))
#lines(log10(xx<- 10^seq((0),log10(2000000),length.out=1000)),MeanGamma2(xx,coef(FitData[[x]])[1],coef(FitData[[x]])[2],coef(FitData[[x]])[3], coef(FitData[[x]])[4]),type='l',col='red',lwd=2)  
#lines(log10(xx<- 10^seq((0),log10(2000000),length.out=1000)),MeanGamma2(xx,coef(FitData[[x]])[1],coef(FitData[[x]])[2],coef(FitData[[x]])[3], coef(FitData[[x]])[4]),type='l',col='red',lwd=2)  
# 

## Coefs 1000:
#   g            k         A         B
#[1,] 1 0.0003231601 0.8604630  4410.095
#[2,] 1 0.0001256127 0.8134885 11647.467  #7961 bp/event; 110 seconds  
#[3,] 1 0.0001804777 0.8330515  7612.104
#[4,] 1 0.0006159702 0.8646847  6668.742
#[5,] 1 0.0005394481 0.8641768  4021.592
#[6,] 1 0.0004140760 0.8433461  4048.605
#[7,] 1 0.0002392780 0.8081242 11142.886
#[8,] 1 0.0089308576 0.9191437  1903.609

#               k         A
#[1,] 0.001647357 0.7689364
#[2,] 0.001577696 0.6502207   # 634 bp/splicing event ; 8.8 sec. half-life @3 kb/min
#[3,] 0.001275278 0.7058163
#[4,] 0.015703446 0.7962496
#[5,] 0.004412396 0.7864589
#[6,] 0.002510646 0.7600609
#[7,] 0.006559696 0.6938854
#[8,] 0.204472603 0.9085582
#

#

# Coefs 100:
#     g            k         A         B
#[1,] 1 0.0003208436 0.8629205  4540.253
#[2,] 1 0.0001183082 0.8231899 12363.069  #8452 bp/ splicing event  ; 117 seconds intron half-life
#[3,] 1 0.0001881989 0.8349544  7247.605
#[4,] 1 0.0005260339 0.8737109  8449.019
#[5,] 1 0.0005077042 0.8682972  4563.735
#[6,] 1 0.0004182649 0.8451115  4094.696
#[7,] 1 0.0002263597 0.8134832 12032.634
#[8,] 1 0.0104906624 0.9204022  2420.488

        #       k         A
#[1,] 0.001736229 0.7686981
#[2,] 0.001581033 0.6523102 # 632 bp/splicing event ; 8.8 sec. half-life @3 kb/min
#[3,] 0.001281080 0.7105821
#[4,] 0.018673272 0.8005563
#[5,] 0.005104690 0.7854219
#[6,] 0.002658559 0.7607866
#[7,] 0.006980512 0.6960977
#[8,] 0.105847535 0.9140437






























    