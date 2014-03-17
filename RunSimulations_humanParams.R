################################################################################################################
################################################################################################################
#### Simulate model!! 
################################################################################################################
################################################################################################################

load(sprintf("%s/ExonPipeline/hg19/Multiple_CellTypes2_splicing100.RData",homedir))
load(sprintf("%s/PlottedData_new3.RData",dataDir)) 


mergedData = mergedData[order(mergedData$UniqueID,mergedData$FeatureCount),]

# Assign each gene the median stability from its group:
S=with(subset(mergedData,isLast==0), tapply(Stability,list(GeneSize,intronCountGroups),median, na.rm=T))
S2 = data.frame(S=as.vector(S),GeneSize=rep(1:4,5),intronCountGroups=  rep(1:5, each=4))
mergedData$GroupStability = S2$S[match(with(mergedData,paste(GeneSize,intronCountGroups )), with(S2,paste(GeneSize,intronCountGroups )) )]
mergedData_forFit = data.frame(Stability=mergedData$GroupStability, Dist2End=mergedData$Dist2End)

column=2
fitData_to_model = coef(FitData[[column]])
fitData_to_model2 = coef(FitData2[[column]])


# Using the fitted parameters
Elong = 1/fitData_to_model[2]
Elong_noRT = 1/fitData_to_model2[2]

readThrough = fitData_to_model[4]

Stability_exp=1
Stability_div = max(S,na.rm=T)*1.20
extraDist = rep(147,nrow(mergedData))

K_splice = 1


# For simplified version: multiplication factor = 1 + 1/Acceptor_div
mult_factor = 2
Acceptor_div = 1/(mult_factor-1) 
K_adjusted=K_splice


myEfficiencies = tapply(1:nrow(mergedData), mergedData$UniqueID, function(x)simulateGene(mergedData_forFit[x ,], extraDist=extraDist,Elong=Elong,readThrough=readThrough, K_splice=K_adjusted,Stability_exp=Stability_exp,Stability_div=Stability_div,Acceptor_div=Acceptor_div))

myEfficiencies_noSplice = tapply(1:nrow(mergedData), mergedData$UniqueID, function(x)simulateGene(mergedData_forFit[x,], extraDist=extraDist,Elong=Elong,readThrough=readThrough, K_splice=K_splice,Acceptor_div=Inf,Stability_exp=Stability_exp,Stability_div=Stability_div))

myEfficiencies_noPause = tapply(1:nrow(mergedData), mergedData$UniqueID, function(x)simulateGene(mergedData[x,], extraDist=extraDist,Elong=Elong,readThrough=readThrough, K_splice=K_adjusted,Stability_div=Inf,Stability_exp=Stability_exp,Acceptor_div=Acceptor_div))

myEfficiencies_plain = tapply(1:nrow(mergedData), mergedData$UniqueID, function(x)simulateGene(mergedData[x,], extraDist=extraDist,Elong=Elong,readThrough=readThrough, K_splice=K_splice,Stability_div=Inf, Acceptor_div=Inf,Stability_exp=Stability_exp))
myEfficiencies_noRT = tapply(1:nrow(mergedData), mergedData$UniqueID, function(x)simulateGene(mergedData[x,], extraDist=extraDist,Elong=Elong,readThrough=0, K_splice=K_splice,Stability_div=Inf, Acceptor_div=Inf,Stability_exp=Stability_exp))

# Get other info
genes.isHK = names(myEfficiencies) %in% subset(mergedData,isHK==1,UniqueID)[,1]
genes.Size = mergedData$GeneSize[match(names(myEfficiencies), mergedData$UniqueID)]
genes.IntronCount = mergedData$MaxFeature[match(names(myEfficiencies), mergedData$UniqueID)]
efficiencyData = data.frame(UniqueID=names( myEfficiencies), isHK=genes.isHK,GeneSize=genes.Size,MaxFeature=genes.IntronCount,intronCountGroups =  splitByVector(genes.IntronCount,c(2.5,4.5,7.5,19.5)),
                            myEfficiencies_noRT,myEfficiencies_plain,myEfficiencies_noSplice,myEfficiencies_noPause,myEfficiencies)


### PLOT better ####
myCols5 = cols=(brewer.pal(9,"YlGnBu"))[5:9]; pch=16
library(reshape2)
eff2 = melt(efficiencyData[,c(!grepl('(RT)|(plain)',dimnames(efficiencyData)[[2]]))],id=c('UniqueID','isHK','GeneSize','MaxFeature','intronCountGroups'))
plainMedians = sapply(c(1,4),function(size)sapply(1:5,function(num)median(subset(efficiencyData,intronCountGroups==num & GeneSize==size)$myEfficiencies_plain, na.rm=T)))


plot.dev(sprintf("GenomeSimulations_numIntrons_by_geneSize_nucElong_147_noAdjust_median%.3f_Extreme_only_3C.pdf",Stability_div),'pdf',height=2,width=5)  
cex=0.6;width=0.4;
par(cex=cex,mai=c(.02,.3,.02,.02),pch=pch)

at = c(c(1:5,7:11), 13 + c(1:5,7:11), 26  + c(1:5,7:11))
x=outer(c(1,1,NA),at)+rep(.4*c(-width,width,NA),30)
boxplot(value ~ intronCountGroups +GeneSize + variable, data= subset(eff2,GeneSize %in% c(1,4)) , at=at,boxwex=width, col=myCols5,pch=20,xaxt='n')
arrows(x[1,],rep(plainMedians,3),x[2,],rep(plainMedians,3),col='maroon',angle=135,code=3,lwd=2,length=width/Inf)
plot.off()


plot.dev("GenomeSimulations_readThrough_comparisons_chromElong_147_noAdjust_2only_2D.pdf",'pdf',width=5,height=1.5)
par(mfcol=c(1,1),cex=0.5,mai=c(0.3,0.2,0.01,0.01))
notch=F; wid=.21;bw=0.1085; offset = 1.4; xlim=c(0.5,9.2); widthFraction = 0.4
cols =c(brewer.pal(9,"YlGnBu")[5:9],brewer.pal(9,"OrRd")[5:9])
for (HK in list(0:1)){
  
  for(sizes in 1:4){
    COL = cols[c(1:5) + HK[1]*5]
    txt = paste(ifelse(length(HK)==1,ifelse(HK,'HK:','Regulated'),'All:'))#, c('Short','Med short','Med long','Long')[sizes])
    x<-subset(efficiencyData,GeneSize==sizes &isHK%in%HK)
    myCats = sort(unique(x$intronCountGroups))
    
    main=''
    # Decide which variable to test
    Form1 = myEfficiencies_noRT ~ intronCountGroups
    Form3 = myEfficiencies_plain ~ intronCountGroups
    
    Xshift = offset*(sizes-1 )
    at = 2*(myCats*wid+Xshift) - (sizes-1)*wid*2
    lwd=1
    b=boxplot(Form1, data=x,main=main,border=COL[myCats],ylim=c(0,1),lwd=1.5,staplewex=0,range=0.00001,outline=F,boxwex=bw,at=at-wid*widthFraction,xlim=xlim,names=NA,yaxt=ifelse(sizes==1,'s','n'),notch=notch,add=sizes-1,xaxt='n')
    
    b=boxplot(Form3, data=x,main=main,col=COL[myCats],border='black',boxwex=bw,at=at+wid*widthFraction,names=NA,notch=notch,add=T,xaxt='n',lwd=lwd,yaxt='n')
  }
  axis(1,at=2*(3*wid+(0:3)*offset)-(0:3)*wid*2, lab = c('Short','Med short','Med long','Long')) 
}
plot.off()
