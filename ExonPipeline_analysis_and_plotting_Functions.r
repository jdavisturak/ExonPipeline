################################################################################################################################
####  JCDT 10/2/2012
####  These functions are to be used after all the files have been created using the functions in ExonPipeline_functions
####
################################################################################################################################

lmText = function(DATA,doLog=T){
	if(doLog){
		myLM <- lm(log10(Value) ~ isLast, data = data.frame(DATA))
	}else		
		myLM <- lm((Value) ~ isLast, data = data.frame(DATA))
	
	values <- summary(myLM)$coefficients[2,c(1,4)] # First value is coefficient, second is p value

	return(sprintf(" %.1e (p = %.1e)",values[1],values[2]))
}


#### Function to generate a simple matrix of the data with 'isLast' column
formatData = function(dat){
    OUT = cbind(Value = dat[,2], isLast = 0)#
    OUT[ grep("_1$",(dat[,1])) ,2] <- 1
    
    OUT
}



##### Function to make ONE set of boxplots
makeBoxplot = function(DATA, mainText ='',doLog=F,NAMES=c('Non-last','last'), addStats=T,at=NULL,labels=NULL,notch=T, ...){
  if (doLog){
    DATA[,1] = log10(DATA[,1])
  }
  yaxt = ifelse(is.null(at),'y','n'); # Don't plot the axis if user supplied a value for 'at'
  main = ifelse(addStats,paste(mainText,lmText(DATA,FALSE),sep="\n"),mainText) # if addStats is set to TRUE, add the lmText() data
  
  boxplot(Value ~ isLast, data = DATA,main=main,yaxt=yaxt,names=NAMES,notch=notch, ...)
 
  # Possibly label the yticks 
  if(!is.null(at)) axis(2,at=at,lab=labels)
}   


##### Function to make ONE set of barplots
makeBarplot = function(DATA, mainText ='', doLog=F, YLIM=NULL, addStats=TRUE,errorcol='black', ...){
#
#  if (doLog){
#    DATA[,1] = log10(DATA[,1])
#  }

  stats = getStats(DATA[,1],DATA[,2])  
  up = stats[,1] + stats[,2]
  down = stats[,1] - stats[,2]
  
  
  if (doLog){ 
    # Take the log10 now
    stats[,1] = log10(stats[,1])
    stats[,3] = log10(stats[,3])

    up = log10(up); print(up)
    stats[,2] = up - stats[,1] # this is so the YLIM will work correctly
    down = log10(down); print(down)
  }

  if(is.null(YLIM)) YLIM = c(0, max(rowSums(stats[,1:2]))*1.2)
  
  if(addStats) mainText=paste(mainText,lmText(DATA,FALSE),sep="\n")

  x=barplot(stats[,1], main=mainText, names=c('Non-last','last'), ylim=YLIM,...)
  arrows(x,stats[,1],x, up ,angle=90,length=0.15,col=errorcol,lwd=2,gap=0);
  arrows(x,stats[,1],x, down,angle=90,length=0.15,col=errorcol,lwd=2,gap=0);
  abline(h=0,lwd=1)
  
}

##### Function to make ONE set of barplots: take LOG first, if applicable
makeBarplot2 = function(DATA, mainText ='', doLog=F, YLIM=NULL, addStats=TRUE,errorcol='black',  names=c('Non-last','last'),...){
  #
    if (doLog){
      DATA[,1] = log10(DATA[,1])
    }
  
  stats = getStats(DATA[,1],DATA[,2])  
  up = stats[,1] + stats[,2]
  down = stats[,1] - stats[,2]
  
  if(is.null(YLIM)) YLIM = c(0, max(rowSums(stats[,1:2]))*1.2)
  
  if(addStats) mainText=paste(mainText,lmText(DATA,FALSE),sep="\n")
  
  x=barplot(stats[,1], main=mainText,names=names, ylim=YLIM,...)
  arrows(x,stats[,1],x, up ,angle=90,length=0.15,col=errorcol,lwd=2,gap=0);
  arrows(x,stats[,1],x, down,angle=90,length=0.15,col=errorcol,lwd=2,gap=0);
  abline(h=0,lwd=1)
  
}



#### Function to split data up by vector 'splitVector': Return a GROUP vector; X is is group(i) if split(i-1) < X <= split(i)

splitByVector=function(dat,splitsVector){
  groups <- rep(NA,length(dat))  # Iniatialize a vector of NAs
  groups[dat <= splitsVector[1]] <- 1  # First group: test if it's SMALLER (or =)
  
  for(x in 1:length(splitsVector))		# All other groups: test if it's BIGGER
    groups[dat > splitsVector[x]] <- x+1  
  groups
}

quantileGroups = function(dat,N=10){
   splitByVector(dat,quantile(dat,na.rm=T,(1:(N-1))/N))
}

getQuantilMedians = function(myData,col1='Dist2End',column='splicing0',N=100){
  myData$Group = quantileGroups(myData[,col1],N)
  data.frame(LengthMedian = tapply(myData[,col1],myData$Group,median,na.rm=T), median = tapply(myData[,column],myData$Group,median,na.rm=T))
}



### Function to convert raw data to something merge-able
getLastData = function(DAT,NAME){ # DAT must be a data frame
  names(DAT)[2] = NAME
  
  # Get Lasts
  DAT = DAT[grep("_1$",DAT[,1]),]
  DAT$UniqueID = as.numeric(sub("(>|^)([0-9]*)_.*","\\2",DAT[,1]))
  DAT[,-1]
}



### Function to keep only the most downstream exons
mostDownstream = function(dat,GeneName='GeneName'){
  
  ## Make a number that is bigger if most downstream: call it 'terminus'
  dat$terminus = dat$end  
  dat$terminus[dat$strand=='-'] = -dat$end[dat$strand=='-']

  # concatenate to ID, so that I can use tapply
  unique_string = paste(dat$terminus,dat$UniqueID,sep="_")

  # For each GeneName, return the ID of the most downstream element
  ID_most_down = tapply(unique_string,dat[,GeneName], function(x){
    id = sub(".*_([0-9]*)","\\1",x)
    terminus = as.numeric(sub("^(.[0-9]*)_.*","\\1",x))
    id[terminus == max(terminus)]  
  })

  # Subset for ONLY those IDs that have the most downstream exon
  
  dat[dat$UniqueID %in% ID_most_down,-which(dimnames(dat)[[2]]=='terminus')]
}



##### Function to merge data.

mergeLastData = function(exons,energyData,donorData,acceptorData,bins,HKlist){
  ### Subset LAST data  
  exons = exons[exons$isLast == 1,] 
  energyData = energyData[grep("_1$",(energyData[,1])),]
  donorData = donorData[grep("_1$",(donorData[,1])),]
  acceptorData = acceptorData[grep("_1$",(acceptorData[,1])),]
  
  ### Get Exon Length
  exons$length = abs(exons$start-exons$end)
  
  ### Keep only the exons whose terminus is most downstream for that gene
  exons= mostDownstream(exons)
  
  ### Merge
  OUT = merge(merge(merge(exons,getLastData(energyData,"Energy"),by='UniqueID'),getLastData(donorData,'Donor'),by='UniqueID'),getLastData(acceptorData,'Acceptor'),by='UniqueID')
  
  ## Add group ID
  OUT$Group = splitByVector(OUT$length,bins) # split into 5 groups
  
  ## Add HK status
  OUT$isHK = toupper(OUT$GeneName) %in% HKlist$Gene.name
  
  OUT   
}

#### Function to get mean, SEM based on group
getStats = function(dat,group,doPrint=TRUE){
	# print out lm fit 
	if (doPrint)
		print(summary(lm(dat ~ group))$coef)


  means = tapply(dat,group,mean,na.rm=T)
  medians = tapply(dat,group, median, na.rm=T)
  SEM = tapply(dat,group,function(x)sd(x,na.rm=T)/sqrt(length(which(!is.na(x)))))  
  data.frame(means,SEM,medians)  
}

############## Function to plot binned data

plotSplits = function(dat,column,main,Group='Group',avg='mean',...){
  require(gplots)
  stats = getStats(dat[,column],dat[,Group])  
  plotCI(stats[,ifelse(avg=='mean',1,3)],uiw=stats$SEM, main=main,...)
}



################################################

##### Function to divide in 2 by a HK status, then plot twice for HK vs Non
plotSplitsHK = function(dat,column,main,Group='Group',names=NULL, checkColumn='isHK', YLIM=NULL, LTY=1,avg='mean',colors=c('red','blue'), lwd=2,sfrac=0.02,type='l',...){
  require(gplots)
  statsHK = getStats(dat[dat[,checkColumn]==1,column],dat[dat[,checkColumn],Group])  
  statsNon = getStats(dat[dat[,checkColumn]==0, column],dat[dat[,checkColumn]==0,Group])
  xaxt=ifelse(is.null(names),'s','n')
   
  avgCol = ifelse(avg=='mean',1,3); # This uses column 1 for mean, 3 for median
  
  if(is.null(YLIM)){
	  YLIM <- c(min(c(min(statsNon[,avgCol]-statsNon[,2]),min(statsHK[,avgCol]-statsHK[,2]))), max(c(max(statsNon[,avgCol]+statsNon[,2]),max(statsHK[,avgCol]+statsHK[,2]))))
  }  
    
  col1=colors[1]; 
  col2=colors[2];
  plotCI(statsHK[,avgCol],uiw=statsHK$SEM, main=main,ylim=YLIM,type='n',lwd=lwd,col=col1,barcol=col1,xaxt=xaxt,gap=0,sfrac=sfrac,labels='',...)
  par(new=T)
  plotCI(statsNon[,avgCol],uiw=statsNon$SEM, ylim=YLIM,lwd=lwd,type='n',col=col2,barcol=col2,xaxt='n',xlab='',ylab='',yaxt='n',gap=0,labels='',sfrac=sfrac)
  
  if(type=='l' | type=='b'){
    lines(statsHK[,avgCol],ylim=YLIM,lwd=lwd,col=col1,xaxt='n',xlab='',ylab='',yaxt='n', lty=LTY)
    lines(statsNon[,avgCol],ylim=YLIM,lwd=lwd,col=col2,xaxt='n',xlab='',ylab='',yaxt='n', lty=LTY)
  }
  
  
  #legend('topright',col=c(col1,col2),lwd=lwd,leg=c("HK","nonHK"))
  if(!is.null(names))
  	axis(1,at=1:length(names),lab=names)

  invisible(list(statsHK = statsHK, statsNon = statsNon))
  
}


plotFig4Like = function(DATA,name,group,main, AVG = 'median'){
  png(sprintf('%s%s.png',settings$CommonName,name),height=1200,width=1200)
  par(mfrow=c(2,2),cex.axis=1.2,cex=1.3,cex.main=1.1,las=3,mar=c(6.1 ,4.1 ,4.1, 2.1)) ; options(warn=-1) -> W
  xlabs = tapply(mergedData$length,mergedData$Group,function(x)paste(range(x),collapse='-'))
  
  plotSplitsHK(DATA,'PredSplice_strength',sprintf('CTS efficiency when modulating splicing rates\nAs a function of %s',main),avg=AVG,Group=group,xlab='',ylab='',names=xlabs,YLIM=c(0.0,1))
  plotSplitsHK(DATA,'PredSplice','',Group=group,add=T,LTY=2,YLIM=c(0.0,1),xlab='Last Exon Length',avg=AVG,ylbias=-1)
  
  plotSplitsHK(DATA,'PredSplice_pause',sprintf('CTS efficiency when modulating pausing rate in last exon\nAs a function of %s',main),avg=AVG,Group=group,xlab='',ylab='',names=xlabs,YLIM=c(0.0,1))
  plotSplitsHK(DATA,'PredSplice','',Group=group,add=T,LTY=2,YLIM=c(0.0,1),xlab='Last Exon Length',avg=AVG,ylbias=-1)
  
  plotSplitsHK(DATA,'PredSplice_runon',sprintf('CTS efficiency when adding Read-through rates\nAs a function of %s',main),avg=AVG,Group=group,xlab='',ylab='',names=xlabs,YLIM=c(0.0,1))
  plotSplitsHK(DATA,'PredSplice','',Group=group,add=T,LTY=2,YLIM=c(0.0,1),xlab='Last Exon Length',avg=AVG,ylbias=-1)
  
  plotSplitsHK(DATA,'PredSplice_ALL',sprintf('CTS efficiency when adding all features\n(converstion Factors %s)\nAs a function of %s',conversionFactor, main),avg=AVG,Group=group,xlab='',ylab='',names=xlabs,YLIM=c(0.0,1))
  plotSplitsHK(DATA,'PredSplice','',Group=group,add=T,LTY=2,YLIM=c(0.0,1),xlab='Last Exon Length',avg=AVG,ylbias=-1)
                   
  options(W); dev.off() 
}


formatAxes = function(){
axis(1,-2:6,sapply(-2:6,function(x) format(10^(x),scientific=F)),las=2)
axis(2,seq(0,1,.1),lab=NA,tcl=-.3)
axis(2,seq(0,1,.5),tcl=-.5)
}


formatAxes_kb = function(){  #format for kb!
axis(1,-2:6,sapply((-2:6)-3,function(x) format(10^(x),scientific=F)),las=2)
axis(2,seq(0,1,.1),lab=NA,tcl=-.3)
axis(2,seq(0,1,.5),tcl=-.5)
}

formatAxes_min = function(Elong){  #format for time to elongate, in minutes.
  axis(1,-2:6,sapply((-2:6),function(x) sprintf("%.1f",10^(x)/Elong)),las=2)
  axis(2,seq(0,1,.1),lab=NA,tcl=-.3)
  axis(2,seq(0,1,.5),tcl=-.5)
}

