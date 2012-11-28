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
makeBoxplot = function(DATA, mainText ='', doLog=F){
  if (doLog){
    DATA[,1] = log10(DATA[,1])
  }
  
  boxplot(Value ~ isLast, data = DATA,main=paste(mainText,lmText(DATA,FALSE),sep="\n"),names=c('Non-last','last'),notch=T)
}   


##### Function to make ONE set of barplots
makeBarplot = function(DATA, mainText ='', doLog=F, ...){

  if (doLog){
    DATA[,1] = log10(DATA[,1])
  }

  stats = getStats(DATA[,1],DATA[,2])  
  YLIM = c(0, max(rowSums(stats))*1.2)
  
  x=barplot(stats[,1],main=paste(mainText,lmText(DATA,FALSE),sep="\n"),names=c('Non-last','last'),...)
  #plotCI(x=x,stats[,1],uiw=stats[,2],add=T,barcol='red',lwd=5,col=NULL,sfrac=5) ## Add error bars
  arrows(x,stats[,1],x,stats[,1]+stats[,2],angle=90,length=0.15,lwd=2);
  arrows(x,stats[,1],x,stats[,1]-stats[,2],angle=90,length=0.15,lwd=2);

}




#### Function to split data up by vector 'splitVector': Return a GROUP vector; X is is group(i) if split(i-1) < X <= split(i)

splitByVector=function(dat,splitsVector){
  groups <- rep(NA,length(dat))  # Iniatialize a vector of NAs
  groups[dat <= splitsVector[1]] <- 1  # First group: test if it's SMALLER (or =)
  
  for(x in 1:length(splitsVector))		# All other groups: test if it's BIGGER
    groups[dat > splitsVector[x]] <- x+1  
  groups
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
  SEM = tapply(dat,group,function(x)sd(x,na.rm=T)/sqrt(length(which(!is.na(x)))))  
  data.frame(means,SEM)  
}

############## Function to plot binned data

plotSplits = function(dat,column,main,Group='Group',...){
  require(gplots)
  stats = getStats(dat[,column],dat[,Group])  
  plotCI(stats$mean,uiw=stats$SEM, main=main,...)
}


################################################

##### Function to divide in 2 by a HK status, then plot twice for HK vs Non
plotSplitsHK = function(dat,column,main,Group='Group',names=NULL, YLIM=NULL, LTY=1,...){
  require(gplots)
  statsHK = getStats(dat[dat$isHK==1,column],dat[dat$isHK==1,Group])  
  statsNon = getStats(dat[dat$isHK==0, column],dat[dat$isHK==0,Group])
  xaxt=ifelse(is.null(names),'s','n') 
  
  if(is.null(YLIM)){
	  YLIM <- c(min(c(min(statsNon[,1]-statsNon[,2]),min(statsHK[,1]-statsHK[,2]))), max(c(rowSums(statsHK),rowSums(statsNon))))
  }  
    
  lwd=2; type='l';col1='red';col2='blue'
  plotCI(statsHK$mean,uiw=statsHK$SEM, main=main,ylim=YLIM,type='p',lwd=lwd,col=col1,barcol=col1,xaxt=xaxt,...)
  par(new=T)
  plotCI(statsNon$mean,uiw=statsNon$SEM, ylim=YLIM,lwd=lwd,type='p',col=col2,barcol=col2,xaxt='n',xlab='',ylab='',yaxt='n')
  if(type=='l' | type=='b'){
    lines(statsHK$mean,uiw=statsHK$SEM,ylim=YLIM,lwd=lwd,col=col1,xaxt='n',xlab='',ylab='',yaxt='n', lty=LTY)
    lines(statsNon$mean,uiw=statsNon$SEM,ylim=YLIM,lwd=lwd,col=col2,xaxt='n',xlab='',ylab='',yaxt='n', lty=LTY)
  }
  
  
  #legend('topright',col=c(col1,col2),lwd=lwd,leg=c("HK","nonHK"))
  if(!is.null(names))
  	axis(1,at=1:length(names),lab=names)

  invisible(list(statsHK = statsHK, statsNon = statsNon))
  
}
