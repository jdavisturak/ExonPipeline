### Do pile-up analysis of bigWig files.
# R 3.0.1
library('rtracklayer')

######
# Read in bigWig file:
#bw = import.bw("myBW.bigWig",asRangedData=FALSE)
#bw=import.bw("/home/RNAseq/ENCODE/wgEncodeSydhTfbsK562Pol2IggmusSig.bigWig",asRangedData=F)


# Import bed file:
# bed = import.bed("myBed.bed",asRangedData=FALSE)
# bed = import.bed("/home/home/jeremy/ExonPipeline/hg19/Human_exons_-73To73_last+.bed",asRangedData=FALSE)

# Overlap bed and bigwig file:

####### USE:

# overlap=findOverlaps(bed,bw)
# RangeScores = cbind(as.matrix(ranges(bw)),score(bw))
# S = computeScores(overlap, bed, RangeScores)


# For each bed file, generate a vector of bigwWig scores
bigWigVector = function(Ranges, bedStart =  min(Ranges[,1]), bedEnd = max(Ranges[,1])){
  ## Get the position in relation to the start (set Start=1)
  indices = as.vector(unlist(apply(Ranges,1,function(x) x[1]-bedStart + 1 + (0:(x[2]-1)))))
    ## Get the scores
  scores = as.vector(unlist(apply(Ranges,1,function(x) rep(x[3],x[2]))))
  data.frame(indices,scores)[indices > 0 & indices <= (bedEnd-bedStart+1) ,T,drop=F]
}


#apply(queryHits(overlap),subjectHits(overlap), bigWigVector)

computeScores = function(overlap, bed, RangeScores) sapply(1:length(bed), function(i){
  # Initialize scores
  scores = numeric(width(ranges(bed)[i]))
  
 #print(i) 
   
  # retreive overlap data
  if (length(which.reads <- subjectHits(overlap)[queryHits(overlap)==i]) > 0){
  
    hitData = bigWigVector(RangeScores[which.reads,T,drop=F], bedStart=start(b<-ranges(bed[i])), bedEnd=end(b))
    scores[hitData[,1]] = hitData[,2]
  }
  scores
  
  }) 

computeScores_Stranded = function(overlap, bed, RangeScores) sapply(1:length(bed), function(i){
  # Initialize scores
  scores = numeric(width(ranges(bed)[i]))
  
 #print(i) 
   
  # retreive overlap data
  if (length(which.reads <- subjectHits(overlap)[queryHits(overlap)==i]) > 0){
  
    hitData = bigWigVector(RangeScores[which.reads,T,drop=F], bedStart=start(b<-ranges(bed[i])), bedEnd=end(b))
    scores[hitData[,1]] = hitData[,2]
  }
  if(as.character(strand(bed)[i]) == '-')
    scores = rev(scores)
  scores
  
  }) 

scoreBW_bed = function(bw, bed, RangeScores){
 overlap = findOverlaps(bed,bw)
 computeScores_Stranded(overlap, bed,RangeScores)
}


make_beds_selection = function(allBeds,goodChrs = paste('chr',c(1:22,'X'),sep='')){
  allBeds2  = allBeds[seqnames(allBeds) %in% goodChrs]
  seqlevels(allBeds2) <- goodChrs
  BigWigSelection(allBeds2)
}






