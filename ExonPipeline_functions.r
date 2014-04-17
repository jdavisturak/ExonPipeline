### Functions to set-up and use my new pipeline.

############################################################################################################
### Function to read in the settings file ############
############################################################################################################
                 
setup = function(Folder=getwd(), chDir=FALSE){

  ## Read in the settings
  settings.txt = read.delim(sprintf('%s/settings.txt',Folder),head=FALSE)
  
  ##  Set up the settings object
  settings = data.frame(t(as.character(settings.txt[,1])),stringsAsFactors=FALSE)
  names(settings)[1:3] = c('genome','CommonName','LatinAbbrev') # these are data the use should supply
  
  ## If the settings.txt file resulted from a call to 'saveData', add the names of the data to this object
  if (length(settings) > 3)
  	names(settings)[-(1:3)] = sub("^# ","",settings.txt[-(1:3),2])
  
  ## Now add some other attributes:
  
  settings$thisDir = Folder
 
  # The name of the genome library loaded by the BSgenome package
  if(settings$LatinAbbrev=='Athaliana'){
  	settings$BSlib = sprintf("BSgenome.%s.TAIR.%s",settings$LatinAbbrev,settings$genome)
  }else{
  	settings$BSlib = sprintf("BSgenome.%s.UCSC.%s",settings$LatinAbbrev,settings$genome)
  }
 
  # Name of file to write later
  settings$NucTestFile = 'exons_forSeqToNuc.txt'
  
  class(settings) <- "ExonPipeline"
  
  ## move to directory
  if (chDir)
    setwd(Folder)


  return(settings)
}



############################################################################################################
### Function to load genome: downloads it if it's needed, otherwise just loads the BSgenome and creates a variable 'genome'
############################################################################################################

loadGenome = function (settings){
  ## First check it it's already installed:
  if(! settings$BSlib %in% installed.packages()[,1]){
    cat(sprintf("Installing %s:\n",settings$BSlib))

    ## install it:
    source("http://bioconductor.org/biocLite.R")
    biocLite(settings$BSlib)
    
  }

  ## Load the genome:
  eval(parse(text=sprintf("library(%s)",settings$BSlib)))
  
  ## Create a new object for the genome that has a common name
  eval(parse(text=sprintf("genome= %s",settings$LatinAbbrev)))
  
  return(genome)
}

############################################################################################################
### This function downloads and reads in the refGene data from refseq (via UCSC's mysql connection)
############################################################################################################

loadRefgene = function(settings){
    refGene = paste(settings$thisDir,'refGene.txt',sep='/')
  # First check if the refGene.txt file exists in this folder
  if (!file.exists(refGene)){
  	cat(sprintf("Downloading refgene:\n"))
    
    ## No file yet: download from UCSC:
    command = sprintf("echo 'use %s; select * from refGene' | mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A > %s/refGene.txt", settings$genome,settings$thisDir)  
    system(command)
  }

  ## Load the file, return
  read.delim(refGene,stringsAsFactors=FALSE)
}


############################################################################################################
### This function adds a UniqID to a refseq file
############################################################################################################

add_UniqueID = function(ref){
  if (!'UniqueID' %in% dimnames(ref)[[2]])
    ref$UniqueID = 1:nrow(ref)
  ref
}


############################################################################################################
### Wrapper function to load or create gene features file
############################################################################################################
loadFeatures= function(settings,refseq,genome){

  # First check if the features.txt file exists in this folder
  if (!file.exists("exons.txt")){
  	
  	# Find exon and intron starts, stops
    cat(sprintf("Parsing refseq file for features:\n"))
    features = retrieveFeatures(settings,refseq,genome)
    write.table(features$exons,  file=paste(settings$thisDir,'exons.txt',sep='/'),  sep="\t",quote=FALSE,row.names=FALSE)
    write.table(features$introns,file=paste(settings$thisDir,'introns.txt',sep='/'),sep="\t",quote=FALSE,row.names=FALSE)
        
   }else{
 	  cat(sprintf("Reading in features\n"))
   	# Read exons, introns in from txt files
  	exons   = read.delim(paste(settings$thisDir,'exons.txt',sep='/'),  stringsAsFactors=FALSE)
  	introns = read.delim(paste(settings$thisDir,'introns.txt',sep='/'),stringsAsFactors=FALSE)   
  	features = list(exons=exons, introns=introns)

   }
   
   features
}



################################################################################################################################################################################
### Function to return the locations of all exons and introns in their own rows from a refseq file (where one GENE is one row) ###
################################################################################################################################################################################

retrieveFeatures = function (settings,refseq,genome){
	require(Biostrings)
	require(GenomicRanges)
	 
	## First get rid of all genes with only 1 exon!
	refseq = refseq[refseq$exonCount > 1,]
	
  	# Initialize output lists
	exon_matrix = matrix(NA,sum(refseq$exonCount),10)
	intron_matrix = matrix(NA,sum(refseq$exonCount-1),10)
	

	exon_counter = 0
	intron_counter = 0

	for( i in 1:nrow(refseq)){
  		if(i%%1000==1)
  			cat(sprintf("Doing refseq %d\r",i))
  		# Retreive the 'exonStarts' and 'exonStops'
  		exonStarts = as.numeric(unlist(strsplit(refseq$exonStarts[i],",")))
  		exonEnds = as.numeric(unlist(strsplit(refseq$exonEnds[i],",")))
  		exonFrames = as.numeric(unlist(strsplit(refseq$exonFrames[i],",")))
  	
  	
    	# Check the strand and chrom, get other info
    	strand = refseq$strand[i]
    	CHR = refseq$chrom[i]
    	exonCount = refseq$exonCount[i]
    	intronCount = exonCount-1
    	uniqueID = refseq$UniqueID[i]
    	geneName = refseq$name2[i]
    	geneLength = tail(exonEnds,1) - exonStarts[1]    	
    	    	
    	## Based on the strand, decide which parts are the real 'starts' and 'ends'

    	if(strand == '-'){
    		# Negative strand: reverse the starts and stops (reverse order, and change definitions)
			  temp = exonStarts
		      exonStarts = rev(exonEnds)
			  exonEnds = rev(temp)
			  exonFrames = rev(exonFrames)
			  
			}
  		# Get the intron starts and stops
   		intronStarts = exonEnds[-length(exonEnds)]
   		intronEnds = exonStarts[-1]    		
		  

	
		## Save this data, along with identifying data
	  	
    	exon_matrix[exon_counter + (1:exonCount),] <- cbind(uniqueID=uniqueID, GeneName=geneName, FeatureCount=1:exonCount, isLast=c(rep(0,exonCount-1),1), CHR, exonStarts, exonEnds, strand,geneLength,exonFrames)
    	exon_counter = exon_counter + exonCount    	
    	
    	intron_matrix[intron_counter + (1:intronCount),] <- cbind(UniqueID=uniqueID, GeneName=geneName, FeatureCount=1:intronCount, isLast=c(rep(0,intronCount-1),1), CHR, intronStarts, intronEnds, strand,upstreamFrame=exonFrames[-length(exonFrames)],downstreamFrame=exonFrames[-1])
    	intron_counter = intron_counter + intronCount    
    	    	
   }
   cat(sprintf("Done                 \n"))
  
   # Add names to the matrix
   dimnames(exon_matrix)[[2]] = c('UniqueID','GeneName','FeatureCount','isLast','chr','start','end','strand','GeneLength','exonFrames')	 
   dimnames(intron_matrix)[[2]] = c('UniqueID','GeneName','FeatureCount','isLast','chr','start','end','strand','upFrame','downFrame')
   exon_matrix = data.frame(exon_matrix,stringsAsFactors=FALSE)
   intron_matrix = data.frame(intron_matrix,stringsAsFactors=FALSE)
    
   # Filter for chromosomes that are in the genome
   exon_matrix = exon_matrix[exon_matrix$chr %in% seqnames(genome),]
   intron_matrix = intron_matrix[intron_matrix$chr %in% seqnames(genome),]

   ## Make the starts/ends numeric   
   exon_matrix$start = as.numeric(exon_matrix$start)
   exon_matrix$end = as.numeric(exon_matrix$end)
   
   intron_matrix$start = as.numeric(intron_matrix$start)
   intron_matrix$end = as.numeric(intron_matrix$end)
    
   
   ## Return the data together in a list
   return ( list(exons=exon_matrix, introns=intron_matrix))
}















############################################################################################################
### Wrapper function to prepare exons and run seqToNuc
############################################################################################################


computeNucleosomeEnergy = function(settings,features,genome,WIDTH = 180){

    ## Get output file name 
    settings$nucTestFileFixed = sub(".txt$",sprintf("_fixed%d.txt",WIDTH),settings$NucTestFile)
    
    ## Export specific sequences 
    if (!file.exists(settings$nucTestFileFixed)){
      firstSeq = 1
      cat(sprintf("Retrieving sequences of exons for nucleosome test (starting with seq %d):\n", firstSeq))
   
      unique_features = removeDuplicateStarts(features$exons)
    
      seqs_for_nucleosomes = retrieveFixedSequences(settings, unique_features, genome, width=WIDTH, overhangUp=37, overhangDown=37, fromStart=TRUE,outputFile=settings$nucTestFileFixed,writeDuring=T,firstSeq=firstSeq)
    }

    ## Run seqToNuc.jar
    settings$nucTestFileFixed_out = sprintf('%s_out.txt',settings$nucTestFileFixed)

    if (!file.exists(settings$nucTestFileFixed_out)){
      cat("Running Seq To Nuc (on fixed exon lengths)\n")

      system(sprintf("java -Xmx1024m -jar ~/bin/seqToNuc.jar %s", settings$nucTestFileFixed))

    }     
    
    
    ## Compute the mean over the entire length.
    settings$nucTestFileFixed_mean = sub(".txt$","_mean.txt",settings$nucTestFileFixed)
    
    if (!file.exists(settings$nucTestFileFixed_mean)){
      cat("Averaging Seq To Nuc output\n")
     
      computeNucMean (settings)
    }
    

    return(settings)
}


############################################################################################################
### Wrapper function to prepare exons and run seqToNuc version 2
### This function start 3bp downstream from the beginning of the exon, so as not to confound the data with the splice site info
############################################################################################################

computeNucleosomeEnergy2 = function(settings,features,genome,WIDTH = 147){

    ## Get output file name 
    settings$nucTestFileFixed2 = sub(".txt$",sprintf("_fixed2_%d.txt",WIDTH),settings$NucTestFile)
    
    ## Export specific sequences 
    if (!file.exists(settings$nucTestFileFixed2)){
      firstSeq = 1
      cat(sprintf("Retrieving sequences of exons for nucleosome test 2 (starting with seq %d):\n", firstSeq))
   
      unique_features = removeDuplicateStarts(features$exons)
    
      if(WIDTH <= 177){# Here, simply take the previous data and select the subset.
        cmd = sprintf("awk 'NR%%2==1{print} NR%%2==0{print(substr($0,4,%d)) }' %s > %s", WIDTH,settings$nucTestFileFixed, settings$nucTestFileFixed2)
        system(cmd) 
             
      }else{
      
        seqs_for_nucleosomes = retrieveFixedSequences(settings, unique_features, genome, width=WIDTH, overhangUp=-3, overhangDown=0, fromStart=TRUE,outputFile=settings$nucTestFileFixed2,writeDuring=T,firstSeq=firstSeq)
      }  
    }

    ### Run BasicSeqToNuc2.jar
    settings$nucTestFileFixed2_out = sprintf('%s_out.txt',settings$nucTestFileFixed2)

    if (!file.exists(settings$nucTestFileFixed2_out)){
      cat("Running Seq To Nuc 2 (on fixed exon lengths)\n")

      system(sprintf("java -Xmx1024m -jar ~/bin/BasicSeqToNuc.jar %s", settings$nucTestFileFixed2))

    }     
    
    
    ## Compute the mean over the entire length.
    settings$nucTestFileFixed2_mean = sub(".txt$","_mean.txt",settings$nucTestFileFixed2)
    
    if (!file.exists(settings$nucTestFileFixed2_mean)){
      cat("Averaging Seq To Nuc output\n")
     
      computeNucMean (settings, settings$nucTestFileFixed2_out, settings$nucTestFileFixed2_mean)
    }
    
    ## Normalize the nucleosome splicing score
    settings$nucNormScore = sub("_mean.txt$","_norm.txt",settings$nucTestFileFixed2_mean)
   
    if (!file.exists(settings$nucNormScore)){
      cat("Normalizing Seq To Nuc output\n")
      normalizeToText(settings$nucTestFileFixed2_mean, settings$nucNormScore)
    }

    return(settings)
}


############################################################################################################
### Wrapper function to prepare exons and run seqToNuc version 2
### 4/2/13 This is constrainted fixed length: limits the length of the sequence to the exon extent, and does not extend past the end of the exon, even if it means using fewer than WIDTH bases.
### This function start 3bp downstream from the beginning of the exon, so as not to confound the data with the splice site info
############################################################################################################

computeNucleosomeEnergy3 = function(settings,features,genome,WIDTH = 147){

    ## Get output file name 
    settings$nucTestFileConstrained = sub(".txt$",sprintf("_constrained_%d.txt",WIDTH),settings$NucTestFile)
    
    ## Export specific sequences 
    if (!file.exists(settings$nucTestFileConstrained)){
      firstSeq = 1
      cat(sprintf("Retrieving sequences of exons for nucleosome test 3 (starting with seq %d):\n", firstSeq))
   
      unique_features = removeDuplicateStarts(features$exons)
    
       seqs_for_nucleosomes = retrieveFixedSequences_min(settings, unique_features, genome, width=WIDTH, overhangUp=-3, overhangDown=0, fromStart=TRUE,outputFile=settings$nucTestFileConstrained,writeDuring=T,firstSeq=firstSeq)
        
    }

    ### Run BasicSeqToNuc2.jar
    settings$nucTestFileConstrained_out = sprintf('%s_out.txt',settings$nucTestFileConstrained)

    if (!file.exists(settings$nucTestFileConstrained_out)){
      cat("Running Seq To Nuc 2 (on fixed exon lengths)\n")

      system(sprintf("java -Xmx1024m -jar ~/bin/BasicSeqToNuc.jar %s", settings$nucTestFileConstrained))

    }     
    
    
    ## Compute the mean over the entire length.
    settings$nucTestFileConstrained_mean = sub(".txt$","_mean.txt",settings$nucTestFileConstrained)
    
    if (!file.exists(settings$nucTestFileConstrained_mean)){
      cat("Averaging Seq To Nuc output\n")
     
      computeNucMean (settings, settings$nucTestFileConstrained_out, settings$nucTestFileConstrained_mean)
    }
    
    ## Normalize the nucleosome splicing score
    settings$nucNormScore_constrained = sub("_mean.txt$","_norm.txt",settings$nucTestFileConstrained_mean)
   
    if (!file.exists(settings$nucNormScore_constrained)){
      cat("Normalizing Seq To Nuc output\n")
      normalizeToText(settings$nucTestFileConstrained_mean, settings$nucNormScore_constrained)
    }

    return(settings)
}



################################################################################################################################################################################
#### Function to compute the mean of the NUC signal.
################################################################################################################################################################################
computeNucMean = function(settings, infile=settings$nucTestFileFixed_out, outfile=settings$nucTestFileFixed_mean ){
    IN = file(infile,open='r')
    OUT = file(outfile,open='w')

    i = 0;
    while(length(HEADER <-readLines(IN,1)) > 0){
         i = i + 1
         cat(sprintf("   Line %i\r",i))
    
         dat = scan(IN,double(),sep=",",nlines=1,quiet=T)       
    
         writeLines(sprintf("%s\t%f",HEADER,mean(dat)), OUT)
    
    } 
    cat("        \n")
    close(IN)
    close(OUT)

}

################################################################################################################################################################################
#### Function to retrieve a FIXED LENGTH of a sequence from a list of features #####
################################################################################################################################################################################

retrieveFixedSequences = function (settings,featureList,genome,width=180,overhangUp=37,overhangDown=37,fromStart=TRUE,firstSeq=1, outputFile=NULL,writeDuring=FALSE,removeBadChars=FALSE){
	# If 'fromStart' is FALSE, then get the sequenced started from the END of the feature 

	
  	## Initialize output list
  	# Length of sequence:
  	seqLength = width -1 + overhangUp + overhangDown
  	
  	# For each feature, add a blank sequence
  	output_seqs = rep(paste(rep('0',seqLength),collapse=''),nrow(featureList))        
    # Make the .fa headers
    headers = paste(">", featureList$UniqueID, '_', featureList$FeatureCount, '_', featureList$isLast, sep='')    
    #browser()
  	  	
  	## Do I start from start or end of each feature?
  	myCol = ifelse(fromStart,'start','end')
  	cat(sprintf ("Collecting features based on column '%s'\n",myCol))
  	#print(names(featureList))
  	
  	## For each feature, return the sequence. NOTE: need to correct for UCSC's 0-based indexing
  	lastSaved = firstSeq-1 
  	for( i in firstSeq:nrow(featureList)){
		
		
		#if(i%%1000==1)
  	 cat(sprintf("Doing seq %d\r",i))
  		
  		
  		if(featureList$strand[i] == '+'){
    		# Positive strand: 'simple'
    		left = featureList[i,myCol] - overhangUp + 1 # Correct for 0-indexing
    		right  = featureList[i,myCol] + overhangDown + width
    		
    		seq = as.character(subseq(genome[[featureList$chr[i]]],start=left,end=right))
	
    	}else{                                                
    		
        # Negative strand: NOTE that left and right are reversed! Because the 'start' in the featureList is the ACTUAL start of the feature
        #		(i.e. to the right on a genome browser track, in Negative strand), but subseq() function requires end >= start
        
        right = featureList[i,myCol] + overhangUp # 'upstream' is in the context of transcription
        left  = featureList[i,myCol] - overhangDown - width + 1 # Correct for 0-indexing
        
        seq = as.character(reverseComplement(subseq(genome[[featureList$chr[i]]],start=left,end=right)))
      }
			
      # save seq to array
    	output_seqs[i] <- seq  
    	 	
    	# If the writeDuring setting is allowed, save last 1,000 rows to file 	
    	if(i%%100==0 & writeDuring){                                 
    			# write some output
  				savePortion(output_seqs[(lastSaved+1):i],headers[(lastSaved+1):i], outputFile)
          lastSaved = i
    		} 	
    	 	
    }
     cat(sprintf("Done                \n"))
   
    ## After looping through all sequences, write the remaining outputs to file
   	if(i != lastSaved){                                 
       # write some output
	     savePortion(output_seqs[(lastSaved+1):i],headers[(lastSaved+1):i], outputFile)      
 		} 
}

################################################################################################################################################################################
#### Function to retrieve a CONSTRAINED FIXED LENGTH of a sequence from a list of features 
## 4/2/12 This revision allows one to constrain the sequence to N characters.
################################################################################################################################################################################

retrieveFixedSequences_min = function (settings,featureList,genome,width=147,overhangUp=37,overhangDown=37,fromStart=TRUE,firstSeq=1, outputFile=NULL,writeDuring=FALSE,removeBadChars=FALSE){
	# If 'fromStart' is FALSE, then get the sequenced started from the END of the feature 

    print(width)

    ### First get the Length of the regions
    featureList$Length = abs(featureList$start-featureList$end)

  	## Initialize output list
  	# Length of sequence:
  	featureList$seqLength = sapply(featureList$Length, function(x) min(c(x,width))) -1 + overhangUp + overhangDown
    # Remove super short ones
    featureList = subset(featureList,seqLength > 10)
	
  	# For each feature, add a blank sequence
  	output_seqs = sapply(featureList$seqLength,function(x)paste(rep('0',x),collapse=''))
    # Make the .fa headers
    headers = paste(">", featureList$UniqueID, '_', featureList$FeatureCount, '_', featureList$isLast, sep='')    
    #browser()
  	  	
  	## Do I start from start or end of each feature?
  	myCol = ifelse(fromStart,'start','end')
  	cat(sprintf ("Collecting features based on column '%s'\n",myCol))
  	#print(names(featureList))
  	
  	## For each feature, return the sequence. NOTE: need to correct for UCSC's 0-based indexing
  	lastSaved = firstSeq-1 
  	for( i in firstSeq:nrow(featureList)){
		
		
		#if(i%%1000==1)
  	 cat(sprintf("Doing seq %d\r",i))
  		
  		
  		if(featureList$strand[i] == '+'){
    		# Positive strand: 'simple'
    		left = featureList[i,myCol] - overhangUp + 1 # Correct for 0-indexing
    		right  = featureList[i,myCol] + overhangDown + featureList$seqLength[i]
    		
    		seq = as.character(subseq(genome[[featureList$chr[i]]],start=left,end=right))
	
    	}else{                                                
    		
        # Negative strand: NOTE that left and right are reversed! Because the 'start' in the featureList is the ACTUAL start of the feature
        #		(i.e. to the right on a genome browser track, in Negative strand), but subseq() function requires end >= start
        
        right = featureList[i,myCol] + overhangUp # 'upstream' is in the context of transcription
        left  = featureList[i,myCol] - overhangDown - featureList$seqLength[i] + 1 # Correct for 0-indexing
        
        seq = try(as.character(reverseComplement(subseq(genome[[featureList$chr[i]]],start=left,end=right))))
        if (inherits(seq,'try-error')){
            print(featureList[i,])
            stop(sprintf("\nSequence %s had error.  Last Saved: %s",i,lastSaved))

        }
      }
			
      # save seq to array
    	output_seqs[i] <- seq  
    	 	
    	# If the writeDuring setting is allowed, save last 1,000 rows to file 	
    	if(i%%100==0 & writeDuring){                                 
    			# write some output
  				savePortion(output_seqs[(lastSaved+1):i],headers[(lastSaved+1):i], outputFile)
          lastSaved = i
    		} 	
    	 	
    }
     cat(sprintf("Done                \n"))
   
    ## After looping through all sequences, write the remaining outputs to file
   	if(i != lastSaved){                                 
       # write some output
	     savePortion(output_seqs[(lastSaved+1):i],headers[(lastSaved+1):i], outputFile)      
 		} 
}

################################################################################################################################################################################
####  ## Subroutine to save a portion of the sequences during retrieval
################################################################################################################################################################################

savePortion = function(seqs,headers, file){
    # get rid of bad characters (no-ACGT)
    badSeqs <- which(regexpr("([^ACGTacgt])",seqs) > 0)
    if(any(badSeqs)){
       headers = headers[-badSeqs]
       seqs = seqs[-badSeqs]      
    }
    if(length(seqs) < 1){
      cat ("All seqs were bad!  Not writing anything\n")
      return (0)
     }      
    # 
    output = matrix(rbind(headers, seqs),nc=1) # reshape the matrix so that it looks like a fasta file
    write.table(output,  file=file,  sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE,append=T)    	
  }







################################################################################################################################################################################
#### Functions to remove duplicate features
# These functions all care about the 'isLast' because an ambiguous last-element gets to go in BOTH 'lasts' and 'non-lasts'
################################################################################################################################################################################

removeDuplicateStarts = function(features){
	features[!duplicated(features[,c('isLast','chr','start')]),]
}
	

removeDuplicateEnds = function(features){
	features[!duplicated(features[,c('isLast','chr','end')]),]
}

removeDuplicateStartANDEnds = function(features){
	features[!duplicated(features[,c('isLast','chr','start','end')]),]
}








                                

############################################################################################################
### Wrapper function to compute splice site strengths
############################################################################################################

computeSpliceSiteStrength = function(settings,features,genome,donorPWM=NULL,acceptorPWM=NULL){

  ####### DONORS #######
  DonorsFile =    "intron_donors_for_Burset.txt"                                
  DonorsScores = sprintf("%s.scores",DonorsFile) 
  DonorsNorm = sprintf("%s.NormScores",DonorsFile)             
  
  ### Create donors if necessary
  if (!file.exists(DonorsFile)){
    firstSeq = 1      
    cat(sprintf("Retrieving sequences of exons for donors (starting with seq %d):\n", firstSeq))

    ## Get rid of duplicate intron starts
    unique_starts = removeDuplicateStarts(features$introns)

    ## Get all donors (writes to file during the loop)
    seqs = retrieveFixedSequences(settings, unique_starts, genome, width=0, overhangUp=3, overhangDown=6, fromStart=TRUE,outputFile=DonorsFile,removeBadChars=T,firstSeq=firstSeq,writeDuring=T)
  }
  
  ### Now run pwm on donors if necessary
  if (!file.exists(DonorsScores)){          
    cat(sprintf("Calculating donor splice scores\n"))
    
    calculatePWM(DonorsFile , outputfile=DonorsScores ,subsetLocations=4:5, subsetString="GT", pwmFile=donorPWM)
  }
  
  ### Now run normalization script if necessary
  if (!file.exists(DonorsNorm)){
	cat(sprintf("Normalizing donor splice scores\n"))
    normalizeToText(DonorsScores, DonorsNorm,LOG=T)
  }  

  ####### Acceptors #######
  AcceptorsFile = "intron_acceptors_for_Burset.txt"                                                                                                                                    
  AcceptorsScores = sprintf("%s.scores",AcceptorsFile) 
  AcceptorsNorm = sprintf("%s.NormScores",AcceptorsFile)             
  
  ### Create Acceptors if necessary  
  if (!file.exists(AcceptorsFile)){      
    firstSeq = 1
    cat(sprintf("Retrieving sequences of exons for acceptors (starting with seq %d):\n", firstSeq))
    
    ## Get rid of duplicate intron ENDS
    unique_ends = removeDuplicateEnds(features$introns)
       
    seqs = retrieveFixedSequences(settings, unique_ends, genome, width=0, overhangUp=14, overhangDown=1, fromStart=FALSE,outputFile=AcceptorsFile,removeBadChars=T,firstSeq=firstSeq,writeDuring=T)        
  }      

  ## Now run pwm on acceptors if necessary
  if (!file.exists(AcceptorsScores)){          
    cat(sprintf("Calculating acceptor splice scores\n"))
    
    calculatePWM(AcceptorsFile , outputfile=AcceptorsScores,subsetLocations=13:14, subsetString="AG", pwmFile=acceptorPWM)
  }
  
  ### Now run normalization script if necessary
  if (!file.exists(AcceptorsNorm)){
	cat(sprintf("Normalizing acceptor splice scores\n"))
    normalizeToText(AcceptorsScores, AcceptorsNorm, LOG=T)
  } 
  
  

  settings$AcceptorScores = AcceptorsScores
  settings$DonorScores = DonorsScores
  settings$AcceptorsNorm = AcceptorsNorm
  settings$DonorsNorm = DonorsNorm
  return(settings)
}
                 

##########################################################################################################
### Function to Read in a PWM and format it for the Biostrings functions
###########################################################################################################
    
readPWM = function(pwmFile){
  PWM.in =read.delim(pwmFile,head=T)
  pwm =  apply(PWM.in,2,function(x)x/sum(x))
  rownames(pwm) = c("A","C","G","T")
  pwm
}

    
############################################################################################################
### Wrapper function to compute splice site strengths
### Uses Gene Yeo's program
############################################################################################################

computeSpliceSiteStrength2 = function(settings,features,genome,donorPWM=NULL,acceptorPWM=NULL){

  ####### DONORS #######
  DonorsFile =    "intron_donors_for_Burset.txt"                                
  DonorsScores = sprintf("%s.scores",DonorsFile) 
  
  ### Create donors if necessary
  if (!file.exists(DonorsFile)){
    firstSeq = 1      
    cat(sprintf("Retrieving sequences of exons for donors (starting with seq %d):\n", firstSeq))

    ## Get rid of duplicate intron starts
    unique_starts = removeDuplicateStarts(features$introns)

    ## Get all donors (writes to file during the loop)
    seqs = retrieveFixedSequences(settings, unique_starts, genome, width=0, overhangUp=3, overhangDown=6, fromStart=TRUE,outputFile=DonorsFile,removeBadChars=T,firstSeq=firstSeq,writeDuring=T)
  }
  
  ### Now run pwm on donors if necessary
  if (!file.exists(DonorsScores)){          
    cat(sprintf("Calculating donor splice scores\n"))
    
  #  calculatePWM(DonorsFile , outputfile=DonorsScores ,subsetLocations=4:5, subsetString="GT", pwmFile=donorPWM)
  }

  ####### Acceptors #######
  AcceptorsFile = "intron_acceptors_for_Yeo.txt"                                                                                                                                    
  AcceptorsScores = sprintf("%s.scores",AcceptorsFile) 
              
  
  ### Create Acceptors if necessary  
  if (!file.exists(AcceptorsFile)){      
    firstSeq = 1
    cat(sprintf("Retrieving sequences of exons for acceptors (starting with seq %d):\n", firstSeq))
    
    ## Get rid of duplicate intron ENDS
    unique_ends = removeDuplicateEnds(features$introns)
       
    seqs = retrieveFixedSequences(settings, unique_ends, genome, width=0, overhangUp=20, overhangDown=3, fromStart=FALSE,outputFile=AcceptorsFile,removeBadChars=T,firstSeq=firstSeq,writeDuring=T)        
  }      

  ## Now run pwm on acceptors if necessary
  if (!file.exists(AcceptorsScores)){          
    cat(sprintf("Calculating acceptor splice scores\n"))
    
    #calculatePWM(AcceptorsFile , outputfile=AcceptorsScores,subsetLocations=13:14, subsetString="AG", pwmFile=acceptorPWM)
  }

  settings$AcceptorScores2 = AcceptorsScores
  settings$DonorScores2 = DonorsScores
  return(settings)
}
                 

##########################################################################################################
### Function to Read in a PWM and format it for the Biostrings functions
###########################################################################################################
    

    
##########################################################################################################
### Function to run a PWM on a fasta file (output of intron retrieval), and write results to a file
###########################################################################################################

calculatePWM = function(inputFile , outputfile = sprintf("%s.scores",inputFile),subsetLocations=NULL, subsetString=NULL, pwmFile=NULL){
# If pwmFile is NULL, generate the pwm from the concensus matrix of the sequences themselves.
# If subsetLocations or subsetString are NULL, don't filter at all: Otherwise filter some bases (e.g. 'AG' or 'GT' intrrons only)

  ## Read in 
  seqs =read.DNAStringSet(inputFile)

  ## Subset only for the 'AG' or 'GT' bases  
  if (!is.null(subsetLocations) & !is.null(subsetString)){
    cat(sprintf("Filtering for bases %s in position %d to %d\n",subsetString,subsetLocations[1],subsetLocations[2]))
    
    myString = DNAString(subsetString)
    which.seqs.good = sapply(seqs,function(x)x[subsetLocations]==myString)
    seqs = seqs[which.seqs.good] 
  }
  
      
  ## Optional way of computing pwm from here:
  if (is.null(pwmFile)){    
    cat(sprintf("Generating Concensus Matrix from Sequence\n"))
   
    pwm = consensusMatrix (seqs)[1:4,]/length(seqs)
    
  }else{
    cat(sprintf("Reading in PWM file %s\n",pwmFile))

    pwm = readPWM(pwmFile)
  }
  
  ## Compute scores
  cat(sprintf("Computing pwm scores\n"))  
  sapply(seqs,function(x)exp(PWMscoreStartingAt(log(pwm),x)))->scores

  ## Write scores to file
  write.table(scores, file=outputfile, quote=FALSE, row.names=T, col.names=FALSE,sep= "\t")

}

    

##########################################################################################################
### Functions to save Data 
###########################################################################################################

saveData = function(settings,file=NULL){
	UseMethod("saveData")
}
 
###### SAVE THE SETTINGS  #########
saveData.ExonPipeline = function(settings){    
	N = names(settings)
	write.table(t(sapply(N,function(x)c(settings[x],paste("#",x)))),quote=FALSE,sep="\t",file='settings.txt',row.names=FALSE,col.names=FALSE)
	
}
  






  
############################################################################################
######  Functions for analyzing Chromatin-associated RNA-seq
############################################################################################


############################################################################################
######  Functions to create a .bed file for the 5' end of all last exons
############################################################################################

makeLastExonBed = function(settings,features){
  ### Use 'start' of exon: write a .bed file that can be used by bedtools to intersect with .bam from alignment
    # Negative strand genes: 'start' bp is actually the 1st bp of the last exon
    # Pos: 'start' bp is actually ONE BASE UPSTREAM OF THE LAST EXON
    
    
  lasts = features$exons[features$exons$isLast == 1,]
  lasts = removeDuplicateStartANDEnds(lasts) 

  # Fix positives to 1-based
  lasts$start[lasts$strand == '+'] = lasts$start[lasts$strand == '+'] + 1 

  ## create .bed for all genes, length 1!
  outputBed = cbind(lasts$chr,lasts$start,lasts$start,lasts$UniqueID,0,lasts$strand)

  ## writeFile
  settings$lastExonBed = paste(settings$thisDir,'last_exons5.bed',sep='/')
  write.table(outputBed, settings$lastExonBed, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
  
  return(settings)
}

make2ndExonBed = function(settings,features){
  ### Use 'start' of exon: write a .bed file that can be used by bedtools to intersect with .bam from alignment
    # Negative strand genes: 'start' bp is actually the 1st bp of the last exon
    # Pos: 'start' bp is actually ONE BASE UPSTREAM OF THE LAST EXON
    
    
  seconds = features$exons[features$exons$isLast == 0 & features$exons$FeatureCount==2,]
  seconds = removeDuplicateStartANDEnds(seconds) 

  # Fix positives to 1-based
  seconds$start[seconds$strand == '+'] = seconds$start[seconds$strand == '+'] + 1 

  ## create .bed for all genes, length 1!
  outputBed = cbind(seconds$chr,seconds$start,seconds$start,seconds$UniqueID,0,seconds$strand)

  ## writeFile
  settings$secondExonBed = paste(settings$thisDir,'second_exons5.bed',sep='/')
  write.table(outputBed, settings$secondExonBed, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
  
  return(settings)
}

### 'All' is actually misleading because it DOES NOT INCLUDE 'first' EXONS
makeAllExonBed = function(settings,features){
  ### Use 'start' of exon: write a .bed file that can be used by bedtools to intersect with .bam from alignment
    # Negative strand genes: 'start' bp is actually the 1st bp of the last exon
    # Pos: 'start' bp is actually ONE BASE UPSTREAM OF THE LAST EXON
    
    
  alls = features$exons[features$exons$FeatureCount > 1,]
  alls = removeDuplicateStartANDEnds(alls) 

  # Fix positives to 1-based
  alls$start[alls$strand == '+'] = alls$start[alls$strand == '+'] + 1 

  ## create .bed for all genes, length 1!
  outputBed = cbind(alls$chr,alls$start,alls$start,paste(alls$UniqueID,alls$FeatureCount-1,sep='_'),0,alls$strand) # adding intron number

  ## writeFile
  settings$allExonBed = paste(settings$thisDir,'all_exons5.bed',sep='/')
  write.table(outputBed, settings$allExonBed, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
  
  return(settings)
}

############################################################################################
######  Function to map chromatin-associated .bam reads against a lastExonBed file (or other similar file)
############################################################################################
mapChromAssocBam = function(settings,chromBamFile,targetFile,lastsFile=settings$lastExonBed){
  #settings$mappedChromReads = paste(sub(".bed$","",settings$lastExonBed),"mapped.bed",sep="_")
  
  if(is.null(chromBamFile))
    stop("No bam-file specified")
  cmd = sprintf("nohup bedtools intersect -abam %s -b %s -wa -wb -bed -split | cut -f 2,3,11-18 > %s & >nohup.log",chromBamFile, lastsFile,targetFile)
  print(cmd)
  system(cmd)

  #return(settings)
}



############################################################################################################################
######  Wrapper for reading in input from running intersectBed command (mapChromAssocBam).
######  This simply returns a list of the input files, as well as concatenates them all into one larger file
############################################################################################################################

get_chromatin_splicing = function(samples,addText){ #"_splicing.txt"
  outData = list()
  
  # Also 'cat' output into one large file:
  combinedOutput = sprintf("Chromatin%s_All.txt",addText)
  combineFiles = !file.exists(combinedOutput)
  
  for(i in 1:length(samples)){
    input = paste(samples[i],addText,".txt",sep="")
    target = paste(samples[i],addText,"_splicing.txt", sep="")
    
    print(target) 
    
    ## Obtain the data for each sample. 
    
    if(!file.exists(target)) {
      # If the data doesn't exist yet, generate it
      outData[[i]] = parseChromOutput(settings,input)  
    }else{
      # Otherwise, read it in.
      outData[[i]] = read.delim(target, stringsAsFactors=F) 
    }
    
    ## If combinedFiles is set to TRUE, append to the end of the file
    if(combineFiles){
      system(sprintf("cat %s >> %s",input ,combinedOutput)	)
    }         
  }
  
  ## Get the aggregate data.
  if(combineFiles){ # We need to generate it first!
    print(combinedOutput);
    combinedChrom = parseChromOutput(settings,combinedOutput)
  }else{
    # Otherwise, read it it
    combinedChrom = read.delim(sub(".txt$","_splicing.txt",combinedOutput),stringsAsFactors=F)   
  }
  outData[[i+1]] = combinedChrom    
  
  outData
}

####################################################################################################################################
### Function to filter and combine the result of 'get_chromatin_splicing'
####################################################################################################################################

mergeChromSplicing = function(splicingList, total_reads_cutoff=10,samples){
  
  ### First filter by introns where the total number of reads is > CUTOFF                            
  chromData = lapply(splicingList[1:(length(splicingList)-1)], function(x) x[x$SplicedReads + x$UnsplicedReads > total_reads_cutoff,])
  
  ### This is a really ugly loop to merge all the data from the 5 samples of last reads
  chromMatrix = chromData[[1]][,c(1,4)]; names(chromMatrix)[2] = samples[1]
  for(x in 2:length(chromData)){
    #a=lapply(2:length(chromData),function(x){
    X <- chromData[[x]]
    names(X)[4] <- samples[x]
    chromMatrix <- merge(chromMatrix,X[,c(1,4)], by=1, all=T)
    NULL
  }
  
  ### Merge the 'combined' data (those where I combined over all samples)
  combinedChrom = chromData[[length(chromData)]]
  chromMatrix <- merge(chromMatrix,combinedChrom[combinedChrom$SplicedReads + combinedChrom$UnsplicedReads > total_reads_cutoff,c(1,4)], by=1, all=T)
  names(chromMatrix)[ncol(chromMatrix)] <- "PooledSplicing"
  
  ## Add names of genes and feature count just in case
  if(length(grep("_",chromMatrix$UniqueID[1])) ==1){
    chromMatrix$FeatureCount = sub("^([0-9]*)_([0-9]*)$","\\2",chromMatrix$UniqueID)
    chromMatrix$UniqueID   = sub("^([0-9]*)_([0-9]*)$","\\1",chromMatrix$UniqueID)
  }
  
  chromMatrix
}


############################################################################################################################
######  Helper functions for reading in and parsing output bed from intersectBed command
############################################################################################################################

###  This function returns gets rid of alignments that don't actually overlap with the last exon, but instead are perfectly adjacent to it
remove_adjacent_alignments = function(aligments){
  BAD = (aligments$strand == '+' & (aligments$readStart+1 == aligments$ExonStart | aligments$readEnd + 1 == aligments$ExonStart)) |
                (aligments$strand == '-' & (aligments$readStart == aligments$ExonStart | aligments$readEnd == aligments$ExonStart ))
  aligments[!BAD,]                
}  

####
getPosSplicers = function(row){
  GAPS = row[4]

  ## For each read (row), loop through the number of segments in the alignment: check for the gapped alignment to coincide with this exon
  while(nchar(GAPS<-sub("^[0-9]*,","",GAPS))){
    gap = as.numeric(sub("^([0-9]*),(.*|$)","\\1",GAPS))

    ## If it does, return success
    if (as.numeric(row[1]) + 1 + gap == as.numeric(row[5])) # start + 1 + gap == exonStart
      return( paste(row[6],1))  # Return "UniqueID I(isSpliced)"
  }
  
  ## If we get here, that means there was no correct splicing that we were querying.  Thus, it should be considered either unspliced or irrelevant
  return (paste(row[6],0))   
}


getNegSplicers = function(row){
  GAPS = row[4]
  LENS = row[3]
  
  ## For negative strand genes, we need to also keep track of the length of each aligned section
  # ALSO, we need to start with the first one!
  
  repeat{
    
    gap = as.numeric(sub("^([0-9]*),(.*|$)","\\1",GAPS))
    len = as.numeric(sub("^([0-9]*),(.*|$)","\\1",LENS))
    #print(gap)      
    if (as.numeric(row[1]) + len + gap == as.numeric(row[5])) # start + 1 + gap + len[current segment] == exonStart
      return( paste(row[6],1))                # Return "UniqueID I(isSpliced)"
      
    GAPS<-sub("^[0-9]*,","",GAPS)  
    LENS<-sub("^[0-9]*,","",LENS)
              
    if(!nchar(GAPS))
      return (paste(row[6],0))
  
  
  }
  # If we get here, that means there was no correct splicing that we were querying.  Thus, it should be considered either unspliced or irrelevant
  return (paste(row[6],0))   
}

############################################################################################
######  Read in and parse output bed from intersectBed command
###         NOTE: bed file is 0-indexed
############################################################################################
parseChromOutput = function(settings,chromFile){
  bed = read.delim(chromFile,head=FALSE,stringsAsFactors=FALSE)
  names(bed) =c("readStart","readEnd","readSections","readGapsFromStart","chr","ExonStart","ExonEnd","UniqueID","Score","strand")
 
  ## Set up initial file
  unsplices = NULL
 
  ### First separate by 'Genome' reads that align with no gaps
  which.Genome = bed$readGapsFromStart == "0,"                                                   
  bedGenome = bed[which.Genome,]
  
  ## Now just make sure that none of the reads start or stop exactly on the splice site.
     
  if (nrow(bedGenome) > 0){          
    ## That that are simply adjacent to the exon boundary get tossed: KEEP THE REST AS UNSPLICED
    unsplices = c(unsplices,remove_adjacent_alignments(bedGenome)$UniqueID)
    cat(sprintf("Tossing %d ambiguous reads.\n", nrow(bedGenome) - length(unsplices) ))
  }                                                             
  
  
  ## Now for those that are NOT Genome-aligned:
  bedSplit = bed[!which.Genome,c(1:4,6,8,10)]
  
  
  ## Positive strand:
  POS = bedSplit[bedSplit$strand=='+',]
  exact_splicing_pos = apply(POS,1,getPosSplicers)
  good.splicing.pos = regexpr("1$",exact_splicing_pos) > 0
    
  ## NEGATIVE strand:
  NEG = bedSplit[bedSplit$strand=='-',] 
  exact_splicing_neg = apply(NEG,1,getNegSplicers)
  good.splicing.neg = regexpr("1$",exact_splicing_neg) > 0
    
  ## Combine and do one final round of cleaning.  
  splices = c(POS$UniqueID[good.splicing.pos], NEG$UniqueID[good.splicing.neg])

  unsplices = c(unsplices, remove_adjacent_alignments(POS[!good.splicing.pos,])$UniqueID, remove_adjacent_alignments(NEG[!good.splicing.neg, ])$UniqueID)
  
  
  ## Combine the splices and unsplices for each gene
  OUT = merge(table(splices),table(unsplices),by=1,all=T) # table() counts each ID's occurence, then merge on ID
  OUT[is.na(OUT)] <-0
  names(OUT) = c("UniqueID","SplicedReads","UnsplicedReads")
  
  OUT$PercentSpliced = OUT$SplicedReads / (OUT$SplicedReads + OUT$UnsplicedReads)
#browser()
  ## Write file  
  write.table(OUT, file=sub("\\.[[:alnum:]]*$", "_splicing.txt", chromFile), row.names=FALSE, sep="\t", quote=FALSE)
  
  OUT
}


########################################################################################################################################
### This function calculates the distance from the 3'end of each intron to the 3' end of the last intron
### It takes all introns and their acceptor scores and puts them into one file, formatted to easily read into the MATLAB model
########################################################################################################################################

formatIntronsForExpectedTimeModel = function(intronsFile='introns.txt', scoreFile='intron_acceptors_for_Burset.txt.scores', outputfile){
    
    ## Read in introns annotation and acceptor score file
    introns = read.delim(intronsFile, stringsAsFactors=F)
    scores = read.delim(scoreFile, stringsAsFactors=F, head=F,col.names=c("intronCode","IntronScore"))
    
    # Get the last introns
    LASTS = introns[introns$isLast == 1, c("UniqueID","end")]
    dimnames(LASTS)[[2]][2] = 'LastEnd'
    
    ## Merge with the LASTS (and the LASTS ends)
    introns.withLastEnd = merge(introns,LASTS, by = "UniqueID")

    # Calculate the distance(in bp) of the 3' end of each intron to the 3' end of the LAST intron
    introns.withLastEnd$distToEnd = abs(introns.withLastEnd$LastEnd - introns.withLastEnd$end)

    ## Merge with scores
    introns.withLastEnd$intronCode = paste(introns.withLastEnd$UniqueID,  introns.withLastEnd$FeatureCount, introns.withLastEnd$isLast , sep='_')
    introns.withScores = merge(introns.withLastEnd,scores, by.x = "intronCode", by.y = 1)

    ## Get rid of genes that don't have ALL their introns listed
    IDs1 = table(introns.withLastEnd$UniqueID)
    IDs2 = table(introns.withScores$UniqueID)
    mergedList = merge(IDs1,IDs2, by = 1) # Merge the counts of each gene, both with and without the scores.
    goodIDs = mergedList[mergedList[,2] == mergedList[,3],1]  # the 'good' ones are those that have the same count in both objects
    length(goodIDs) #[1] 19874 in Mouse
    
    ## Use only those that are goodNames
    intronsToWrite = introns.withScores[introns.withScores$UniqueID %in% goodIDs,]
    
    ## WRITE ONE LINE FOR EACH GENE:
    intronsToWrite = intronsToWrite[order(intronsToWrite$UniqueID, intronsToWrite$FeatureCount),]
    dataByGene = aggregate(intronsToWrite[,c("distToEnd","IntronScore")],by=list(UniqueID=intronsToWrite$'UniqueID'),paste,collapse=",") 
    write.table(dataByGene, file = outputfile, sep="\t",quote=F,row.names=F,col.names=F)
    
}

########################################################################################################################################
### This function runs a matlab model to calculate the expected time to splice an entire gene
########################################################################################################################################
calculateExpectedTime_MATLAB = function(inFile,elongation_rate=3000,splicingFactor=5){

	# Write a tiny matlab script
	write.table(sprintf("run /home/jeremy/startup.m\n[IDs ReadTimes] = run_genomic_expected_time('%s/%s',%f,%f);",getwd(),inFile,elongation_rate,splicingFactor),file='tmp.matlab.m',quote=F,row.names=F, col.names=F)
	# Call that script!
	system("matlab -nojvm < tmp.matlab.m")			

}

########################################################################################################################################
### This function runs a matlab model to calculate the % splicing of all genes with Runon data
########################################################################################################################################
calculateSplicingRunon_MATLAB = function(inFile,elongation_rate=3000,k_unpause=30,conversionString ='0.5,5,30,150',VERSION=1){

	# Write a tiny matlab script
	if(VERSION==1)
    write.table(sprintf("run /home/jeremy/startup.m\n[IDs ReadTimes] = run_genomic_withRunon('%s/%s',%f,%f,'%s');",getwd(),inFile,elongation_rate,k_unpause,conversionString),file='tmp.matlab.m',quote=F,row.names=F, col.names=F)
  if(VERSION==2)
    write.table(sprintf("run /home/jeremy/startup.m\n[IDs ReadTimes] = run_genomic_withRunon2('%s/%s',%f,%f,'%s');",getwd(),inFile,elongation_rate,k_unpause,conversionString),file='tmp.matlab.m',quote=F,row.names=F, col.names=F)    

	# Call that script
	system("matlab -nojvm < tmp.matlab.m")
  # read the file back in
  cat("Finished MATLAB process\n")
  if(VERSION==1)
    cat(sprintf("Writing to %s\n", OUT<-sprintf('%s_splicingRatesRunon_%.2f_%.2f_%s.txt',sub('.txt$','',inFile), elongation_rate, k_unpause, conversionString)))
  if(VERSION==2)
    cat(sprintf("Writing to %s\n", OUT<-sprintf('%s_splicingRatesRunon2_%.2f_%.2f_%s.txt',sub('.txt$','',inFile), elongation_rate, k_unpause, conversionString)))
  read.table(OUT,col.names=c("UniqueID","PredSplice","PredSplice_strength","PredSplice_pause","PredSplice_runon","PredSplice_ALL")) 
}

################################################################################################################################################
## This function normalizes a vector to the unit normal
################################################################################################################################################
normalize = function(vec){
	(vec - mean(vec, na.rm=T)) / sd(vec,na.rm=T)
}

################################################################################################################################################
## Read in, normalize and write: for a file where the data is in the 2nd column (default)
################################################################################################################################################
normalizeToText = function(fileIn, fileOut, myCol=2, header=FALSE, LOG=F, ...){
	dat = read.delim(fileIn, header=header, ...)
	
	if (LOG){
		dat[,myCol] = normalize(log10(dat[,myCol]))
	}else
		dat[,myCol] = normalize(dat[,myCol])
	write.table(dat, file=fileOut, row.names=F, col.names=F, quote=F, sep='\t')
}


###################################################################################################################################################
## Get first, last, middle exons for a genome
###################################################################################################################################################
getFirstMiddleLast = function(exons){
  ## First go back and figure out the # exons per gene
  # For each isoform, get the isLast numExons
  myLastExons = exons[exons$isLast==1,]
  myFirstExons = exons[exons$FeatureCount==1,]
  myMiddleExonLengths = data.frame(UniqueID=myLastExons$UniqueID, FeatureCount=round(myLastExons$FeatureCount/2))

  myMiddleExons = merge(exons, myMiddleExonLengths, by=c('UniqueID','FeatureCount')) 
  #myMiddleExons = subset(myMiddleExonLengths, FeatureCount != MiddleCount & MiddleCount > 1)

  list(Firsts = myFirstExons, Middles = myMiddleExons, Lasts = myLastExons)

}

###################################################################################################################################################
## Take some GRO-seq data from Karmel Allison's database and remove genes that are too close to other genes
###################################################################################################################################################
makeGROdata = function(inputFile="/home/home/jeremy/RNAseq/Glass/refseq_and_subsequent_transcripts_with_tag_counts.txt", refGene,outputFile = "GRO_Data.RData", geneDistance=1000){ 
  
  ## Read in GRO data
  GRO = read.delim(inputFile,stringsAsFactors=F)
  
  ## Calculate the Read-Through (hereafter referred to as 'Runon') and get rid of genes with no measurable read-through.
  GRO$Runon = abs(GRO$post_transcript_end - GRO$post_transcript_start)
  GRO = GRO[!is.na(GRO$Runon),]
  
  ## Annotate by merging with refGene
  GRO2 = merge(GRO,refGene[,c("name","UniqueID","strand")],by=1)
  
  ## Create .bed file, overlap with refseq.bed
  refGene$TSS = refGene$txStart
  refGene$TSS[refGene$strand=='-'] <- refGene$txEnd[refGene$strand=='-']
  
  options(scipen=10) # make it so I don't write scientific numbers here
  
  ## Omit any GRO-seq data when the Read-Through ends within 1Kb of a gene.
  write.table(cbind(refGene$chr, refGene$TSS-geneDistance, refGene$TSS +geneDistance , refGene$name, 0, refGene$strand),file='temp_refseq.bed',sep="\t",row.names=F,col.names=F, quote=F)
  write.table(cbind(GRO2$chr, GRO2$post_transcript_start-1, GRO2$post_transcript_end+1, GRO2$UniqueID, 0, GRO2$strand.y),file='temp_GRO.bed',sep="\t",row.names=F,col.names=F, quote=F)
  system("bedtools intersect -a temp_GRO.bed -b temp_refseq.bed -v -s -wa > temp_GRO_no_TSS_interuptions.bed")
  options(scipen=0)
  
  ## Read data back in
  good_GRO =read.delim("temp_GRO_no_TSS_interuptions.bed", head=F)
  # Filter for good data
  GRO3 <- GRO2[GRO2$UniqueID %in% good_GRO[,4],]
  
  ## Save
  save(GRO3,file=outputFile)
  
  ## Erase temp files
  unlink("temp_GRO_no_TSS_interuptions.bed")
  unlink("temp_refseq.bed")
  unlink("temp_GRO.bed")

  return(GRO3)
}

###################################################################################################################################################
## Take some GRO-seq data from Karmel Allison's database and remove genes that are too close to other genes
###################################################################################################################################################
makeGROdata2 = function(inputFile="/home/RNAseq/Glass/post_gene_transcripts.txt", refGene,outputFile = "GRO_Data_redone.RData", geneDistance=1000){ 
  
  ## Read in GRO data
  GRO = read.delim(inputFile,stringsAsFactors=F)
  
  ## Calculate the Read-Through (hereafter referred to as 'Runon') and get rid of genes with no measurable read-through.
  #GRO$Runon = abs(GRO$post_gene_end - GRO$gene_ends)
  #GRO = GRO[!is.na(GRO$Runon),]
  
  ## Annotate by merging with refGene
  GRO2 = merge(GRO,refGene[,c("name","UniqueID","strand")],by=1)
  
  ## Create .bed file, overlap with refseq.bed
  refGene$TSS = refGene$txStart
  refGene$TSS[refGene$strand=='-'] <- refGene$txEnd[refGene$strand=='-']
  
  options(scipen=10) # make it so I don't write scientific numbers here
  
  ## Omit any GRO-seq data when the Read-Through ends within 1Kb of a gene.
  write.table(cbind(refGene$chr, refGene$TSS-geneDistance, refGene$TSS +geneDistance , refGene$name, 0, refGene$strand),file='temp_refseq.bed',sep="\t",row.names=F,col.names=F, quote=F)
  write.table(cbind(GRO2$chr, GRO2$post_gene_start-1, GRO2$post_gene_end+1, GRO2$UniqueID, 0, GRO2$strand),file='temp_GRO.bed',sep="\t",row.names=F,col.names=F, quote=F)
  system("bedtools intersect -a temp_GRO.bed -b temp_refseq.bed -v -s -wa > temp_GRO_no_TSS_interuptions.bed")
  options(scipen=0)
  
  ## Read data back in
  good_GRO =read.delim("temp_GRO_no_TSS_interuptions.bed", head=F)
  
  # Filter for good data
  GRO3 <- GRO2[GRO2$UniqueID %in% good_GRO[,4],]
  
  GRO3$Runon = ifelse(GRO3$strand=='+',   abs(GRO3$post_gene_end - GRO3$gene_ends),  abs(GRO3$post_gene_start - GRO3$gene_start))
  
  ## Save
  save(GRO3,file=outputFile)
  
  ## Erase temp files
  unlink("temp_GRO_no_TSS_interuptions.bed")
  unlink("temp_refseq.bed")
  unlink("temp_GRO.bed")

  return(GRO3)
}


