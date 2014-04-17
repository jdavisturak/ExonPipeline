setwd("/home/jeremy/ExonPipeline/hg19")

rm(list=ls())
library(gplots)           
library(RColorBrewer)
library(Hmisc)
source("/home/jeremy/ExonPipeline/ExonPipeline_functions.r")
source("/home/jeremy/ExonPipeline/ExonPipeline_analysis_and_plotting_Functions.r")
source("/home/home/jeremy/Code/useful_R/useful_R.r")

settings = setup(getwd())
refseq<-add_UniqueID(loadRefgene(settings))
features = loadFeatures(genome=genome<-loadGenome(settings<-setup( getwd())),settings=settings,refseq)

### Now merge with Veloso data
elongData = read.delim("SuppTable1_Veloso.txt")
refseq2 = refseq
refseq2$name2 = toupper(refseq$name2)
elongFeatures = merge(subset(refseq2,!grepl('hap',chrom)), elongData,by.x='name2',by.y="Gene.Symbol")
elongFeatures$start = ifelse(elongFeatures$strand=='+',elongFeatures$Start..analysis.region.,elongFeatures$End..analysis.region.)
elongFeatures$end = ifelse(elongFeatures$strand=='-',elongFeatures$Start..analysis.region.,elongFeatures$End..analysis.region.)
elongFeatures$Length = abs(elongFeatures$end-elongFeatures$start)
elongFeatures$isLast = 0
elongFeatures$FeatureCount=0

computeNucleosomeEnergy_Veloso = function(settings,unique_features,genome,WIDTH = Inf){

    ## Get output file name 
    settings$nucTestFile_Veloso = "regions_forSeqToNuc_Veloso.txt"
    
    ## Export specific sequences 
    if (!file.exists(settings$nucTestFile_Veloso)){
      firstSeq = 1
      cat(sprintf("Retrieving sequences of exons for nucleosome test Veloso (starting with seq %d):\n", firstSeq))
       
       seqs_for_nucleosomes = retrieveFixedSequences_min(settings, unique_features, genome, width=WIDTH, overhangUp=0, overhangDown=0, fromStart=TRUE,outputFile=settings$nucTestFile_Veloso,writeDuring=T,firstSeq=firstSeq)
        
    }

    ### Run BasicSeqToNuc2.jar
    settings$nucTestFile_Veloso_out = sprintf('%s_out.txt',settings$nucTestFile_Veloso)

    if (!file.exists(settings$nucTestFile_Veloso_out)){
      cat("Running Seq To Nuc 2 (on fixed exon lengths)\n")

      system(sprintf("java -Xmx1024m -jar ~/bin/BasicSeqToNuc.jar %s", settings$nucTestFile_Veloso))

    }     
    
    
    ## Compute the mean over the entire length.
    settings$nucTestFile_Veloso_mean = sub(".txt$","_mean.txt",settings$nucTestFile_Veloso)
    
    if (!file.exists(settings$nucTestFile_Veloso_mean)){
      cat("Averaging Seq To Nuc output\n")
     
      computeNucMean (settings, settings$nucTestFile_Veloso_out, settings$nucTestFile_Veloso_mean)
    }
    
    ## Normalize the nucleosome splicing score
    settings$nucNormScore_constrained = sub("_mean.txt$","_norm.txt",settings$nucTestFile_Veloso_mean)
   
    if (!file.exists(settings$nucNormScore_constrained)){
      cat("Normalizing Seq To Nuc output\n")
      normalizeToText(settings$nucTestFile_Veloso_mean, settings$nucNormScore_constrained)
    }

    return(settings)
}


seqs_for_nucleosomes = retrieveFixedSequences_min(settings, elongFeatures, genome, width=Inf, overhangUp=0, overhangDown=0, fromStart=TRUE,outputFile=settings$nucTestFile_Veloso,writeDuring=T,firstSeq=3001)
    
settings = computeNucleosomeEnergy_Veloso(settings,elongFeatures,genome,WIDTH = Inf)

