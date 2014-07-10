WDs = c("/home/jeremy/ExonPipeline/hg19", "/home/jeremy/ExonPipeline/mm10")
for (WD in WDs){

    setwd(WD)

    rm(list=ls())
    library(gplots)           
    library(RColorBrewer)
    library(Hmisc)
    load("PlottedData_new3.RData") 
    source("/home/jeremy/ExonPipeline/ExonPipeline_functions.r")
    source("/home/jeremy/ExonPipeline/ExonPipeline_analysis_and_plotting_Functions.r")
    source("/home/home/jeremy/Code/useful_R/useful_R.r")
    
    settings = setup(getwd())
    refseq<-add_UniqueID(loadRefgene(settings))
    genome<-loadGenome(settings<-setup( getwd()))
    features = loadFeatures(genome=genome,settings=settings,refseq)

    if (WD=='/home/jeremy/ExonPipeline/hg19'){
        ### Now merge with Veloso data
        elongData = read.delim("SuppTable1_Veloso.txt")
        refseq2 = refseq
        refseq2$name2 = toupper(refseq$name2)
        elongFeatures = merge(subset(refseq2,!grepl('hap',chrom)), elongData,by.x='name2',by.y="Gene.Symbol")
        elongFeatures$start = ifelse(elongFeatures$strand=='+',elongFeatures$Start..analysis.region.,elongFeatures$End..analysis.region.) # my 'start' and 'end' are always in the direction of txn: start > end if strand is Neg
        elongFeatures$end = ifelse(elongFeatures$strand=='-',elongFeatures$Start..analysis.region.,elongFeatures$End..analysis.region.)
        elongFeatures$Length = abs(elongFeatures$end-elongFeatures$start)
        elongFeatures$isLast = 0
        elongFeatures$FeatureCount=0



        settingsV = computeNucleosomeEnergy_additional(settings,elongFeatures,genome,"regions_forSeqToNuc_Veloso.txt",WIDTH = Inf)
    }

    ## Ok, now take the regions used for computing splicing (see methods in ExonPipeline_simulate_model)
    # First, just get the entire gene region of interest (start with the second exon, until the end of the gene)
    elongFeatures = subset(mergedData,FeatureCount==1)
    elongFeatures$start <- elongFeatures$end # mergedData keeps track of INTRONS: start with the end of intron 1!
    elongFeatures$end = ifelse(elongFeatures$strand=='-',elongFeatures$txStart,elongFeatures$txEnd)
    elongFeatures$Length = abs(elongFeatures$end-elongFeatures$start)
    elongFeatures$isLast = 0
    elongFeatures$FeatureCount=0

    settings = computeNucleosomeEnergy_additional(settings,elongFeatures,genome, "mergedData_splicingSimpleRegions_forSeqToNuc.txt",WIDTH = Inf)
}








