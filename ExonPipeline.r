### Pipeline to do my analysis of exons in different species.

# Going to use Biostrings package to help me!  Also will require some mysql.

source("/home/jeremy/ExonPipeline/ExonPipeline_functions.r")

### Setup tasks:

# Set up file names and extensions
settings = setup( getwd())

# Load genome
genome = loadGenome(settings)

# Load refseq data
refseq = add_UniqueID(loadRefgene(settings))

# Load the features (exon, intron info)
features = loadFeatures(settings,refseq,genome)  # 	This excludes intron-less genes


### Compute Nucleosome Energy
## Old way
#settings = computeNucleosomeEnergy(settings,features,genome,WIDTH = 180)

## New way
settings = computeNucleosomeEnergy2(settings,features,genome,WIDTH = 147)

### Compute splice site strength
## Simply using 
settings = computeSpliceSiteStrength(settings,features,genome)

saveData(settings)#,file='settings.RData')
                                               
#source("/home/jeremy/ExonPipeline/ExonPipeline_plottingOnly.r")
#source("/home/jeremy/ExonPipeline/ExonPipeline_plottingOnly2.r")
source("/home/jeremy/ExonPipeline/ExonPipeline_plottingOnly3.r")

#settings = computeSpliceSiteStrength2(settings,features,genome)
#saveData(settings)#,file='settings.RData')


