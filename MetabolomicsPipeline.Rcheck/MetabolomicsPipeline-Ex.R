pkgname <- "MetabolomicsPipeline"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('MetabolomicsPipeline')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("allSigSubpath")
### * allSigSubpath

flush(stderr()); flush(stdout())

### Name: allSigSubpath
### Title: Table of Significant Subpathways
### Aliases: allSigSubpath

### ** Examples

# Load data
data("demoDataSmall", package = "MetabolomicsPipeline")
dat <- demoDataSmall

# Runsubpathay analysis
sub_analysis <- subpathwayAnalysis(dat,
    treat_var = "GROUP_NAME",
    block_var = "TIME1",
    strat_var = NULL,
    Assay = "normalized"
)

#############################################################################
### Results Plots ###########################################################
#############################################################################

# significant subpathways by model type
subpathByModel(sub_analysis)

# Percentage of signficant subpathways within superpathways
subpathWithinSuperpath(sub_analysis)

metWithinSub(sub_analysis, subpathway = "Aminosugar Metabolism")

# All signifiicant subpathways
allSigSubpath(sub_analysis)




cleanEx()
nameEx("createMetSe")
### * createMetSe

flush(stderr()); flush(stdout())

### Name: createMetSe
### Title: Create a SummarizedExperiment from Metabolomics Data
### Aliases: createMetSe

### ** Examples

#Get demo sample metadata
data("demoSampleMeta", package = "MetabolomicsPipeline")

#Get demo chemical annotation file
data("demoChemAnno", package = "MetabolomicsPipeline")

# Get demo peak data
data("demoPeak", package = "MetabolomicsPipeline")

dat <- createMetSe(chemical_annotation = demoChemAnno,
                    sample_metadata = demoSampleMeta,
                    peak_data = demoPeak,
                    chemical_id = "CHEM_ID",
                    sample_names = "PARENT_SAMPLE_NAME")




cleanEx()
nameEx("logTransformation")
### * logTransformation

flush(stderr()); flush(stdout())

### Name: logTransformation
### Title: Log Transformation of Metabolite Data
### Aliases: logTransformation

### ** Examples

library(SummarizedExperiment)
data("demoDataSmall", package = "MetabolomicsPipeline")

# Median standardization
demoDataSmall <- medianStandardization(met_se = demoDataSmall,
 assay = "peak")

# Minimum value imputation
demoDataSmall <- minValImpute(met_se = demoDataSmall, assay = "median_std")

# Log transformation
demoDataSmall <- logTransformation(met_se = demoDataSmall,
 assay = "min_impute")

# Access log-transformed data
assay(demoDataSmall, "normalized")[1:5, 1:5]




cleanEx()
nameEx("medianStandardization")
### * medianStandardization

flush(stderr()); flush(stdout())

### Name: medianStandardization
### Title: Median standardization for metabolite data
### Aliases: medianStandardization

### ** Examples

library(SummarizedExperiment)
data("demoDataSmall", package = "MetabolomicsPipeline")

# Median standardization
peak_med <- medianStandardization(met_se = demoDataSmall, assay = "peak")

# Access the median standardized data within the SummarizedExperiment
assay(peak_med, "median_std")[1:5, 1:5]




cleanEx()
nameEx("metEstHeatmap")
### * metEstHeatmap

flush(stderr()); flush(stdout())

### Name: metEstHeatmap
### Title: Metabolite Pairwise Estimate Interactive Heatmap.
### Aliases: metEstHeatmap

### ** Examples

# Load data
data("demoDataSmall", package = "MetabolomicsPipeline")
dat <- demoDataSmall

# Run pairwise analysis
strat_pairwise <- metabolitePairwise(dat,
    form = "GROUP_NAME*TIME1",
    strat_var = "Gender"
)

#############################################################################
## Create Estimate Heatmap ##################################################
#############################################################################

metEstHeatmap(strat_pairwise$Female, dat,
               interactive = FALSE,
                SUB_PATHWAY = "SUB_PATHWAY",
               CHEMICAL_NAME = "CHEMICAL_NAME",
               plotlyTitle = "Metabolite log fold change",
               main = "Log fold change heatmap", show_rownames = FALSE
)


##############################################################################
## Create P-value Heatmap ####################################################
##############################################################################
# Female
metPHeatmap(strat_pairwise$Female, dat,
             interactive = FALSE, show_rownames = FALSE,
             plotlyTitle = "P-Value Heatmap",
             main = "Pvalue Heatmap"
)





cleanEx()
nameEx("metPHeatmap")
### * metPHeatmap

flush(stderr()); flush(stdout())

### Name: metPHeatmap
### Title: Metabolite Pairwise P-Value Interactive Heatmap.
### Aliases: metPHeatmap

### ** Examples

# Load data
data("demoDataSmall", package = "MetabolomicsPipeline")
dat <- demoDataSmall

# Run pairwise analysis
strat_pairwise <- metabolitePairwise(dat,
    form = "GROUP_NAME*TIME1",
    strat_var = "Gender"
)

#############################################################################
## Create Estimate Heatmap ##################################################
############################################################################

metEstHeatmap(strat_pairwise$Female, dat,
               interactive = FALSE,
               CHEM_ID = "CHEM_ID", SUB_PATHWAY = "SUB_PATHWAY",
               CHEMICAL_NAME = "CHEMICAL_NAME",
               plotlyTitle = "Estimate Heatmap", 
               main = "Log fold change heatmap", show_rownames = FALSE
)


##############################################################################
## Create P-value Heatmap ####################################################
##############################################################################
# Female
metPHeatmap(strat_pairwise$Female, dat,
             interactive = FALSE, show_rownames = FALSE,
             plotlyTitle = "P-Value Heatmap",
             main = "Pvalue Heatmap"
)






cleanEx()
nameEx("metWithinSub")
### * metWithinSub

flush(stderr()); flush(stdout())

### Name: metWithinSub
### Title: Metabolites within Subpathway Table
### Aliases: metWithinSub

### ** Examples

data("demoDataSmall", package = "MetabolomicsPipeline")
dat <- demoDataSmall

# Runsubpathay analysis
sub_analysis <- subpathwayAnalysis(dat,
    treat_var = "GROUP_NAME",
    block_var = "TIME1",
    strat_var = NULL,
    Assay = "normalized"
)

#############################################################################
### Results Plots ###########################################################
#############################################################################

# significant subpathways by model type
subpathByModel(sub_analysis)

# Percentage of signficant subpathways within superpathways
subpathWithinSuperpath(sub_analysis)

metWithinSub(sub_analysis, subpathway = "Aminosugar Metabolism")




cleanEx()
nameEx("metaboliteHeatmap")
### * metaboliteHeatmap

flush(stderr()); flush(stdout())

### Name: metaboliteHeatmap
### Title: Create metabolite heatmap
### Aliases: metaboliteHeatmap

### ** Examples

# load data
data("demoDat", package = "MetabolomicsPipeline")
dat <- demoDat

# Heatmap with one group
treat_heatmap <- metaboliteHeatmap(dat,
    top_mets = 50,
    group_vars = "GROUP_NAME",
    strat_var = NULL,
    caption = "Heatmap Arranged By Group",
    Assay = "normalized",
    GROUP_NAME
)




cleanEx()
nameEx("metabolitePairwise")
### * metabolitePairwise

flush(stderr()); flush(stdout())

### Name: metabolitePairwise
### Title: Metabolite Pairwise Comparisons.
### Aliases: metabolitePairwise

### ** Examples

# Load data
data("demoDataSmall", package = "MetabolomicsPipeline")
dat <- demoDataSmall

# Run pairwise analysis
strat_pairwise <- metabolitePairwise(dat,
    form = "GROUP_NAME*TIME1",
    strat_var = "Gender"
)

#############################################################################
## Create Estimate Heatmap #####################################################
################################################################################

metEstHeatmap(strat_pairwise$Female, dat,
               interactive = FALSE,
               CHEM_ID = "CHEM_ID", SUB_PATHWAY = "SUB_PATHWAY",
               CHEMICAL_NAME = "CHEMICAL_NAME",
               main = "Log fold change heatmap", show_rownames = FALSE
)


################################################################################
## Create P-value Heatmap ######################################################
################################################################################
# Female
metPHeatmap(strat_pairwise$Female, dat,
             interactive = FALSE, show_rownames = FALSE,
             main = "Pvalue Heatmap"
)




cleanEx()
nameEx("metabolitePca")
### * metabolitePca

flush(stderr()); flush(stdout())

### Name: metabolitePca
### Title: Metabolite PCA
### Aliases: metabolitePca

### ** Examples


# load data
data("demoDat", package = "MetabolomicsPipeline")
dat <- demoDat

# Define PCA label from metadata
meta_var <- "Gender"

# Run PCA
pca <- metabolitePca(dat,
    meta_var = meta_var
)


# Show PCA
pca




cleanEx()
nameEx("minValImpute")
### * minValImpute

flush(stderr()); flush(stdout())

### Name: minValImpute
### Title: Minimum Value Imputation for Metabolite Data
### Aliases: minValImpute

### ** Examples

library(SummarizedExperiment)
data("demoDataSmall", package = "MetabolomicsPipeline")

# Median standardization
demoDataSmall <- medianStandardization(met_se = demoDataSmall,
 assay = "peak")

# Minimum value imputation
demoDataSmall <- minValImpute(met_se = demoDataSmall, assay = "median_std")

# Access the imputed data
assay(demoDataSmall, "min_impute")[1:5, 1:5]




cleanEx()
nameEx("subWithinSuperResults")
### * subWithinSuperResults

flush(stderr()); flush(stdout())

### Name: subWithinSuperResults
### Title: Display Subpathways Within a Superpathway
### Aliases: subWithinSuperResults

### ** Examples

# 
data("demoDataSmall", package = "MetabolomicsPipeline")
dat <- demoDataSmall

# Runsubpathay analysis
sub_analysis <- subpathwayAnalysis(dat,
    treat_var = "GROUP_NAME",
    block_var = "TIME1",
    strat_var = NULL,
    Assay = "normalized"
)
#' # Display subpathways within a specific superpathway
subWithinSuperResults(sub_analysis, Superpath = "Amino Acid")




cleanEx()
nameEx("subpathByModel")
### * subpathByModel

flush(stderr()); flush(stdout())

### Name: subpathByModel
### Title: Subpathway model type table
### Aliases: subpathByModel

### ** Examples

# Load data
data("demoDataSmall", package = "MetabolomicsPipeline")
dat <- demoDataSmall

# Runsubpathay analysis
sub_analysis <- subpathwayAnalysis(dat,
    treat_var = "GROUP_NAME",
    block_var = "TIME1",
    strat_var = NULL,
    Assay = "normalized"
)

#############################################################################
### Results Plots ###########################################################
#############################################################################

# significant subpathways by model type
subpathByModel(sub_analysis)

# Percentage of signficant subpathways within superpathways
subpathWithinSuperpath(sub_analysis)

metWithinSub(sub_analysis, subpathway = "Aminosugar Metabolism")




cleanEx()
nameEx("subpathWithinSuperpath")
### * subpathWithinSuperpath

flush(stderr()); flush(stdout())

### Name: subpathWithinSuperpath
### Title: Proportion of the Significant Subpathways Within Superpathways
### Aliases: subpathWithinSuperpath

### ** Examples


# Load data
data("demoDataSmall", package = "MetabolomicsPipeline")
dat <- demoDataSmall

# Runsubpathay analysis
sub_analysis <- subpathwayAnalysis(dat,
    treat_var = "GROUP_NAME",
    block_var = "TIME1",
    strat_var = NULL,
    Assay = "normalized"
)

#############################################################################
### Results Plots ###########################################################
#############################################################################

# significant subpathways by model type
subpathByModel(sub_analysis)

# Percentage of signficant subpathways within superpathways
subpathWithinSuperpath(sub_analysis)

metWithinSub(sub_analysis, subpathway = "Aminosugar Metabolism")




cleanEx()
nameEx("subpathwayAnalysis")
### * subpathwayAnalysis

flush(stderr()); flush(stdout())

### Name: subpathwayAnalysis
### Title: Subpathway Analysis
### Aliases: subpathwayAnalysis

### ** Examples

# Load data
data("demoDataSmall", package = "MetabolomicsPipeline")
dat <- demoDataSmall

# Runsubpathay analysis
sub_analysis <- subpathwayAnalysis(dat,
    treat_var = "GROUP_NAME",
    block_var = "TIME1",
    strat_var = NULL,
    Assay = "normalized"
)

#############################################################################
### Results Plots ###########################################################
#############################################################################

# significant subpathways by model type
subpathByModel(sub_analysis)

# Percentage of signficant subpathways within superpathways
subpathWithinSuperpath(sub_analysis)

metWithinSub(sub_analysis, subpathway = "Aminosugar Metabolism")

# All signifiicant subpathways
allSigSubpath(sub_analysis)




cleanEx()
nameEx("subpathwayBoxplots")
### * subpathwayBoxplots

flush(stderr()); flush(stdout())

### Name: subpathwayBoxplots
### Title: Subpathway Boxplots
### Aliases: subpathwayBoxplots

### ** Examples

# load data
data("demoDat", package = "MetabolomicsPipeline")
dat <- demoDat

#############################################################################
### BoxPlots ###############################################################
############################################################################

subpathwayBoxplots(dat,
    subpathway = "Lactoyl Amino Acid", block_var = TIME1,
    treat_var = GROUP_NAME, Assay = "normalized",
    CHEMICAL_NAME = "CHEMICAL_NAME",
    SUB_PATHWAY = "SUB_PATHWAY", Gender == "Female"
)


############################################################################
## Line plots ##############################################################
############################################################################

# Set up data
dat$TIME1 <- as.numeric(factor(dat$TIME1,
    levels = c("PreSymp", "Onset", "End")
))
# Create line plots
subpathwayLineplots(dat,
                    subpathway = "Lactoyl Amino Acid",
                    block_var = TIME1, treat_var = GROUP_NAME,
                    Assay = "normalized",
                    CHEMICAL_NAME = "CHEMICAL_NAME",
                    SUB_PATHWAY="SUB_PATHWAY",Gender == "Female")






cleanEx()
nameEx("subpathwayLineplots")
### * subpathwayLineplots

flush(stderr()); flush(stdout())

### Name: subpathwayLineplots
### Title: Subpathway Lineplots
### Aliases: subpathwayLineplots

### ** Examples

data("demoDat", package = "MetabolomicsPipeline")
dat <- demoDat

#############################################################################
### BoxPlots ###############################################################
############################################################################

subpathwayBoxplots(dat,
    subpathway = "Lactoyl Amino Acid", block_var = TIME1,
    treat_var = GROUP_NAME, Assay = "normalized",
    CHEMICAL_NAME = "CHEMICAL_NAME",
    SUB_PATHWAY = "SUB_PATHWAY", Gender == "Female"
)


############################################################################
## Line plots ##############################################################
############################################################################

# Set up data
dat$TIME1 <- as.numeric(factor(dat$TIME1,
    levels = c("PreSymp", "Onset", "End")
))

# Create line plots
subpathwayLineplots(dat,
    subpathway = "Lactoyl Amino Acid",
    block_var = TIME1, treat_var = GROUP_NAME,
    Assay = "normalized",
    CHEMICAL_NAME = "CHEMICAL_NAME",
    SUB_PATHWAY = "SUB_PATHWAY", Gender == "Female"
)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
