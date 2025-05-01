pkgname <- "MetabolomicsPipeline"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('MetabolomicsPipeline')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("all_sig_subpath")
### * all_sig_subpath

flush(stderr()); flush(stdout())

### Name: all_sig_subpath
### Title: Table of Significant Subpathways
### Aliases: all_sig_subpath

### ** Examples

# Load data
data("demoDataSmall", package = "MetabolomicsPipeline")
dat <- demoDataSmall

# Runsubpathay analysis
sub_analysis <- subpathway_analysis(dat,
    treat_var = "GROUP_NAME",
    block_var = "TIME1",
    strat_var = NULL,
    Assay = "normalized"
)

#############################################################################
### Results Plots ###########################################################
#############################################################################

# significant subpathways by model type
subpath_by_model(sub_analysis)

# Percentage of signficant subpathways within superpathways
subpath_within_superpath(sub_analysis)

met_within_sub(sub_analysis, subpathway = "Aminosugar Metabolism")

# All signifiicant subpathways
all_sig_subpath(sub_analysis)




cleanEx()
nameEx("create_met_se")
### * create_met_se

flush(stderr()); flush(stdout())

### Name: create_met_se
### Title: Create a SummarizedExperiment from Metabolomics Data
### Aliases: create_met_se

### ** Examples

#Get demo sample metadata
data("demoSampleMeta", package = "MetabolomicsPipeline")

#Get demo chemical annotation file
data("demoChemAnno", package = "MetabolomicsPipeline")

# Get demo peak data
data("demoPeak", package = "MetabolomicsPipeline")

dat <- create_met_se(chemical_annotation = demoChemAnno,
                    sample_metadata = demoSampleMeta,
                    peak_data = demoPeak,
                    chemical_id = "CHEM_ID",
                    sample_names = "PARENT_SAMPLE_NAME")




cleanEx()
nameEx("log_transformation")
### * log_transformation

flush(stderr()); flush(stdout())

### Name: log_transformation
### Title: Log Transformation of Metabolite Data
### Aliases: log_transformation

### ** Examples

library(SummarizedExperiment)
data("demoDataSmall", package = "MetabolomicsPipeline")

# Median standardization
demoDataSmall <- median_standardization(met_se = demoDataSmall,
 assay = "peak")

# Minimum value imputation
demoDataSmall <- min_val_impute(met_se = demoDataSmall, assay = "median_std")

# Log transformation
demoDataSmall <- log_transformation(met_se = demoDataSmall,
 assay = "min_impute")

# Access log-transformed data
assay(demoDataSmall, "normalized")[1:5, 1:5]




cleanEx()
nameEx("median_standardization")
### * median_standardization

flush(stderr()); flush(stdout())

### Name: median_standardization
### Title: Median standardization for metabolite data
### Aliases: median_standardization

### ** Examples

library(SummarizedExperiment)
data("demoDataSmall", package = "MetabolomicsPipeline")

# Median standardization
peak_med <- median_standardization(met_se = demoDataSmall, assay = "peak")

# Access the median standardized data within the SummarizedExperiment
assay(peak_med, "median_std")[1:5, 1:5]




cleanEx()
nameEx("met_est_heatmap")
### * met_est_heatmap

flush(stderr()); flush(stdout())

### Name: met_est_heatmap
### Title: Metabolite Pairwise Estimate Interactive Heatmap.
### Aliases: met_est_heatmap

### ** Examples

# Load data
data("demoDataSmall", package = "MetabolomicsPipeline")
dat <- demoDataSmall

# Run pairwise analysis
strat_pairwise <- metabolite_pairwise(dat,
    form = "GROUP_NAME*TIME1",
    strat_var = "Gender"
)

#############################################################################
## Create Estimate Heatmap ##################################################
#############################################################################

met_est_heatmap(strat_pairwise$Female, dat,
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
met_p_heatmap(strat_pairwise$Female, dat,
             interactive = FALSE, show_rownames = FALSE,
             plotlyTitle = "P-Value Heatmap",
             main = "Pvalue Heatmap"
)





cleanEx()
nameEx("met_p_heatmap")
### * met_p_heatmap

flush(stderr()); flush(stdout())

### Name: met_p_heatmap
### Title: Metabolite Pairwise P-Value Interactive Heatmap.
### Aliases: met_p_heatmap

### ** Examples

# Load data
data("demoDataSmall", package = "MetabolomicsPipeline")
dat <- demoDataSmall

# Run pairwise analysis
strat_pairwise <- metabolite_pairwise(dat,
    form = "GROUP_NAME*TIME1",
    strat_var = "Gender"
)

#############################################################################
## Create Estimate Heatmap ##################################################
############################################################################

met_est_heatmap(strat_pairwise$Female, dat,
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
met_p_heatmap(strat_pairwise$Female, dat,
             interactive = FALSE, show_rownames = FALSE,
             plotlyTitle = "P-Value Heatmap",
             main = "Pvalue Heatmap"
)






cleanEx()
nameEx("met_within_sub")
### * met_within_sub

flush(stderr()); flush(stdout())

### Name: met_within_sub
### Title: Metabolites within Subpathway Table
### Aliases: met_within_sub

### ** Examples

data("demoDataSmall", package = "MetabolomicsPipeline")
dat <- demoDataSmall

# Runsubpathay analysis
sub_analysis <- subpathway_analysis(dat,
    treat_var = "GROUP_NAME",
    block_var = "TIME1",
    strat_var = NULL,
    Assay = "normalized"
)

#############################################################################
### Results Plots ###########################################################
#############################################################################

# significant subpathways by model type
subpath_by_model(sub_analysis)

# Percentage of signficant subpathways within superpathways
subpath_within_superpath(sub_analysis)

met_within_sub(sub_analysis, subpathway = "Aminosugar Metabolism")




cleanEx()
nameEx("metabolite_heatmap")
### * metabolite_heatmap

flush(stderr()); flush(stdout())

### Name: metabolite_heatmap
### Title: Create metabolite heatmap
### Aliases: metabolite_heatmap

### ** Examples

# load data
data("demoDat", package = "MetabolomicsPipeline")
dat <- demoDat

# Heatmap with one group
treat_heatmap <- metabolite_heatmap(dat,
    top_mets = 50,
    group_vars = "GROUP_NAME",
    strat_var = NULL,
    caption = "Heatmap Arranged By Group",
    Assay = "normalized",
    GROUP_NAME
)




cleanEx()
nameEx("metabolite_pairwise")
### * metabolite_pairwise

flush(stderr()); flush(stdout())

### Name: metabolite_pairwise
### Title: Metabolite Pairwise Comparisons.
### Aliases: metabolite_pairwise

### ** Examples

# Load data
data("demoDataSmall", package = "MetabolomicsPipeline")
dat <- demoDataSmall

# Run pairwise analysis
strat_pairwise <- metabolite_pairwise(dat,
    form = "GROUP_NAME*TIME1",
    strat_var = "Gender"
)

#############################################################################
## Create Estimate Heatmap #####################################################
################################################################################

met_est_heatmap(strat_pairwise$Female, dat,
               interactive = FALSE,
               CHEM_ID = "CHEM_ID", SUB_PATHWAY = "SUB_PATHWAY",
               CHEMICAL_NAME = "CHEMICAL_NAME",
               main = "Log fold change heatmap", show_rownames = FALSE
)


################################################################################
## Create P-value Heatmap ######################################################
################################################################################
# Female
met_p_heatmap(strat_pairwise$Female, dat,
             interactive = FALSE, show_rownames = FALSE,
             main = "Pvalue Heatmap"
)




cleanEx()
nameEx("metabolite_pca")
### * metabolite_pca

flush(stderr()); flush(stdout())

### Name: metabolite_pca
### Title: Metabolite PCA
### Aliases: metabolite_pca

### ** Examples


# load data
data("demoDat", package = "MetabolomicsPipeline")
dat <- demoDat

# Define PCA label from metadata
meta_var <- "Gender"

# Run PCA
pca <- metabolite_pca(dat,
    meta_var = meta_var
)


# Show PCA
pca




cleanEx()
nameEx("min_val_impute")
### * min_val_impute

flush(stderr()); flush(stdout())

### Name: min_val_impute
### Title: Minimum Value Imputation for Metabolite Data
### Aliases: min_val_impute

### ** Examples

library(SummarizedExperiment)
data("demoDataSmall", package = "MetabolomicsPipeline")

# Median standardization
demoDataSmall <- median_standardization(met_se = demoDataSmall,
 assay = "peak")

# Minimum value imputation
demoDataSmall <- min_val_impute(met_se = demoDataSmall, assay = "median_std")

# Access the imputed data
assay(demoDataSmall, "min_impute")[1:5, 1:5]




cleanEx()
nameEx("subpath_by_model")
### * subpath_by_model

flush(stderr()); flush(stdout())

### Name: subpath_by_model
### Title: Subpathway model type table
### Aliases: subpath_by_model

### ** Examples

# Load data
data("demoDataSmall", package = "MetabolomicsPipeline")
dat <- demoDataSmall

# Runsubpathay analysis
sub_analysis <- subpathway_analysis(dat,
    treat_var = "GROUP_NAME",
    block_var = "TIME1",
    strat_var = NULL,
    Assay = "normalized"
)

#############################################################################
### Results Plots ###########################################################
#############################################################################

# significant subpathways by model type
subpath_by_model(sub_analysis)

# Percentage of signficant subpathways within superpathways
subpath_within_superpath(sub_analysis)

met_within_sub(sub_analysis, subpathway = "Aminosugar Metabolism")




cleanEx()
nameEx("subpath_within_superpath")
### * subpath_within_superpath

flush(stderr()); flush(stdout())

### Name: subpath_within_superpath
### Title: Proportion of the Significant Subpathways Within Superpathways
### Aliases: subpath_within_superpath

### ** Examples


# Load data
data("demoDataSmall", package = "MetabolomicsPipeline")
dat <- demoDataSmall

# Runsubpathay analysis
sub_analysis <- subpathway_analysis(dat,
    treat_var = "GROUP_NAME",
    block_var = "TIME1",
    strat_var = NULL,
    Assay = "normalized"
)

#############################################################################
### Results Plots ###########################################################
#############################################################################

# significant subpathways by model type
subpath_by_model(sub_analysis)

# Percentage of signficant subpathways within superpathways
subpath_within_superpath(sub_analysis)

met_within_sub(sub_analysis, subpathway = "Aminosugar Metabolism")




cleanEx()
nameEx("subpathway_analysis")
### * subpathway_analysis

flush(stderr()); flush(stdout())

### Name: subpathway_analysis
### Title: Subpathway Analysis
### Aliases: subpathway_analysis

### ** Examples

# Load data
data("demoDataSmall", package = "MetabolomicsPipeline")
dat <- demoDataSmall

# Runsubpathay analysis
sub_analysis <- subpathway_analysis(dat,
    treat_var = "GROUP_NAME",
    block_var = "TIME1",
    strat_var = NULL,
    Assay = "normalized"
)

#############################################################################
### Results Plots ###########################################################
#############################################################################

# significant subpathways by model type
subpath_by_model(sub_analysis)

# Percentage of signficant subpathways within superpathways
subpath_within_superpath(sub_analysis)

met_within_sub(sub_analysis, subpathway = "Aminosugar Metabolism")

# All signifiicant subpathways
all_sig_subpath(sub_analysis)




cleanEx()
nameEx("subpathway_boxplots")
### * subpathway_boxplots

flush(stderr()); flush(stdout())

### Name: subpathway_boxplots
### Title: Subpathway Boxplots
### Aliases: subpathway_boxplots

### ** Examples

# load data
data("demoDat", package = "MetabolomicsPipeline")
dat <- demoDat

#############################################################################
### BoxPlots ###############################################################
############################################################################

subpathway_boxplots(dat,
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
subpathway_lineplots(dat,
                    subpathway = "Lactoyl Amino Acid",
                    block_var = TIME1, treat_var = GROUP_NAME,
                    Assay = "normalized",
                    CHEMICAL_NAME = "CHEMICAL_NAME",
                    SUB_PATHWAY="SUB_PATHWAY",Gender == "Female")






cleanEx()
nameEx("subpathway_lineplots")
### * subpathway_lineplots

flush(stderr()); flush(stdout())

### Name: subpathway_lineplots
### Title: Subpathway Lineplots
### Aliases: subpathway_lineplots

### ** Examples

data("demoDat", package = "MetabolomicsPipeline")
dat <- demoDat

#############################################################################
### BoxPlots ###############################################################
############################################################################

subpathway_boxplots(dat,
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
subpathway_lineplots(dat,
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
