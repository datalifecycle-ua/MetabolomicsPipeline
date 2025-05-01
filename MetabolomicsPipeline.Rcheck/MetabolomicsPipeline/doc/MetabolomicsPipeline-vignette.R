## ----eval=FALSE---------------------------------------------------------------
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#     install.packages("BiocManager")
# }
# BiocManager::install("MetabolomicsPipeline")

## ----eval=FALSE---------------------------------------------------------------
# if (!requireNamespace("devtools", quietly=TRUE))
#     install.packages("devtools")
# 
#  devtools::install_github("datalifecycle-ua/MetabolomicsPipeline", build_vignettes = TRUE)

## ----setup, warning=FALSE, message=FALSE , include = TRUE---------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>",
    fig.width = 10,
    fig.height = 10,
    warning = FALSE
)

# Tables
library(table1)

# Load Metabolomics Pipeline
library(MetabolomicsPipeline)

library(SummarizedExperiment)

# Figures
library(ggplot2)

## ----AnalysisDataCreation-----------------------------------------------------
################################################################################
### Load Data ##################################################################
################################################################################

# if sample metadata, chemical annotation and peak data are stored in .xlsx
# you can use.

# dat <- load_met_excel(
  # path,
  # raw_sheet = "Peak Area Data",
  # chemical_sheet = "Chemical Annotation",
  # sample_meta = "Sample Meta Data",
  # normalized_peak = "Log Transformed Data",
  # sample_names = "PARENT_SAMPLE_NAME",
  # chemicalID = "CHEM_ID"
# )

# Get demo sample metadata
data("demoSampleMeta", package = "MetabolomicsPipeline")

# Get demo chemical annotation file
data("demoChemAnno", package = "MetabolomicsPipeline")

# Get demo peak data
data("demoPeak", package = "MetabolomicsPipeline")

dat <- create_met_se(chemical_annotation = demoChemAnno,
                     sample_metadata = demoSampleMeta,
                     peak_data = demoPeak,
                     chemical_id = "CHEM_ID",
                     sample_names = "PARENT_SAMPLE_NAME")


################################################################################
### Create Table 1 #############################################################
################################################################################
# Create table 1
tbl1 <- table1(~ GROUP_NAME + TIME1 | Gender,
  data = colData(dat)
)

# Display table 1
tbl1

## ----Data Processing----------------------------------------------------------
# Median standardization
dat <- median_standardization(dat, assay = "peak")

# Min value imputation
dat <- min_val_impute(dat, assay = "median_std") 

# log transformation
dat <- log_transformation(dat, assay = "min_impute")

## ----ExploratoryAnalysis_PCA--------------------------------------------------
###############################################################################
### Run PCA ###################################################################
###############################################################################

# Define PCA label from metadata
meta_var <- "Gender"

# Run PCA
pca <- metabolite_pca(dat,
    meta_var = meta_var
)


# Show PCA
pca

## ----ExploratoryAnalysis_heatmap----------------------------------------------
################################################################################
### Run Heatmaps ###############################################################
################################################################################

# Heatmap with one group
metabolite_heatmap(dat,
    top_mets = 50,
    group_vars = "GROUP_NAME",
    strat_var = NULL,
    caption = "Heatmap Arranged By Group",
    Assay = "normalized",
    sample_names = "PARENT_SAMPLE_NAME",
    GROUP_NAME
)




# Heatmap with two groups
metabolite_heatmap(dat,
    top_mets = 50,
    group_vars = c("GROUP_NAME", "TIME1"),
    strat_var = NULL,
    caption = "Heatmap Arranged By Group and TIME",
    Assay = "normalized",
    sample_names = "PARENT_SAMPLE_NAME",
    GROUP_NAME, desc(TIME1)
)


# Heatmap with 2 group and stratified
metabolite_heatmap(dat,
    top_mets = 50,
    group_vars = c("GROUP_NAME", "TIME1"),
    strat_var = "Gender",
    caption = "Heatmap Arranged By Group and TIME",
    Assay = "normalized",
    sample_names = "PARENT_SAMPLE_NAME",
    GROUP_NAME, desc(TIME1)
)



## ----SubpathwayAnalysis-------------------------------------------------------
################################################################################
## Stratified Analysis #########################################################
################################################################################

# Stratified Analysis
stratified <- subpathway_analysis(dat,
    treat_var = "GROUP_NAME",
    block_var = "TIME1",
    strat_var = "Gender",
    Assay = "normalized"
)


################################################################################
### Results Plots ##############################################################
################################################################################

# 1. significant subpathways by model type
subpath_by_model(stratified)

# 2. Percentage of signficant subpathways within superpathways
subpath_within_superpath(stratified)

# 3. Metabolites within subpathway
tables <- met_within_sub(stratified,
    subpathway = "Partially Characterized Molecules"
)

### Females
tables[[1]]

### Males
tables[[2]]

## ----PairwiseAnalysis---------------------------------------------------------
################################################################################
#### Run Pairwise Comparisons ##################################################
################################################################################

strat_pairwise <- metabolite_pairwise(dat,
    form = "GROUP_NAME*TIME1",
    strat_var = "Gender"
)


###############################################################################
## Create Estimate Heatmap #####################################################
################################################################################

met_est_heatmap(strat_pairwise$Female, dat,
    interactive = FALSE,
     SUB_PATHWAY = "SUB_PATHWAY",
    CHEMICAL_NAME = "CHEMICAL_NAME",
    plotlyTitle = "Estimate Heatmap",
    main = "Log fold change heatmap", show_rownames = FALSE
)


################################################################################
## Create P-value Heatmap ######################################################
################################################################################
# Female
met_p_heatmap(strat_pairwise$Female, dat,
    interactive = FALSE, show_rownames = FALSE,
    plotlyTitle = "Pvalue Heatmap",
    main = "Pvalue Heatmap"
)

## ----BoxPlotsAndLinePlots-----------------------------------------------------
################################################################################
### BoxPlots ###################################################################
################################################################################

subpathway_boxplots(dat,
    subpathway = "Lactoyl Amino Acid", block_var = TIME1,
    treat_var = GROUP_NAME, Assay = "normalized",
    CHEMICAL_NAME = "CHEMICAL_NAME",
   SUB_PATHWAY = "SUB_PATHWAY", Gender == "Female"
)


################################################################################
## Line plots ##################################################################
################################################################################

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
) +
    xlab("Time")

## ----Session Info-------------------------------------------------------------
sessionInfo()

