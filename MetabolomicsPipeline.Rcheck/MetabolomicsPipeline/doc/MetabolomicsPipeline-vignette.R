## ----eval=FALSE---------------------------------------------------------------
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#     install.packages("BiocManager")
# }
# BiocManager::install("MetabolomicsPipeline", force = TRUE)

## ----eval=FALSE---------------------------------------------------------------
# if (!requireNamespace("devtools", quietly=TRUE))
#     install.packages("devtools")
# 
#  devtools::install_github("datalifecycle-ua/MetabolomicsPipeline", build_vignettes = TRUE)

## ----setup--------------------------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>",
    fig.width = 10,
    fig.height = 10
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

dat <- createMetSe(chemical_annotation = demoChemAnno,
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
dat <- medianStandardization(dat, assay = "peak")

# Min value imputation
dat <- minValImpute(dat, assay = "median_std") 

# log transformation
dat <- logTransformation(dat, assay = "min_impute")

## ----ExploratoryAnalysis_PCA--------------------------------------------------
###############################################################################
### Run PCA ###################################################################
###############################################################################

# Define PCA label from metadata
meta_var <- "Gender"

# Run PCA
pca <- metabolitePca(dat,
    meta_var = meta_var
)


# Show PCA
pca

## ----ExploratoryAnalysis_heatmap----------------------------------------------
################################################################################
### Run Heatmaps ###############################################################
################################################################################

# Heatmap with one group
metaboliteHeatmap(dat,
    top_mets = 50,
    group_vars = "GROUP_NAME",
    strat_var = NULL,
    caption = "Heatmap Arranged By Group",
    Assay = "normalized",
    GROUP_NAME
)




# Heatmap with two groups
metaboliteHeatmap(dat,
    top_mets = 50,
    group_vars = c("GROUP_NAME", "TIME1"),
    strat_var = NULL,
    caption = "Heatmap Arranged By Group and TIME",
    Assay = "normalized",
    GROUP_NAME, desc(TIME1)
)


# Heatmap with 2 group and stratified
metaboliteHeatmap(dat,
    top_mets = 50,
    group_vars = c("GROUP_NAME", "TIME1"),
    strat_var = "Gender",
    caption = "Heatmap Arranged By Group and TIME",
    Assay = "normalized",
    GROUP_NAME, desc(TIME1)
)



## ----SubpathwayAnalysis-------------------------------------------------------
################################################################################
## Stratified Analysis #########################################################
################################################################################

# Stratified Analysis
stratified <- subpathwayAnalysis(dat,
    treat_var = "GROUP_NAME",
    block_var = "TIME1",
    strat_var = "Gender",
    Assay = "normalized"
)


################################################################################
### Results Plots ##############################################################
################################################################################

# 1. significant subpathways by model type
subpathByModel(stratified)

# 2. Percentage of signficant subpathways within superpathways
subpathWithinSuperpath(stratified)

# Show significant subpathways within superpathways
subWithinSuperResults(stratified, Superpath = "Amino Acid")

# 3. Metabolites within subpathway
tables <- metWithinSub(stratified,
    subpathway = "Lysine Metabolism"
)


### Females
tables$Female

### Males
tables$Male

## ----PairwiseAnalysis---------------------------------------------------------
################################################################################
#### Run Pairwise Comparisons ##################################################
################################################################################

strat_pairwise <- metabolitePairwise(dat,
    form = "GROUP_NAME*TIME1",
    strat_var = "Gender"
)


###############################################################################
## Create Estimate Heatmap #####################################################
################################################################################

metEstHeatmap(strat_pairwise$Female, dat,
              diff_cutoff = .7, pv_cutoff = 0.05,
    interactive = TRUE,
     SUB_PATHWAY = "SUB_PATHWAY",
    CHEMICAL_NAME = "CHEMICAL_NAME",
    plotlyTitle = "Estimate Heatmap",
    main = "Log fold change heatmap", show_rownames = FALSE
)


################################################################################
## Create P-value Heatmap ######################################################
################################################################################
# Female
metPHeatmap(strat_pairwise$Female, dat,
            pv_cutoff = 0.05,
    interactive = TRUE, show_rownames = FALSE,
    plotlyTitle = "Pvalue Heatmap",
    main = "Pvalue Heatmap"
)

## ----BoxPlotsAndLinePlots-----------------------------------------------------
################################################################################
### BoxPlots ###################################################################
################################################################################

subpathwayBoxplots(dat,
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
subpathwayLineplots(dat,
    subpathway = "Lactoyl Amino Acid",
    block_var = TIME1, treat_var = GROUP_NAME,
    Assay = "normalized",
    CHEMICAL_NAME = "CHEMICAL_NAME",
   SUB_PATHWAY = "SUB_PATHWAY", Gender == "Female"
) +
    xlab("Time")

## ----Session Info-------------------------------------------------------------
sessionInfo()

