## ----eval=FALSE---------------------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE)) {
#      install.packages("BiocManager")
#  }
#  BiocManager::install("MetabolomicsPipeline")

## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>",
    fig.width = 10,
    fig.height = 10,
    warning = FALSE
)

# Data
library(readxl)
library(dplyr)

# Tables
library(table1)

# Plots
library(ggplot2)
library(ggplotify)

# Load Metabolomics Pipeline
library(MetabolomicsPipeline)

library(SummarizedExperiment)

## ----AnalysisDataCreation-----------------------------------------------------
################################################################################
### Load Data ##################################################################
################################################################################

# Load Metabolon data from the Excel file
# dat <- loadMetabolon(path = "../data/Metabolon.xlsx")
data("demoDat", package = "MetabolomicsPipeline")

dat <- demoDat


################################################################################
### Create Table 1 #############################################################
################################################################################
# Create table 1
tbl1 <- table1(~ GROUP_NAME + TIME1 | Gender,
    data = colData(dat)
)

# Display table 1
tbl1

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
treat_heatmap <- metabolite_heatmap(dat,
    top_mets = 50,
    group_vars = "GROUP_NAME",
    strat_var = NULL,
    caption = "Heatmap Arranged By Group",
    Assay = "normalized",
    GROUP_NAME
)


as.ggplot(treat_heatmap)


# Heatmap with two groups
treatandtime <- metabolite_heatmap(dat,
    top_mets = 50,
    group_vars = c("GROUP_NAME", "TIME1"),
    strat_var = NULL,
    caption = "Heatmap Arranged By Group and TIME",
    Assay = "normalized",
    GROUP_NAME, desc(TIME1)
)


as.ggplot(treatandtime)


# Heatmap with 2 group and stratified
strat_heat <- metabolite_heatmap(dat,
    top_mets = 50,
    group_vars = c("GROUP_NAME", "TIME1"),
    strat_var = "Gender",
    caption = "Heatmap Arranged By Group and TIME",
    Assay = "normalized",
    GROUP_NAME, desc(TIME1)
)


## Female
as.ggplot(strat_heat[[1]])

# Male
as.ggplot(strat_heat[[2]])

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

strat_pairwise <- metabolite_pairwise(dat, form = "GROUP_NAME*TIME1", strat_var = "Gender")


###############################################################################
## Create Estimate Heatmap #####################################################
################################################################################

met_est_heatmap(strat_pairwise$Female, dat,
    interactive = FALSE,
     CHEM_ID="CHEM_ID",SUB_PATHWAY="SUB_PATHWAY",
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

## ----BoxPlotsAndLinePlots-----------------------------------------------------
################################################################################
### BoxPlots ###################################################################
################################################################################

subpathway_boxplots(dat,
    subpathway = "Lactoyl Amino Acid", block_var = TIME1,
    treat_var = GROUP_NAME, Assay = "normalized",
    CHEMICAL_NAME = "CHEMICAL_NAME",
    CHEM_ID="CHEM_ID",SUB_PATHWAY="SUB_PATHWAY",Gender == "Female")


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
    CHEM_ID="CHEM_ID",SUB_PATHWAY="SUB_PATHWAY",Gender == "Female"
) +
    xlab("Time") 

## ----Session Info-------------------------------------------------------------
sessionInfo()

