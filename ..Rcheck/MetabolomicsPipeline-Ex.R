pkgname <- "MetabolomicsPipeline"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('MetabolomicsPipeline')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("metabolite_heatmap")
### * metabolite_heatmap

flush(stderr()); flush(stdout())

### Name: metabolite_heatmap
### Title: Create metabolite heatmap
### Aliases: metabolite_heatmap

### ** Examples

# load data
dat = MetabolomicsPipeline::demoDat

# Heatmap with one group
treat_heatmap <- metabolite_heatmap(dat,top_mets = 50,
                                   group_vars = "GROUP_NAME",
                                   strat_var = NULL,
                                   caption = "Heatmap Arranged By Group",
                                   Assay = "normalized",
                                   GROUP_NAME)






cleanEx()
nameEx("metabolite_pairwise")
### * metabolite_pairwise

flush(stderr()); flush(stdout())

### Name: metabolite_pairwise
### Title: Metabolite Pairwise Comparisons.
### Aliases: metabolite_pairwise

### ** Examples

# Load data
dat = MetabolomicsPipeline::demoDat

# Run pairwise analysis
strat_pairwise = metabolite_pairwise(dat,form = "GROUP_NAME*TIME1",strat_var = "Gender")





cleanEx()
nameEx("metabolite_pca")
### * metabolite_pca

flush(stderr()); flush(stdout())

### Name: metabolite_pca
### Title: Metabolite PCA
### Aliases: metabolite_pca

### ** Examples


# load data
dat = MetabolomicsPipeline::demoDat

# Define PCA label from metadata
 meta_var = "Gender"

# Run PCA
pca <- metabolite_pca( dat,
                      meta_var = meta_var)


# Show PCA
pca





cleanEx()
nameEx("subpathway_analysis")
### * subpathway_analysis

flush(stderr()); flush(stdout())

### Name: subpathway_analysis
### Title: Subpathway Analysis
### Aliases: subpathway_analysis

### ** Examples


# Load data
dat = MetabolomicsPipeline::demoDat

# Run subpathway analysis
stratified = subpathway_analysis(dat,
  treat_var = "GROUP_NAME",
  block_var = "TIME1",
  strat_var = NULL,
  Assay = "normalized")
  
################################################################################
### Results Plots ##############################################################
################################################################################

# significant subpathways by model type
subpath_by_model(stratified)

# Percentage of signficant subpathways within superpathways
subpath_within_superpath(stratified)





cleanEx()
nameEx("subpathway_boxplots")
### * subpathway_boxplots

flush(stderr()); flush(stdout())

### Name: subpathway_boxplots
### Title: Subpathway Boxplots
### Aliases: subpathway_boxplots

### ** Examples

# load data
dat = MetabolomicsPipeline::demoDat

################################################################################
### BoxPlots ###################################################################
################################################################################

subpathway_boxplots(dat, subpathway = "Lactoyl Amino Acid", block_var = TIME1,
                   treat_var = GROUP_NAME, Assay = "normalized",Gender =="Female")


################################################################################
## Line plots ##################################################################
################################################################################

# Set up data
dat$TIME1 <- as.numeric(factor(dat$TIME1,
                              levels = c("PreSymp","Onset","End")))
# Create line plots 
subpathway_lineplots(dat, subpathway = "Lactoyl Amino Acid",
                    block_var = TIME1,treat_var = GROUP_NAME, Assay = "normalized",Gender=="Female" ) 






cleanEx()
nameEx("subpathway_lineplots")
### * subpathway_lineplots

flush(stderr()); flush(stdout())

### Name: subpathway_lineplots
### Title: Subpathway Lineplots
### Aliases: subpathway_lineplots

### ** Examples

#' # load data
dat = MetabolomicsPipeline::demoDat

################################################################################
### BoxPlots ###################################################################
################################################################################

subpathway_boxplots(dat, subpathway = "Lactoyl Amino Acid", block_var = TIME1,
                   treat_var = GROUP_NAME, Assay = "normalized",Gender =="Female")


################################################################################
## Line plots ##################################################################
################################################################################

# Set up data
dat$TIME1 <- as.numeric(factor(dat$TIME1,
                              levels = c("PreSymp","Onset","End")))
# Create line plots 
subpathway_lineplots(dat, subpathway = "Lactoyl Amino Acid",
                    block_var = TIME1,treat_var = GROUP_NAME, Assay = "normalized",Gender=="Female" )







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
