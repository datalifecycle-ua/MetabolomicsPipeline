# Load required packages for analysis
source("Code/Setup.R")

# Load Data from AnalysisDataCreation.R
dat <- readRDS("data/dat.Rds")


###############################################################################
### Run PCA ###################################################################
###############################################################################

# Define PCA label from metadata
meta_var = "Gender"

# Run PCA
pca <- metabolite_pca( dat,
                       meta_var = meta_var)


# Show heatmap
pca

################################################################################
### Run Heatmaps ###############################################################
################################################################################

# Heatmap with one group
treat_heatmap <- metabolite_heatmap(dat,top_mets = 50,
                                    group_vars = "GROUP_NAME",
                                    strat_var = NULL,
                                    caption = "Heatmap Arranged By Group",
                                    GROUP_NAME)


as.ggplot(treat_heatmap)


# Heatmap with two groups
treatandtime <-  metabolite_heatmap(dat,top_mets = 50,
                                    group_vars = c("GROUP_NAME","TIME1"),
                                    strat_var = NULL,
                                    caption = "Heatmap Arranged By Group and TIME",
                                    GROUP_NAME, desc(TIME1))


as.ggplot(treatandtime)


# Heatmap with 2 group and stratified
strat_heat <- metabolite_heatmap(dat,top_mets = 50,
                                 group_vars = c("GROUP_NAME","TIME1"),
                                 strat_var = "Gender",
                                 caption = "Heatmap Arranged By Group and TIME",
                                 GROUP_NAME, desc(TIME1))


# Heatmap with 2 group and stratified
strat_heat <- metabolite_heatmap(dat,top_mets = 50,
                                 group_vars = c("GROUP_NAME","TIME1"),
                                 strat_var = "Gender",
                                 caption = "Heatmap Arranged By Group and TIME",
                                 GROUP_NAME, desc(TIME1))


## Female
as.ggplot(strat_heat[[1]])

# Male
as.ggplot(strat_heat[[2]])
