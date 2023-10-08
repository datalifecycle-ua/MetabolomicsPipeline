################################################################################
### Load Data ##################################################################
################################################################################

# MetPipe data
load("data/demo.rda")

###############################################################################
### Run PCA ###################################################################
###############################################################################

# Define PCA label from metadata
meta_var = "Gender"

# Run PCA
metabolite_pca( demo,
               meta_var = meta_var)

################################################################################
### Run Heatmaps ###############################################################
################################################################################

# Heatmap with one group
metabolite_heatmap(demo,top_mets = 50,
                   group_vars = "GROUP_NAME",
                   strat_var = NULL,
                   caption = "Heatmap Arranged By Group",
                   GROUP_NAME)


# Heatmap with two groups
metabolite_heatmap(demo,top_mets = 50,
                   group_vars = c("GROUP_NAME","TIME1"),
                   strat_var = NULL,
                   caption = "Heatmap Arranged By Group and TIME",
                   GROUP_NAME, desc(TIME1))


# Heatmap with 3 groups
 metabolite_heatmap(demo,top_mets = 50,
                   group_vars = c("GROUP_NAME","TIME1"),
                   strat_var = "Gender",
                   caption = "Heatmap Arranged By Group and TIME",
                   GROUP_NAME, desc(TIME1))


