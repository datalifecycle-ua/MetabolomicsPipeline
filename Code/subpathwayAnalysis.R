################################################################################
### Load Data  #################################################################
################################################################################

# MetPipe data
load("data/demo.rda")

################################################################################
## Stratified Analysis #########################################################
################################################################################

# Non stratified results
non_stratified <- subpathway_analysis(demo,
                                      treat_var = "GROUP_NAME",
                                      block_var = "TIME1")


# Stratified Analysis
stratified = subpathway_analysis(demo,
                                     treat_var = "GROUP_NAME",
                                     block_var = "TIME1",
                                     strat_var = "Gender")


################################################################################
### Results Plots ##############################################################
################################################################################

# significant subpathways by model type
subpath_by_model(stratified)

# Percentage of signficant subpathways within superpathways
subpath_within_superpath(stratified)

# Metabolites within subpathway
met_within_sub(stratified, subpathway = "Partially Characterized Molecules")

