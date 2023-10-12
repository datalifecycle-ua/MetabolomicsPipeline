################################################################################
## Stratified Analysis #########################################################
################################################################################

# Stratified Analysis
stratified = subpathway_analysis(dat,
                                 treat_var = "GROUP_NAME",
                                 block_var = "TIME1",
                                 strat_var = "Gender")


################################################################################
### Results Plots ##############################################################
################################################################################

# 1. significant subpathways by model type
subpath_by_model(stratified)

# 2. Percentage of signficant subpathways within superpathways
subpath_within_superpath(stratified)

# 3. Metabolites within subpathway
tables <- met_within_sub(stratified, subpathway = "Partially Characterized Molecules")

### Females
tables[[1]]

### Males
tables[[2]]