## ----Enter_stratifiedVars----
# 1. Get stratified variable and stratified values
strat_var <- "Gender" #<1>

values <- unique(analysis_data[,strat_var]) #<2>

# 2 Get heatmap Variables
heatmap_variables <- c("GROUP_NAME", #<3>
                       "TIME1") #<3>
