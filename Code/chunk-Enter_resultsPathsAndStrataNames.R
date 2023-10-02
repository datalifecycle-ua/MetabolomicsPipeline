## ----Enter_resultsPathsAndStrataNames----
# 1. Read in Results from the analysis step. 
results_path <- c("../Data/Results/Stratified_Male_subpathway_results.csv", #<1>
                  "../Data/Results/Stratified_Female_subpathway_results.csv") # <1>


# 2. Define the names of the strata
strata <- c("Male", "Female") #<2>
