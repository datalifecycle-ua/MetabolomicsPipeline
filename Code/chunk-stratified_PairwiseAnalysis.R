## ----stratified_PairwiseAnalysis----
# 1. Read in data from the analysis above
analysis_data <- read.csv(path3, check.names = F) #<1>

# 2. Make the variable names charactors
names(analysis_data) <- as.character(names(analysis_data)) #<2>


stratas <- unique(analysis_data[,strata_var])
for (i in 1:length(stratas)) {
  
  # 3. Filter analysis data
  dat <- analysis_data[analysis_data[,strata_var]==stratas[i],] #<3>
  
  
  # 4. Run the metabolite_pairwise function for the strata
  strat_results <- metabolite_pairwise(form = form, #<4>
                                       data=dat, #<4>
                                       metabolites = metabolites) #<4>
  
  # 5. Save results
  write.csv(strat_results, file = paste0(pair_strat_results_path,"/Pairwise_", #<5>
                                         stratas[i],".csv"), row.names = F) #<5>
  
  }
