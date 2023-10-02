## ----runStratifiedAnalysis----
# For each value perform the four steps listed above.
for(i in uni_vals){
  
  print(paste0("Running the ", i, " strata."))
  
  # 1. Subset analysis data
  strat_data = analysis_data[analysis_data[,stratified_var]==i,]
  
  
  # 2. Run subpathway function on the subsetted data
  strat_path_results <- subpathway_analysis(strat_data,
                                        chem_data = chem_data,
                                  block_var = blockvar,
                                  treat_var = treatvar)
  
  # 3 Save stratified data
  write.csv(strat_data, paste0(statified_data_path,"/analysis_data_",i,".csv"), row.names = F)
  
  # 4 Save results
  write.csv(strat_path_results, paste0(stratified_results_path,"/Stratified_",i,"_subpathway_results.csv"), row.names = F)
  
}
