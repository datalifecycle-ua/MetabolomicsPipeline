median_standardization <- function(peak_data){
  # 1. Initialize a new peak_data_med matrix
  peak_data_std <-  peak_data 
  
  # Summarize median of all metabolites
  peak_data_med <- peak_data_std %>% 
    select(-PARENT_SAMPLE_NAME) %>% 
    summarise_all(median, na.rm = T) 
  
  
  # 3. Divide each value for each metabolite by the median value of that metabolite
  for(i in colnames(peak_data_med)){  
    peak_data_std[,i] <- peak_data_std[,i]/peak_data_med[,i]  
  }  
  
  return(peak_data_norm)
}
