#' Median standardization for metabolite data
#' 
#' This function standardizes the metabolites by the median of the metabolite.
#' 
#' @param peak_data Peak data with metabolites in the columns. The data also 
#' must include the "Parent_SAMPLE_NAME". 
#' 
#' @returns Median standardized peak data. 
#' 
#' @import dplyr
#' @export



median_standardization <- function(peak_data){
  # 1. Initialize a new peak_data_med matrix
  peak_data_std <-  peak_data 
  
  # Summarize median of all metabolites
  peak_data_med <- peak_data_std %>% 
    dplyr::select(-PARENT_SAMPLE_NAME) %>% 
    dplyr::summarise_all(median, na.rm = T) 
  
  
  # 3. Divide each value for each metabolite by the median value of that metabolite
  for(i in colnames(peak_data_med)){  
    peak_data_std[,i] <- peak_data_std[,i]/peak_data_med[,i]  
  }  
  
  return(peak_data_norm)
}
