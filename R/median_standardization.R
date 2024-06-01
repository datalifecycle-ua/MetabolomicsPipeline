#' Median standardization for metabolite data
#' 
#' This function standardizes the metabolites by the median of the metabolite.
#' 
#' @param peak_data Peak data with metabolites in the columns. The data also 
#' must include the "PARENT_SAMPLE_NAME". 
#' 
#' @returns Median standardized peak data. 
#' 
#' @examples
#' data("demoDataSmall", package = "MetabolomicsPipeline")
#' peak <- assay(demoDataSmall, "peak")
#' 
#' # Median standardization
#' peak_med <- median_standardization(peak_data = peak)
#' 
#' # Min value imputation
#' peakImpute <- min_val_impute(peak_data = peak_med)
#' 
#' # log transformation
#' peak_log <- log_transformation(peak_data = peakImpute)
#' 

#' 
#' 
#' @importFrom dplyr summarise_all
#' @export



median_standardization <- function(peak_data){
  # 1. Initialize a new peak_data_med matrix
  peak_data_std <-  as.data.frame(t(peak_data)) 
  
  # Summarize median of all metabolites
  peak_data_med <- peak_data_std %>% 
    #dplyr::select(-PARENT_SAMPLE_NAME) %>% 
    dplyr::summarise_all(median, na.rm = T) 
  
  
  # 3. Divide each value for each metabolite by the median value of that metabolite
  for(i in colnames(peak_data_med)){  
    peak_data_std[,i] <- as.matrix(peak_data_std[,i])/as.numeric(peak_data_med[,i]) 
  }  
  
  return(peak_data_std)
}
