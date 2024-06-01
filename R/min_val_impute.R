#' Minimum Value Imputation
#' 
#' Imputes the minimum value for each metabolite
#' 
#' @param peak_data Peak data matrix with metabolites in the columns.
#' 
#' @returns Metabolite imputed peak data.
#' 
#' @examples
#' data("demoDataSmall", package = "MetabolomicsPipeline")
#' peak <- SummarizedExperiment::assay(demoDataSmall, "peak")
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
#' @importFrom dplyr select
#' @importFrom dplyr summarise_all
#' @export


min_val_impute <- function(peak_data){
  
  # 1. Initialize the new peak_data_imputed matrix
  peak_data_imputed <- as.data.frame(t(peak_data)) 
  
  
  # 2. Find the minimum value for each metabolite and compute 1/5 of that value
  peak_data_mins <- peak_data_imputed %>%  
   # dplyr::select(-PARENT_SAMPLE_NAME) %>%  
    dplyr::summarise_all(min, na.rm = T)  
  
  
  # 3. Impute the value
  for(i in colnames(peak_data_mins)){  
    if(length(peak_data_imputed[,i][is.na(peak_data_imputed[,i])]) > 0){ 
      peak_data_imputed[which(is.na(peak_data_imputed[,i])),i] <- as.numeric(peak_data_mins[i]) 
    } 
  }
  
  return(peak_data_imputed)
}