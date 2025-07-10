#' Log Transformation of Metabolite Data
#'
#' This function performs a natural logarithm (log)
#' transformation on metabolite data stored within a 
#' `SummarizedExperiment` object. 
#' All numeric metabolite values are log-transformed,
#' and the resulting data are added as a new assay named `"normalized"`.
#'
#' @param met_se A `SummarizedExperiment` object containing metabolite data.
#' @param assay A character string specifying the assay 
#' name within `met_se` to log-transform.
#' Default is `"min_impute"`.
#'
#' @return The input `SummarizedExperiment` object with a
#' new assay `"normalized"` containing the log-transformed data.
#'
#' @examples
#' library(SummarizedExperiment)
#' data("demoDataSmall", package = "MetabolomicsPipeline")
#'
#' # Median standardization
#' demoDataSmall <- median_standardization(met_se = demoDataSmall,
#'  assay = "peak")
#'
#' # Minimum value imputation
#' demoDataSmall <- min_val_impute(met_se = demoDataSmall, assay = "median_std")
#'
#' # Log transformation
#' demoDataSmall <- log_transformation(met_se = demoDataSmall,
#'  assay = "min_impute")
#'
#' # Access log-transformed data
#' assay(demoDataSmall, "normalized")[1:5, 1:5]
#'
#' @importFrom dplyr mutate_if
#' @importFrom magrittr %>%
#' @importFrom SummarizedExperiment assay
#' @export
#'
#'


logTransformation <- function(met_se, assay = "min_impute") {

  # Get data for log transformation
  assay_dat = as.data.frame(t(assay(met_se, assay)))

  # 1. Log transform all of the values
  peak_data_log <- assay_dat %>%
    dplyr::mutate_if(is.numeric, log)

  assay(met_se, "normalized") = t(peak_data_log)
  
  # Return log tranformed data
  return(met_se)
}