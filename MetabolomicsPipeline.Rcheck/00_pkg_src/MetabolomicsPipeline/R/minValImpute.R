#' Minimum Value Imputation for Metabolite Data
#'
#' This function imputes missing values in metabolite
#' data stored within a `SummarizedExperiment` object.
#' For each metabolite, missing values are replaced
#' with the minimum observed value for that metabolite.
#' The imputed data are added to the `SummarizedExperiment`
#' object as a new assay named `"min_impute"`.
#'
#' @param met_se A `SummarizedExperiment` object containing metabolite data.
#' @param assay A character string specifying
#' the assay name within `met_se` to perform minimum value imputation on.
#' Default is `"median_std"`.
#'
#' @return The input `SummarizedExperiment` object with a
#' new assay `"min_impute"` containing the minimum value-imputed data.
#'
#' @examples
#' library(SummarizedExperiment)
#' data("demoDataSmall", package = "MetabolomicsPipeline")
#'
#' # Median standardization
#' demoDataSmall <- medianStandardization(met_se = demoDataSmall,
#'  assay = "peak")
#'
#' # Minimum value imputation
#' demoDataSmall <- minValImpute(met_se = demoDataSmall, assay = "median_std")
#'
#' # Access the imputed data
#' assay(demoDataSmall, "min_impute")[1:5, 1:5]
#'
#' @importFrom dplyr summarise_all
#' @importFrom SummarizedExperiment assay
#' @export

minValImpute <- function(met_se, assay = "median_std") {
  # 1. Initialize the new peak_data_imputed matrix
  peak_data_imputed <- as.data.frame(t(assay(met_se, assay)))


  # 2.Find the minimum value for each metabolite and compute
  peak_data_imputed <- apply(peak_data_imputed,2, function(col) {
    if (is.numeric(col)) {
      min_val <- min(col, na.rm = TRUE)
      col[is.na(col)] <- min_val
    }
    return(col)
  })
  
  # Add min_imput assay element
  assay(met_se, "min_impute") = t(peak_data_imputed)

return(met_se)
}