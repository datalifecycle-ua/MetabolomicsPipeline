#' Median standardization for metabolite data
#'
#' This function performs median standardization of metabolite
#'  data stored within a `SummarizedExperiment` object. 
#' For each metabolite, the values are divided by the
#'  median value of that metabolite across samples.
#' The standardized data are returned as a data frame
#'  and also added to the original `SummarizedExperiment`
#' object as a new assay named `"median_std"`.
#'
#' @param met_se A `SummarizedExperiment` object containing metabolite data.
#' @param assay A character string specifying the assay name
#'  within `met_se` to be median standardized. 
#' Default is `"peak"`.
#'
#' @return A data frame of median standardized metabolite data. 
#' The `SummarizedExperiment` object is also updated
#'  internally with a new assay `"median_std"`.
#'
#' @examples
#' library(SummarizedExperiment)
#' data("demoDataSmall", package = "MetabolomicsPipeline")
#'
#' # Median standardization
#' peak_med <- median_standardization(met_se = demoDataSmall, assay = "peak")
#'
#' # Access the median standardized data within the SummarizedExperiment
#' assay(demoDataSmall, "median_std")[1:5, 1:5]
#'
#' @importFrom dplyr summarise_all
#' @importFrom stats median
#' @importFrom SummarizedExperiment assay 
#' @export


median_standardization <- function(met_se, assay = "peak") {
  # 1. Initialize a new peak_data_med matrix
  peak_data_std <- as.data.frame(t(assay(met_se, assay)))

  # Summarize median of all metabolites
  peak_data_med <- peak_data_std %>%
    dplyr::summarise_all(median, na.rm = TRUE)


  # 3. Divide each value for each metabolite by the median value of that
  # metabolite
  for (i in colnames(peak_data_med)) {
    peak_data_std[, i] <- as.matrix(
      peak_data_std[, i]
    ) / as.numeric(peak_data_med[, i])
  }

  # Add Median standardized data to met_se
  assay(met_se, "median_std") = t(peak_data_std)


    return(met_se)
}