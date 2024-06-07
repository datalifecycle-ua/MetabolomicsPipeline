#' Log transformation of metabolite data
#'
#' This function log transforms each metabolite in the Metabolon data.
#'
#' @param peak_data A matrix of peak data with metabolites in the columns
#'
#' @return log transformed peak data
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
#' @importFrom dplyr mutate_if
#' @importFrom magrittr `%>%`
#' @export
#'
#'
#'


log_transformation <- function(peak_data) {
    # 1. Log transform all of the values
    peak_data_log <- peak_data %>%
        dplyr::mutate_if(is.numeric, log)

    # Return log tranformed data
    return(peak_data_log)
}
