% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/min_val_impute.R
\name{min_val_impute}
\alias{min_val_impute}
\title{Minimum Value Imputation}
\usage{
min_val_impute(peak_data)
}
\arguments{
\item{peak_data}{Peak data matrix with metabolites in the columns.}
}
\value{
Metabolite imputed peak data.
}
\description{
Imputes the minimum value for each metabolite
}
\examples{
data("demoDataSmall", package = "MetabolomicsPipeline")
peak <- SummarizedExperiment::assay(demoDataSmall, "peak")

# Median standardization
peak_med <- median_standardization(peak_data = peak)

# Min value imputation
peakImpute <- min_val_impute(peak_data = peak_med)

# log transformation
peak_log <- log_transformation(peak_data = peakImpute)

}
