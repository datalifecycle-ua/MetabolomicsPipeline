#' Create a SummarizedExperiment from Targeted Metabolomics Data
#'
#' This function constructs a `SummarizedExperiment`
#' object from targeted metabolomics data.
#' It uses peak intensity data as the assay, 
#' sample metadata as `colData`, 
#' and chemical annotationsas `rowData`.
#' All inputs must be aligned via identifiers
#' in the first columns of the sample metadata and chemical annotation tables.
#'
#' @param chemical_annotation A `data.frame` where each row represents a chemical. 
#'   The first column must contain chemical identifiers matching the column names of `peak_data`.
#' @param sample_metadata A `data.frame` where each row represents a biological sample. 
#'   The first column must contain sample IDs matching the row names of `peak_data`.
#' @param peak_data A `data.frame` with samples as rows and chemicals as columns. 
#'   Sample names must be the first column. 
#'
#' @return A `SummarizedExperiment` object containing:
#'   - assays: the peak data
#'   - rowData: the chemical annotation
#'   - colData: the sample metadata
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @export
#'
#' @examples
#' # Assuming chemical_annotation_df, sample_metadata_df, and peak_data_df are defined
#' se <- createMetabolomicsSE(chemical_annotation_df, sample_metadata_df, peak_data_df)

create_met_se <- function(chemical_annotation,
                          sample_metadata,
                          peak_data) {
  
  # Check that all inputs are data.frames
  if (!is.data.frame(chemical_annotation)) {
    stop("`chemical_annotation` must be a data.frame.")
  }
  if (!is.data.frame(sample_metadata)) {
    stop("`sample_metadata` must be a data.frame.")
  }
  if (!is.data.frame(peak_data)) {
    stop("`peak_data` must be a data.frame.")
  }

  # Convert first column to rownames
  peak_data <- column_to_rownames(peak_data,
                                  var = colnames(peak_data)[1])

  chemical_annotation <- chemical_annotation %>%
    column_to_rownames(colnames(chemical_annotation)[1])

  sample_metadata <- sample_metadata %>%
    column_to_rownames(colnames(sample_metadata)[1])

  # Check if identifiers match dimensions
  if (!all(rownames(chemical_annotation) %in% colnames(peak_data))) {
    stop("Not all chemical identifiers in
        chemical_annotation are present in peak_data columns.")
  }
  if (!all(rownames(sample_metadata) %in% rownames(peak_data))) {
    stop(
      paste(
        "Not all sample IDs in sample_metadata are present
        in the rownames of peak_data.",
        "Please ensure first column of the peak data corresponds
        the first column of the sample_metadata",
        sep = " "
      )
    )
  }

  # Align row and column names
  raw_data <- peak_data[rownames(sample_metadata),
                        rownames(chemical_annotation)]



  # Create SummarizedExperiment
  exp <- SummarizedExperiment::SummarizedExperiment(
    assays = list(peak = t(raw_data)),
    rowData = chemical_annotation,
    colData = sample_metadata
  )

  return(exp)
}