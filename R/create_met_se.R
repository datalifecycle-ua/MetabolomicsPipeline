#' Create a SummarizedExperiment from Targeted Metabolomics Data
#'
#' This function constructs a `SummarizedExperiment`
#' object from targeted metabolomics data.
#' It uses peak intensity data, 
#' sample metadata as `colData`, 
#' and chemical annotations as `rowData`.
#' All inputs must be aligned via identifiers (sample_names and chemicalID)
#'
#' @param chemical_annotation A `data.frame` where each row represents a chemical. 
#' 
#' @param sample_metadata A `data.frame` where each row represents a biological sample. 
#' 
#' @param peak_data A `data.frame` with samples as rows and chemicals as columns. 
#'   Sample names must be the first column.
#' 
#' @param sample_names Column name in the meta data containing the sample names.
#'  This must correspond to the row names of the raw peak data in the excel
#'  file.
#'
#' @param chemicalID Column name in the meta data containing the sample names.
#'  This must correspond to the column names of the raw peak data. 
#'
#' @return A `SummarizedExperiment` object containing:
#'   - assays: the peak data
#'   - rowData: the chemical annotation
#'   - colData: the sample metadata
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom dplyr select
#' @importFrom tibble column_to_rownames
#' @export
#'
#' @examples
#' # Assuming chemical_annotation_df, sample_metadata_df, and peak_data_df are defined
#' se <- createMetabolomicsSE(chemical_annotation_df, sample_metadata_df, peak_data_df)

create_met_se <- function(chemical_annotation,
                          sample_metadata,
                          peak_data,
                          sample_names = "PARENT_SAMPLE_NAME",
                          chemical_id = "CHEM_ID") {
  
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

  # Step 1: Move specified column to rownames
  peak_data <- peak_data %>%
    column_to_rownames(var = sample_names)

  # Step 2: Select the column dynamically
  peak_data <- peak_data %>%
    select(all_of(as.character(chemical_annotation[[chemical_id]])))

  chemical_annotation <- chemical_annotation %>%
    column_to_rownames(chemical_id)

  sample_metadata <- sample_metadata %>%
    column_to_rownames(sample_names)

  # Check if identifiers match dimensions
  if (!all(rownames(chemical_annotation) %in% colnames(peak_data))) {
    stop("Not all chemical identifiers in
        chemical_annotation (chemcial_id) are present in peak_data columns.")
  }
  if (!all(rownames(sample_metadata) %in% rownames(peak_data))) {
    stop(
      paste(
        "Not all sample IDs in sample_metadata are present
        in peak_data.",
        "Please ensure sample_names are a column in peak_data",
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
