#' Load Metabolomic Data as SummarizedExperiment
#'
#' Automatically load metabolomic data from excel file
#'
#' @param path Path to excel file with peak data, chemical
#'  annotations, sample meta data, and (optionally) the normalized peak counts
#'
#' @param raw_sheet Sheet name for the raw peak data.
#'
#' @param chemical_sheet Sheet name for chemical annotation.
#'
#' @param sample_meta Sheet name for sample meta data.
#'
#' @param normalized_peak Sheet name for the normalized peak data.
#'  If you are not adding the normalized data from the excel file then set
#'  normalized_peak=NA.
#'
#' @param sample_names Column name in the meta data containing the sample names.
#'  This must correspond to the row names of the raw peak data in the excel
#'  file.
#'
#' @param chemicalID Column name in the meta data containing the sample names.
#'  This must correspond to the column names of the raw peak data.
#'
#' @returns A SummarizedExperiment containing metabolomics expirement data.
#'
#' @details
#' The metabolomics experiment data are stored in a SummarizedExperiment.
#'
#' @seealso [SummarizedExperiment::SummarizedExperiment]
#'
#'
#'
#' @import readxl
#' @import methods
#' @import tibble
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
#' @export
#'

load_met_excel <- function(path,
    raw_sheet = "Peak Area Data",
    chemical_sheet = "Chemical Annotation",
    sample_meta = "Sample Meta Data",
    normalized_peak = "Log Transformed Data",
    sample_names = "PARENT_SAMPLE_NAME",
    chemicalID = "CHEM_ID"
) {
    # Get attributes
    raw_data <- readxl::read_excel(path, sheet = "Peak Area Data")
    chem_data <- readxl::read_excel(path, sheet = "Chemical Annotation")
    meta <- readxl::read_excel(path, sheet = "Sample Meta Data")

    if (!is.na(normalized_peak)) {
        normalized <- readxl::read_excel(path, sheet = "Log Transformed Data")
    }
    # Make row names align for summarized expirement
    meta <- tibble::column_to_rownames(meta, sample_names)
    chem_data <- tibble::column_to_rownames(chem_data, chemicalID)
    raw_data <- tibble::column_to_rownames(raw_data, sample_names)

    if (!is.na(normalized_peak)) {
        normalized <- tibble::column_to_rownames(normalized, sample_names)
    }

    # Align row and column names
    raw_data <- raw_data[rownames(meta), rownames(chem_data)]

    if (!is.na(normalized_peak)) {
        normalized <- normalized[rownames(meta), rownames(chem_data)]
    }

    # Load summarized expirement
    if (!is.na(normalized_peak)) {
        exp <- SummarizedExperiment::SummarizedExperiment(
            assays = list(
                peak = t(raw_data),
                normalized = t(normalized)
            ),
            colData = meta, rowData = chem_data
        )
    }

    if (is.na(normalized_peak)) {
        exp <- SummarizedExperiment::SummarizedExperiment(
            assays = list(peak = t(raw_data)),
            colData = meta, rowData = chem_data
        )
    }


    return(exp)
}
