#' Subpathway Lineplots
#'
#' Create line plots for each metabolite within a subpathway.
#'
#' @param data SummarizedExperiment with metabolomics experiment data.
#'
#' @param subpathway Character value of the subpathway of interest. This is case
#' sensitive and must be in the chemical annotation file.
#'
#' @param block_var This the the name of the variable in the meta data that is
#' used for the X axis of the line plots. We recommend using the "block_var"
#' variable from the subpathway analyis.
#'
#' @param treat_var This is a grouping variable. As a recommendation the
#' treatment groups should be used in the groupBy argument as this will provide
#' a different color for each of the treatments making it easier to identify.
#'
#' @param Assay Name of the assay to be used for the pairwise analysis
#'  (default='normalized')
#'
#'
#' @param SUB_PATHWAY Column name in the chemical annotation worksheet which
#' contains the subpathway information.
#'
#' @param CHEMICAL_NAME Column name in the chemical annotation worksheet which
#' contains the chemical name.
#'
#' @param ... Additional arguments to filter the analysis data by.
#'
#'
#' @returns Line plots stratified by metabolite.
#'
#'
#' @examples
#' data("demoDat", package = "MetabolomicsPipeline")
#' dat <- demoDat
#'
#' #############################################################################
#' ### BoxPlots ###############################################################
#' ############################################################################
#'
#' subpathwayBoxplots(dat,
#'     subpathway = "Lactoyl Amino Acid", block_var = TIME1,
#'     treat_var = GROUP_NAME, Assay = "normalized",
#'     CHEMICAL_NAME = "CHEMICAL_NAME",
#'     SUB_PATHWAY = "SUB_PATHWAY", Gender == "Female"
#' )
#'
#'
#' ############################################################################
#' ## Line plots ##############################################################
#' ############################################################################
#'
#' # Set up data
#' dat$TIME1 <- as.numeric(factor(dat$TIME1,
#'     levels = c("PreSymp", "Onset", "End")
#' ))
#'
#' # Create line plots
#' subpathwayLineplots(dat,
#'     subpathway = "Lactoyl Amino Acid",
#'     block_var = TIME1, treat_var = GROUP_NAME,
#'     Assay = "normalized",
#'     CHEMICAL_NAME = "CHEMICAL_NAME",
#'     SUB_PATHWAY = "SUB_PATHWAY", Gender == "Female"
#' )
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment rowData
#'
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr rename
#' @importFrom dplyr filter
#' @importFrom dplyr select
#'
#'
#' @export
#'
#'


subpathwayLineplots <- function(data, subpathway, block_var, treat_var,
                                Assay = "normalized",
                                CHEMICAL_NAME = "CHEMICAL_NAME",
                                SUB_PATHWAY = "SUB_PATHWAY", ...) {
  
    # Create chem_id variable name
     CHEM_ID = "rownames"
  
    # Create analysis data
    analysis <- SummarizedExperiment::colData(data) %>%
        merge(t(SummarizedExperiment::assay(data, Assay)), by = "row.names")
    #%>%
     #   dplyr::rename(PARENT_SAMPLE_NAME = Row.names)


    # Create chem data
    chem <- SummarizedExperiment::rowData(data) %>%
        as.data.frame() %>%
        tibble::rownames_to_column(CHEM_ID)


    return(
        analysis %>%
            dplyr::filter(...) %>%
            dplyr::select(
                treat = {{ treat_var }}, X = {{ block_var }},
                as.character(chem[
                    which(chem[, SUB_PATHWAY] == subpathway),
                    CHEM_ID
                ])
            ) %>%
            tidyr::pivot_longer(cols = -c(treat, X)) %>%
            merge(chem, by.x = "name", by.y = CHEM_ID) %>%
            ggplot2::ggplot(aes(x = X, y = value, color = treat)) +
            ggplot2::geom_jitter() +
            ggplot2::geom_smooth(method = "lm") +
            ggplot2::facet_wrap(as.formula(paste0("~", CHEMICAL_NAME))) +
            ggplot2::theme_bw()
    )
}
