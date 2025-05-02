#' Metabolite Pairwise P-Value Interactive Heatmap.
#'
#' Produce an interactive heatmap of the p-values produced in
#' \code{\link{metabolite_pairwise}}.
#'
#' @param results_data Results data frame of the pairwise comparisons produced
#'  by \code{\link{metabolite_pairwise}}.
#'
#' @param data A SummarizedExperiment containing metabolomics experiment data.
#'
#' @param interactive boolean (T/F) for whether or not the plot should be
#'     interactive. Use interactive=T to produce an interactive plot using
#'     plotly. Use interactive=F to produce a static heatmap using pheatmap.
#'
#'
#' @param SUB_PATHWAY Column name in the chemical annotation worksheet which
#' contains the subpathway information.
#'
#' @param CHEMICAL_NAME Column name in the chemical annotation worksheet which
#' contains the chemical name.
#' 
#' @param plotlyTitle Title for the interactive heatmap.
#'
#' @param ... Additional arguments that can be passed to pheatmap.
#'
#' @details
#' For the metabolites which had a significant overall p-value (which tested if
#' the treatment group means were equal under the null hypothesis), we will
#' produce a heatmap of the p-values.
#'
#'
#' @returns An interactive heatmap of pairwise p-values.
#'
#'
#' @examples
#' # Load data
#' data("demoDataSmall", package = "MetabolomicsPipeline")
#' dat <- demoDataSmall
#'
#' # Run pairwise analysis
#' strat_pairwise <- metabolite_pairwise(dat,
#'     form = "GROUP_NAME*TIME1",
#'     strat_var = "Gender"
#' )
#' 
#' #############################################################################
#' ## Create Estimate Heatmap ##################################################
#' ############################################################################
#'
#' met_est_heatmap(strat_pairwise$Female, dat,
#'                interactive = FALSE,
#'                CHEM_ID = "CHEM_ID", SUB_PATHWAY = "SUB_PATHWAY",
#'                CHEMICAL_NAME = "CHEMICAL_NAME",
#'                plotlyTitle = "Estimate Heatmap", 
#'                main = "Log fold change heatmap", show_rownames = FALSE
#' )
#'
#'
#'##############################################################################
#'## Create P-value Heatmap ####################################################
#'##############################################################################
#' # Female
#' met_p_heatmap(strat_pairwise$Female, dat,
#'              interactive = FALSE, show_rownames = FALSE,
#'              plotlyTitle = "P-Value Heatmap",
#'              main = "Pvalue Heatmap"
#')
#'
#'
#'
#' @importFrom dplyr filter
#' @importFrom dplyr arrange
#' @importFrom dplyr mutate
#' @importFrom tibble rownames_to_column
#' @importFrom reshape2 melt
#' @importFrom plotly plot_ly
#' @importFrom SummarizedExperiment rowData
#' @importFrom stats as.formula
#'
#'
#' @export
#'


met_p_heatmap <- function(results_data, data, interactive = FALSE,
                         SUB_PATHWAY = "SUB_PATHWAY",
                        CHEMICAL_NAME = "CHEMICAL_NAME",
                        plotlyTitle = "P-Value Heatmap" ,...) {
  
  CHEM_ID = "CHEM_ID"
    # 2. Merge the chemical annotation fill with the results from the pairwise
    #    comparisons.
    dat <- rowData(data) %>%
        as.data.frame() %>%
        tibble::rownames_to_column(CHEM_ID) %>%
        dplyr::select(dplyr::all_of(c(SUB_PATHWAY, CHEMICAL_NAME, CHEM_ID))) %>%
        merge(results_data, by.x = CHEM_ID, by.y = "metabolite") %>%
        dplyr::filter(Overall_pval < 0.05) %>%
        dplyr::arrange(!!as.name(SUB_PATHWAY)) %>%
        dplyr::select(
            dplyr::all_of(c(CHEM_ID, SUB_PATHWAY, CHEMICAL_NAME)),
            all_of(names(results_data)[grepl("PVALS", names(results_data))])
        ) %>%
        reshape2::melt(
            id.vars = c(CHEM_ID, SUB_PATHWAY, CHEMICAL_NAME),
            variable.name = "Contrast", value.name = "P_value"
        ) %>%
        dplyr::mutate(
            Contrast = gsub("_PVALS", "", Contrast),
            P_value = ifelse(P_value < 0.05, round(P_value, 3), NA)
        ) %>%
        dplyr::arrange(!!as.name(SUB_PATHWAY))



    # Produce static heatmap
    if (interactive == FALSE) {
        # Create matrix for heatmap
        matr <- dat %>%
            reshape2::dcast(as.formula(paste0(CHEMICAL_NAME, "~", "Contrast")),
                value.var = "P_value"
            )

        rownames(matr) <- matr[, 1]

        matr <- matr[, -1]


        # Create row annotation
        # rowAnno <- dat %>%
        #   dplyr::select(CHEMICAL_NAME, SUB_PATHWAY)

        # Create heatmap
        pal <- grDevices::colorRampPalette(
            rev(RColorBrewer::brewer.pal(10, "RdBu"))
        )(256)

        p <- pheatmap::pheatmap(matr,
            cluster_rows = FALSE, cluster_cols = FALSE,
            color = pal, ...
        )
    }



    # 3. Produce interactive Heatmap
    if (interactive == TRUE) {
        p <- dat %>%
            plotly::plot_ly(
                type = "heatmap",
                x = ~Contrast,
                y = as.formula(paste0("~", CHEMICAL_NAME)),
                z = ~P_value,
                text = as.formula(paste0("~", SUB_PATHWAY)),
                hovertemplate = paste(
                    "<b>Metabolite: %{y}</b><br><br>",
                    "Subpathway: %{text}<br>",
                    "Contrast: %{x}<br>",
                    "P-Value: %{z}<br>",
                    "<extra></extra>"
                ),
                colorbar = list(title = "<b>P-value</b>")
            ) %>%
            plotly::layout(
                title = paste0("<b>",plotlyTitle,"</b>"),
                xaxis = list(title = "<b>Contrasts</b>"),
                yaxis = list(title = "")
            )
    }

    # Return heatmap
    return(p)
}
