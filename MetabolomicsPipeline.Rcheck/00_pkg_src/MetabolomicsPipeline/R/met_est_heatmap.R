#' Metabolite Pairwise Estimate Interactive Heatmap.
#'
#' Produce an interactive heatmap of the estimates produced in
#'  \code{\link{metabolite_pairwise}}.
#'
#' @param results_data Results data frame of the pairwise comparisons produced
#'  by \code{\link{metabolite_pairwise}}.
#'
#' @param data A SummarizedExperiment containing the Metabolon experiment data.
#'
#' @param interactive boolean (T/F) for whether or not the plot should be
#'     interactive. Use interactive=T to produce an interactive plot using
#'     plotly. Use interactive=F to produce a static heatmap using pheatmap.
#'
#' 
#' @param CHEM_ID Column name in the chemical annotation worksheet that contains
#' the chemical ID. 
#' 
#' @param SUB_PATHWAY Column name in the chemical annotation worksheet which 
#' contains the subpathway information.
#' 
#' @param CHEMICAL_NAME Column name in the chemical annotation worksheet which
#' contains the chemical name. 
#' 
#' @param ... Additional arguments that can be passed to pheatmap.
#'
#' @details

#' This function will produce a heatmap of the log fold changes for the
#' metabolites with a significant overall p-value (which tested if the treatment
#'  group means were equal under the null hypothesis). The heatmap colors will
#'  only show if the log fold-change is greater than log(2) or less than
#'  log(.5). Therefore, this heatmap will only focus on comparisons with a
#'  fold change of two or greater.
#'
#'
#' @returns An interactive heatmap of pairwise estimates.
#'
#' @importFrom dplyr filter
#' @importFrom dplyr arrange
#' @importFrom dplyr mutate
#' @importFrom tibble rownames_to_column
#' @importFrom reshape2 melt
#' @importFrom plotly plot_ly
#' @importFrom plotly layout
#' @importFrom SummarizedExperiment rowData
#'
#' @export
#'


met_est_heatmap <- function(results_data, data, interactive = FALSE,
                            CHEM_ID="CHEM_ID",SUB_PATHWAY="SUB_PATHWAY",
                            CHEMICAL_NAME = "CHEMICAL_NAME",...) {
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
            dplyr::all_of(names(results_data)[grepl("ESTS",
                                                    names(results_data))])
        ) %>%
        reshape2::melt(
            id.vars = c(CHEM_ID, SUB_PATHWAY, CHEMICAL_NAME),
            variable.name = "Contrast", value.name = "logFoldChange"
        ) %>%
        dplyr::mutate(
            Contrast = gsub("_ESTS", "", Contrast),
            logFoldChange = ifelse(
              logFoldChange < log(0.5) | logFoldChange > log(2),
                round(logFoldChange, 3), NA
            )
        ) %>%
        dplyr::arrange(!!as.name(SUB_PATHWAY))



    # Produce static heatmap
    if (interactive == FALSE) {
        # Create matrix for heatmap
        matr <- dat %>%
            reshape2::dcast(as.formula(paste0(CHEMICAL_NAME,"~","Contrast")),
                            value.var = "logFoldChange")

        rownames(matr) <- matr[, 1]

        matr <- matr[, -1]

        # Create heatmap
        pal <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(10, "RdBu")))(256)

        p <- pheatmap::pheatmap(matr,
            cluster_rows = F, cluster_cols = F,
            color = pal, ...
        )
    }



    # 3. Produce interactive Heatmap
    if (interactive == TRUE) {
        p <- dat %>%
            plotly::plot_ly(
                type = "heatmap",
                x = ~Contrast,
                y = as.formula(paste0("~",CHEMICAL_NAME)),
                z = ~logFoldChange,
                text = as.formula(paste0("~",SUB_PATHWAY)),
                hovertemplate = paste(
                    "<b>Metabolite: %{y}</b><br><br>",
                    "Sub-pathway: %{text}<br>",
                    "Contrast: %{x}<br>",
                    "Log Fold-Change: %{z}<br>",
                    "<extra></extra>"
                ),
                colorbar = list(title = "<b>Log Fold-Change</b>")
            ) %>%
            plotly::layout(
                title = "<b>Log Fold-Change Heatmap</b>",
                xaxis = list(title = "<b>Contrasts</b>"),
                yaxis = list(title = "")
            )
    }

    # Return heatmap
    return(p)
}
