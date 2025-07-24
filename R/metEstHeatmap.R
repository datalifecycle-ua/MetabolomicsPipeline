#' Metabolite Pairwise Estimate Interactive Heatmap.
#'
#' Produce an interactive heatmap of the estimates produced in
#'  \code{\link{metabolitePairwise}}.
#'
#' @param results_data Results data frame of the pairwise comparisons produced
#'  by \code{\link{metabolitePairwise}}.
#'
#' @param data A SummarizedExperiment containing the metabolomics experiment data.
#'
#' @param diff_cutoff Numeric value for the difference cutoff. Must be larger
#' default is .7.
#' 
#' @param pv_cutoff Numeric value for the p-value cutoff. Default is 0.05.
#'
#' @param interactive boolean (TRUE/FALSE) for whether or not the plot should be
#'     interactive. Use interactive=TRUE to produce an interactive plot using
#'     plotly. Use interactive=FALSE to produce a static heatmap using pheatmap.
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
#'
#'
#'
#' @details
#' This function will produce a heatmap of the log fold changes for the
#' metabolites with a significant overall p-value (which tested if the treatment
#'  group means were equal under the null hypothesis). The heatmap colors will
#'  only show if the log fold-change is greater than log(2) or less than
#'  log(.5). Therefore, this heatmap will only focus on comparisons with a
#'  fold change of two or greater.
#'
#' @examples
#' # Load data
#' data("demoDataSmall", package = "MetabolomicsPipeline")
#' dat <- demoDataSmall
#'
#' # Run pairwise analysis
#' strat_pairwise <- metabolitePairwise(dat,
#'     form = "GROUP_NAME*TIME1",
#'     strat_var = "Gender"
#' )
#' 
#' #############################################################################
#' ## Create Estimate Heatmap ##################################################
#' #############################################################################
#'
#' metEstHeatmap(strat_pairwise$Female, dat,
#'                interactive = FALSE,
#'                 SUB_PATHWAY = "SUB_PATHWAY",
#'                CHEMICAL_NAME = "CHEMICAL_NAME",
#'                plotlyTitle = "Metabolite log fold change",
#'                main = "Log fold change heatmap", show_rownames = FALSE
#' )
#'
#'
#'##############################################################################
#'## Create P-value Heatmap ####################################################
#'##############################################################################
#' # Female
#' metPHeatmap(strat_pairwise$Female, dat,
#'              interactive = FALSE, show_rownames = FALSE,
#'              plotlyTitle = "P-Value Heatmap",
#'              main = "Pvalue Heatmap"
#' )
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
metEstHeatmap <- function(results_data, data, diff_cutoff = .7,
                          pv_cutoff = 0.05,
                          interactive = FALSE,
                            SUB_PATHWAY = "SUB_PATHWAY",
                            CHEMICAL_NAME = "CHEMICAL_NAME",
                            plotlyTitle = "Metabolite log fold change", ...) {
  
    # Check if FC_cutoff is larger than 1
    if (diff_cutoff <= 0) {
        stop("diff cut off must be greater than 0")
    }
    # Check if pv_cutoff is >0<1
    if (pv_cutoff <= 0 || pv_cutoff >= 1) {
        stop("pv_cutoff must be between 0 and 1.")
    }
  
    # check if data is SummarizedExperiment
    if (!inherits(data, "SummarizedExperiment")) {
        stop("data must be a SummarizedExperiment object.")
    }
  
    CHEM_ID = "row.names"
    
    # Create data frame with all of the estimates and pvalues
    pvals <- results_data %>%
      as.data.frame() %>%
      dplyr::select(metabolite,Overall_pval,
                    dplyr::all_of(ends_with("PVALS"))) %>%
      tidyr::pivot_longer(-c(metabolite,Overall_pval),
                             names_to = "Contrast", values_to = "pval") %>%
      dplyr::mutate(Contrast = gsub("_PVALS", "", Contrast)) 
    
    
    ests <- results_data %>%
      as.data.frame() %>%
      dplyr::select(metabolite,Overall_pval, all_of(ends_with("ESTS"))) %>%
      tidyr::pivot_longer(-c(metabolite,Overall_pval), names_to = "Contrast",
                          values_to = "Difference") %>%
      dplyr::mutate(Contrast = gsub("_ESTS", "", Contrast))
    
    est_pvals <- merge(ests, pvals, by = c("metabolite", "Contrast","Overall_pval")) 
    
    # 2. Merge the chemical annotation fill with the results from the pairwise
    #    comparisons.
    dat <- est_pvals %>%
      merge(
        SummarizedExperiment::rowData(data),
        by.x = "metabolite", by.y = CHEM_ID) %>%
      as.data.frame() %>%
      dplyr::filter(
        Overall_pval <= pv_cutoff & 
        pval <= pv_cutoff &
          abs(Difference) >= diff_cutoff)



    # Produce static heatmap
    if (interactive == FALSE) {
        # Create matrix for heatmap
        matr <- dat %>%
            reshape2::dcast(as.formula(paste0(CHEMICAL_NAME, "~", "Contrast")),
                value.var = "Difference"
            )

        rownames(matr) <- matr[, 1]

        matr <- matr[, -1]

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
        mutate(
          text = paste0(
            "<br><b>Metabolite: ",!!as.name(CHEMICAL_NAME), "</b><br><br>",
            "<b>Subpathway: </b>", !!as.name(SUB_PATHWAY), "<br>",
            "<b>Contrast: </b>", Contrast, "<br>",
            "<b>Difference: </b>", round(Difference,3), "<br>",
            "<b>P-value:</b> ",round(pval,3)  ,"</b>")
        ) %>%
        plotly::plot_ly(
          type = "heatmap",
          x = ~Contrast,
          y = as.formula(paste0("~", CHEMICAL_NAME)),
          z = ~Difference,
          text = ~text,
          hoverinfo = 'text',
          colorbar = list(title = "<b>Difference</b>")
        ) %>%
        plotly::layout(
          title = paste0("<b>",plotlyTitle,"</b>"),
          xaxis = alist(title = "<b>Contrasts</b>"),
          yaxis = list(title = "")
        )
    }

    # Return heatmap
    return(p)
}
