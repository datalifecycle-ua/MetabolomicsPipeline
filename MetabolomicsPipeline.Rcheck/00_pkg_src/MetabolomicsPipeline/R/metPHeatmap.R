#' Metabolite Pairwise P-Value Interactive Heatmap.
#'
#' Produce an interactive heatmap of the p-values produced in
#' \code{\link{metabolitePairwise}}.
#'
#' @param results_data Results data frame of the pairwise comparisons produced
#'  by \code{\link{metabolitePairwise}}.
#'
#' @param data A SummarizedExperiment containing metabolomics experiment data.
#' 
#' @param pv_cutoff Numeric value for the p-value cutoff. Default is 0.05.
#'
#' @param interactive boolean (TRUE/FALSE) for whether or not the plot should be
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
#' strat_pairwise <- metabolitePairwise(dat,
#'     form = "GROUP_NAME*TIME1",
#'     strat_var = "Gender"
#' )
#' 
#' #############################################################################
#' ## Create Estimate Heatmap ##################################################
#' ############################################################################
#'
#' metEstHeatmap(strat_pairwise$Female, dat,
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
#' metPHeatmap(strat_pairwise$Female, dat,
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


metPHeatmap <- function(results_data, data, pv_cutoff = 0.05,
                        interactive = FALSE,
                         SUB_PATHWAY = "SUB_PATHWAY",
                        CHEMICAL_NAME = "CHEMICAL_NAME",
                        plotlyTitle = "P-Value Heatmap" ,...) {
  
  # check pv_cutoff is between 0 and 1
  if (pv_cutoff <= 0 || pv_cutoff >= 1) {
    stop("pv_cutoff must be between 0 and 1.")
  }

    # Check if data is a SummarizedExperiment
    if (!inherits(data, "SummarizedExperiment")) {
        stop("data must be a SummarizedExperiment object.")
    }

    # Check if results_data is a data frame, or DFrame
    if (!is.data.frame(results_data) && !inherits(results_data, "DFrame")) {
        stop("results_data must be a data frame or DFrame.")
    }
  
  # Get chemical ID
  CHEM_ID = "CHEM_ID"
    
    # Create data frame with all of the p-values
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
      by.x = "metabolite", by.y = "row.names") %>%
    as.data.frame() %>%
    dplyr::filter(
      Overall_pval <= pv_cutoff & 
        pval <= pv_cutoff)
  
  


    # Produce static heatmap
    if (interactive == FALSE) {
        # Create matrix for heatmap
        matr <- dat %>%
            reshape2::dcast(as.formula(paste0(CHEMICAL_NAME, "~", "Contrast")),
                value.var = "pval"
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
            z = ~pval,
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
