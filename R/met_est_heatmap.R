#' Metabolite Pairwise Estimate Interactive Heatmap.
#' 
#' Produce an interactive heatmap of the estimates produced in \code{\link{metabolite_pairwise}}.
#' 
#' @param results_data Results data frame of the pairwise comparisons produced by \code{\link{metabolite_pairwise}}.
#' 
#' @param data A SummarizedExperiment containing the Metabolon experiment data.   
#' 
#' @details
#' This function will produce a heatmap of the log fold changes for the metabolites
#'  with a significant overall p-value (which tested if the treatment group means
#'   were equal under the null hypothesis). The heatmap colors will only show if
#'    the log fold-change is greater than log(2) or less than log(.5). Therefore,
#'     this heatmap will only focus on comparisons with a fold change of two or greater.
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


met_est_heatmap <- function(results_data, data){
  
  #1. filter results
  results_data <- results_data  %>% 
    dplyr::filter(Overall_pval < 0.05) 
  
  # get chemical annotation file
  chem <- SummarizedExperiment::rowData(data) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "CHEM_ID")
  
  # 2. Merge the chemical annotation fill with the results from the pairwise
  #    comparisons.
  data <- chem %>% 
    dplyr::select(SUB_PATHWAY,CHEMICAL_NAME,CHEM_ID) %>% 
    merge(results_data, by.x = "CHEM_ID",by.y = "metabolite") %>% 
    dplyr::arrange(SUB_PATHWAY)  
  
  
  # 3. Produce Heatmap
  p <- data %>% 
    dplyr::select(CHEM_ID,SUB_PATHWAY,CHEMICAL_NAME,
                  all_of(names(data)[grepl("EST",names(data))])) %>%
    reshape2::melt(id.vars = c("CHEM_ID","SUB_PATHWAY","CHEMICAL_NAME"),
                   variable.name = "Contrast", 
         value.name = "logFoldChange") %>%
    dplyr::mutate(Contrast = gsub("_ESTS","",Contrast),
      logFoldChange = ifelse(logFoldChange < log(0.5) | logFoldChange > log(2),
                                   round(logFoldChange,3), NA)) %>% 
    plotly::plot_ly(
      type = "heatmap",
      x= ~Contrast,
      y = ~CHEMICAL_NAME,
      z = ~logFoldChange,
      text = ~SUB_PATHWAY,
      hovertemplate = paste("<b>Metabolite: %{y}</b><br><br>",
                            "Sub-pathway: %{text}<br>",
                            "Contrast: %{x}<br>",
                            "Log Fold-Change: %{z}<br>",
                            "<extra></extra>"),
      colorbar = list(title ="<b>Log Fold-Change</b>")) %>%
      plotly::layout(title = "<b>Log Fold-Change Heatmap</b>",
           xaxis = list(title="<b>Contrasts</b>"),
           yaxis = list(title = "")) 
  
  return(p)
  
  
}
