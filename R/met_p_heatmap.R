#' Metabolite Pairwise P-Value Interactive Heatmap.
#' 
#' Produce an interactive heatmap of the pvalues produced in \code{\link{metabolite_pairwise}}.
#' 
#' @param results_data Results data frame of the pairwise comparisons produced by \code{\link{metabolite_pairwise}}.
#' 
#' @param chem_data The chemical annotation data provided by Metabolome
#' 
#' @details
#' For the metabolites which had a significant overall p-value (which tested if the
#' treatment group means were equal under the null hypothesis), we will produce a 
#' heatmap of the p-values.
#'
#'
#' @returns An interactive heatmap of pairwise pvalues. 
#' 
#' #' @details
#' For the metabolites which had a significant overall p-value (which tested if the
#' treatment group means were equal under the null hypothesis), we will produce a 
#' heatmap of the pairwise p-values. 
#' 
#' @import dplyr
#' @import reshape2
#' @import plotly
#' 
#' @export
#' 


met_p_heatmap <- function(results_data, chem_data){

  # 2. Merge the chemical annotation fill with the results from the pairwise comparisons.
  data <- chem_data %>% 
    dplyr::select(SUB_PATHWAY,CHEMICAL_NAME,CHEM_ID) %>% 
    merge(results_data, by.x = "CHEM_ID",by.y = "metabolite") %>% 
    dplyr::filter(Overall_pval < 0.05) %>%
    dplyr::arrange(SUB_PATHWAY)  
  
  
  # 3. Produce Heatmap
  p = data %>% 
    dplyr::select(CHEM_ID,SUB_PATHWAY,CHEMICAL_NAME, all_of(names(data)[grepl("PVALS",names(data))])) %>%  
    reshape2::melt(id.vars = c("CHEM_ID","SUB_PATHWAY","CHEMICAL_NAME"), variable.name = "Contrast",  
         value.name = "P_value") %>%  
    dplyr::mutate(Contrast = gsub("_PVALS","",Contrast), 
           P_value = round(P_value,3)) %>% 
    dplyr::arrange(SUB_PATHWAY) %>% 
    plotly::plot_ly( 
      type = "heatmap", 
      x= ~Contrast, 
      y = ~CHEMICAL_NAME, 
      z = ~P_value, 
      text = ~SUB_PATHWAY, 
      hovertemplate = paste("<b>Metabolite: %{y}</b><br><br>", 
                            "Sub-pathway: %{text}<br>", 
                            "Contrast: %{x}<br>", 
                            "P-Value: %{z}<br>", 
                            "<extra></extra>"), 
      colorbar = list(title ="<b>P-value</b>")) %>% 
    layout(title = "<b>P-Value Heatmap</b>", 
           xaxis = list(title="<b>Contrasts</b>"), 
           yaxis = list(title = "")) 
  
  
  # Return heatmap
  return(p)
}
