#' Create metabolite heatmap
#' 
#' This function creates metabolite heatmaps
#' 
#' @param .analysis_data The analysis data containing all of the metabolite and 
#' metadata information. 
#' 
#' @param heatmap_variables A vector of variable names that are NOT metabolites to include
#' in the heatmap.
#' 
#' @param metabolites a vector of metabolite names to include in the heatmap.
#' 
#' @param ... Additional arguments that can be passed into the arrange function. 
#' This parameter will order the columns of the heatmap. 
#' 
#' @return Heatmap plot with the columns arranged as specified in ...
#' 
#' 
#' 
#' @import dplyr
#' @import pheatmap
#' @import grDevices
#' 
#' @export
#' 
#' 



# This function creates heatmap data
metabolite_heatmap<- function(.analysis_data,heatmap_variables, metabolites,...){
  
  heatmap_data =.analysis_data %>%
   dplyr::arrange(...)
  
  heatmap_meta_data <- heatmap_data %>%
    dplyr::column_to_rownames("PARENT_SAMPLE_NAME") %>%
    dplyr::select(all_of(heatmap_variables))
  
  
  heatmap_data2 <- heatmap_data %>%
    dplyr::select(PARENT_SAMPLE_NAME,all_of(metabolites)) %>%
    dplyr::column_to_rownames("PARENT_SAMPLE_NAME") %>%
    as.matrix() %>% t()
  
  
  # Create heatamp
  
  # Heat map colors 
  palette <- grDevices::colorRampPalette(rev(brewer.pal(10, "RdBu")))(256) 
  
  
  # Create heatmap 
  map <- pheatmap::pheatmap(mat = heatmap_data2,cluster_cols = F, cluster_rows = F, color = palette, 
                  annotation_col = heatmap_meta_data, show_rownames = F, border_color = NA, show_colnames = F) 
  
  
  # return plot
  
  return(map)
  
}


