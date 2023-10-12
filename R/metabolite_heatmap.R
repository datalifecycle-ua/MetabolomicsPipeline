#' Create metabolite heatmap
#' 
#' Create heatmaps which are arranged by the experimental conditions. 
#' 
#' @param MetPipe Experiment MetPipe object 
#' 
#' @param top_mets Number of metabolites to include in the heatmap. Metabolites
#' are chosen based on the highest variability.
#' 
#' @param group_vars Vector of variables to annotate heatmap with. Columns will
#' be grouped by these variables.
#' 
#' @param caption A title for the heatmap. If strat_var is used, the title will
#' automatically include the stratum with the tile.
#' 
#' @param strat_var Variable to stratify the heatmap by. 
#' 
#' @param ...  Additional arguments can be passed into the arrange function. 
#' This parameter will order the columns of the heatmap.  
#' 
#' @return A gtable class with all of the information to build the heatmap. To view
#' the heatmap use ggplotify::as.ggplot().
#' 
#' 
#' 
#' @import dplyr
#' @import pheatmap
#' @import grDevices
#' @import tibble
#' @import RColorBrewer
#' @import tidyr
#' 
#' @export
#' 
#' 

# This function creates heatmap data
metabolite_heatmap<- function(MetPipe, top_mets=50, group_vars, strat_var = NULL,
                              caption=NULL,...){
  
  ## Get top metabolites
  select_variables <- MetPipe@standardized_peak %>% 
    select(-PARENT_SAMPLE_NAME) %>%
    summarise(across(everything(),\(x) mean(x,na.rm = T))) %>%
    tidyr::pivot_longer(cols = everything()) %>%
    arrange(desc(value)) %>%
    slice(c(1:top_mets))
  

############ Non Stratified case ##############################################  
  if(is.null(strat_var)){
    
      heatmap_data = MetPipe@analysis %>%
        dplyr::select(all_of(c("PARENT_SAMPLE_NAME",
                                  group_vars,select_variables$name))) %>%
        dplyr::arrange(...)
      
      
      rownames(heatmap_data) = NULL
      
      heatmap_meta_data <- heatmap_data %>%
        tibble::column_to_rownames("PARENT_SAMPLE_NAME") %>%
        dplyr::select(dplyr::all_of(group_vars))
      
      
      heatmap_data2 <- heatmap_data %>%
        dplyr::select(PARENT_SAMPLE_NAME,dplyr::all_of(select_variables$name)) %>%
        tibble::column_to_rownames("PARENT_SAMPLE_NAME") %>%
        as.matrix() %>% t()
      
      
      # Create heatamp
      
      # Heat map colors 
      palette <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(10, "RdBu")))(256) 
      
      
      # Create heatmap 
      map <- pheatmap::pheatmap(mat = heatmap_data2,cluster_cols = F, cluster_rows = F, color = palette, 
                      annotation_col = heatmap_meta_data, show_rownames = F, border_color = NA, show_colnames = F,
                      main = caption, silent = T) 
      
      
      # return plot
      
      return(map$gtable)
  }
  
  
  
  #################### Stratified case #########################################
  
  if(!is.null(strat_var)){
    
    # Get stratas
    strats = split(MetPipe@analysis,MetPipe@analysis[[strat_var]])
    
    
    tabs <- lapply(names(strats),FUN = function(X){
      
      # Get heatmap data
      heatmap_data <- strats[[X]] %>%
        dplyr::select(all_of(c("PARENT_SAMPLE_NAME",
                               group_vars, select_variables$name))) %>%
        dplyr::arrange(...)
      
      rownames(heatmap_data) = NULL
      
      heatmap_meta_data <- heatmap_data %>%
        tibble::column_to_rownames("PARENT_SAMPLE_NAME") %>%
        dplyr::select(dplyr::all_of(group_vars))
      
      
      heatmap_data2 <- heatmap_data %>%
        dplyr::select(PARENT_SAMPLE_NAME,dplyr::all_of(select_variables$name)) %>%
        tibble::column_to_rownames("PARENT_SAMPLE_NAME") %>%
        as.matrix() %>% t()
      
      
      # Create heatamp
      
      # Heat map colors 
      palette <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(10, "RdBu")))(256) 
      
      
      # Create heatmap 
      
      map <- pheatmap::pheatmap(mat = heatmap_data2,cluster_cols = F, cluster_rows = F, color = palette, 
                                annotation_col = heatmap_meta_data, show_rownames = F, border_color = NA, show_colnames = F,
                                main = paste0(caption," (",X,")"),silent = T)
      
      return( map$gtable)
      
    })
    
    
    # Return list
    return(tabs)
  
  }
  
}






