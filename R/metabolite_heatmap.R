#' Create metabolite heatmap
#' 
#' This function creates metabolite heatmaps
#' 
#' @param MetPipe Expirement MetPipe object 
#' 
#' @param top_mets Number of metabolites to include in heatmap. Metabolites
#' are choose based on highest variability.
#' 
#' @param group_vars Vector of variables to annotate heatmap with. Columns will
#' be grouped by these variables
#' 
#' @param caption A title for the heatmap.
#' 
#' @param strat_var Variable to stratify the heatmap by. 
#' 
#' @param ... Additional arguments that can be passed into the arrange function. 
#' This parameter will order the columns of the heatmap. 
#' 
#' @return Heatmap plot with the columns arranged as specified in ...
#' if a strata is specified then a list is returned, if strata=NULL
#' then a single heatmap is returned.
#' 
#' 
#' 
#' @import dplyr
#' @import pheatmap
#' @import grDevices
#' @import tibble
#' @import RColorBrewer
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
    pivot_longer(cols = everything()) %>%
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
                      main = caption) 
      
      
      # return plot
      
      return(map)
  }
  
  
  
  #################### Stratified case #########################################
  
  if(!is.null(strat_var)){
    
    # Get stratas
    strats = unique(MetPipe@analysis[,strat_var])
    
    
    # Create heatmaps list
    heats = list()
    
    # Run for loop
    for (i in 1:length(strats)) {
      
      heatmap_data = MetPipe@analysis[MetPipe@analysis[,strat_var]==strats[i],]%>%
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
                                main = paste0(caption," (",strats[i],")"))
    
      heats[[strats[i]]] = map
      
    }
    
    # Return list
    return(heats)
  
  }
  
}




