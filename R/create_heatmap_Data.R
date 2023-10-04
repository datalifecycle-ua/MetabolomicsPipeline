#' Create metadata and matricies for metabolite heatmaps
#' 
#' This function creates the required matrices for the metabolite heatmaps. 
#' 
#' @param .analysis_data The analysis data containing all of the metabolite and 
#' metadata information. 
#' 
#' @param heatmap_variables A vector of variable names that are NOT metabolites.
#' 
#' @param ... Additional arguments that can be passed into the arrange function. 
#' This parameter will order the columns of the heatmap data. 
#' 
#' @return A list of matricies including the heatmap variable (meta data for heatmap)
#' and the values for the heatmap.
#' 
#' 
#' @export
#' 



# This function creates heatmap data
create_heatmap_Data <- function(.analysis_data,heatmap_variables,...){
  
  heatmap_data =.analysis_data %>%
   arrange(...)
  
  heatmap_meta_data <- heatmap_data %>%
    column_to_rownames("PARENT_SAMPLE_NAME") %>%
    select(all_of(heatmap_variables))
  
  
  heatmap_data2 <- heatmap_data %>%
    select(-heatmap_variables[!(heatmap_variables %in% "PARENT_SAMPLE_NAME")]) %>%
    column_to_rownames("PARENT_SAMPLE_NAME") %>%
    as.matrix() %>% t()
  
  return(list( heatmap_variables = heatmap_meta_data, heatmap_data_vals = heatmap_data2))
  
}


