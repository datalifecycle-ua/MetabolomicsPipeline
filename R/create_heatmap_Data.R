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


