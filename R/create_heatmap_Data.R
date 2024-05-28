#' Create metadata and matricies for metabolite heatmaps
#' 
#' This function creates the required matrices for the metabolite heatmaps. 
#' 
#' @param data A SummarizedExperiment containing Metabolon data. 
#' 
#' @param heatmap_variables A vector of variable names that are NOT metabolites.
#' 
#' @param Assay Name of assay data to be used for heatmaps. Default="normalized".
#' 
#' @param ... Additional arguments that can be passed into the arrange function. 
#' This parameter will order the columns of the heatmap data. 
#' 
#' @return A list of matrices including the heatmap variable (meta data for heatmap)
#' and the values for the heatmap.
#' 
#' 
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment assay
#' 
#' @export
#' 



# This function creates heatmap data
create_heatmap_Data <- function(data,heatmap_variables,Assay="normalized",...){
  
  # Create analysis data
  meta <- SummarizedExperiment::colData(data)
  ass <- t(SummarizedExperiment::assay(data,Assay))
  
  analysis_data <- merge(meta,ass,by="row.names")
  
  
  heatmap_data =analysis_data %>%
   arrange(...)
  
  heatmap_meta_data <- heatmap_data %>%
    column_to_rownames("Row.names") %>%
    select(all_of(heatmap_variables))
  
  
  heatmap_data2 <- heatmap_data %>%
    select(-heatmap_variables[!(heatmap_variables %in% "Row.names")]) %>%
    column_to_rownames("Row.names") %>%
    as.matrix() %>% t()
  
  return(list( heatmap_variables = heatmap_meta_data, heatmap_data_vals = heatmap_data2))
  
}


