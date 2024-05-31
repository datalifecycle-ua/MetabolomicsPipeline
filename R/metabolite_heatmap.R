#' Create metabolite heatmap
#' 
#' Create heatmaps which are arranged by the experimental conditions. 
#' 
#' @param data A SummarizedExperiment containing the Metabolon experiment data. 
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
#' @param Assay Which assay data to use for the heatmap (default="normalized").
#' 
#' @param ...  Additional arguments can be passed into the arrange function. 
#' This parameter will order the columns of the heatmap.  
#' 
#' @return A gtable class with all of the information to build the heatmap. To view
#' the heatmap use ggplotify::as.ggplot().
#' 
#' 
#' @examples
#' # load data
#' dat = MetabolomicsPipeline::demoDat
#' 
#' # Heatmap with one group
#' treat_heatmap <- metabolite_heatmap(dat,top_mets = 50,
#'                                    group_vars = "GROUP_NAME",
#'                                    strat_var = NULL,
#'                                    caption = "Heatmap Arranged By Group",
#'                                    Assay = "normalized",
#'                                    GROUP_NAME)
#' 
#' as.ggplot(treat_heatmap)
#'
#' 
#' @importFrom dplyr rename
#' @importFrom dplyr all_of
#' @importFrom pheatmap pheatmap 
#' @importFrom grDevices colorRampPalette 
#' @importFrom tibble column_to_rownames 
#' @import RColorBrewer
#' 
#' @export
#' 
#' 

# This function creates heatmap data
metabolite_heatmap<- function(data, top_mets=50, group_vars, strat_var = NULL,
                              caption=NULL,Assay="normalized",...){
  
  ## Get top metabolites
  select_variables <- SummarizedExperiment::assay(data,Assay) %>% 
    apply(1, mean)
  
  select_variables = select_variables[order(select_variables, decreasing = T)][1:top_mets]
  

  

############ Non Stratified case ##############################################  
  if(is.null(strat_var)){
    
      
    # Create analysis data
    analysis <- SummarizedExperiment::colData(data) %>%
      merge(t(SummarizedExperiment::assay(data,Assay)), by="row.names") %>%
      dplyr::rename(PARENT_SAMPLE_NAME=Row.names)
    
    heatmap_data = analysis%>%
        dplyr::select(all_of(c("PARENT_SAMPLE_NAME",
                                  group_vars,names(select_variables)))) %>%
        dplyr::arrange(...)
      
      
      rownames(heatmap_data) = NULL
      
      heatmap_meta_data <- heatmap_data %>%
        tibble::column_to_rownames("PARENT_SAMPLE_NAME") %>%
        dplyr::select(dplyr::all_of(group_vars))
      
      
      heatmap_data2 <- heatmap_data %>%
        dplyr::select(PARENT_SAMPLE_NAME,dplyr::all_of(names(select_variables))) %>%
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
    
    # Create analysis data
    analysis <- SummarizedExperiment::colData(data) %>%
      merge(t(SummarizedExperiment::assay(data,Assay)), by="row.names") %>%
      dplyr::rename(PARENT_SAMPLE_NAME=Row.names)
    
    heatmap_data = analysis%>%
      dplyr::select(all_of(c("PARENT_SAMPLE_NAME",
                             group_vars,names(select_variables)))) %>%
      dplyr::arrange(...)
    
    
    
    # Get stratas
    strats = split(analysis,analysis[[strat_var]])
    
    
    tabs <- lapply(names(strats),FUN = function(X){
      
      # Get heatmap data
      heatmap_data <- strats[[X]] %>%
        dplyr::select(all_of(c("PARENT_SAMPLE_NAME",
                               group_vars, names(select_variables)))) %>%
        dplyr::arrange(...)
      
      rownames(heatmap_data) = NULL
      
      heatmap_meta_data <- heatmap_data %>%
        tibble::column_to_rownames("PARENT_SAMPLE_NAME") %>%
        dplyr::select(dplyr::all_of(group_vars))
      
      
      heatmap_data2 <- heatmap_data %>%
        dplyr::select(PARENT_SAMPLE_NAME,dplyr::all_of(names(select_variables))) %>%
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






