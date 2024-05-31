#' Metabolite PCA
#' 
#' Computes and plots the first two components of the PCA from the metabolite data.
#' 
#' @param data SummarizedExperiment with Metabolon experiment data.
#' 
#' @param Assay Name of the assay to be used for the pairwise analysis (default='normalized')
#
#' @param meta_var A metadata variable to color code the PCA plot by.  
#' 
#' @returns A PCA plot of the first two principal components, colored by the metadata
#' variable. 
#' 
#' @examples
#' 
#' # load data
#' dat = MetabolomicsPipeline::demoDat
#' 
#' # Define PCA label from metadata
#'  meta_var = "Gender"
#' 
#' # Run PCA
#' pca <- metabolite_pca( dat,
#'                       meta_var = meta_var)
#'
#'
#' # Show PCA
#' pca
#' 
#' 
#' @importFrom  FactoMineR PCA
#' @importFrom  factoextra fviz_pca_ind
#' 
#' @export


metabolite_pca <- function(data, Assay ="normalized",meta_var){
  
  # Create analsysis data
  # Create analysis data
  analysis <- SummarizedExperiment::colData(data) %>%
    merge(t(SummarizedExperiment::assay(data,Assay)), by="row.names") %>%
    dplyr::rename(PARENT_SAMPLE_NAME=Row.names)
  
  
  # Run PCA of the pca_dat matrix containing only the metabolites.   
  res.pca <- FactoMineR::PCA(analysis[,rownames(data)], 
                 graph = F)
  
  
  # Create figure 
  pca_plot <- factoextra::fviz_pca_ind(res.pca, 
                 label = "none",
                 habillage = as.factor(analysis[,meta_var])) 
  
  return(pca_plot)
}
