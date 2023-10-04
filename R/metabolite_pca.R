#' Metabolite PCA
#' 
#' Computes and plots the first two components of the PCA from metabolite data.
#' 
#' @param analysis_data -  analysis data with metabolites in the columns
#' 
#' @param A vector of metabolites to include in the PCA. This is primarily to 
#' to remove columns that do not include metabolite peak data.
#' 
#' @param meta_var A metadata variable to color code the PCA plot by.  
#' 
#' @returns A PCA plot of the first two principal components, colored by the metadata
#' variable. 
#' 
#' @import FactoMineR
#' @import factoextra
#' 
#' @export


metabolite_pca <- function(analysis_data,metabolites,meta_var){
  
  # Create PCA data containing only metabolite data
  pca_dat <- analysis_data[,metabolites] #<4>
  
  
  # Run PCA of the pca_dat matrix containing only the metabolites.   
  res.pca <- FactoMineR::PCA(pca_dat, 
                 graph = F)
  
  
  # Create figure 
  pca_plot <- factoextra::fviz_pca_ind(res.pca, 
                 label = "none",
                 habillage = as.factor(analysis_data[,meta_var])) 
  
  return(pca_plot)
}