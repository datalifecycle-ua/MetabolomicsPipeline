#' Metabolite PCA
#' 
#' Computes and plots the first two components of the PCA from the metabolite data.
#' 
#' @param analysis_data Analysis data, which contains the analysis variables and metabolites in the columns
#' 
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


metabolite_pca <- function(MetPipe,meta_var){
  
  # Get metabolite
  mets <- intersect(names(MetPipe@standardized_peak),
                    MetPipe@chemical_annotation$CHEM_ID)
  
  # Create PCA data containing only metabolite data
  pca_dat <- MetPipe@analysis[,mets]
  
  
  # Run PCA of the pca_dat matrix containing only the metabolites.   
  res.pca <- FactoMineR::PCA(pca_dat, 
                 graph = F)
  
  
  # Create figure 
  pca_plot <- factoextra::fviz_pca_ind(res.pca, 
                 label = "none",
                 habillage = as.factor(MetPipe@analysis[,meta_var])) 
  
  return(pca_plot)
}
