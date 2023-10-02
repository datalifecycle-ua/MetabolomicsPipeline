metabolite_pca <- function(analysis_data,non_metabolite){
  
  # Create PCA data containing only metabolite data
  pca_dat <- analysis_data[,!names(analysis_data) %in% non_metabolite] #<4>
  
  
  # Run PCA of the pca_dat matrix containing only the metabolites.   
  res.pca <- PCA(pca_dat, 
                 graph = F)
  
  
  # Create figure 
  pca_plot <- fviz_pca_ind(res.pca, 
                 label = "none",
                 habillage = as.factor(analysis_data[,meta_var])) 
  
  return(pca_plot)
}