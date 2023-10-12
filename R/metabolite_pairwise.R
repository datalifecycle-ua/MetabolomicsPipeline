#' Metabolite Pairwise Comparisons.
#' 
#' Computes the pairwise comparison estimates and p-values for each metabolite. 
#' 
#' 
#' @details This function will analyze each metabolite individually. For each metabolite, the metabolite_pairwise function
#' will first test whether the model explained a significant proportion
#' of the variance in the metabolite using an F-test. Since we will be looking at 
#' multiple comparisons for the metabolite, it is good practice to first look at the
#' overall p-value from the F-test before looking at the pairwise comparisons.
#' The metabolite_pairwise function then looks at all pairwise comparisons utilizing
#' the [emmeans](https://cran.r-project.org/web/packages/emmeans/index.html) package. 
#' The metabolite_pairwise function returns a data frame with the metabolite overall
#' p-value, log fold change for each group, and the p-value
#' for each comparison.
#' 
#' @param  form: This is a character string the resembles the right hand side of a 
#' simple linear regression model in R. For example form = "Group1 + Group2". 
#'
#' @param  Metpipe: The analysis data we will use for the pairwise comparisons. The 
#' data must be organized in a MetPipe object.  
#'
#' @param  strat_var: A variable in the analysis data to stratify the model by. If this
#' is specified, a list of results will be returned.  
#' 
#' @return The overall F-test p-value, and the estimate and pvalue for each pairwise comparison.
#' 
#' @import emmeans
#' 
#' @export
#' 
#' 




metabolite_pairwise <- function(MetPipe,form,strat_var=NULL){
  
  if(is.null(strat_var)){
    
    # Get metabolites
    mets <- intersect(names(MetPipe@analysis),MetPipe@chemical_annotation$CHEM_ID)
    
    # Get pairwise comparisons
    res <- apply(MetPipe@analysis[,mets], MARGIN=2, FUN = function(X){
                                                            pairwise(X,form = form,data = MetPipe@analysis)})
    
    results <- do.call(rbind,res)
    
    # Add metabolite names
    results$metabolite = rownames(results)
    
    return(results)
  }
  
  if(!is.null(strat_var)){
    
    # Get Metabolites
    mets <- intersect(names(MetPipe@analysis),MetPipe@chemical_annotation$CHEM_ID)
    
    # Split data
    data <- split(MetPipe@analysis,f=MetPipe@analysis[,strat_var])
    
    # Get results
    results <- lapply(data, function(X){
      
      # Get pairwise comparisons
      res <- apply(X[,mets], MARGIN=2, FUN = function(m){
        
        pairwise(m,form = form,data = X )
        
      })
      
      res2 <- do.call(rbind,res)
      
      # Add metabolite names
      res2$metabolite = rownames(res2)
      
      return(res2)
    })
   
  }

}
