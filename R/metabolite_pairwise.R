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
#' @param  data SummarizedExperiment with Metabolon experiment data. 
#' 
#' @param  form This is a character string the resembles the right hand side of a 
#' simple linear regression model in R. For example form = "Group1 + Group2". 
#'
#'  
#' 
#' @param Assay Name of the assay to be used for the pairwise analysis (default='normalized')
#'
#' @param  strat_var A variable in the analysis data to stratify the model by. If this
#' is specified, a list of results will be returned. 
#' 
#' @param mets Chemical ID for the metabolites of interest. If NULL then the pairwise analysis
#' is completed for all metabololites.  
#' 
#' @return The overall F-test p-value, and the estimate and pvalue for each pairwise comparison.
#' 
#' 
#' @examples
#' # Load data
#' data("demoDat", package = "MetabolomicsPipeline")
#' dat = demoDat
#' 
#' # Run pairwise analysis
#' strat_pairwise = metabolite_pairwise(dat,form = "GROUP_NAME*TIME1",strat_var = "Gender")
#' 
#' 
#' @export
#' 
#' 




metabolite_pairwise <- function(data,form,Assay="normalized",strat_var=NULL,mets=NULL){
  
  # Create analysis data
  analysis <- SummarizedExperiment::colData(data) %>%
    merge(t(SummarizedExperiment::assay(data,Assay)), by="row.names") %>%
    dplyr::rename(PARENT_SAMPLE_NAME=Row.names)
  
  
  # Get metabolites
  if(is.null(mets)){
    mets <- rownames(SummarizedExperiment::rowData(data))
  }
  
  
  if(is.null(strat_var)){
    
    # Get metabolites
    
    
    # Get pairwise comparisons
    res <- apply(analysis[,mets], MARGIN=2, FUN = function(X){
                                                            pairwise(X,form = form,data = analysis)})
    
    results <- do.call(rbind,res)
    
    # Add metabolite names
    results$metabolite = rownames(results)
    
    return(results)
  }
  
  if(!is.null(strat_var)){
    
   
    
    # Split data
    data <- split(analysis,f= analysis[,strat_var])
    
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
