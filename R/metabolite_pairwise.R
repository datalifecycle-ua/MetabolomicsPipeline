#' Metabolite Pairwise Comparisons.
#' 
#' Computes the pairwise comparison estimates and pvalues for each metabolite. 
#' 
#' 
#' @details This function will analyze each metabolite individually. For each metabolite, the the metabolite_pairwise function
#' will first test if any of the experimental groups explained a significant proportion
#' of the variance in the metabolite using an F-test. Since we will be looking at 
#' multiple comparisons for the metabolite it is good practice to first look at the
#' over-all p-value from the F-test before looking at the pairwise comparisons.
#' The metabolite_pairwise function then looks at all pairwise comparisons utilizing
#' the [emmeans](https://cran.r-project.org/web/packages/emmeans/index.html) package. 
#' The metabolite_pairwise function returns a data frame with the metabolite, overall
#' p-value, estimated difference in means for all experimental group combinations, and the p-value
#' which test if the difference the groups is significantly different. 
#' 
#' @param  form: This is a character string the resembles the right hand side of a 
#' simple linear regression model in R. For example form = "Group1 + Group2". 
#'
#' @param  data: The analysis data we will use for the pairwise comparisons. 
#'
#' @param  Metabolites: A list of metabolites that we want to be analyzed. 
#' 
#' @return The overall F-test pvalue, and the estimate and pvalue for each pairwise comparison.
#' 
#' @import emmeans
#' 
#' @export
#' 
#' 




metabolite_pairwise <- function(form,data,metabolites){
  
  # Create progress bar
  pb = txtProgressBar(min = 0, max = length(metabolites), initial = 0) 
  
  for (i in 1:length(metabolites)) {
    
    # update progress bar
    setTxtProgressBar(pb,i)
    
    # Get gene
    outcome <- data[,metabolites[i]]
    
    # Define mode
    model <- as.formula(
      paste("outcome", 
            paste(form,collapse = " + "),
            sep = "~"
      )
    )
    
    # Run model
    mod <- lm(model, data = data)
    
    # run anova
    tryCatch(
      {anov <- anova(mod,lm(outcome~1))
      },
     warning = function(c){
       print(paste0(c," for CHEM_ID ",metabolites[i] ))
     },
     finally = {anov <- anova(mod,lm(outcome~1))}
    )
    
    # Get the F statistic pvalue
    overall <- anov$`Pr(>F)`[2]
    
    
    # Create model for the pairwise comparisons
    means_mod <- as.formula(
      paste("pairwise", 
            paste(form,collapse = " + "),
            sep = "~"
      )
    )
    
    # Find pairwise comparisons
    pairs <- emmeans::emmeans(mod, means_mod,adjust='none')$contrasts
    
    
    # Create results table
    if(i ==1 ){
      cols <- c("Metabolite","Overall_pval",paste0(summary(pairs)[,1],"_ESTS"),paste0(summary(pairs)[,1],"_PVALS"))
      
      results <- data.frame(matrix(nrow = length(metabolites),ncol = length(cols)))
      
      colnames(results) <- cols
    }
    
    # Put results in the results data frame
    results[i,"Metabolite"] = metabolites[i]
    
    results[i,"Overall_pval"] = overall
    
    results[i,paste0(summary(pairs)[,1],"_ESTS")] <- summary(pairs)[,"estimate"]
    
    results[i,paste0(summary(pairs)[,1],"_PVALS")] <- summary(pairs)[,"p.value"]
    
    # Order column names
    new_cols <- c("Metabolite","Overall_pval",names(results)[-1:-2][order(names(results)[-1:-2])])
    
    results <- results[,new_cols]
    
    
  }
  close(pb)
  return(results)
}
