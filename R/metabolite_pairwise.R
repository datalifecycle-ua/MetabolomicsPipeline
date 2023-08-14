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
    pairs <- emmeans(mod, means_mod,adjust='none')$contrasts
    
    
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
