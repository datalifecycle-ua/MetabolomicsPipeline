subpathway_analysis <- function(data, chem_data, treat_var, block_var){
  
  if(!is.null(block_var)){
    
    # Create results dataframe
    path_data <- data.frame(CHEM_ID = chem_data$CHEM_ID,
                              sub_pathway = chem_data$SUB_PATHWAY,
                             chem_name = chem_data$CHEMICAL_NAME,
                             interaction_pval = NA,
                             interaction_fisher = NA,
                             parallel_pval = NA,
                             parallel_fisher = NA,
                             single_pval = NA,
                             single_fisher = NA,
                             model = NA) %>%
      mutate(sub_pathway = ifelse(is.na(sub_pathway), "Unknown Metabolite", sub_pathway))
    
    
    #For each pathways run the models for each metabolite within the pathway.
    for(pathway in unique(path_data$sub_pathway)){
      var_names <- path_data$CHEM_ID[which(path_data$sub_pathway == pathway)]
      
      
      # Calculate p values for each linear model
      for(var in var_names){
        outcome <- data[[as.character(var)]]
        
        if(is.null(outcome)) next
        
        # This Generalizes Model formulation
        inter_formula = as.formula(
          paste("outcome", 
                paste(c(treat_var,block_var,paste0(treat_var,"*",block_var)),
                      collapse = " + "),
                sep = "~"
          )
        )
        
        int_mod <- lm(inter_formula, data = data)
        path_data$interaction_pval[which(path_data$CHEM_ID == var)] <- anova(int_mod)[paste0(treat_var,":",block_var),5]
        
        # Parallel model, this generalizes the model formulation
        par_formula = as.formula(
          paste("outcome", 
                paste(c(treat_var, block_var),collapse = " + "),
                sep = "~"
          )
        )
        par_mod <- lm(par_formula, data = data)
        path_data$parallel_pval[which(path_data$CHEM_ID == var)] <- anova(par_mod)[block_var,5]
        
        # Treat model
        treat_formula <- as.formula(
          paste("outcome", 
                paste(c(treat_var),collapse = " + "),
                sep = "~"
          )
        )
        treat_mod <- lm(treat_formula, data = data)
        path_data$single_pval[which(path_data$CHEM_ID == var)] <- anova(treat_mod)[treat_var,5]
      }
      
      # Compute combined fisher probabilities for each model.
      if(length(var_names) == 1){
        path_data$interaction_fisher[which(path_data$CHEM_ID == var)] <- anova(int_mod)[paste0(treat_var,":",block_var),5]
        path_data$parallel_fisher[which(path_data$CHEM_ID == var)] <- anova(par_mod)[block_var,5]
        path_data$single_fisher[which(path_data$CHEM_ID == var)] <- anova(treat_mod)[treat_var,5]
        
        path_data$model[which(path_data$sub_pathway == pathway)] <- case_when(anova(int_mod)[paste0(treat_var,":",block_var),5] < 0.05 ~ "Interaction",
                                                                                anova(par_mod)[2,5] < 0.05 ~ "Parallel", 
                                                                                anova(treat_mod)[1,5] < 0.05 ~ "Single",
                                                                                TRUE ~ "None")
      }else{
        interaction_fisher <- -2*sum(log(path_data$interaction_pval[which(path_data$sub_pathway == pathway)]), na.rm = T)
        interaction_p <- round(pchisq(interaction_fisher, 2*length(var_names), lower.tail = F), 2)
        path_data$interaction_fisher[which(path_data$sub_pathway == pathway)] <- interaction_p
        
        parallel_fisher <- -2*sum(log(path_data$parallel_pval[which(path_data$sub_pathway == pathway)]), na.rm = T)
        parallel_p <- round(pchisq(parallel_fisher, 2*length(var_names), lower.tail = F), 2)
        path_data$parallel_fisher[which(path_data$sub_pathway == pathway)] <- parallel_p
        
        single_fisher <- -2*sum(log(path_data$single_pval[which(path_data$sub_pathway == pathway)]), na.rm = T)
        single_p <- round(pchisq(single_fisher, 2*length(var_names), lower.tail = F), 2)
        path_data$single_fisher[which(path_data$sub_pathway == pathway)] <- single_p
        
        path_data$model[which(path_data$sub_pathway == pathway)] <- case_when(interaction_p < 0.05 ~ "Interaction",
                                                                                parallel_p < 0.05 ~ "Parallel", 
                                                                                single_p < 0.05 ~ "Single",
                                                                                TRUE ~ "None")
        
      }
      
    }
  }
  
  # Case when treatment only = True
  if(is.null(block_var)==T){
    
    # Create path data dataframe. 
    path_data <- data.frame(sub_pathway = chem_data$SUB_PATHWAY,
                            chem_name = chem_data$CHEMICAL_NAME,
                            single_pval = NA,
                            single_fisher = NA,
                            model = NA) %>%
      mutate(sub_pathway = ifelse(is.na(sub_pathway), "Unknown Metabolite", sub_pathway))
    
    for(pathway in unique(path_data$sub_pathway)){
      var_names <- path_data$CHEM_ID[which(path_data$sub_pathway == pathway)]
      
      
      # Calculate p values for each linear model
      for(var in var_names){
        outcome <- data[[as.character(var)]]
        
        if(is.null(outcome)) next
        
        
        # Treat model
        treat_formula <- as.formula(
          paste("outcome", 
                paste(c(treat_var),collapse = " + "),
                sep = "~"
          )
        )
        treat_mod <- lm(treat_formula, data = data)
        path_data$single_pval[which(path_data$CHEM_ID == var)] <- anova(treat_mod)[treat_var,5]
      }
      
      if(length(var_names) == 1){
        path_data$single_fisher[which(path_data$CHEM_ID== var)] <- anova(treat_mod)[treat_var,5]
        
        path_data$model[which(path_data$sub_pathway == pathway)] <- case_when(anova(treat_mod)[1,5] < 0.05 ~ "Single",
                                                                                TRUE ~ "None")
      }else{
        
        single_fisher <- -2*sum(log(path_data$single_pval[which(path_data$sub_pathway == pathway)]), na.rm = T)
        single_p <- round(pchisq(single_fisher, 2*length(var_names), lower.tail = F), 2)
        path_data$single_fisher[which(path_data$sub_pathway == pathway)] <- single_p
        
        path_data$model[which(path_data$sub_pathway == pathway)] <- case_when(single_p < 0.05 ~ "Single",
                                                                                TRUE ~ "None")
        
      }
      
    }
  }
  
  # Return path data
  return(path_data)
}
