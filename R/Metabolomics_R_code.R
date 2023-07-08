ancova_function <- function(.data, chem_data, meta_analysis_variables,
                            additive_vars=meta_analysis_variables,
                            treatment.only =F,
                            interaction_vars = NULL,
                            treat_var =NULL,...){
  
  # Filter for any pairwise compairisons
  .data = .data %>%
    filter(...)
  
  if(treatment.only==F){
    chem_sub_pathway <- data.frame(chem_name = names(.data)[-which(names(.data) %in% meta_analysis_variables)]) %>%
      left_join(chem_data, by = c("chem_name" = "CHEMICAL_NAME"))
    
    pathh_data <- data.frame(sub_pathway = chem_sub_pathway$SUB_PATHWAY,
                             chem_name = chem_sub_pathway$chem_name,
                             interaction_pval = NA,
                             interaction_fisher = NA,
                             parallel_pval = NA,
                             parallel_fisher = NA,
                             single_pval = NA,
                             single_fisher = NA,
                             model = NA) %>%
      mutate(sub_pathway = ifelse(is.na(sub_pathway), "Unknown Metabolite", sub_pathway))
    
    for(pathway in unique(pathh_data$sub_pathway)){
      var_names <- pathh_data$chem_name[which(pathh_data$sub_pathway == pathway)]
      
      
      # Calculate p values for each linear model
      for(var in var_names){
        outcome <- .data[[var]]
        
        if(is.null(outcome)) next
        
        # This Generalizes Model formulation
        inter_formula = as.formula(
          paste("outcome", 
                paste(c(additive_vars,interaction_vars),collapse = " + "),
                sep = "~"
          )
        )
        
        int_mod <- lm(inter_formula, data = .data)
        pathh_data$interaction_pval[which(pathh_data$chem_name == var)] <- anova(int_mod)[3,5]
        
        # Parallel model, this generalizes the model formulation
        par_formula = as.formula(
          paste("outcome", 
                paste(c(additive_vars),collapse = " + "),
                sep = "~"
          )
        )
        par_mod <- lm(par_formula, data = .data)
        pathh_data$parallel_pval[which(pathh_data$chem_name == var)] <- anova(par_mod)[2,5]
        
        # Treat model
        treat_formula <- as.formula(
          paste("outcome", 
                paste(c(treat_var),collapse = " + "),
                sep = "~"
          )
        )
        treat_mod <- lm(treat_formula, data = .data)
        pathh_data$single_pval[which(pathh_data$chem_name == var)] <- anova(treat_mod)[1,5]
      }
      
      if(length(var_names) == 1){
        pathh_data$interaction_fisher[which(pathh_data$chem_name == var)] <- anova(int_mod)[3,5]
        pathh_data$parallel_fisher[which(pathh_data$chem_name == var)] <- anova(par_mod)[2,5]
        pathh_data$single_fisher[which(pathh_data$chem_name == var)] <- anova(treat_mod)[1,5]
        
        pathh_data$model[which(pathh_data$sub_pathway == pathway)] <- case_when(anova(int_mod)[3,5] < 0.05 ~ "Interaction",
                                                                                anova(par_mod)[2,5] < 0.05 ~ "Parallel", 
                                                                                anova(treat_mod)[1,5] < 0.05 ~ "Single",
                                                                                TRUE ~ "None")
      }else{
        interaction_fisher <- -2*sum(log(pathh_data$interaction_pval[which(pathh_data$sub_pathway == pathway)]), na.rm = T)
        interaction_p <- round(pchisq(interaction_fisher, 2*length(var_names), lower.tail = F), 2)
        pathh_data$interaction_fisher[which(pathh_data$sub_pathway == pathway)] <- interaction_p
        
        parallel_fisher <- -2*sum(log(pathh_data$parallel_pval[which(pathh_data$sub_pathway == pathway)]), na.rm = T)
        parallel_p <- round(pchisq(parallel_fisher, 2*length(var_names), lower.tail = F), 2)
        pathh_data$parallel_fisher[which(pathh_data$sub_pathway == pathway)] <- parallel_p
        
        single_fisher <- -2*sum(log(pathh_data$single_pval[which(pathh_data$sub_pathway == pathway)]), na.rm = T)
        single_p <- round(pchisq(single_fisher, 2*length(var_names), lower.tail = F), 2)
        pathh_data$single_fisher[which(pathh_data$sub_pathway == pathway)] <- single_p
        
        pathh_data$model[which(pathh_data$sub_pathway == pathway)] <- case_when(interaction_p < 0.05 ~ "Interaction",
                                                                                parallel_p < 0.05 ~ "Parallel", 
                                                                                single_p < 0.05 ~ "Single",
                                                                                TRUE ~ "None")
        
      }
      
    }
  }
  
  # Case when treatment only = True
  if(treatment.only==T){
    
    
    chem_sub_pathway <- data.frame(chem_name = names(.data)[-which(names(.data) %in% meta_analysis_variables)]) %>%
      left_join(chem_data, by = c("chem_name" = "CHEMICAL_NAME"))
    
    pathh_data <- data.frame(sub_pathway = chem_sub_pathway$SUB_PATHWAY,
                             chem_name = chem_sub_pathway$chem_name,
                             single_pval = NA,
                             single_fisher = NA,
                             model = NA) %>%
      mutate(sub_pathway = ifelse(is.na(sub_pathway), "Unknown Metabolite", sub_pathway))
    
    for(pathway in unique(pathh_data$sub_pathway)){
      var_names <- pathh_data$chem_name[which(pathh_data$sub_pathway == pathway)]
      
      
      # Calculate p values for each linear model
      for(var in var_names){
        outcome <- .data[[var]]
        
        if(is.null(outcome)) next
        
        
        # Treat model
        treat_formula <- as.formula(
          paste("outcome", 
                paste(c(treat_var),collapse = " + "),
                sep = "~"
          )
        )
        treat_mod <- lm(treat_formula, data = .data)
        pathh_data$single_pval[which(pathh_data$chem_name == var)] <- anova(treat_mod)[1,5]
      }
      
      if(length(var_names) == 1){
        pathh_data$single_fisher[which(pathh_data$chem_name == var)] <- anova(treat_mod)[1,5]
        
        pathh_data$model[which(pathh_data$sub_pathway == pathway)] <- case_when(anova(treat_mod)[1,5] < 0.05 ~ "Single",
                                                                                TRUE ~ "None")
      }else{
        
        single_fisher <- -2*sum(log(pathh_data$single_pval[which(pathh_data$sub_pathway == pathway)]), na.rm = T)
        single_p <- round(pchisq(single_fisher, 2*length(var_names), lower.tail = F), 2)
        pathh_data$single_fisher[which(pathh_data$sub_pathway == pathway)] <- single_p
        
        pathh_data$model[which(pathh_data$sub_pathway == pathway)] <- case_when(single_p < 0.05 ~ "Single",
                                                                                TRUE ~ "None")
        
      }
      
    }
  }
  
  # Return path data
  return(pathh_data)
}


# This function creates heatmap data
create_heatmap <- function(.peak_data_scale, .chem_data, .meta_data ,meta_analysis_variables,...){
  
  heatmap_data =.peak_data_scale %>%
    pivot_longer(cols = -PARENT_SAMPLE_NAME,
                 names_to = "CHEM_ID") %>%
    left_join(select(.chem_data, CHEM_ID, CHEMICAL_NAME) %>%
                mutate(CHEM_ID= as.character(CHEM_ID)),
              by = "CHEM_ID") %>%
    left_join(.meta_data[,meta_analysis_variables],
              by = "PARENT_SAMPLE_NAME") %>%
    pivot_wider(id_cols = meta_analysis_variables,
                names_from = CHEMICAL_NAME, values_from = value) %>%
    arrange(...)
  
  heatmap_meta_data <- heatmap_data %>%
    column_to_rownames("PARENT_SAMPLE_NAME") %>%
    select(all_of(meta_analysis_variables[!(meta_analysis_variables %in% "PARENT_SAMPLE_NAME")]))
  
  
  heatmap_data2 <- heatmap_data %>%
    select(-meta_analysis_variables[!(meta_analysis_variables %in% "PARENT_SAMPLE_NAME")]) %>%
    column_to_rownames("PARENT_SAMPLE_NAME") %>%
    as.matrix() %>% t()
  
  return(list(heatmap_data_all = heatmap_data, heatmap_meta_data = heatmap_meta_data, heatmap_data_vals = heatmap_data2))
  
}






# Filtering heatmap
filter_heatmap <- function(.heatmap_data,drop_vars,...){
  
  
  
  # Heatmap_all
  heatmap_data_all = .heatmap_data$heatmap_data_all %>%
    filter(...) %>%
    select(-drop_vars)
  
  heatmap_meta_data = .heatmap_data$heatmap_meta_data %>%
    filter(...) %>%
    select(-drop_vars)
  
  heatmap_data_vals = .heatmap_data$heatmap_data_vals[,row.names(heatmap_meta_data)]
  
  return(list(heatmap_data_all=heatmap_data_all,
              heatmap_meta_data= heatmap_meta_data,
              heatmap_data_vals= heatmap_data_vals))
  
}



# Create table for pairwise compairisons
sub_path_pairwise_comparison   <- function(path, model){
  
  
  
  # Get pairwise files 
  pairwise_files <- list.files(path)
  
  
  
  # Get pairs
  pairs <- sub("PairwiseCompairison_","",sub(".csv","",pairwise_files))
  
  
  # Get model name
  dat <- read.csv(paste0(path,pairwise_files[1]), check.names = F)
  
  mod = names(dat %>%
          select(sub_pathway,{{model}}))[2]
  
  # Create initial table
  strat1_path_comparisons <- select(read.csv(paste0(path,pairwise_files[1]), check.names = F), sub_pathway, comp = {{model}}) %>%
    distinct() 
  
  colnames(strat1_path_comparisons)[which(colnames(strat1_path_comparisons)=="comp")] = paste0(pairs[1],"_",mod)
  
  
  # Run for loop
  for(i in 2:length(pairwise_files)){
    
    strat1_path_comparisons <- strat1_path_comparisons %>%
      left_join(select(read.csv(paste0(path,pairwise_files[i]), check.names = F), sub_pathway, comp = interaction_fisher) %>%
                  distinct(), 
                by = "sub_pathway")
    
    colnames(strat1_path_comparisons)[which(colnames(strat1_path_comparisons)=="comp")] = paste0(pairs[i], "_",mod)
    
  }
  
  
  return(strat1_path_comparisons)
}


# Table for all compairisons
sub_path_pairwise_comparison_all   <- function(path, ...){
  
  # Get pairwise files 
  pairwise_files <- list.files(path = path)
  
  # Get pairs
  pairs <- sub("PairwiseCompairison_","",sub(".csv","",pairwise_files))
  
  
  # Create first data compairison table
  comptable <- read.csv(paste0(path,pairwise_files[1])) %>% 
    select(sub_pathway, ends_with("_fisher")) %>% distinct() %>%
    mutate(model = pairs[1])
  
  for (i in 2:length(pairs)) {
    comptable <- comptable %>%
      rbind(
        read.csv(paste0(path,pairwise_files[i])) %>% 
          select(sub_pathway, ends_with("_fisher")) %>% distinct() %>%
          mutate(model = pairs[i]))
  }
  
  

  
  comptable <- comptable %>%
    select(model, sub_pathway, everything()) %>%
    arrange(model, sub_pathway) %>%
    filter(...) %>%
    flextable() %>%
    set_header_labels(model = "Pairwise Model", sub_pathway="Sub Pathway", interaction_fisher="Interaction", parallel_fisher="Parallel", single_fisher = "Single" ) %>%
    add_header_row(values = c(NA,"Fisher Combined Probabilities"), colwidths = c(2,3)) %>%
    theme_vanilla() %>%
    set_table_properties(layout = "autofit")
  
  return(comptable)
  
}







subpathway_boxplots <- function(analysis_data,chem_data,subpathway,X, groupBy,...){
  
  analysis_data %>%
    filter(...) %>%
    select(group = {{groupBy}}, X = {{X}}, 
           chem_data$CHEMICAL_NAME[which(chem_data$SUB_PATHWAY == subpathway)]) %>%
    pivot_longer(cols = -c(group, X)) %>%
    mutate(X = factor(X)) %>%
    ggplot(aes(x = X, y = value, color = group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position=position_jitterdodge(jitter.width = 0.2)) +
    facet_wrap(~name) +
    theme_bw()
  
  
}



subpathway_lineplots <- function(analysis_data,chem_data,subpathway,X, groupBy,...){
  
  analysis_data %>%
    select(treat = {{groupBy}}, X = {{X}},
           intersect(chem_data$CHEMICAL_NAME[which(chem_data$SUB_PATHWAY == 
                                                     subpathway)],
                     names(analysis_data))) %>%
    pivot_longer(cols = -c(treat, X)) %>%
    ggplot(aes(x = X, y = value, color = treat)) +
    geom_jitter() +
    geom_smooth(method = "lm") +
    facet_wrap( ~ name) +
    theme_bw()
  
}






