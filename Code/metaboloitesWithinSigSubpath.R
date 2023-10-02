## ----metaboloitesWithinSigSubpath----
# Create for loop
for(i in 1:length(strata)){
  # 1. Read in table data and results from overall analysis. 
  path_data <- read.csv(results_path[i], check.names = F) #<1>


 
  # 2. Structure levels
  if(sum(grepl("interaction",names(path_data)))==0){ #<2>
    levs <- c("Single","None") #<2>
    
    cases <- expression( case_when(model == "Single" ~ single_fisher)) #<2>
  }
  if(!sum(grepl("interaction",names(path_data)))==0){ #<2>
    levs =  c("Interaction", "Parallel", "Single", "None") #<2>
    
    cases <- expression( case_when(model == "Interaction" ~ interaction_fisher, #<2>
                                        model == "Parallel" ~ parallel_fisher, #<2>
                                        model == "Single" ~ single_fisher)) #<2>
  }


   # 3. Create table data
  table_data <- path_data %>% #<3>
    mutate(model = factor(model, levels = levs)) %>% #<3>
    select(sub_pathway, ends_with("_fisher"), model) %>% #<3>
    distinct() #<3>



  # 4. Create top pathways table 
  top_pathway_table <- table_data %>% # <4>
    filter(model != "None") %>% # <4>
    mutate(pval = eval(cases)) %>% # <4>
    select(sub_pathway, model, pval) %>% # <4>
    arrange(model, pval) %>% # <4>
    filter(model == mod) %>% # <4>
    group_by(model) %>% slice(1:2) %>% ungroup() %>% # <4>
    left_join(select(path_data, sub_pathway, chem_name, ends_with("_pval")), # <4>
              by = "sub_pathway") %>% # <4>
    distinct() # <4>


  # 5. list the top two subpathways. 
  pathways <- unique(top_pathway_table$sub_pathway) #<5>
  
    # Create pathways for look
    for(j in 1:length(pathways)){
      
      
      # 6. Create table for specific model type
      top_paths_table <- top_pathway_table %>% #<6>
            filter(sub_pathway == pathways[j]) %>% #<6>
            select(chem_name, all_of(paste0(tolower(mod),"_pval"))) %>% #<6>
            arrange(all_of(paste0(tolower(mod),"p_val"))) %>% #<6>
            knitr::kable(format = "pipe",col.names = c("Chemical Name","P-value"), caption = paste0("Metabolites Within Pathways of Significant ",mod, " Model(", #<6>
                                                   strata[i],")."))  #<6>

    
      # 7. Display table 
      print(top_paths_table) #<7>
      
      
      # 8. Save table 
    save_kable(top_pathway_table, 
                  file = paste0(tablePath,"/Metabolite_model_pvalues_",gsub("[^A-Za-z0-9 ]","",pathways[j]),"_",strata[i],".png")) #<8>
      
    }
}  
