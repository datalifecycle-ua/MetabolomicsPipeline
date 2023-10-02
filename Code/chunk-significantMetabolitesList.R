## ----significantMetabolitesList----
# Create for loop
for (i in 1:length(strata)){
  
  
  # 1. Read Path data
  path_data <- read.csv(results_path[i]) #<1>
  
  
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

  

  # 3. Create pathway table
  pathway_table <- path_data %>% # <3>
    mutate(model = factor(model, levels = levs)) %>% # <3>
    select(sub_pathway, ends_with("_fisher"), model) %>% # <3>
    distinct() %>%  # <3>
    filter(model != "None") %>% # <3>
    mutate(pval = round(eval(cases),3)) %>% # <3>
    select(sub_pathway, model, pval) %>% # <3>
    arrange(model, pval) # <3>
  
  
  # 4. Format table
  path_tab = pathway_table %>% #<4>
    knitr::kable(format = 'pipe',col.names = c("Subpathway", "Model Type","P-value"), #<4>
                 caption = paste0("Significant Subpathways (", #<4>
                                                   strata[i],").") ) #<4>
  
  
  # 5. Display table
    print(path_tab) #<5>

  
  
  # 6. Save table
  save_kable(path_tab ,file  = paste0(tablePath,"/SigSubpathwayTable_",strata[i],".png")) #<6>
}
