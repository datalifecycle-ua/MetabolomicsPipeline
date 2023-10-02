## ----sigSubpathsWithinSuperPath----
for(i in 1:length(strata)){
  
  # 1. Read results for the strata
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
  
  
    # 4. Formulate the Super-pathway results table.
  superPath <- table_data %>% #<4>
                  filter(model != "None") %>% #<4>
                  mutate(pval = eval(cases)) %>% #<4>
                  select(sub_pathway, pval) %>% #<4>
                  right_join(chem_data %>% select(SUPER_PATHWAY, SUB_PATHWAY) %>% distinct(), #<4>
                            by = c("sub_pathway" = "SUB_PATHWAY"))  %>% #<4>
                  filter(!is.na(sub_pathway)) %>% #<4>
                  mutate(sig = ifelse(is.na(pval), 0, 1)) %>% #<4>
                  group_by(SUPER_PATHWAY) %>% #<4>
                  summarise(percent_significant = round(mean(sig) * 100, 2), #<4>
                            number_significant = sum(sig), #<4>
                            pathway_count = n()) %>% #<4>
                  ungroup() %>% arrange(-percent_significant) %>%
                  transmute(SUPER_PATHWAY, percent_significant = paste0(number_significant, " / ", pathway_count,  #<4>
                                                                        " (", percent_significant, "%)")) %>% #<4>
                  knitr::kable(format = 'pipe',col.names = c("Super Pathway", "Percent Significant"), #<4>
               caption = paste0("Proportion of significant <br> subpathways within super-pathways (", #<4>
                                                   strata[i],").") ) #<4>
              
  
  
  # 5. Display table
  print(superPath) #<5>
  cat("\n") #<5>
  
  # 6. Save table
  save_kable(superPath, file = paste0(tablePath,"/SigSuperPathwayPecentages_",strata[i],".png")) #<6>
}     
