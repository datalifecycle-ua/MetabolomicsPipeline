## ----sigMetabolitesByModelType----
for (i in 1:length(strata)) {
  
    # 1. Read in results data
    path_data <- read.csv(results_path[i], check.names = F) #<1>
    
    
    # 2. Structure levels
    if(sum(grepl("interaction",names(path_data)))==0){ #<2>
      levs <- c("Single","None") #<2>
    }
    if(!sum(grepl("interaction",names(path_data)))==0){ #<2>
      levs =  c("Interaction", "Parallel", "Single", "None") #<2>
    } #<2>
    
    
    # 3. Create table data
    table_data <- path_data %>% #<3>
      mutate(model = factor(model, levels = levs)) %>%  #<3>
      select(sub_pathway, ends_with("_fisher"), model) %>% #<3>
      distinct() #<3>
    
    
    # 4. Create table
     sig_subPaths <- count(table_data, model) %>% #<4>
                        arrange((model)) %>% #<4>
                        kable(format = "pipe", col.names = c("Model Type","Count"),
                              caption = paste0("Sigificant Pathways by Model (", #<4>
                                                   strata[i],")."))
     
     
     # 5. Display table
     print(sig_subPaths) #<5>
    
     
     # 6. Save table
     save_kable(sig_subPaths,file=paste0( #<6>
      tablePath ,"/NumberOfSigPathwaysByModelType_",strata[i],".png")) #<6>
}     
