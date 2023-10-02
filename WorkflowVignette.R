
# Set global chunk obptions
knitr::opts_chunk$set(fig.width=12, fig.height=12, warning = F)

# Data 
library(tidyverse)
library(openxlsx)
library(kableExtra)
library(table1)
library(flextable)
library(reshape2)


# Heat Map
library(pheatmap)
library(RColorBrewer)
library(rstatix)
library(plotly)

# PCA
library(factoextra)
library(FactoMineR)


## Combining plots
library(ggpubr)
library(grid)
library(gridExtra)

## Image package
library(magick)

## Pairwise comparisons
library(emmeans)

# Metabolomics Pipeline source files
files.sources = list.files("../R")
sapply(paste0("../R/",files.sources), source)




img = image_read("../Workflow.png")

image_trim(img)



## # 1 Provide a path to Metabolon .xlsx file.
## metabolon_path <- "../Data/Metabolon/UNAZ-0501-22VW_ DATA TABLES.xlsx" #<1>



## # 2. Load raw peak data
## peak_data <- read.xlsx(metabolon_path, sheet = "Peak Area Data") # <2>



## # 1. Select the first five metabolites for the box plot.
## metabolites <- names(peak_data)[2:6] #<1>
## 
## 
## # 2. Create the boxplot data
## plot_data <- peak_data %>% #<2>
##   reshape2::melt(id.vars="PARENT_SAMPLE_NAME") %>% #<2>
##   filter(variable %in% metabolites) #<2>
## 
## 
## # 3. Plot the box plots for each of the five metabolites.
## ggplot(plot_data,aes(x=variable,y=value)) + #<3>
##   geom_boxplot() + #<3>
##   theme_bw() #<3>



## # 1. Initialize a new peak_data_med matrix
## peak_data_norm <-  peak_data #<1>
## 
## 
## # 2. Create a matrix containing the median value for each metabolite
## peak_data_med <- peak_data_norm %>% #<2>
##   select(-PARENT_SAMPLE_NAME) %>% #<2>
##   summarise_all(median, na.rm = T) #<2>
## 
## 
## # 3. Divide each value for each metabolite by the median value of that metabolite
## for(i in colnames(peak_data_med)){  #<3>
##   peak_data_norm[,i] <- peak_data_norm[,i]/peak_data_med[,i]  #<3>
## }  #<3>
## 



## # 1. Select the first five metabolites for the box plot.
## metabolites <- names(peak_data)[2:6] # <1>
## 
## 
## # 2. Create the boxplot data
## plot_data <- peak_data_norm %>%  # <2>
##   reshape2::melt(id.vars="PARENT_SAMPLE_NAME") %>% # <2>
##   filter(variable %in% metabolites) # <2>
## 
## 
## # 3. Plot the box plots for each of the five metabolites.
## ggplot(plot_data,aes(x=variable,y=value)) + # <3>
##   geom_boxplot() + # <3>
##   theme_bw() # <3>



## # 1. Initialize the new peak_data_imputed matrix
## peak_data_imputed <- peak_data_norm # <1>
## 
## 
## # 2. Find the minimum value for each metabolite and compute 1/5 of that value
## peak_data_mins <- peak_data_imputed %>%  # <2>
##   select(-PARENT_SAMPLE_NAME) %>%  # <2>
##   summarise_all(min, na.rm = T)  # <2>
## 
## 
## # 3. Impute the value
## for(i in colnames(peak_data_mins)){  # <3>
##   if(length(peak_data_imputed[,i][is.na(peak_data_imputed[,i])]) > 0){ # <3>
##     peak_data_imputed[which(is.na(peak_data_imputed[,i])),i] <- as.numeric(peak_data_mins[i]) # <3>
##   } # <3>
## } # <3>



## # 1. Log transform all of the values
## peak_data_log <- peak_data_imputed %>% #<1>
##   mutate_if(is.numeric, log) # <1>



## # 1. Save log transformed data.
## write.csv(peak_data_log, file = "../Data/Processed/Log_transformed_data.csv", #<1>
##           row.names = F) #<1>



# 1. Metabolon excel file 
met_excel <- "../Data/Metabolon/UNAZ-0501-22VW_ DATA TABLES.xlsx" # <1>



# 2. Read in sample metadata
meta_data <- read.xlsx(met_excel,sheet = "Sample Meta Data") # <2>


# 3. Read Log transformed data
peak_data_log <- read.xlsx(met_excel,sheet = "Log Transformed Data") #<3>



# If you have additional meta data set additional_meta=T
additional_meta=T #<1>

# 1. Provide path to additional metadata
additional_metaPath <- "../Data/Metabolon/AdditionalVars.xlsx" #<2>



# 2. Merge additional vars to the meta data
if(!is.null(additional_meta)){ #<3>
  
  # Read in additional meta data #<3>
  meta_data_additional <- read.xlsx(additional_metaPath) #<3>
  
  meta_data <- meta_data_additional %>% #<3>
    left_join(meta_data,"PARENT_SAMPLE_NAME") #<3>
} #<3>



# 1. Select metadata variables
metadata_variables <- c("PARENT_SAMPLE_NAME", # <1>
                        "GROUP_NAME", #<1>
                        "TIME1", # <1>
                        "Gender") # <1>



# 2. Create analysis data
analysis_data <- meta_data %>% #<2>
  select(all_of(metadata_variables)) %>% #<2>
  left_join(peak_data_log,"PARENT_SAMPLE_NAME") #<2>



# 1. Choose the meta data variable name.
var1 = "GROUP_NAME" #<1>

var2 = "TIME1" # <1>


# 2. Choose a single variable to stratify the table by.
stratified_var = "Gender" # <2>


# 3. Assign the variable a label name
label(analysis_data[,var1]) <- "Treatment Group" #<3>

label(analysis_data[,var2]) <- "Time" #<3>


# 4. Assign a stratified variable a label name for the table. 
label(analysis_data[,stratified_var]) = "Gender" #<4>



# 1. Creates table1
tbl1 <- table1(~ TIME1 + GROUP_NAME| Gender # <1>
               , data = analysis_data) #<2>



# 2. Displays table 1
tbl1 # <2>



# 3. saves table 1 
t1flex(tbl1) %>% #<3>
  save_as_docx(path = paste0("../Outputs/Tables/","table1.docx")) #<3>



# Save analysis data
write.csv(analysis_data,"../Data/Processed/analysis_data.csv", row.names = F) #<1>



# 1. Read analysis dataset
analysis_data_path <- "../Data/Processed/analysis_data.csv" #<1>



out <- "../Outputs/Plots" #<1>



# 1. Create the non-metabolite vector. 
non_metabolite <- c("PARENT_SAMPLE_NAME", #<1>
                             "GROUP_NAME", #<1>
                             "TIME1", #<1>
                             "Gender") #<1>

# 2. Define labels
meta_var = "Gender" # <2>




# 3. Read in analysis data
analysis_data <- read.csv(analysis_data_path, check.names = F) #<3>


# 4. Create PCA data containing only metabolite data
pca_dat <- analysis_data[,!names(analysis_data) %in% non_metabolite] #<4>


# 5. Run PCA of the pca_dat matrix containing only the metabolites.   
res.pca <- PCA(pca_dat, #<5>
               graph = F)#<5>


# 6. Create figure 
fviz_pca_ind(res.pca, #<6>
             label = "none",#<6>
             habillage = as.factor(analysis_data[,meta_var])) #<6>


# 7. Save figure
ggsave(paste0(out,"/PCA.pdf"), width = 10, height = 10)#<7>



# 1. Define grouping variables
heatmap_variables <- c("GROUP_NAME") #<1>



# 2. Find the top 50 metabolites
select_variables <- analysis_data %>% # <2>
  select(-all_of(c("PARENT_SAMPLE_NAME",heatmap_variables))) %>%# <2>
  summarise(across(everything(),\(x) mean(x,na.rm = T))) %>%# <2>
  pivot_longer(cols = everything()) %>%# <2>
  arrange(desc(value)) %>%# <2>
  slice(c(1:50))# <2>


# 3. Filter to the top 50 metabolites
analysis_data_top50 <- analysis_data %>% #<3>
  select(all_of(c("PARENT_SAMPLE_NAME", heatmap_variables)), select_variables$name) # <3>



#  Create heatmap data
heatmap_data <- create_heatmap_Data(analysis_data_top50, #<1>
                                    heatmap_variables = heatmap_variables , #<1>
                               GROUP_NAME) #<2>



# Heat map colors 
palette <- colorRampPalette(rev(brewer.pal(10, "RdBu")))(256) # <1>


# Values for heatmap
vals = heatmap_data$heatmap_data_vals #<2>

meta = heatmap_data$heatmap_variables #<2>


# Create heatmap 
pheatmap(vals, cluster_cols = F, cluster_rows = F, color = palette, #<3>
         annotation_col = meta, show_rownames = F, border_color = NA, show_colnames = F) #<3>

# Save heatmap
pheatmap(vals, cluster_cols = F, cluster_rows = F, color = palette, #<4>
         annotation_col = meta, show_rownames = F, show_colnames = F,border_color = NA,filename = paste0(out,"/HeatmapTop50Metabolites.pdf")) #<4>




# 1. Define meta analysis variables
heatmap_variables <- c("GROUP_NAME", #<1>
                       "TIME1") #<1>



# 1. Find the top 50 metabolites
select_variables <- analysis_data %>% #<2>
  select(-all_of(c("PARENT_SAMPLE_NAME",heatmap_variables))) %>% #<2>
  summarise(across(everything(),\(x) mean(x,na.rm = T))) %>% #<2>
  pivot_longer(cols = everything()) %>% #<2>
  arrange(desc(value)) %>% #<2>
  slice(c(1:50)) #<2>


# 2. Filter to the top 50 metabolites
analysis_data_top50 <- analysis_data %>%
  select(all_of(c("PARENT_SAMPLE_NAME", heatmap_variables)), select_variables$name) 



# 4. Create heatmap data
heatmap_data <- create_heatmap_Data(analysis_data_top50, #<1>
                                    heatmap_variables = heatmap_variables , #<1>
                               GROUP_NAME, desc(TIME1)) # #<2>



# Heat map colors 
palette <- colorRampPalette(rev(brewer.pal(10, "RdBu")))(256) # <1>


# Values for heatmap
vals = heatmap_data$heatmap_data_vals #<2>

meta = heatmap_data$heatmap_variables #<2>


# Create heatmap 
pheatmap(vals, cluster_cols = F, cluster_rows = F, color = palette, #<3>
         annotation_col = meta, show_rownames = F, border_color = NA, show_colnames = F) #<3>

# Save heatmap
pheatmap(vals, cluster_cols = F, cluster_rows = F, color = palette, #<4>
         annotation_col = meta, show_rownames = F, show_colnames = F,border_color = NA,filename = paste0(out,"/HeatmapTop50MetabolitesTreatAndBlock.pdf")) #<4>




# 1. Get stratified variable and stratified values
strat_var <- "Gender" #<1>

values <- unique(analysis_data[,strat_var]) #<2>

# 2 Get heatmap Variables
heatmap_variables <- c("GROUP_NAME", #<3>
                       "TIME1") #<3>




# Find the top 50 metabolites
select_variables <- analysis_data %>% #<1>
  select(-all_of(c("PARENT_SAMPLE_NAME",heatmap_variables))) %>% #<1>
  summarise(across(everything(),\(x) mean(x,na.rm = T))) %>% #<1>
  pivot_longer(cols = everything()) %>% #<1>
  arrange(desc(value)) %>% #<1>
  slice(c(1:50)) #<1>


# Create for loop
maps <- list()
for (i in 1:length(values)) {
  
    # Filter to the top 50 metabolites
    analysis_data_top50_i <- analysis_data[analysis_data[,strat_var]==values[i],] %>% #<2>
    select(all_of(c("PARENT_SAMPLE_NAME", heatmap_variables)), select_variables$name) #<2>
    
    
    # Create heatmap data
    heatmap_data <- create_heatmap_Data(analysis_data_top50_i, #<3>
                                      heatmap_variables = heatmap_variables , #<3>
                                 GROUP_NAME, desc(TIME1)) # <4>
    
    # Heat map colors 
    palette <- colorRampPalette(rev(brewer.pal(10, "RdBu")))(256) #<5>
    
    
    # Values for heatmap
    vals = heatmap_data$heatmap_data_vals #<6>
    
    meta = heatmap_data$heatmap_variables #<6>
    
    
    #  Create heatmap
    maps[[values[i]]] = pheatmap(vals, cluster_cols = F, cluster_rows = F, color = palette, #<7>
           annotation_col = meta, show_rownames = F, border_color = NA, show_colnames = F, # <7>
           main = values[i], silent = T)
    
}


# 8. Save heatmaps for all levels
grids <- grid.arrange(grobs= list(maps$Female[[4]], # <8>
                         maps$Male[[4]])) #<8>

ggsave(filename = paste0(out,"/HeatmapTop50MetabolitesTreamentandBlockStratified.pdf"), grids) #<8>




# 1. Path to Metabolon .xlsx file.
met_excel <- "../Data/Metabolon/UNAZ-0501-22VW_ DATA TABLES.xlsx" #<1>

analysis_dat_path <- "../Data/Processed/analysis_data.csv" #<2>




# 1. Analysis Data
analysis_data <- read.csv(analysis_dat_path, check.names = F) #<1>


# 2. Read chemical annotations
chem_data <- read.xlsx(met_excel,sheet = "Chemical Annotation") #<2>



# 1. Run the subpathway function
path_dat <- subpathway_analysis(analysis_data,#<1>
                            chem_data = chem_data,#<1>
                                  block_var = "TIME1",#<1>
                                  treat_var = "GROUP_NAME") #<1>



# 1. Save results
write.csv(path_dat,file = "../Data/Results/NonStratified_Full_subpathway_results.csv", row.names = F)#<1>



# Enter block variable name
blockvar = "TIME1" #<1>

# Enter treatment variable
treatvar = "GROUP_NAME" #<1>


# 1. Name of the variable to stratify by.
stratified_var = "Gender" # <2>

# 3. Find unique values of this variable
uni_vals <- unique(analysis_data[,stratified_var]) #<3>



# Enter Statified Variable path (Do not include / at the end)
statified_data_path = "../Data/Processed" #<1>

# Enter Statified Results path
stratified_results_path = "../Data/Results" #<2>



# For each value perform the four steps listed above.
for(i in uni_vals){
  
  print(paste0("Running the ", i, " strata."))
  
  # 1. Subset analysis data
  strat_data = analysis_data[analysis_data[,stratified_var]==i,]
  
  
  # 2. Run subpathway function on the subsetted data
  strat_path_results <- subpathway_analysis(strat_data,
                                        chem_data = chem_data,
                                  block_var = blockvar,
                                  treat_var = treatvar)
  
  # 3 Save stratified data
  write.csv(strat_data, paste0(statified_data_path,"/analysis_data_",i,".csv"), row.names = F)
  
  # 4 Save results
  write.csv(strat_path_results, paste0(stratified_results_path,"/Stratified_",i,"_subpathway_results.csv"), row.names = F)
  
}



# 1. Read in Results from the analysis step. 
results_path <- c("../Data/Results/Stratified_Male_subpathway_results.csv", #<1>
                  "../Data/Results/Stratified_Female_subpathway_results.csv") # <1>


# 2. Define the names of the strata
strata <- c("Male", "Female") #<2>



# Enter path to save tables
tablePath = "../Outputs/Tables" #<1>


#| tbl-cap: "Significant Metabolites by Model"

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



# Enter Model type (Interaction, Parallel, Single)
mod = "Interaction" #<1>



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



path3 = "../Data/Processed/analysis_data.csv" #<1>



pair_strat_results_path = "../Data/Results/Pairwise_Comparisons" #<1>



# 1. Define the variable in the analysis data that we want to stratify 
strata_var <- "Gender" #<1>


# 2. Define the form variable
form <- "GROUP_NAME*TIME1" #<2>



# Enter Non-metabolite variables
non_metabolite <- c("PARENT_SAMPLE_NAME", "GROUP_NAME", #<3>
                                           "TIME1", "Gender") #<3>

metabolites <- names(analysis_data)[!names(analysis_data) %in% non_metabolite] #<4>



# 1. Read in data from the analysis above
analysis_data <- read.csv(path3, check.names = F) #<1>

# 2. Make the variable names charactors
names(analysis_data) <- as.character(names(analysis_data)) #<2>


stratas <- unique(analysis_data[,strata_var])
for (i in 1:length(stratas)) {
  
  # 3. Filter analysis data
  dat <- analysis_data[analysis_data[,strata_var]==stratas[i],] #<3>
  
  
  # 4. Run the metabolite_pairwise function for the strata
  strat_results <- metabolite_pairwise(form = form, #<4>
                                       data=dat, #<4>
                                       metabolites = metabolites) #<4>
  
  # 5. Save results
  write.csv(strat_results, file = paste0(pair_strat_results_path,"/Pairwise_", #<5>
                                         stratas[i],".csv"), row.names = F) #<5>
  
  }



# Enter pairwise results path
pair_results_path = "../Data/Results/Pairwise_Comparisons/Pairwise_Female.csv" #<1>

# Enter path to metabolome .xlsx file
met_excel <- "../Data/Metabolon/UNAZ-0501-22VW_ DATA TABLES.xlsx" #<2>



# 1. Read in pairwise results
results_pair <-read.csv(pair_results_path, #<1>
                         check.names = F) %>% #<1>
  filter(Overall_pval < 0.05) #<1>

# Read chemical annotation file
chem_data <- read.xlsx(met_excel,sheet = "Chemical Annotation") #<1>

 

# 2. Merge the chemical annotation fill with the results from the pairwise comparisons.
data <- chem_data %>% #<2>
  select(SUB_PATHWAY,CHEMICAL_NAME,CHEM_ID) %>% #<2>
  merge(results_pair, by.x = "CHEM_ID",by.y = "Metabolite") %>% #<2>
  arrange(SUB_PATHWAY)  #<2>


# 3. Produce Heatmap
p = data %>% # <3>
      select(CHEM_ID,SUB_PATHWAY,CHEMICAL_NAME, all_of(names(data)[grepl("EST",names(data))])) %>%# <3>
      melt(id.vars = c("CHEM_ID","SUB_PATHWAY","CHEMICAL_NAME"), variable.name = "Contrast", # <3>
           value.name = "logFoldChange") %>%# <3>
      mutate(Contrast = gsub("_ESTS","",Contrast),# <3>
             logFoldChange =  ifelse(logFoldChange < log(0.5) | logFoldChange > log(2),# <3>
                                     round(logFoldChange,3), NA)) %>% # <3>
      plot_ly(# <3>
        type = "heatmap",# <3>
        x= ~Contrast,# <3>
        y = ~CHEMICAL_NAME,# <3>
        z = ~logFoldChange,# <3>
        text = ~SUB_PATHWAY,# <3>
        hovertemplate = paste("<b>Metabolite: %{y}</b><br><br>",# <3>
                              "Sub-pathway: %{text}<br>",# <3>
                              "Contrast: %{x}<br>",# <3>
                              "Log Fold-Change: %{z}<br>",# <3>
                              "<extra></extra>"),# <3>
        colorbar = list(title ="<b>Log Fold-Change</b>")) %>%# <3>
      layout(title = "<b>Log Fold-Change Heatmap</b>",# <3>
             xaxis = list(title="<b>Contrasts</b>"),# <3>
             yaxis = list(title = "")) # <3>

# 4. Display heatmap
p #<4>



# 1. Read in pairwise results
results_pair <-read.csv(pair_results_path, #<1>
                         check.names = F) %>%  #<1>
  filter(Overall_pval < 0.05)  #<1>

# Read chemical annotation file
chem_data <- read.xlsx(met_excel,sheet = "Chemical Annotation") #<1>

 

# 2. Merge the chemical annotation fill with the results from the pairwise comparisons.
data <- chem_data %>% #<2>
  select(SUB_PATHWAY,CHEMICAL_NAME,CHEM_ID) %>% #<2>
  merge(results_pair, by.x = "CHEM_ID",by.y = "Metabolite") %>% #<2>
  arrange(SUB_PATHWAY)  #<2>


# 3. Produce Heatmap
p = data %>% #<3>
      select(CHEM_ID,SUB_PATHWAY,CHEMICAL_NAME, all_of(names(data)[grepl("PVALS",names(data))])) %>%  #<3>
      melt(id.vars = c("CHEM_ID","SUB_PATHWAY","CHEMICAL_NAME"), variable.name = "Contrast",  #<3>
           value.name = "P_value") %>%  #<3>
      mutate(Contrast = gsub("_PVALS","",Contrast), #<3>
             P_value = round(P_value,3)) %>% #<3>
  arrange(SUB_PATHWAY) %>% #<3>
      plot_ly( #<3>
        type = "heatmap", #<3>
        x= ~Contrast, #<3>
        y = ~CHEMICAL_NAME, #<3>
        z = ~P_value, #<3>
        text = ~SUB_PATHWAY, #<3>
        hovertemplate = paste("<b>Metabolite: %{y}</b><br><br>", #<3>
                              "Sub-pathway: %{text}<br>", #<3>
                              "Contrast: %{x}<br>", #<3>
                              "P-Value: %{z}<br>", #<3>
                              "<extra></extra>"), #<3>
        colorbar = list(title ="<b>P-value</b>")) %>% #<3>
      layout(title = "<b>P-Value Heatmap</b>", #<3>
             xaxis = list(title="<b>Contrasts</b>"), #<3>
             yaxis = list(title = "")) #<3>

# 4. Display heatmap
p



# Enter pairwise results path
analysis = "../Data/Processed/analysis_data.csv" #<1>

# Enter chemical annotation path
met_excel <- "../Data/Metabolon/UNAZ-0501-22VW_ DATA TABLES.xlsx" #<2>



# 1. Specify subpathway of interest. 
subpath = "Gamma-glutamyl Amino Acid" #<1>



plot_save = "../Outputs" #<1>




# 1. Read in analysis data
analysis_data <- read.csv(analysis, check.names = F) #<1>

# 2. Reach in chemical annotation data
# Read chemical annotation file
chem_data <- read.xlsx(met_excel,sheet = "Chemical Annotation") # <2>


# 3. Run the subpathway_boxplots function
box_1 <- subpathway_boxplots(analysis_data =  analysis_data, #<3>
                chem_data = chem_data, #<3>
                subpathway = subpath, #<3>
                X = TIME1, #<3>
                groupBy = GROUP_NAME, #<3>
                # Additional analysis data filtering options. 
                Gender == "Male", GROUP_NAME != "Control") #<3>

# 4. Customize labels for the boxplots
box_1 + #<4>
  labs(x = "Time", y="Scaled Intensity", color="Treatment Group", #<4>
       title = "Metabolite boxplots for the Gamma-glutamyl Amino Acid Subpathway") #<4>

# 5. Save boxplots. You can change the file type for example change pdf to png. 
ggsave(filename = paste0(plot_save,"/Plots/Metabolite_boxplots_",subpath,".pdf")) #<5>





# 1. Read in analysis data
analysis_data <- read.csv(analysis, check.names = F) #<1>

# 2. Reach in chemical annotation data
# Read chemical annotation file
chem_data <- read.xlsx(met_excel,sheet = "Chemical Annotation") #<2>


# 3. Prep data for line plots
analysis_data <- analysis_data %>% #<3>
  filter(Gender=="Male") %>% #<3>
  mutate(TIME1 = case_when(TIME1== "PreSymp" ~ 1, #<3>
                           TIME1 =="Onset" ~ 2, #<3>
                       TIME1 == "End" ~ 3)) #<3>




# 4. Run subpathway_line plots function.
line_plots <- subpathway_lineplots(analysis_data = analysis_data, #<4>
                     chem_data = chem_data, #<4>
                     subpathway = subpath, #<4>
                     X = TIME1, #<4>
                     groupBy = GROUP_NAME) #<4>


# 5. Edit asthetics of plots
line_plots + #<5>
  scale_x_continuous(breaks = c(1, 2, 3), labels = c("PreSymp", "Onset", "End")) +#<5>
  labs(x= "Time",y="Scaled Intensity",color = "Treatment Group")#<5>


# 6. Save plots 
ggsave(filename = paste0(plot_save ,"/LinePlots_",subpath,".pdf")) #<6>


