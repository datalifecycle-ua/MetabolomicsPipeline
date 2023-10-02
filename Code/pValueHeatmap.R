## ----pValueHeatmap----
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
