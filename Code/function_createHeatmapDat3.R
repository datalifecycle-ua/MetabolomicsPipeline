## ----function_createHeatmapDat3----
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
