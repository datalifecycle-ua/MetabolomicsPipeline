## ----generateHeatmap2----
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
