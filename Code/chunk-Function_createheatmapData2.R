## ----Function_createheatmapData2----
# 4. Create heatmap data
heatmap_data <- create_heatmap_Data(analysis_data_top50, #<1>
                                    heatmap_variables = heatmap_variables , #<1>
                               GROUP_NAME, desc(TIME1)) # #<2>
