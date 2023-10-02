## ----mergeAnalysisData----
# 2. Create analysis data
analysis_data <- meta_data %>% #<2>
  select(all_of(metadata_variables)) %>% #<2>
  left_join(peak_data_log,"PARENT_SAMPLE_NAME") #<2>
