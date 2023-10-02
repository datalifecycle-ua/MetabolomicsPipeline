## ----filterTop50Metabolites----
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
