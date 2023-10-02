## ----merge_additionalVars----
# 2. Merge additional vars to the meta data
if(!is.null(additional_meta)){ #<3>
  
  # Read in additional meta data #<3>
  meta_data_additional <- read.xlsx(additional_metaPath) #<3>
  
  meta_data <- meta_data_additional %>% #<3>
    left_join(meta_data,"PARENT_SAMPLE_NAME") #<3>
} #<3>
