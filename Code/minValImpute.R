## ----minValImpute----
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
