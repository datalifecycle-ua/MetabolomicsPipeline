## ----medianStandardization----
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
