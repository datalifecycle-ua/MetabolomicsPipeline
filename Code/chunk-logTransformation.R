## ----logTransformation----
## # 1. Log transform all of the values
## peak_data_log <- peak_data_imputed %>% #<1>
##   mutate_if(is.numeric, log) # <1>
