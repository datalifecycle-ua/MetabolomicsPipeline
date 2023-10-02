## ----beforeMedianStandardization----
## # 1. Select the first five metabolites for the box plot.
## metabolites <- names(peak_data)[2:6] #<1>
## 
## 
## # 2. Create the boxplot data
## plot_data <- peak_data %>% #<2>
##   reshape2::melt(id.vars="PARENT_SAMPLE_NAME") %>% #<2>
##   filter(variable %in% metabolites) #<2>
## 
## 
## # 3. Plot the box plots for each of the five metabolites.
## ggplot(plot_data,aes(x=variable,y=value)) + #<3>
##   geom_boxplot() + #<3>
##   theme_bw() #<3>
