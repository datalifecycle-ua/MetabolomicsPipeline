subpathway_boxplots <- function(analysis_data,chem_data,subpathway,X, groupBy,...){
  
  analysis_data %>%
    filter(...) %>%
    select(group = {{groupBy}}, X = {{X}}, 
           as.character(chem_data$CHEM_ID[which(chem_data$SUB_PATHWAY == subpathway)])) %>%
    pivot_longer(cols = -c(group, X)) %>%
    mutate(X = factor(X)) %>%
    merge(chem_data,by.x = "name",by.y = "CHEM_ID") %>%
    ggplot(aes(x = X, y = value, color = group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position=position_jitterdodge(jitter.width = 0.2)) +
    facet_wrap(~CHEMICAL_NAME) +
    theme_bw()
  
  
}
