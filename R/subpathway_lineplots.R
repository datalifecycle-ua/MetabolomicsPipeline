subpathway_lineplots <- function(analysis_data,chem_data,subpathway,X, groupBy){
  
  analysis_data %>%
    select(treat = {{groupBy}}, X = {{X}},
           as.character(chem_data$CHEM_ID[which(chem_data$SUB_PATHWAY == subpathway)])) %>%
    pivot_longer(cols = -c(treat, X)) %>%
    merge(chem_data,by.x = "name",by.y = "CHEM_ID") %>%
    ggplot(aes(x = X, y = value, color = treat)) +
    geom_jitter() +
    geom_smooth(method = "lm") +
    facet_wrap( ~ CHEMICAL_NAME) +
    theme_bw()
  
}
