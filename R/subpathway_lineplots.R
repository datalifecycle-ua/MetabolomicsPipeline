#' Subpathway Lineplots
#' 
#' Create line plots for each metabolite within a subpathway. 
#' 
#' @param MetPipe MetPipe data object
#' 
#' @param subpathway Character value of the subpathway of interest. This is case 
#' sensitive and must be in the chemical annotation file.
#' 
#' @param X This the the name of the variable in the meta data that is used for the 
#' X axis of the line plots
#' 
#' @param groupBy This is a grouping variable. As a recommendation the treatment
#' groups should be used in the groupBy argument as this will provide a different color
#' for each of the treatments making it easier to identify.
#' 
#' @param ... Additional arguments to filter the analysis data by. 
#' 
#' 
#' @returns Line plots stratified by metabolite.
#' 
#' 
#' @import ggplot2
#' @import tidyr
#' @import dplyr
#' 
#' @export
#' 
#' 


subpathway_lineplots <- function(MetPipe,subpathway,X, groupBy,...){
  
  return(
  MetPipe@analysis %>%
    dplyr::filter(...) %>%
    dplyr::select(treat = {{groupBy}}, X = {{X}},
           as.character(MetPipe@chemical_annotation$CHEM_ID[which(MetPipe@chemical_annotation$SUB_PATHWAY == subpathway)])) %>%
    tidyr::pivot_longer(cols = -c(treat, X)) %>%
    merge(MetPipe@chemical_annotation,by.x = "name",by.y = "CHEM_ID") %>%
    ggplot2::ggplot(aes(x = X, y = value, color = treat)) +
    ggplot2::geom_jitter() +
    ggplot2::geom_smooth(method = "lm") +
    ggplot2::facet_wrap( ~ CHEMICAL_NAME) +
    ggplot2::theme_bw()
  )
}
