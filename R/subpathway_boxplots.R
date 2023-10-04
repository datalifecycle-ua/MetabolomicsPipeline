#' Subpathway Boxplots
#' 
#' Creates boxplots for each metabolite within a specified subpathway. 
#' 
#' @param analysis_data The analysis data containing sample meta data and metabolite peak data.
#' 
#' @param chem_data Chemical Annotation data provided by Metabolon.
#' 
#' @param subpathway Character value of the subpathway of interest. This is case 
#' sensitive and must be in the chem_data.
#' 
#' @param X This the the name of the variable in the meta data that is used for the 
#' X axis of the box plots.
#' 
#' @param groupBy This is a grouping variable. As a recommendation the treatment
#' groups should be used in the groupBy argument as this will provide a different color
#' for each of the treatments making it easier to identify.
#' 
#' 
#' @returns Boxplots stratified by metabolites.
#' 
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' 
#' @export

subpathway_boxplots <- function(analysis_data,chem_data,subpathway,X, groupBy,...){
  
  return(
  analysis_data %>%
    dplyr::filter(...) %>%
    dplyr::select(group = {{groupBy}}, X = {{X}}, 
           as.character(chem_data$CHEM_ID[which(chem_data$SUB_PATHWAY == subpathway)])) %>%
    tidyr::pivot_longer(cols = -c(group, X)) %>%
    dplyr::mutate(X = factor(X)) %>%
    merge(chem_data,by.x = "name",by.y = "CHEM_ID") %>%
    ggplot2::ggplot(aes(x = X, y = value, color = group)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_point(position=position_jitterdodge(jitter.width = 0.2)) +
    ggplot2::facet_wrap(~CHEMICAL_NAME) +
    ggplot2::theme_bw()
  )
}
