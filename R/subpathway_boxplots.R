#' Subpathway Boxplots
#' 
#' Creates boxplots for each metabolite within a specified subpathway. 
#' 
#' .
#' @param MetPipe MetPipe data object
#' @param subpathway Character value of the subpathway of interest. This is case 
#' sensitive and must be in the chemical annotation file.
#' 
#' @param block_var This the the name of the variable in the meta data that is used for the 
#' X axis of the box plots. We recommend using the "block_var" from the subpathway 
#' analysis.  
#' 
#' @param treat_var This is a grouping variable. As a recommendation the treatment
#' groups should be used in the treat_var argument as this will provide a different color
#' for each of the treatments making it easier to identify.
#' 
#' @param ... Additional arguments to filter the analysis data by.
#' 
#' 
#' @returns Boxplots stratified by metabolites.
#' 
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' 
#' @export

subpathway_boxplots <- function(MetPipe,subpathway,block_var, treat_var,...){
  
  return(
  MetPipe@analysis %>%
    dplyr::filter(...) %>%
    dplyr::select(group = {{treat_var}}, X = {{block_var}}, 
           as.character(MetPipe@chemical_annotation$CHEM_ID[which(MetPipe@chemical_annotation$SUB_PATHWAY == subpathway)])) %>%
    tidyr::pivot_longer(cols = -c(group, X)) %>%
    dplyr::mutate(X = factor(X)) %>%
    merge(MetPipe@chemical_annotation,by.x = "name",by.y = "CHEM_ID") %>%
    ggplot2::ggplot(aes(x = X, y = value, color = group)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_point(position=position_jitterdodge(jitter.width = 0.2)) +
    ggplot2::facet_wrap(~CHEMICAL_NAME) +
    ggplot2::theme_bw()
  )
}
