#' Subpathway Boxplots
#' 
#' Creates boxplots for each metabolite within a specified subpathway. 
#' 
#' .
#' @param data SummarizedExperiment with Metabolon experiment data.
#' 
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
#' @param Assay Name of the assay to be used for the pairwise analysis (default='normalized')
#' 
#' @param ... Additional arguments to filter the analysis data by.
#' 
#' 
#' @examples
#' # load data
#' dat = MetabolomicsPipeline::demoDat
#' 
#' ################################################################################
#' ### BoxPlots ###################################################################
#' ################################################################################
#'
#' subpathway_boxplots(dat, subpathway = "Lactoyl Amino Acid", block_var = TIME1,
#'                    treat_var = GROUP_NAME, Assay = "normalized",Gender =="Female")
#'
#'
#' ################################################################################
#' ## Line plots ##################################################################
#' ################################################################################
#'
#' # Set up data
#' dat$TIME1 <- as.numeric(factor(dat$TIME1,
#'                               levels = c("PreSymp","Onset","End")))

#'# Create line plots 
#'subpathway_lineplots(dat, subpathway = "Lactoyl Amino Acid",
#'                     block_var = TIME1,treat_var = GROUP_NAME, Assay = "normalized",Gender=="Female" ) +
#'  xlab("Time")
#' 
#' 
#' 
#' @returns Boxplots stratified by metabolites.
#' 
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_boxplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 theme_bw
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' 
#' @export

subpathway_boxplots <- function(data,subpathway,block_var, treat_var,Assay="normalized",...){
  
  # Create analysis data
  analysis <- SummarizedExperiment::colData(data) %>%
    merge(t(SummarizedExperiment::assay(data,Assay)), by="row.names") %>%
    dplyr::rename(PARENT_SAMPLE_NAME=Row.names)
  
  
  # Create chem data
  chem <- SummarizedExperiment::rowData(data) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("CHEM_ID")
  
  return(
    analysis %>%
    dplyr::filter(...) %>%
    dplyr::select(group = {{treat_var}}, X = {{block_var}}, 
           as.character(chem$CHEM_ID[which(chem$SUB_PATHWAY == subpathway)])) %>%
    tidyr::pivot_longer(cols = -c(group, X)) %>%
    dplyr::mutate(X = factor(X)) %>%
    merge(chem,by.x = "name",by.y = "CHEM_ID") %>%
    ggplot2::ggplot(aes(x = X, y = value, color = group)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_point(position=position_jitterdodge(jitter.width = 0.2)) +
    ggplot2::facet_wrap(~CHEMICAL_NAME) +
    ggplot2::theme_bw()
  )
}
