#' Display Subpathways Within a Superpathway
#'
#' This function displays subpathway model results based on a specified 
#' superpathway.
#'
#' @param subpath_results Results from a subpathway analysis, either as a data frame or a list of data frames
#' @param Superpath A character string specifying the name of the superpathway to filter by.
#'
#' @return A `knitr::kable` formatted HTML table (if `subpath_results` is a data frame), or
#'         a list of such tables (if `subpath_results` is a list of data frames), showing
#'         subpathways and models within the specified superpathway.
#'
#' @import dplyr
#' @importFrom stringr str_to_title
#' @importFrom knitr kable
#' @importFrom kableExtra kable_paper
#'
#' @examples
#' # 
#' data("demoDataSmall", package = "MetabolomicsPipeline")
#' dat <- demoDataSmall
#'
#' # Runsubpathay analysis
#' sub_analysis <- subpathwayAnalysis(dat,
#'     treat_var = "GROUP_NAME",
#'     block_var = "TIME1",
#'     strat_var = NULL,
#'     Assay = "normalized"
#' )
#' #' # Display subpathways within a specific superpathway
#' subWithinSuperResults(sub_analysis, Superpath = "Amino Acid")
#'
#' @export

subWithinSuperResults <- function(subpath_results, Superpath) {
  
  # Check if subpath_results is a data frame or a list
  if (inherits(subpath_results, "data.frame") |
      inherits(subpath_results, "DFrame")) {
    
    # Filter the results for the specified Superpathway
   table <- subpath_results %>%
     as.data.frame() %>%
      dplyr::filter(super_pathway == Superpath) %>%
      dplyr::select(sub_pathway,model) %>%
      distinct() %>%
      knitr::kable(
        digits = 3,
        col.names = c("Subpathway", "Model"),
        caption = paste0("Subpathways within ",
                         stringr::str_to_title(Superpath)," Superpathway")
      ) %>%
      kableExtra::kable_paper(full_width = FALSE, html_font = "Cambria")
    
    return(table)
  }
  
  if (inherits(subpath_results, "list")) {
    tables <- lapply(names(subpath_results), function(name) {
      stratum <- subpath_results[[name]]
      
      # Filter the results for the specified Superpathway
      stratum %>%
        as.data.frame() %>%
        dplyr::filter(super_pathway == Superpath) %>%
        dplyr::select(sub_pathway, model) %>%
        distinct() 
    })
    
    table <- Reduce(function(x, y) dplyr::full_join(x, y, by = "sub_pathway"),
                    tables)
    
    tab_res <- table %>%
      knitr::kable(
        digits = 3,
        col.names = c("Subpathway", paste0(names(subpath_results), " Model")),
        caption = paste0("Subpathways within ",
                         stringr::str_to_title(Superpath), " Superpathway")
      ) %>%
      kableExtra::kable_paper(full_width = FALSE, html_font = "Cambria")
    
    return(tab_res)
  }
  
}
