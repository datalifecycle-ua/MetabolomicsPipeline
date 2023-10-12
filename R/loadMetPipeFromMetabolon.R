#' Load Metabolomic Pipeline Data
#' 
#' Automatically load metabolomic data from Metabolon
#' 
#' @param metabolon Path to Metabolon .xlsx file
#' 
#' @returns A MetPipe data object which contains 5 data slots.
#' * raw_peak: The raw peak data
#' * standardized_peak: Standardized peak data from Metabolon
#' * meta: Sample metadata
#' * chemical_annotation: Chemical annotation file
#' * analysis: Analysis data used for the downstream analysis
#' 
#' @import readxl
#' @import methods
#'
#' @export
#' 

loadMetPipeFromMetabolon <- function(metabolon_path){
  # Get attributes
  raw_data <- readxl::read_excel(metabolon_path, sheet = "Peak Area Data")
  chem_data <- readxl::read_excel(metabolon_path, sheet = "Chemical Annotation")
  meta <- readxl::read_excel(metabolon_path, sheet = "Sample Meta Data")
  standardized <- readxl::read_excel(metabolon_path, sheet = "Log Transformed Data")
  
  # Create the MetPipe object using the constructor
  met_pipe <- createMetPipe(as.data.frame(raw_data), as.data.frame(standardized), as.data.frame(meta), as.data.frame(chem_data))
  
  # Return the MetPipe object
  return(met_pipe)
}




