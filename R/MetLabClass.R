#' Define metlab class
#' 
#' @slot raw_peak Raw Peak Data
#' 
#' @slot standardized_peak Standardized Peak data
#' 
#' @slot meta Sample Metadata
#' 
#' @slot chemical_annotation Chemical Annotation file
#' 
#' @slot analysis Analysis data file
#' 
#' @name MetPipe
#' 
#' @rdname MetPipe
#' 
#' @export
#' 
#' @docType methods
#' 

# Define the class and slots outside the function
MetLab <- methods::setClass("MetPipe", slots = list(
  raw_peak = "data.frame",
  standardized_peak = "data.frame",
  meta = "data.frame",
  chemical_annotation = "data.frame",
  analysis = "data.frame"
))