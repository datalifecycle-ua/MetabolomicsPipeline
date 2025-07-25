% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createHeatmapData.R
\name{createHeatmapData}
\alias{createHeatmapData}
\title{Create metadata and matrices for metabolite heatmaps}
\usage{
createHeatmapData(data, heatmap_variables, Assay = "normalized", ...)
}
\arguments{
\item{data}{A SummarizedExperiment containing metabolomics data.}

\item{heatmap_variables}{A vector of variable names that are NOT metabolites.}

\item{Assay}{Name of assay data to be used for heatmaps.
Default="normalized".}

\item{...}{Additional arguments that can be passed into the arrange function.
This parameter will order the columns of the heatmap data.}
}
\value{
A list of matrices including the heatmap variable
(meta data for heatmap)
and the values for the heatmap.
}
\description{
This function creates the required matrices for the metabolite heatmaps.
}
