% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/met_est_heatmap.R
\name{met_est_heatmap}
\alias{met_est_heatmap}
\title{Metabolite Pairwise Estimate Interactive Heatmap.}
\usage{
met_est_heatmap(
  results_data,
  data,
  interactive = FALSE,
  SUB_PATHWAY = "SUB_PATHWAY",
  CHEMICAL_NAME = "CHEMICAL_NAME",
  plotlyTitle = "Metabolite log fold change",
  ...
)
}
\arguments{
\item{results_data}{Results data frame of the pairwise comparisons produced
by \code{\link{metabolite_pairwise}}.}

\item{data}{A SummarizedExperiment containing the metabolomics experiment data.}

\item{interactive}{boolean (T/F) for whether or not the plot should be
interactive. Use interactive=T to produce an interactive plot using
plotly. Use interactive=F to produce a static heatmap using pheatmap.}

\item{SUB_PATHWAY}{Column name in the chemical annotation worksheet which
contains the subpathway information.}

\item{CHEMICAL_NAME}{Column name in the chemical annotation worksheet which
contains the chemical name.}

\item{plotlyTitle}{Title for the interactive heatmap.}

\item{...}{Additional arguments that can be passed to pheatmap.}
}
\value{
An interactive heatmap of pairwise estimates.
}
\description{
Produce an interactive heatmap of the estimates produced in
\code{\link{metabolite_pairwise}}.
}
\details{
This function will produce a heatmap of the log fold changes for the
metabolites with a significant overall p-value (which tested if the treatment
group means were equal under the null hypothesis). The heatmap colors will
only show if the log fold-change is greater than log(2) or less than
log(.5). Therefore, this heatmap will only focus on comparisons with a
fold change of two or greater.
}
\examples{
# Load data
data("demoDataSmall", package = "MetabolomicsPipeline")
dat <- demoDataSmall

# Run pairwise analysis
strat_pairwise <- metabolite_pairwise(dat,
    form = "GROUP_NAME*TIME1",
    strat_var = "Gender"
)

#############################################################################
## Create Estimate Heatmap ##################################################
#############################################################################

met_est_heatmap(strat_pairwise$Female, dat,
               interactive = FALSE,
                SUB_PATHWAY = "SUB_PATHWAY",
               CHEMICAL_NAME = "CHEMICAL_NAME",
               plotlyTitle = "Metabolite log fold change",
               main = "Log fold change heatmap", show_rownames = FALSE
)


##############################################################################
## Create P-value Heatmap ####################################################
##############################################################################
# Female
met_p_heatmap(strat_pairwise$Female, dat,
             interactive = FALSE, show_rownames = FALSE,
             plotlyTitle = "P-Value Heatmap",
             main = "Pvalue Heatmap"
)


}
