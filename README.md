
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MetabolomicsPipeline

<!-- badges: start -->
<!-- badges: end -->

The purpose of the MetabolomicsPipeline is to provide tools for
additional analysis to complement the metabolomic analysis done by
Metabolon.

<img src="Workflow.png"/>

Each step of the workflow is found in the Code folder.

1.  Normalization and Standardization
    (NormalizationAndStandardization.R)

- This step is turned off by default since this step is completed by
  Metabolon

- More detail in the Workflow vignett.

2.  Analysis Data Creation (AnalysisDataCreation.R)

3.  Exploratory Analysis (ExploratoryAnalysis.R)

4.  SubpathwayAnalysis (subpathwayAnalysis.R)

5.  Pairwise Comparisons (PairwiseAnalysis)

6.  Boxplots and Lineplots (BoxPlotsAndLinePlots.R)

## Installation of Package

You can install the development version of MetabolomicsPipeline from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("JoelParkerUofA/MetabolomicsPipeline")
```

## Getting Started

## Summary

MetabolomicsPipeline can be utilized as a package or a pipline. If you
are using the pipeline option, you must install the MetabolomicsPipeline
package prior to cloning the github project.
