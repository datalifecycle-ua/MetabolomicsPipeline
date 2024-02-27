
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MetabolomicsPipeline

<!-- badges: start -->
<!-- badges: end -->

The purpose of the MetabolomicsPipeline is to provide tools for
additional analyses to complement the metabolomic analysis done by
Metabolon. Below is a demonstration of a potential workflow using the
MetabolomicsPipeline Package. We demonstrate this workflow in the
“Workflow” vignette.

<img src="Workflow.png"/>

Each step of the workflow is found in the “Code” folder.

1.  Normalization and Standardization
    (NormalizationAndStandardization.R)

- This step is optional as it is already provided by Metabolon.

2.  Analysis Data Creation (AnalysisDataCreation.R)

3.  Exploratory Analysis (ExploratoryAnalysis.R)

4.  SubpathwayAnalysis (subpathwayAnalysis.R)

5.  Pairwise Comparisons (PairwiseAnalysis)

6.  Boxplots and Lineplots (BoxPlotsAndLinePlots.R)

## Download MetabolomicsPipeline

### Check R version

You will need to have the R-version \> 4.3.2 to install the development
version of MetabolomicsPipeline from
[GitHub](https://github.com/JoelParkerUofA/MetabolomicsPipeline). To see
the R version you currently have installed you can use the command

``` r
R.version
```

Check to make sure the R version is at least R 4.3.2. If this
requirement is not met, you will need to update to a newer version of R.
Instructions for how to do this can be found at the [R project
website](https://www.r-project.org/).

### Download pipeline

You can download the pipeline one of two different ways:

1.) **Download pipeline as a zip file:** Navigate to the top of the page
and click the green button that says “\<code\>” and then click “download
zip”. This will download a zip file containing all the files necessary
for the pipeline. You will need to unzip this file to the destination of
your choice.

2.) **Clone pipeline:** If you are familiar with git and github, the
“MetabolomicsPipeline” can be cloned using:

    git clone https://github.com/JoelParkerUofA/MetabolomicsPipeline.git

You only need to utilize one of the two options above to download the
pipeline. In the downloaded/cloned folder, open the “Metabolomics
Pipeline.Rproj” to open the R project.

### Update and install packages

To aid in the reproducibility of the results, it is best practice to
ensure you are using the same package versions that were used to create
the pipeline. To do this, you will need to synchronize your package
versions with the versions used in this pipeline by using:

``` r
renv::restore()
```

Once all of the packages are up to date, you can install the
“MetabolomicsPipeline” package using:

``` r
install.packages("devtools")
devtools::install_github("JoelParkerUofA/MetabolomicsPipeline", upgrade="never")
```

This package will give you access to the many functions used throughout
this pipeline. After this step, you will be ready to use the
MetabolomicsPipeline.

## Getting Started

We demonstrate a workflow using the MetabolomicsPipeline in the
“Workflow” vignette. In this vignette, we use data which consists of 86
samples (42 males, 44 females), three treatment groups, and the samples
were taken at three different time points. We walk through each analysis
step and demonstrate how to use the MetabolomicPipeline package on
metabolomic data from Metabolon. We include this data in the “data”
folder. While the vignette is an excellent starting place, you can also
run each analysis step from the .R files in the “Code” folder.
