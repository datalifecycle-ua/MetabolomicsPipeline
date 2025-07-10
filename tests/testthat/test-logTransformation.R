library(SummarizedExperiment)

# Get demo sample metadata
data("demoSampleMeta", package = "MetabolomicsPipeline")

# Get demo chemical annotation file
data("demoChemAnno", package = "MetabolomicsPipeline")

# Get demo peak data
data("demoPeak", package = "MetabolomicsPipeline")

dat <- createMetSe(chemical_annotation = demoChemAnno,
                     sample_metadata = demoSampleMeta,
                     peak_data = demoPeak,
                     chemical_id = "CHEM_ID",
                     sample_names = "PARENT_SAMPLE_NAME")


# Median standardization
dat <- medianStandardization(met_se = dat,
                              assay = "peak")

# Minimum value imputation
dat <- minValImpute(met_se = dat, assay = "median_std")

# Log transformation
dat <- logTransformation(met_se = dat,
                          assay = "min_impute")



# Get peak data
peak <- assay(dat, "peak")

# Get normalized data
normalized <- assay(dat, "normalized")

# Get peak median for each metabolite
peak_med <- apply(peak, 1, function(x){median(x, na.rm = TRUE)})


test_that("Testing normalization for Metabolite 35", {
  expect_equal(log(peak["35",] / peak_med["35"]), normalized["35", ])
})




