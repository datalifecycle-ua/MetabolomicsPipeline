data("demoDat")
library(SummarizedExperiment)
# Sample metadata
sample_metadata <- colData(demoDat)

sample_metadata <- sample_metadata %>%
  as.data.frame() %>%
  tibble::rownames_to_column("PARENT_SAMPLE_NAME")

# Chemical Annotation
chemical_annotation <- rowData(demoDat)

chemical_annotation <- chemical_annotation %>%
  as.data.frame() %>%
  tibble::rownames_to_column("CHEM_ID")

# Peak data
peak_data <- assay(demoDat, "peak")

peak_data <- peak_data %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("PARENT_SAMPLE_NAME")


# Create test set.
test_set <- create_met_se(chemical_annotation = chemical_annotation,
                          sample_metadata = sample_metadata,
                          peak_data = peak_data)


names(colData(demoDat)) <- make.names(names(colData(demoDat)))

test_that("Check data loaded properly (colData)", {
  expect_identical(colData(demoDat), colData(test_set))
})

test_that("Check data loaded properly (rowData)", {
  expect_identical(rowData(demoDat), rowData(test_set))
})

test_that("Check data loaded properly (peak)", {
  expect_identical(assay(demoDat, "peak"), assay(test_set, "peak"))
})
