
# Setup ---------------------------
library(testthat)
library(tidyFlowCore)

## generate simulated data
simulations <- simulate_cytometry_data()

## generate simulated flowFrame (100 cells, 10 features each)
test_flowframe <- simulations$flowframe
# colnames(flowCore::exprs(test_flowframe)) <-
#   toupper(flowCore::featureNames(test_flowframe))

test_flowset <- simulations$flowset


## read real flow cytometry data from the flowCore package
real_flowframe <-
  system.file("extdata", "0877408774.B08", package="flowCore") |>
  flowCore::read.FCS()

fcs_file_location <- system.file("extdata",package="flowCore")
real_flowset <-
  paste(
    fcs_file_location,
    dir(fcs_file_location),
    sep="/"
  )[1:3] |>
  flowCore::read.flowSet()

# tof_find_panel_info ----------------------------------------------------------

test_that("tof_find_panel_info produces the correct columns", {
  expect_equal(
    test_flowframe |>
      tof_find_panel_info() |>
      colnames(),
    c("channels", "antigens", ".flowCore_featureNames", ".flowCore_colnames")
  )
})

test_that("tof_find_panel_info produces the correct number of rows", {
  expect_equal(
    test_flowframe |>
      tof_find_panel_info() |>
      nrow(),
    10L
  )
})

# as_tof_tbl -------------------------------------------------------------------

## flowFrame -----------------------

### Correct number of rows and columns
test_that("as_tof_tbl produces the correct number of rows and columns", {
expect_equal(
  test_flowframe |>
    as_tof_tbl(.name_method = "tidyFlowCore") |>
    nrow(),
  100L
)

  expect_equal(
    test_flowframe |>
      as_tof_tbl(.name_method = "tidyFlowCore") |>
      ncol(),
    10L
  )
})

### Columns named correctly

test_that("as_tof_tbl produces the correct column names", {

  expect_equal(
    test_flowframe |>
      as_tof_tbl(.name_method = "colnames") |>
      colnames(),
    paste0("feature_", seq_len(10))
  )

  expect_equal(
    test_flowframe |>
      as_tof_tbl(.name_method = "featureNames") |>
      colnames(),
    paste0("feature_", seq_len(10))
  )

  expect_equal(
    test_flowframe |>
      as_tof_tbl(.name_method = "tidyFlowCore") |>
      colnames(),
    paste(
      paste0("FEATURE_", seq_len(10)),
      paste0("feature_", seq_len(10)),
      sep = "|"
      )
  )

})



## flowSet -------------------------

### Correct number of rows and columns
test_that(desc = "as_tof_tbl produces the correct number of rows and columns", {
  expect_equal(
    test_flowset |>
      as_tof_tbl(.name_method = "tidyFlowCore") |>
      nrow(),
    500L
  )

  expect_equal(
    test_flowset |>
      as_tof_tbl(.name_method = "tidyFlowCore") |>
      ncol(),
    10L
  )
})

### Columns named correctly

test_that(desc = "as_tof_tbl produces the correct column names", {

  expect_equal(
    test_flowset |>
      as_tof_tbl(.name_method = "colnames") |>
      colnames(),
    paste0("feature_", seq_len(10))
  )

  expect_equal(
    test_flowset |>
      as_tof_tbl(.name_method = "featureNames") |>
      colnames(),
    paste0("feature_", seq_len(10))
  )

  expect_equal(
    test_flowset |>
      as_tof_tbl(.name_method = "tidyFlowCore") |>
      colnames(),
    paste(
      paste0("FEATURE_", seq_len(10)),
      paste0("feature_", seq_len(10)),
      sep = "|"
    )
  )
})

test_that(desc = "as_tof_tbl can merge metadata columns successfully", {
  expect_contains(
    test_flowset |>
      as_tof_tbl(include_metadata = TRUE) |>
      colnames(),
    c("patient", "cell_type")
  )
})

test_that(desc = "as_tof_tbl can merge .fcs file identifier column successfully", {
  expect_contains(
    test_flowset |>
      as_tof_tbl(include_tidyFlowCore_identifier = TRUE) |>
      colnames(),
    c(".tidyFlowCore_identifier")
  )
})


# tof_get_panel and tof_set_panel ----------------------------------------------

test_that("tof_get_panel returns the correct columns", {
  expect_equal(
    test_flowframe |>
      as_tof_tbl() |>
      tof_get_panel() |>
      colnames(),
    c("channels", "antigens", ".flowCore_featureNames", ".flowCore_colnames")
  )
})

test_that("tof_get_panel can alter a panel", {
  expect_equal(
    test_flowframe |>
      as_tof_tbl() |>
      tof_set_panel(
        panel = data.frame(channels = 'test', antigens = 'test')
      ) |>
      tof_get_panel(),
    data.frame(channels = 'test', antigens = 'test')
  )
})

# tof_tbl tidyverse methods ----------------------------------------------------

my_tof_tbl <-
  test_flowframe |>
  as_tof_tbl(.name_method = "colnames")
my_tof_tbl$group <- sample(c("a", "b"), size = nrow(my_tof_tbl), replace = TRUE)

test_that('nesting a tof_tbl produces a nested tibble', {
  expect_equal(
    my_tof_tbl |>
      tidyr::nest(data = -group) |>
      nrow(),
    2L)
})

test_that('unnesting a nested tof_tbl produces a tof_tbl', {
  expect_equal(
    my_tof_tbl |>
      tidyr::nest(data = -group) |>
      tidyr::unnest(data) |>
      nrow(),
    nrow(my_tof_tbl)
  )
})


test_that('grouping a tof_tbl produces a grouped tof_tbl', {
  expect_contains(
    my_tof_tbl |>
      dplyr::group_by(group) |>
      class(),
    "grouped_tof_tbl"
  )
})




