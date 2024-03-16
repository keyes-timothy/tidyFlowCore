library(tidyFlowCore)
library(dplyr)
library(tidyr)
library(testthat)

# preliminaries ----------------------------------------------------------------
simulations <- simulate_cytometry_data()

## generate simulated flowFrame (100 cells, 10 features each)
test_flowframe <-
  simulations$flowframe |>
  mutate(
    random_group =
      sample(
        c("a", "b"),
        size = nrow(simulations$flowframe),
        replace = TRUE
      )
  )
test_flowset <- simulations$flowset



# mutate ----------------------------------------------------------------------


## flowFrame

test_that("mutate produces a new column in a flowFrame", {
  expect_contains(
    test_flowframe |>
      mutate(new_feature = feature_1 + feature_2) |>
      flowCore::featureNames(),
    "new_feature"
  )
})

## flowSet

test_that("mutate produces a new column in a flowSet", {
  expect_contains(
    test_flowset |>
      mutate(new_feature = feature_1 + feature_2) |>
      (\(.x) .x[[1]])() |>
      flowCore::featureNames(),
    "new_feature"
  )
})


# select -----------------------------------------------------------------------

## flowFrame
test_that("select chooses specified columns in a flowFrame", {
  expect_equal(
    test_flowframe |>
      select(feature_1, feature_2) |>
      flowCore::featureNames() |>
      as.character(),
    c("feature_1", "feature_2")
  )

  expect_equal(
    test_flowframe |>
      select(any_of(c("feature_1", "feature_2"))) |>
      flowCore::featureNames() |>
      as.character(),
    c("feature_1", "feature_2")
  )
})

## flowSet

test_that("select chooses specified columns in a flowSet", {
  expect_equal(
    test_flowset |>
      select(feature_1, feature_2) |>
      (\(.x) .x[[1]])() |>
      flowCore::featureNames() |>
      as.character(),
    c("feature_1", "feature_2")
  )

  expect_equal(
    test_flowset |>
      select(any_of(c("feature_1", "feature_2"))) |>
      (\(.x) .x[[1]])() |>
      flowCore::featureNames() |>
      as.character(),
    c("feature_1", "feature_2")
  )
})


# group_by ---------------------------------------------------------------------

test_that(
  "Grouping a flowFrame returns a flowSet with the correct number of
          constituent flowFrames.",
  {
    expect_equal(
      test_flowframe |>
        group_by(random_group) |>
        flowCore::pData() |>
        nrow(),
      2L
    )
  })

test_that("Grouping a flowFrame returns a flowSet with the correct columns", {
  expect_equal(
    test_flowframe |>
      group_by(random_group) |>
      colnames(),
    colnames(dplyr::select(test_flowframe, -random_group))
  )
})

# summarize --------------------------------------------------------------------

test_that(
  "flowFrame summarize/summarise produce a tibble with the requested summaries",
  {
    expect_equal(
      test_flowframe |>
        summarize(feature_1 = mean(feature_1)) |>
        nrow(),
      1L
    )
    expect_equal(
      test_flowframe |>
        summarize(feature_1 = mean(feature_1)) |>
        ncol(),
      1L
    )
    expect_equal(
      test_flowframe |>
        summarize(feature_1 = mean(feature_1)) |>
        colnames(),
      "feature_1"
    )
  })


test_that(
  "flowFrame summarize/summarise produce a tibble with grouped summaries",
  {
    expect_equal(
      test_flowframe |>
        summarize(feature_1 = mean(feature_1), .by = random_group) |>
        nrow(),
      2L
    )
    expect_equal(
      test_flowframe |>
        summarize(feature_1 = mean(feature_1), .by = random_group) |>
        ncol(),
      2L
    )
    expect_equal(
      test_flowframe |>
        summarize(feature_1 = mean(feature_1), .by = random_group) |>
        colnames(),
      c("random_group", "feature_1")
    )
  })


test_that(
  "flowSet summarize/summarise produce a tibble with the requested summaries",
  {
    expect_equal(
      test_flowset |>
        summarize(feature_1 = mean(feature_1)) |>
        nrow(),
      5L
    )
    expect_equal(
      test_flowset |>
        summarize(feature_1 = mean(feature_1)) |>
        ncol(),
      4L
    )
    expect_equal(
      test_flowset |>
        summarize(feature_1 = mean(feature_1)) |>
        colnames(),
      c("feature_1", ".flowframe_identifier", "patient", "cell_type")
    )
  })

# ungroup ----------------------------------------------------------------------


num_cell_types <-
  test_flowset |>
  flowCore::pData() |>
  dplyr::count(cell_type) |>
  nrow()

test_that(
  "Ungrouping a flowSet entirely returns a flowFrame with the correct number of
  rows",
  {
    expect_equal(
      test_flowset |>
        ungroup(patient, cell_type) |>
        nrow() |>
        as.numeric(),
      test_flowset |>
        as_tof_tbl() |>
        nrow()
    )
  }
)

test_that(
  "Trying to ungroup a flowSet with columns that are not in pData produces an
  error",
  {
    expect_error(
      test_flowset |>
        ungroup(starts_with("feature"))
    )
  }
)

test_that(
  "Ungrouping a flowSet partially (with some remaining groups) returns a
  flowSet with the correct number of experiments.",
  {
    expect_equal(
      test_flowset |>
        ungroup(patient) |>
        length(),
      num_cell_types
    )
  }
)


test_that(
  "Ungrouping a flowSet with no cols specified returns a
  flowFrame with the correct number of rows.",
  {
    expect_equal(
        test_flowset |>
          ungroup() |>
          nrow() |>
          as.numeric(),
      nrow(as_tof_tbl(test_flowset))
    )
  }
)

test_that("Ungrouping a flowSet returns a flowFrame with the correct columns", {
  expect_contains(
      test_flowset |>
        ungroup() |>
        flowCore::colnames(),
    c(
      flowCore::colnames(dplyr::select(test_flowframe, -random_group)),
      "patient",
      "cell_type",
      ".tidyFlowCore_name"
    )
  )
})


# count ------------------------------------------------------------------------


## flowFrame
test_that("flowFrame count produces a tibble with the correct results", {
  expect_equal(
    test_flowframe |>
      dplyr::count() |>
      nrow(),
    1L
  )
  expect_equal(
    test_flowframe |>
      dplyr::count(random_group) |>
      nrow(),
    2L
  )

  expect_equal(
    test_flowframe |>
      dplyr::count() |>
      ncol(),
    1L
  )
  expect_equal(
    test_flowframe |>
      dplyr::count(random_group) |>
      ncol(),
    2L
  )

  expect_equal(
    test_flowframe |>
      dplyr::count() |>
      colnames(),
    "n"
  )
  expect_contains(
    test_flowframe |>
      dplyr::count(random_group) |>
      colnames(),
    c("n", "random_group")
  )

})

## flowSet

test_that("flowSet count produces a tibble with the correct results", {
  expect_equal(
    test_flowset |>
      dplyr::count() |>
      nrow(),
    5L
  )
  expect_equal(
    test_flowset |>
      dplyr::count(cell_type) |>
      nrow(),
    2L
  )

  expect_equal(
    test_flowset |>
      dplyr::count() |>
      ncol(),
    2L
  )
  expect_equal(
    test_flowset |>
      dplyr::count(cell_type) |>
      ncol(),
    2L
  )

  expect_equal(
    test_flowset |>
      dplyr::count() |>
      colnames(),
    c(".tidyFlowCore_identifier", "n")
  )
  expect_contains(
    test_flowset |>
      dplyr::count(cell_type) |>
      colnames(),
    c("n", "cell_type")
  )

})

# pull -------------------------------------------------------------------------

test_that("flowFrame pull takes a single numeric vector", {
  expect_equal(
    test_flowframe |>
      pull(feature_1) |>
      as.numeric(),
    test_flowframe@exprs[,'feature_1']
  )
})

test_that("flowSet pull takes a single numeric vector", {
  expect_equal(
    test_flowset |>
      pull(feature_1) |>
      as.numeric(),
    flowCore::fsApply(test_flowset, exprs)[, 'feature_1']
  )
})

# rename  ----------------------------------------------------------------------

test_that("flowFrame rename renames a column", {
  expect_contains(
    test_flowframe |>
      rename(new_feature = feature_1) |>
      colnames(),
    "new_feature"
  )
})

test_that("flowSet rename renames a column", {
  expect_contains(
    test_flowset |>
      rename(new_feature = feature_1) |>
      (\(.x) .x[[1]])() |>
      colnames(),
    "new_feature"
  )
})

# rename_with  -----------------------------------------------------------------
test_that("flowFrame rename_with rename all columns", {
  expect_contains(
    test_flowframe |>
      rename_with(.fn = toupper) |>
      flowCore::colnames(),
    paste0("FEATURE_", seq_len(ncol(test_flowframe)-1))
  )
})

test_that("flowSet renames all columns", {
  expect_contains(
    test_flowset |>
      rename_with(.fn = toupper) |>
      (\(.x) .x[[1]])() |>
      flowCore::colnames(),
    paste0("FEATURE_", seq_len(ncol(test_flowframe)-1))
  )
})

# slice_*  ---------------------------------------------------------------------

## slice

test_that("flowFrame slice chooses only specified rows", {
  expect_equal(
    test_flowframe |>
      slice(1) |>
      nrow() |>
      as.numeric(),
    1L
  )
})

test_that("flowFrame slice chooses only specified rows", {
  expect_equal(
    test_flowset |>
      slice(1) |>
      fsApply(nrow) |>
      as.numeric() |>
      sum(),
    5L
  )
})

## slice_sample
test_that("flowFrame slice_sample chooses the correct number of rows", {
  expect_equal(
    test_flowframe |>
      slice_sample(n = 100) |>
      nrow() |>
      as.numeric(),
    100L
  )
})

test_that("flowFrame slice_sample chooses the correct number of rows", {
  expect_equal(
    test_flowset |>
      slice_sample(n = 100) |>
      fsApply(nrow) |>
      as.numeric() |>
      sum(),
    500L
  )
})

## slice_max

test_that("flowFrame slice_max chooses the correct number of rows", {
  expect_equal(
    test_flowframe |>
      slice_max(order_by = feature_1, n = 100) |>
      nrow() |>
      as.numeric(),
    100L
  )
})

test_that("flowFrame slice_max chooses the correct number of rows", {
  expect_equal(
    test_flowset |>
      slice_max(order_by = feature_1, n = 100) |>
      fsApply(nrow) |>
      as.numeric() |>
      sum(),
    500L
  )
})





