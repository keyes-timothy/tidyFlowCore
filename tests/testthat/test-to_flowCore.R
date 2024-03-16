
# Setup ---------------------------
library(testthat)
library(tidyFlowCore)

## generate simulated data
simulations <- simulate_cytometry_data()

## generate simulated flowFrame (100 cells, 10 features each)
test_flowframe <- simulations$flowframe

test_flowframe_tbl <-
  test_flowframe |>
  as_tof_tbl(.name_method = "colnames")

test_flowset <- simulations$flowset

test_flowset_tbl <-
  test_flowset |>
  as_tof_tbl(.name_method = "colnames", include_metadata = TRUE)


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


# as_flowFrame -----------------------------------------------------------------

## without character vectors
test_that("as_flowFrame produces the correct number of rows and columns", {
  expect_equal(
    test_flowframe_tbl |>
      as_flowFrame() |>
      nrow(),
    nrow(test_flowframe)
  )
  expect_equal(
    test_flowframe_tbl |>
      as_flowFrame() |>
      ncol(),
    ncol(test_flowframe)
  )
})

## with character vectors

test_that("as_flowFrame produces the correct number of rows and columns", {
  expect_equal(
    test_flowset_tbl |>
      as_flowFrame() |>
      nrow() |>
      as.numeric(),
    test_flowset |>
      flowCore::fsApply(nrow) |>
      sum()
  )
  expect_equal(
    test_flowset_tbl |>
      as_flowFrame() |>
      ncol() |>
      as.numeric(),
    (test_flowset |>
      flowCore::fsApply(ncol) |>
      mean()) + 2L
  )
})

test_that("as_flowFrame produces the correct number character codes", {
  expect_contains(
    test_flowset_tbl |>
      as_flowFrame() |>
      flowCore::keyword("character_codes") |>
      (\(.x) .x[[1]])() |>
      names(),
    c("patient", "cell_type")
  )

  expect_contains(
    test_flowset_tbl |>
      as_flowFrame() |>
      flowCore::keyword() |>
      names(),
    c("patient_codes", "cell_type_codes")
  )

})


# as_flowSet ---------------------------------------------

no_grouping_flowframe <-
  test_flowset_tbl |>
  as_flowSet()

patient_flowset <-
  test_flowset_tbl |>
  as_flowSet(group_cols = patient)

cell_type_flowset <-
  test_flowset_tbl |>
  as_flowSet(group_cols = cell_type)

combined_flowset <-
  test_flowset_tbl |>
  as_flowSet(group_cols = c(patient, cell_type))

test_that(
  "as_flowSet produces a flowFrame when no group_cols are used, but produces
  a flowSet otherwise.",
  {
    expect_s4_class(no_grouping_flowframe, "flowFrame")

    expect_s4_class(patient_flowset, "flowSet")

    expect_s4_class(cell_type_flowset, "flowSet")

    expect_s4_class(combined_flowset, "flowSet")
  }
)

test_that(
  "as_flowSet generates the correct metadata in flowCore::pData()
          from group_cols",
  {
    expect_contains(
      colnames(flowCore::pData(patient_flowset)),
      c("name", ".tidyFlowCore_unique_identifier", "patient")
    )

    expect_contains(
      colnames(flowCore::pData(cell_type_flowset)),
      c("name", ".tidyFlowCore_unique_identifier", "cell_type")
    )

    expect_contains(
      colnames(flowCore::pData(combined_flowset)),
      c("name", ".tidyFlowCore_unique_identifier", "cell_type", "patient")
    )
  }
)

test_that(
  "as_flowSet converts character columns to codes in constituent flowFrames",
  {
    expect_type(
      patient_flowset[[1]]@exprs[, "cell_type"],
      "double"
    )

    expect_type(
      cell_type_flowset[[1]]@exprs[, "patient"],
      "double"
    )

    expect_equal(
      patient_flowset[[1]] |>
        flowCore::keyword("character_codes") |>
        (\(.x) .x[[1]]$cell_type$code)() |>
        sort(),
      patient_flowset[[1]]@exprs[, "cell_type"] |>
        unique() |>
        sort()
      )

    expect_equal(
      cell_type_flowset[[1]] |>
        flowCore::keyword("character_codes") |>
        (\(.x) .x[[1]]$patient$code)() |>
        sort(),
      cell_type_flowset[[1]]@exprs[, "patient"] |>
        unique() |>
        sort()
    )
  }
)


