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


# nest -------------------------------------------------------------------------

test_that(
  "Nesting a flowFrame returns a flowSet with the correct number of
          constituent flowFrames.",
  {
    expect_equal(
      test_flowframe |>
        tidyr::nest(.by = random_group) |>
        flowCore::pData() |>
        nrow(),
      2L
    )
  })

test_that("Nesting a flowFrame returns a flowSet with the correct columns", {
  expect_equal(
    test_flowframe |>
      tidyr::nest(.by = random_group) |>
      colnames(),
    colnames(dplyr::select(test_flowframe, -random_group))
  )
})



# unnest -----------------------------------------------------------------------


num_cell_types <-
  test_flowset |>
  flowCore::pData() |>
  dplyr::count(cell_type) |>
  nrow()

test_that(
  "Unnesting a flowSet entirely returns a flowFrame with the correct number of
  rows",
  {
    expect_equal(
      test_flowset |>
        tidyr::unnest(cols = c(patient, cell_type)) |>
        nrow() |>
        as.numeric(),
      test_flowset |>
        as_tof_tbl() |>
        nrow()
    )
  }
)

test_that(
  "Trying to unnest a flowSet with columns that are not in pData produces an
  error",
  {
    expect_error(
      test_flowset |>
        tidyr::unnest(cols = starts_with("feature"))
    )
  }
)

test_that(
  "Unnesting a flowSet partially (with some remaining groups) returns a
  flowSet with the correct number of experiments.",
  {
    expect_equal(
      test_flowset |>
        tidyr::unnest(cols = patient) |>
        length(),
      num_cell_types
    )
  }
)


test_that(
  "Unnesting a flowSet with no cols specified returns a
  flowFrame with the correct number of rows.",
  {
    expect_equal(
      suppressWarnings(
        test_flowset |>
          tidyr::unnest() |>
          nrow() |>
          as.numeric()
      ),
      nrow(as_tof_tbl(test_flowset))
    )
  }
)

test_that("Unnesting a flowSet returns a flowFrame with the correct columns", {
  expect_contains(
    suppressWarnings(
      test_flowset |>
      tidyr::unnest() |>
      flowCore::colnames()
    ),
    c(
      flowCore::colnames(dplyr::select(test_flowframe, -random_group)),
      "patient",
      "cell_type",
      ".tidyFlowCore_name"
    )
  )
})


