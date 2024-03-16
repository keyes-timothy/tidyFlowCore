library(tidyFlowCore)
library(dplyr)
library(ggplot2)
library(testthat)

# preliminaries ----------------------------------------------------------------
simulations <- simulate_cytometry_data()

## generate simulated flowFrame (100 cells, 10 features each)
test_flowframe <- simulations$flowframe
test_flowset <- simulations$flowset

# ggplot -----------------------------------------------------------------------


## flowFrame

flowframe_plot <-
  test_flowframe |>
  ggplot2::ggplot(aes(x = feature_1, y = feature_2)) +
  ggplot2::geom_point()

test_that("ggplot creates a ggplot object when applied to a flowFrame", {
  expect_s3_class(
    flowframe_plot,
    "ggplot"
  )
})

## flowSet


flowset_plot <-
  test_flowset |>
  ggplot2::ggplot(aes(x = feature_1, y = feature_2)) +
  ggplot2::geom_point()

flowset_plot_with_metadata <-
  test_flowset |>
  ggplot2::ggplot(aes(x = feature_1, y = feature_2, color = patient)) +
  ggplot2::geom_point()

test_that("ggplot creates a ggplot object when applied to a flowSet", {
  expect_s3_class(
    flowset_plot,
    "ggplot"
  )

  expect_s3_class(
    flowset_plot_with_metadata,
    "ggplot"
  )
})
