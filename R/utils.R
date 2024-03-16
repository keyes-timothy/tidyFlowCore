# utils.R
# This file contains utility and helper functions for the tidyFlowCore package.


#' Simulate Cytometry Data for FlowSet and FlowFrame Analysis
#'
#' @param num_cells An integer indicating the number of cells to simulate.
#' @param num_features An integer indicating how many features to simulate.
#' @param num_flowframes An integer indicating how many flowFrames to simulate
#'                       for the simulated flowSet.
#'
#' @return A list containing two entries: a flowFrame and a flowSet.
#'
#' @importFrom Biobase AnnotatedDataFrame
#'
#' @importFrom flowCore flowFrame
#' @importFrom flowCore flowSet
#' @importFrom flowCore pData
#' @importFrom flowCore sampleNames
#'
#' @importFrom stats runif
#'
#' @export
#'
#' @examples
#' simulate_cytometry_data()
simulate_cytometry_data <-
  function(
    num_cells = 100,
    num_features = 10,
    num_flowframes = 5
  ) {

    flowFrames <- list()
    ### generate `num_flowframes` matrices with random data and convert them to flowFrame objects
    for(i in seq_len(num_flowframes)) {
      exprs <-
        matrix(
          stats::runif(num_cells * num_features, min = 1, max = 100),
          ncol = num_features
        )
      colnames(exprs) <- paste0("feature_", seq_len(num_features))
      feature_names <- toupper(colnames(exprs))
      parameters <-
        Biobase::AnnotatedDataFrame(
          data =
            data.frame(
              name = colnames(exprs),
              desc = feature_names,
              range = 100,
              minRange = 0,
              maxRange = 100
            )
        )
      flowFrames[[i]] <- flowCore::flowFrame(exprs = exprs, parameters = parameters)
    }
    ### create a flowSet from the list of flowFrame objects
    test_flowset <- flowCore::flowSet(flowFrames)

    ### add metadata to the flowSet
    patient_ids <-
      rep(
        seq_len(num_flowframes %/% 2 + min(num_flowframes %% 2, 1)),
        each = 2,
        length.out = num_flowframes
      )
    cell_types <- rep(c("a", "b"), length.out = num_flowframes)
    new_metadata <- data.frame(
      name = flowCore::sampleNames(test_flowset),
      patient = paste0("patient_", patient_ids),
      cell_type = paste0("population_", cell_types)
    )
    rownames(new_metadata) <- flowCore::sampleNames(test_flowset)
    flowCore::pData(test_flowset) <- new_metadata

    ## generate simulated flowFrame: 100 cells and 10 features
    test_flowframe <- test_flowset[[1]]

    result <-
      list(
        "flowframe" = test_flowframe,
        "flowset" = test_flowset
      )

    return(result)
  }
