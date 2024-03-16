# tof_tbl.R
# This file contains functions for creating and manipulating a behind-the-scenes
# flowCore-tibble abstraction (called a "tof_tbl") for use in tidyFlowCore.


data(metal_masterlist, envir = environment())

# tof_find_panel_info ----------------------------------------------------------

#' Use tidyFlowCore's opinionated heuristic for extracting a high-dimensional cytometry panel's channel-antigen pairs
#' from a flowFrame (read from a .fcs file.)
#'
#' Using the character vectors obtained from the `name` and `desc` columns of
#' the parameters of the data of a flowFrame, infer the cytometry panel used
#' to collect the data and return it as a tidy tibble.
#'
#' @param input_flowFrame A flowFrame (just read from an .fcs file) from which
#' a high-dimensional cytometry panel should be extracted
#'
#' @return A tibble with 4 columns (`channels`, `antigens`, `.flowCore_featureNames`
#' and `.flowCore_colnames`). The first two columns correspond to the
#' channels and antigens of the high-dimensional cytometry panel used during
#' data acquisition, respectively. The last two channels represent the featureNames
#' and colnames attributes used to represent each channel in the input flowFrame.
#'
#' @importFrom stringr str_detect
#' @importFrom stringr str_extract
#' @importFrom stringr str_remove
#' @importFrom stringr str_remove_all
#' @importFrom stringr str_c
#'
#' @importFrom dplyr if_else
#' @importFrom dplyr tibble
#'
#'
tof_find_panel_info <- function(input_flowFrame) {

  # extract "desc" and "names" information from the flowFrame
  data_desc <- as.character(input_flowFrame@parameters@data$desc)
  data_names <- as.character(input_flowFrame@parameters@data$name)

  # find channels
  # start by manually looking for CyTOF metals, as they have irregular naming
  # conventions across many mass cytometers
  channels <-
    dplyr::if_else(
      stringr::str_detect(
        data_names,
        pattern = stringr::str_c(metal_masterlist, collapse = "|")
      ),
      stringr::str_c(
        stringr::str_extract(
          data_names,
          pattern = stringr::str_c(metal_masterlist, collapse = "|")
        ) |>
          stringr::str_extract("[:alpha:]+"),
        stringr::str_extract(
          data_names,
          pattern = stringr::str_c(metal_masterlist, collapse = "|")
        ) |>
          stringr::str_extract("[:digit:]+")
      ),
      stringr::str_c(
        stringr::str_extract(
          data_desc,
          pattern = stringr::str_c(metal_masterlist, collapse = "|")
        ) |>
          stringr::str_extract("[:alpha:]+"),
        stringr::str_extract(
          data_desc,
          pattern = stringr::str_c(metal_masterlist, collapse = "|")
        ) |>
          stringr::str_extract("[:digit:]+")
      )
    )


  # if no metal could be detected, just use whatever was in the names
  # slot
  channels <-
    dplyr::if_else(
      is.na(channels),
      data_names,
      channels
    )

  # find antigens --------------------------------------------------------------

  # first, look in the description slot and remove any metal patterns. What
  # remains (minus any punctuation) is a candidate antigen name.
  antigens <-
    dplyr::if_else(
      stringr::str_detect(data_desc, pattern = stringr::str_c(metal_masterlist, collapse = "|")),
      stringr::str_remove(data_desc, pattern = stringr::str_c(metal_masterlist, collapse = "|")),
      data_desc
    ) |>
    stringr::str_remove("^[:punct:]|[:punct:]$") |>
    stringr::str_remove_all("\\(|\\)|Di")

  # if a given antigen name is empty after the first round of candidates is
  # explored, check the names slot. Remove any metal patterns (and punctuation)
  # and what remains should be the antigen name.
  antigens <-
    dplyr::if_else(
      antigens == "" | is.na(antigens),
      stringr::str_remove(data_names, pattern = stringr::str_c(metal_masterlist, collapse = "|")),
      antigens
    ) |>
    stringr::str_remove("^[:punct:]|[:punct:]$") |>
    stringr::str_remove_all("\\(|\\)|Di")

  # if the antigen name of any given channel is still empty (or NA), just put
  # the word "empty"
  antigens <-
    dplyr::if_else((antigens == "" | is.na(antigens)), "empty", antigens)

  # return result
  result <-
    dplyr::tibble(
      channels = channels,
      antigens = antigens,
      .flowCore_featureNames =
        as.character(flowCore::featureNames(input_flowFrame)),
      .flowCore_colnames =
        as.character(flowCore::colnames(input_flowFrame))
    )

  return(result)
}


# new_tof_tbl ------------------------------------------------------------------

#' Constructor for a tof_tibble.
#'
#' @param x A data.frame or tibble containing single-cell mass cytometry data
#' such that rows are cells and columns are CyTOF measurements.
#'
#' @param panel A data.frame or tibble containing information about the panel
#' for the mass cytometry data in x.
#'
#' @return A `tof_tbl`, a tibble extension that tracks a few other attributes
#' that are useful for CyTOF data analysis.
#'
#' @family tof_tbl utilities
#'
#' @importFrom dplyr tibble
#' @importFrom tibble new_tibble
#'
new_tof_tibble <- function(x = dplyr::tibble(), panel = dplyr::tibble()) {

  stopifnot(inherits(x, "tbl_df"))
  stopifnot(inherits(panel, "tbl_df"))

  if("grouped_df" %in% class(x)) {
    subclasses <- c("grouped_tof_tbl", "grouped_df", "tof_tbl")
  } else {
    subclasses <- "tof_tbl"
  }

  tibble::new_tibble(
    x,
    panel = panel,
    nrow = nrow(x),
    class = subclasses
  )
}

# as_tof_tbl -------------------------------------------------------------------

#' Coerce flowFrames or flowSets into tibbles.
#'
#' @param flow_data A flowFrame or flowSet
#'
#' @param .name_method A string indicating how tidyFlowCore should extract column
#' names from `flow_data`. Available options are "tidyFlowCore" (the default), which
#' uses tidyFlowCore's internal heuristic to name columns; "featureNames", which
#' uses \code{\link[flowCore]{featureNames}} to name the columns; and "colnames",
#' which uses \code{\link[base]{colnames}} to name the columns. Note that,
#' in most cases, "featureNames" and "colnames" will give identical results.
#'
#' @param sep A string indicating which symbol should be used to separate
#' antigen names and channel names in the columns of the output tof_tbl
#' when .name_method = 'tidyFlowCore'.
#'
#' @param ... Optional method-specific arguments.
#'
#' @export
#'
#' @return A cytometry-specialized tibble called a `tof_tbl`.
#'
#' @examples
#'
#' input_file <- system.file("extdata", "0877408774.B08", package="flowCore")
#'
#' input_flowframe <- flowCore::read.FCS(input_file)
#'
#' tof_tibble <- as_tof_tbl(input_flowframe)
#'
as_tof_tbl <-
  function(
    flow_data,
    .name_method = c("tidyFlowCore", "featureNames", "colnames"),
    sep = "|",
    ...
  ) {
    UseMethod("as_tof_tbl")
  }

#' Convert an object into a tibble-flowCore abstraction (a `tof_tbl`)
#'
#' @param flow_data A FlowSet
#'
#' @param .name_method A string indicating how tidyFlowCore should extract column
#' names for the output tof_tbl from `flow_data`.
#' Available options are "tidyFlowCore" (the default), which
#' uses tidyFlowCore's internal heuristic to name columns; "featureNames", which
#' uses \code{\link[flowCore]{featureNames}} to name the columns; and "colnames",
#' which uses \code{\link[base]{colnames}} to name the columns.
#'
#' @param sep A string to use to separate the antigen name and its associated
#' channel name in the column names of the output tibble. Defaults to "|".
#'
#' @param include_metadata A boolean value indicating if the metadata for each
#' .fcs file read by flowCore (stored in \code{\link[Biobase]{pData}})
#' should be merged into the final result. Defaults to FALSE.
#'
#' @param include_tidyFlowCore_identifier A boolean value indicating if
#' tidyFlowCore's internal identifier for each flowFrame in the flowSet
#' should be included in the output tof_tbl result. Defaults to FALSE.
#'
#' @param ... Currently unused.
#'
#' @export
#'
#' @importFrom flowCore colnames
#' @importFrom flowCore fsApply
#' @importFrom flowCore featureNames
#' @importFrom flowCore exprs
#' @importFrom flowCore identifier
#' @importFrom flowCore pData
#'
#' @importFrom dplyr as_tibble
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate
#' @importFrom dplyr rename
#'
#' @importFrom rlang arg_match
#'
#' @return A cytometry-specialized tibble called a `tof_tbl`.
#'
#'
as_tof_tbl.flowSet <-
  function(
    flow_data,
    .name_method = c("tidyFlowCore", "featureNames", "colnames"),
    sep = "|",
    ...,
    include_metadata = FALSE,
    include_tidyFlowCore_identifier = FALSE
  ) {

    # check .name_method argument
    .name_method <- rlang::arg_match(.name_method)

    # check if flowset is empty
    if (length(flow_data) < 1) {
      stop("This flowSet is empty.")
    }
    # get panel information
    panel_info <-
      flow_data[[1]] |>
      tof_find_panel_info()

    # get protein measurements for all cells
    flowset_exprs <-
      flow_data |>
      flowCore::fsApply(FUN = flowCore::exprs) |>
      dplyr::as_tibble()

    if (include_metadata | include_tidyFlowCore_identifier) {
      # get metadata from each flowFrame's identifier and flowCore::pData
      flowset_identifiers <-
        flow_data |>
        flowCore::fsApply(FUN = flowCore::identifier) |>
        as.character()

      flowset_num_cells <-
        flow_data |>
        flowCore::fsApply(FUN = flowCore::nrow) |>
        as.integer()

      flowset_metadata <-
        flow_data |>
        flowCore::pData()

      # experimental
      ## adds a unique tidyflowcore identifier to each row of the metadata
      ## matrix
      if (!(".tidyFlowCore_unique_identifier" %in% colnames(flowset_metadata))) {
        flowset_metadata <-
          flowset_metadata |>
          dplyr::mutate(
            .tidyFlowCore_unique_identifier = as.character(.data$name)
          )
      }

      ## TODO: check
      result_metadata <-
        data.frame(
          .tidyFlowCore_unique_identifier =
            rep(flowset_identifiers, times = flowset_num_cells)
        ) |>
        # location of key mismatch
        #dplyr::left_join(flowset_metadata, by = "name") |>
        ## TODO: This section is experimental.
        ##Test with repeated groupings and ungroupings
        dplyr::left_join(
          flowset_metadata,
          by = c(".tidyFlowCore_unique_identifier")
        ) |>
        dplyr::select(-"name")

      if (!include_tidyFlowCore_identifier) {
        result_metadata <-
          result_metadata |>
          dplyr::select(-".tidyFlowCore_unique_identifier")
      } else {
        if (!include_metadata) {
          result_metadata <-
            result_metadata |>
            dplyr::select(".tidyFlowCore_unique_identifier")
        }
        result_metadata <-
          result_metadata |>
          dplyr::rename(
            .tidyFlowCore_identifier = ".tidyFlowCore_unique_identifier"
          )

      }

    }

    if (.name_method == "tidyFlowCore") {
      col_names <-
        base::paste(panel_info$antigens, panel_info$channels, sep = sep)

      # prevent repeating names twice when antigen and channel are identical
      repeat_indices <-
        which(panel_info$channels == panel_info$antigens)
      col_names[repeat_indices] <- panel_info$antigens[repeat_indices]

    } else if (.name_method == "featureNames") {
      col_names <-
        panel_info$.flowCore_featureNames

    } else if (.name_method == "colnames") {
      col_names <-
        panel_info$.flowCore_colnames


    } else {
      stop("Not a valid .name_method.")
    }

    colnames(flowset_exprs) <- col_names

    result <- new_tof_tibble(x = flowset_exprs, panel = panel_info)

    if (include_metadata | include_tidyFlowCore_identifier) {
      result <- dplyr::bind_cols(result, result_metadata)
    }

    return(result)
  }

#' @export
#'
#' @importFrom dplyr as_tibble
#'
#' @importFrom flowCore colnames
#' @importFrom flowCore exprs
#' @importFrom flowCore featureNames
#'
#' @importFrom rlang arg_match
#'
#' @importFrom stats setNames
#'
as_tof_tbl.flowFrame <-
  function(
    flow_data,
    .name_method = c("tidyFlowCore", "featureNames", "colnames"),
    sep = "|",
    ...
  ) {

    # check .name_method argument
    .name_method <- rlang::arg_match(.name_method)

    panel_info <-
      flow_data |>
      tof_find_panel_info()

    if (.name_method == "tidyFlowCore") {
      col_names <-
        base::paste(panel_info$antigens, panel_info$channels, sep = sep)

      # prevent repeating names twice when antigen and metal are identical
      repeat_indices <-
        which(panel_info$channels == panel_info$antigens)
      col_names[repeat_indices] <- panel_info$antigens[repeat_indices]

    } else if (.name_method == "featureNames") {
      col_names <- flowCore::featureNames(flow_data)

    } else if (.name_method == "colnames") {
      col_names <- flowCore::colnames(flow_data)

    } else {
      stop("Not a valid .name_method.")
    }

    flowframe_exprs <-
      stats::setNames(
        object = dplyr::as_tibble(flowCore::exprs(flow_data)),
        nm = col_names
      )

    result <-
      new_tof_tibble(
        x = flowframe_exprs,
        panel = panel_info
      )

    return(result)

  }

#' @export
as_tof_tbl.data.frame <-
  function(
    flow_data,
    .name_method = c("tidyFlowCore", "featureNames", "colnames"),
    sep = "|",
    ...
  ) {
    result <- new_tof_tibble(
      x = dplyr::as_tibble(flow_data),
      panel = dplyr::tibble()
    )
    return(result)
  }



# tof_get_panel ----------------------------------------------------------------

#' Get panel information from a tof_tibble
#'
#' @param tof_tibble A `tof_tbl`.
#'
#' @return A tibble containing information about the CyTOF panel
#' that was used during data acquisition for the data contained
#' in `tof_tibble`.
#'
#' @family tof_tbl utilities
#'
#' @examples
#' NULL
#'
#'
tof_get_panel <- function(tof_tibble) {
  panel <-
    tof_tibble |>
    attr(which = "panel")

  return(panel)
}


# tof_set_panel ----------------------------------------------------------------

#' Set panel information from a tof_tbl
#'
#' @param tof_tibble A `tof_tbl`.
#'
#' @param panel A data.frame containing two columns (`channels` and `antigens`) representing
#' the information about a panel
#'
#' @return A `tof_tbl` containing information about the CyTOF panel
#' that was used during data acquisition for the data contained
#' in the input `tof_tibble`. Two columns are required: "metals" and "antigens".
#'
#' @family tof_tbl utilities
#'
#'
#' @examples
#' NULL
#'
#'
tof_set_panel <- function(tof_tibble, panel) {
  attr(tof_tibble, which = "panel") <- panel
  return(tof_tibble)
}

# tof_tbl tidyverse methods ----------------------------------------------------

## tidyr methods

#' @export
#'
#' @importFrom tidyr nest
nest.tof_tbl <- function(.data, ..., .names_sep = NULL) {
  panel <- tof_get_panel(.data)
  return(new_tof_tibble(x = NextMethod(), panel = panel))
}

#' @export
#'
#' @importFrom tidyr unnest
unnest.tof_tbl <- function(data, ...) {
  start_panel <- tof_get_panel(data)
  return(new_tof_tibble(x = NextMethod(), panel = start_panel))
}

#' @export
#'
#' @importFrom tidyr pivot_longer
pivot_longer.tof_tbl <- function(data, ...) {
  panel <- tof_get_panel(data)
  return(new_tof_tibble(x = NextMethod(), panel = panel))
}

#' @export
#'
#' @importFrom tidyr pivot_wider
pivot_wider.tof_tbl <- function(data, ...) {
  panel <- tof_get_panel(data)
  return(new_tof_tibble(x = NextMethod(), panel = panel))
}

## dplyr methods

#' @export
#'
#' @importFrom dplyr group_by
#' @importFrom dplyr group_by_drop_default
group_by.tof_tbl <-
  function(.data, ..., .add = FALSE, .drop = dplyr::group_by_drop_default(.data)) {
    panel <- tof_get_panel(.data)
    return(new_tof_tibble(x = NextMethod(), panel = panel))
  }

#' @export
#'
#' @importFrom dplyr mutate
#'
mutate.tof_tbl <- function(.data, ...) {
  panel <- tof_get_panel(.data)
  return(new_tof_tibble(x = NextMethod(), panel = panel))
}

#' @export
#'
#' @importFrom dplyr slice_sample
slice_sample.tof_tbl <-
  function(.data, ..., n, prop, weight_by = NULL, replace = FALSE) {
    panel <- tof_get_panel(.data)
    return(new_tof_tibble(x = NextMethod(), panel = panel))
  }

# grouped_tof_tbl methods ------------------------------------------------------

## tidyr methods

#' @export
#'
#' @importFrom tidyr nest
nest.grouped_tof_tbl <- nest.tof_tbl

#' @export
#'
#' @importFrom dplyr ungroup
ungroup.grouped_tof_tbl <- function(x, ...) {
  panel <- tof_get_panel(x)
  return(new_tof_tibble(x = NextMethod(), panel = panel))
}
