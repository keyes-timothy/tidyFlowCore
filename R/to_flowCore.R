# to_flowCore.R
# This file contains methods for converting data frames and tibbles into
# flowCore native data structures (flowFrame's and flowSet's).

#' Coerce an object into a \code{\link[flowCore]{flowFrame}}
#'
#' @param x An object.
#'
#' @param ... Method-specific arguments
#'
#' @export
#'
#' @return A \code{\link[flowCore]{flowFrame}}
#'
#' @examples
#' NULL
#'
as_flowFrame <- function(x, ...) {
  UseMethod("as_flowFrame")
}

#' @export
as_flowFrame.data.frame <- function(x, ...) {
  result <- as_flowFrame.tof_tbl(x, ...)
  return(result)
}


#' @export
as_flowSet.data.frame <- function(x, group_cols, ...) {
  result <- as_flowFrame.tof_tbl(x, group_cols, ...)
  return(result)
}

#' Coerce a data.frame, tbl_df, or tof_tbl into a \code{\link[flowCore]{flowFrame}}
#'
#' @param x A data.frame, tbl_df, or tof_tbl.
#'
#' @param ... Unused.
#'
#' @return A \code{\link[flowCore]{flowFrame}}. Note that all non-numeric
#' columns in `x` will be removed.
#'
#' @rdname as_flowFrame
#'
#' @export
#'
#' @importFrom dplyr across
#' @importFrom dplyr everything
#' @importFrom dplyr rename_with
#' @importFrom dplyr select
#' @importFrom dplyr distinct
#' @importFrom dplyr where
#'
#' @importFrom flowCore flowFrame
#'
#' @importFrom purrr map
#'
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr pivot_wider
#'
#' @examples
#' NULL
as_flowFrame.tof_tbl <- function(x, ...) {
  tof_tibble <-
    x |>
    dplyr::select(dplyr::where(is.numeric))

  # find all columns in x that are character vectors
  tof_tibble_characters <-
    x |>
    dplyr::select(where(\(.x) !is.numeric(.x)))

  ## make a list of codes that can be used to convert all character vectors
  ## into numeric factors
  character_codes <-
    purrr::map(
      .x = tof_tibble_characters,
      .f = \(.x) data.frame(
        original_value = unique(.x),
        code = as.numeric(as.factor(unique(.x)))
      )
    )

  tof_tibble_characters <-
    tof_tibble_characters |>
    dplyr::mutate(
      dplyr::across(dplyr::everything(), \(.x) as.numeric(as.factor(.x)))
    )

  tof_tibble <-
    dplyr::bind_cols(tof_tibble, tof_tibble_characters)

  maxes_and_mins <-
    tof_tibble |>
    dplyr::summarize(
      dplyr::across(
        dplyr::everything(),
        .fns =
          list(max = ~ max(.x, na.rm = TRUE), min = ~ min(.x, na.rm = TRUE)),
        # use the many underscores because it's unlikely this will come up
        # in column names on their own
        .names = "{.col}_____{.fn}"
      )
    ) |>
    tidyr::pivot_longer(
      cols = dplyr::everything(),
      names_to = c("antigen", "value_type"),
      values_to = "value",
      names_sep = "_____"
    )  |>
    tidyr::pivot_wider(
      names_from = "value_type",
      values_from = "value"
    )

  # extract the names of all columns to used in the flowFrame
  data_cols <- maxes_and_mins$antigen

  # make the AnnotatedDataFrame flowCore needs
  parameters <-
    make_flowcore_annotated_data_frame(maxes_and_mins = maxes_and_mins)

  # assemble flowFrame
  result <-
    tof_tibble |>
    dplyr::rename_with(
      .fn = stringr::str_replace,
      pattern = "\\|",
      replacement = "_"
    ) |>
    as.matrix() |>
    flowCore::flowFrame(
      parameters = parameters
    )

  ## store codes for each character vector column in its own keyword slot
  if (length(character_codes) > 0) {
    for (i in seq_along(character_codes)) {
      keyword_name <- paste0(names(character_codes)[[i]], "_codes")
      flowCore::keyword(result)[[keyword_name]] <- character_codes[[i]]
    }
  }

  ## and store all the character vector codes together in a single keyword slot
  flowCore::keyword(result)[["character_codes"]] <- character_codes

  return(result)
}



#' Coerce an object into a \code{\link[flowCore]{flowSet}}
#'
#' @param x An object.
#'
#' @param ... Method-specific arguments
#'#'
#' @return A \code{\link[flowCore]{flowSet}}
#'
#' @export
#'
#' @examples
#' NULL
#'
as_flowSet <- function(x, ...) {
  UseMethod("as_flowSet")
}


#' Coerce a tof_tbl into a \code{\link[flowCore]{flowSet}}
#'
#' @param x A tof_tbl.
#'
#' @param group_cols Unquoted names of the columns in `x` that should
#' be used to group cells into separate \code{\link[flowCore]{flowFrame}}s.
#' Supports tidyselect helpers. Defaults to
#' NULL (all cells are written into a single \code{\link[flowCore]{flowFrame}}).
#' Note that the metadata column name "name" is a special value in the
#' \code{\link[flowCore]{flowSet}}) class, so if any of `group_cols` refers to
#' a column named "name," an error will be thrown.
#'
#' @param ... Unused.
#'
#' @return A \code{\link[flowCore]{flowSet}} in which cells are grouped into
#' constituent \code{\link[flowCore]{flowFrame}}s based on the values in
#' `group_cols`. If no `group_cols` are specified, a
#' \code{\link[flowCore]{flowFrame}} will be returned instead.
#' Note that all non-numeric columns in will be removed.
#'
#' @rdname as_flowSet
#'
#' @importFrom dplyr across
#' @importFrom dplyr everything
#' @importFrom dplyr mutate
#' @importFrom dplyr rename_with
#' @importFrom dplyr select
#' @importFrom dplyr summarize
#'
#' @importFrom flowCore flowFrame
#' @importFrom flowCore flowSet
#' @importFrom flowCore phenoData
#' @importFrom flowCore sampleNames
#' @importFrom flowCore keyword
#'
#' @importFrom purrr map
#'
#' @importFrom stringr str_replace
#'
#' @importFrom tidyr nest
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr pivot_wider
#'
#' @export
#'
#' @examples
#' NULL
as_flowSet.tof_tbl <- function(x, group_cols, ...) {

  if (missing(group_cols)) {
    result <- as_flowFrame(x)

  } else {
    tof_tibble <-
      x |>
      dplyr::select({{ group_cols }}, where(is.numeric))

    # throw an error if `group_cols` includes a column named "name", which
    # flowSets do not like.
    if (any(colnames(tof_tibble) == "name")) {
      stop(
        "x cannot have a column named `name`
           because this is a special column name in flowSets.
           Rename this column and try again. See help file for details."
      )
    }

    ## Converts all character vectors into numeric codes that can be
    ## stored in the constituent flowFrames of the flowSet
    tof_tibble_characters <-
      x |>
      dplyr::select(-{{ group_cols }}) |>
      dplyr::select(where(\(.x) !is.numeric(.x)))

    character_codes <-
      purrr::map(
        .x = tof_tibble_characters,
        .f = \(.x) data.frame(
          original_value = unique(.x),
          code = as.numeric(as.factor(unique(.x)))
        )
      )

    tof_tibble_characters <-
      tof_tibble_characters |>
      dplyr::mutate(
        dplyr::across(dplyr::everything(), \(.x) as.numeric(as.factor(.x)))
      )

    tof_tibble <-
      dplyr::bind_cols(tof_tibble, tof_tibble_characters)

    maxes_and_mins <-
      tof_tibble |>
      dplyr::summarize(
        dplyr::across(
          -{{group_cols}},
          .fns =
            list(
              max = ~ max(as.numeric(.x), na.rm = TRUE),
              min = ~ min(as.numeric(.x), na.rm = TRUE)
            ),
          # use the many underscores because it's unlikely this will come up
          # in column names on their own
          .names = "{.col}_____{.fn}"
        )
      ) |>
      tidyr::pivot_longer(
        cols = dplyr::everything(),
        names_to = c("antigen", "value_type"),
        values_to = "value",
        names_sep = "_____"
      )  |>
      tidyr::pivot_wider(
        names_from = "value_type",
        values_from = "value"
      )

    # extract the names of all non-grouping columns to be saved to the .fcs file
    data_cols <- maxes_and_mins$antigen

    tof_tibble <-
      tof_tibble |>
      tidyr::nest(.by = {{ group_cols }})

    # make the AnnotatedDataFrame flowCore needs
    parameters <-
      make_flowcore_annotated_data_frame(maxes_and_mins = maxes_and_mins)

    tof_tibble <-
      tof_tibble |>
      dplyr::mutate(
        flowFrames =
          purrr::map(
            .x = .data$data,
            ~ flowCore::flowFrame(
              exprs =
                # have to change any instances of "|" in column names to another
                # separator, as "|" has special meaning as an .fcs file delimiter
                as.matrix(
                  dplyr::rename_with(
                    .x,
                    stringr::str_replace,
                    pattern = "\\|",
                    replacement = "_"
                  )
                ),
              parameters = parameters
            )
          )
      ) |>
      dplyr::select(
        {{ group_cols }},
        "flowFrames"
      )

    metadata_frame <-
      tof_tibble |>
      dplyr::select({{ group_cols }}) |>
      as.data.frame()

    # store group_cols metadata in a data.frame for the flowSet
    row.names(metadata_frame) <-
      paste0('tidyFlowCore_id_', seq_len(nrow(metadata_frame)))
    metadata_frame$.tidyFlowCore_unique_identifier <- row.names(metadata_frame)

    result <- flowCore::flowSet(tof_tibble$flowFrames)

    # add the grouping columns' information to the flowSet as metadata
    # access this using flowCore::pData()
    flowCore::sampleNames(result) <-
      paste0('tidyFlowCore_id_', seq_along(result))
    flowCore::pData(result)$name <- flowCore::sampleNames(result)
    flowCore::pData(result)$.tidyFlowCore_unique_identifier <-
      flowCore::sampleNames(result)

    result_pData <-
      flowCore::pData(result) |>
      dplyr::left_join(metadata_frame, by = ".tidyFlowCore_unique_identifier")
    row.names(result_pData) <- result_pData$.tidyFlowCore_unique_identifier
    flowCore::pData(result) <- result_pData

    if (length(character_codes) > 0) {
      for (j in seq_along(result)) {
        for (i in seq_along(character_codes)) {
          keyword_name <- paste0(names(character_codes)[[i]], "_codes")
          flowCore::keyword(result[[j]])[[keyword_name]] <- character_codes[[i]]
        }
        flowCore::keyword(result[[j]])[["character_codes"]] <- character_codes
      }
    }

  }

  return(result)
}




#' Make the AnnotatedDataFrame needed for the flowFrame class
#'
#' @param maxes_and_mins a data.frame containing information about the max
#' and min values of each channel to be saved in the flowFrame.
#'
#' @return An AnnotatedDataFrame.
#'
#' @importFrom dplyr transmute
#'
#' @importFrom methods new
#'
#' @importFrom stringr str_replace
#' @importFrom stringr str_c
#'
#' @examples
#' NULL
make_flowcore_annotated_data_frame <- function(maxes_and_mins) {
  fcs_varMetadata <-
    data.frame(
      labelDescription =
        c(
          "Name of Parameter",
          "Description of Parameter",
          "Range of Parameter",
          "Minimum Parameter Value after Transformation",
          "Maximum Parameter Value after Transformation"
        )
    )

  fcs_data <-
    maxes_and_mins |>
    dplyr::transmute(
      # have to change any instances of "|" in column names to another
      # separator, as "|" has special meaning as an .fcs file delimiter
      name = stringr::str_replace(.data$antigen, "\\|", "_"),
      desc = stringr::str_replace(.data$antigen, "\\|", "_"),
      range = max - min,
      minRange = min,
      maxRange = max
    ) |>
    as.data.frame()

  row.names(fcs_data) <- stringr::str_c("$", "P", seq_len(nrow(fcs_data)))

  # make the AnnotatedDataFrame
  parameters <-
    methods::new(
      "AnnotatedDataFrame",
      data = fcs_data,
      varMetadata = fcs_varMetadata
    )

  return(parameters)
}

