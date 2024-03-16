
#' A character vector of CyTOF metal name patterns supported by tidyFlowCore
#'
#' A character vector used by `tof_find_panel_info` to detect and
#' parse which CyTOF metals correspond to each channel in an input .fcs file.
#'
#' @usage data(metal_masterlist)
#'
#' @format A character vector in which each entry is a pattern that tidyFlowCore searches
#' for in every CyTOF channel in input .fcs files. These patterns are an amalgamate
#' of example .fcs files sampled from the studies linked below.
#'
#' @returns None
#'
#' @source \url{https://github.com/kara-davis-lab/DDPR}
#' \url{https://cytobank.org/nolanlab/reports/Levine2015.html}
#' \url{https://cytobank.org/nolanlab/reports/Spitzer2015.html}
#' \url{https://cytobank.org/nolanlab/reports/Spitzer2017.html}
#' \url{https://community.cytobank.org/cytobank/projects/609}
"metal_masterlist"
