#' @title Trimming specification object
#'
#' @description
#' Specify trimming parameters for \code{transportIP} and \code{transportTADA}.
#'
#' @param value Threshold for trimming, which should be between 0 and 0.5.
#' 
#' @details
#' Here, trimming refers to excluding observations with extreme weights from the analysis. This package only supports Crump trimming (cite), which trims based on propensity scores rather than weights. Thus, trimming is not supported for weighting components that do not have a well-defined propensity score (\code{final} in \code{transportIP}; \code{participation} and \code{final} in \code{transportTADA}).
#' 
#'
#' @return
#' A \code{trim} object that can be processed by \code{transportIP} and \code{transportTADA} for weight trimming.
#' 
#' @export
trim <- function(value) {
  trimObj <- list(value = value)
  class(trimObj) <- "trim"
  return(trimObj)
}

trimInd <- function(weightGLM, trim) {
  score <- weightGLM$fitted.values
  return(as.integer(score >= trim$value & score <= 1 - trim$value))
}

is.trim <- function(trim) {
  if (!all(names(trim) %in% c("value"))) return(F)
  return(inherits(trim, "trim") &
           is.numeric(trim$value))
}