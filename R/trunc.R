#' @title Truncation specification object
#'
#' @description
#' Specify truncation parameters for \code{transportIP} and \code{transportTADA}.
#'
#' @param type One of \code{"quantile"} or \code{"raw"}. If \code{"quantile"}, the truncation threshold is the quantile of the weights specified by \code{value}. If \code{"raw"}, the truncation threshold is \code{value}.
#' @param value Value of threshold as described above. When \code{type = "quantile"}, this argument should be between 0 and 1.
#'
#' @details
#' Here, truncation refers to the procedure of setting weights above a threshold equal to the threshold. Thus, observations with extreme weights are included in the analysis, albeit with reduced weights.
#'
#' @return
#' A \code{trunc} object that can be processed by \code{transportIP} and \code{transportTADA} for weight truncation.
#' 
#' @export
trunc <- function(type = c("quantile", "raw"),
                  value) {
  type <- match.arg(type, c("quantile", "raw"))
  truncObj <- list(type = type,
                   value = value)
  class(truncObj) <- "trunc"
  return(truncObj)
}

truncWeights <- function(weights,
                         trunc) {
  if (!inherits(trunc, "trunc")) stop("Please provide trunc object.")
  
  if (trunc$type == "raw") limit <- trunc$value
  else if (trunc$type == "quantile") limit <- stats::quantile(weights, probs = trunc$value)
  
  truncatedWeights <- pmin(weights, limit)
  
  return(truncatedWeights)
}

is.trunc <- function(trunc) {
  if (!all(names(trunc) %in% c("quantile", "raw"))) return(F)
  return(inherits(trunc, "trunc") &
           trunc$type %in% c("quantile", "raw") &
           is.numeric(trunc$value))
}