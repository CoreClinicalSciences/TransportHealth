#' @title
#' Transportability analysis using interpolated g-computation
#'
#' @description
#' Estimates the average treatment effect (ATE) using interpolated g-computation in a generalizability or transportability analysis. In particular, the estimators should be unbiased for the ATE in the superpopulation or the target population, respectively.
#'
#' @param link Defaults to \code{"identity"}, which corresponds to absolute treatment effects for continuous responses. The \code{"log"} option accommodates relative treatment effects such as relative risk, odds ratio and hazard ratio.
#' @param effectModifiers Vector of strings indicating effect modifiers to adjust for
#' @param mainTreatmentEffect Estimate of ATE in original study
#' @param mainSE Estimate of standard error of estimator of ATE in original study
#' @param subgroupTreatmentEffects Estimates of subgroup ATEs in original study. Please provide subgroup ATEs in the order of effect modifiers listed in \code{effectModifiers}, and provide the ATE of the subgroup whose proportion is provided in \code{summaryAggregateData} first in each pair
#' @param subgroupSEs Estimates of standard errors of subgroup ATEs in original study. Please provide SEs in the order of effect modifiers listed in \code{effectModifiers}, and provide the SE of the subgroup whose proportion is provided in \code{summaryAggregateData} first in each pair
#' @param corrStructure Correlation structure of dichotomized effect modifiers. If target IPD is provided, this will be estimated from the target data, if user input is omitted. If target aggregate data is provided, this will be specified by the user and default to an independent correlation structure if left unspecified.
#' @param studySampleSize Sample size of original study
#' @param aggregateStudyData Vector of proportions of dichotomized effect modifiers in study data. Please provide proportions of only one category for each effect modifier. This category should correspond to the the first ATE and SE provided for each effect modifier.
#' @param targetData May be IPD or aggregate. If aggregate, provide proportions of only one category of dichotomized effect modifiers in a named vector (not a data frame)
#'
#' @details
#' This function transports the ATE estimate from the original study data to the target data by utilizing subgroup effect estimates in the same way as network meta-interpolation (cite Ofir's paper). As standard errors are transported manually, no bootstrapping is done, unlike other methods supported by \code{TransportHealthR}.
#' 
#' @return
#' A \code{transportInterpolated} object with the following components:
#' * \code{effect}: Transported ATE estimate
#' * \code{link}: Denotes the link function used. Absolute treatment effects (i.e. for continuous outcomes) correspond to \code{"identity"}, while relative treatment effects correspond to \code{"log"}
#' * \code{var}: Transported variance estimate of effect estimate
#' * \code{effectModifiers}: Vector of strings indicating effect modifiers adjusted for
#' * \code{mainTreatmentEffect}: Estimate of ATE in original study
#' * \code{mainSE}: Standard error of estimator of ATE in original study
#' * \code{subgroupTreatmentEffects}: Estimates of subgroup ATEs in original study, as provided
#' * \code{subgroupSEs}: Estimates of standard errors of subgroup ATEs in original study, as provided
#' * \code{corrStructure}: Correlation structure of effect modifiers used in analysis
#' * \code{studySampleSize}: Sample size of original study
#' * \code{aggregateStudyData}: Aggregate-level study data, as provided
#' * \code{targetData}: Target data, as provided
#' 
#' @export
#' 
#' @md
transportInterpolated <- function (link = c("identity", "log"),
                           effectModifiers,
                           mainTreatmentEffect,
                           mainSE,
                           subgroupTreatmentEffects,
                           subgroupSEs,
                           corrStructure = NULL,
                           studySampleSize,
                           aggregateStudyData,
                           targetData) {
  if (length(subgroupTreatmentEffects) != length(effectModifiers) * 2 |
      length(subgroupSEs) != length(effectModifiers) * 2 |
      length(effectModifiers) != length(aggregateStudyData))
    stop("Incongruent lengths of list of effect modifiers and list of subgroup effects/SEs. Make sure that all effect modifiers are dichotomized. Provide treatment effects for each marginal subgroup.")
  
  link <- match.arg(link, c("identity", "log"))
  
  # Obtain treatment effects and SEs on the linear predictor scale
  if (link == "log") {
    lpMainTreatmentEffect  <- log(mainTreatmentEffect)
    lpMainSE <- mainSE / lpMainTreatmentEffect
    lpSubgroupTreatmentEffects <- log(subgroupTreatmentEffects)
    lpSubgroupSEs <- subgroupSEs / lpSubgroupTreatmentEffects
  } else {
    lpMainTreatmentEffect <- mainTreatmentEffect
    lpMainSE <- mainSE
    lpSubgroupTreatmentEffects <- subgroupTreatmentEffects
    lpSubgroupSEs <- subgroupSEs
  }
  
  m <- length(effectModifiers)
  
  if (is.data.frame(targetData)) {
    # Calculate aggregate statistics of target data if IPD
    if (!all(effectModifiers %in% names(targetData)))
      stop("Please provide effect modifier names as in target data.")
    numericTargetData <- apply(as.matrix(targetData[, effectModifiers]), 2, as.numeric)
    emTargetProps <- apply(targetData[, effectModifiers], 2, function(x) mean(as.numeric(x)))
    if (is.null(corrStructure)) corrStructure <- stats::cor(numericTargetData)
  } else {
    # Extract aggregate statistics of target data if already aggregate
    if (!all(effectModifiers %in% names(targetData)))
      stop("Please provide effect modifier names as in target data.")
    emTargetProps <- targetData[effectModifiers]
    if (is.null(corrStructure)) {
      warning("Correlation structure not provided while aggregate target data is provided. Defaulting to independent correlation structure.")
      corrStructure <- diag(m)
    }
  }
  
  emStudyVars <- (aggregateStudyData * (1 - aggregateStudyData)) / studySampleSize
  
  # Construct enriched data matrix for treatment effect
  enrichedTEMatrix <- matrix(double((2*m + 1) * (m + 1)), ncol = m + 1)
  enrichedTEMatrix[, 1] <- 1 # Intercept column
  enrichedTEMatrix[1, (2:(m+1))] <- aggregateStudyData # Overall TE row
  for (j in 1:m) enrichedTEMatrix[2 * j, 1 + j] <- 1 # Filling in 1s at corresponding rows of subgroup treatment effect
  for (i in 1:m) { # BLUP
    enrichedTEMatrix[2*i, - (i + 1)] <- corrStructure[i, -i] * sqrt(emStudyVars[i]) / sqrt(emStudyVars[-i]) * (1 - aggregateStudyData[-i]) + aggregateStudyData[-i]
    enrichedTEMatrix[2*i + 1, - (i + 1)] <- corrStructure[i, -i] * sqrt(emStudyVars[i]) / sqrt(emStudyVars[-i]) * (0 - aggregateStudyData[-i]) + aggregateStudyData[-i]
  }
  
  # Construct enriched data matrix for SE
  enrichedSEMatrix <- matrix(double((2*m + 1) * ((m^2 + 3*m + 2) / 2)), nrow = 2*m + 1)
  enrichedSEMatrix[, 1:(m + 1)] <- enrichedTEMatrix^2 # Square columns
  enrichedSEMatrix[, (m + 2):(2*m + 1)] <- 2 * enrichedTEMatrix[, -1] # Double columns
  for (j in 1:(m-1)) {
    for (i in ((j+1):m)) { # Entrywise products of columns
      enrichedSEMatrix[, 2*m + 1 + (m * (m - 1))/2 - ((m - j) * (m - j + 1))/2 + (i - j)] <- 2 * enrichedTEMatrix[, j] * enrichedTEMatrix[, i]
    }
  }
  
  # Construct corresponding EM configuration for target data
  emTargetVarsAllTerms <- double((m^2 + 3*m + 2) / 2)
  emTargetVarsAllTerms[1:(m + 1)] <- c(1, emTargetProps)^2
  emTargetVarsAllTerms[(m + 2):(2*m + 1)] <- 2 * emTargetProps
  for (j in 1:(m-1)) {
    for (i in ((j+1):m)) {
      emTargetVarsAllTerms[2*m + 1 + (m * (m - 1))/2 - ((m - j) * (m - j + 1))/2 + (i - j)] <- 2 * emTargetProps[j] * emTargetProps[i]
    }
  }
  
  # Obtain coefficient estimates for treatment effect
  betaHat <- solve(t(enrichedTEMatrix) %*% enrichedTEMatrix) %*% t(enrichedTEMatrix) %*% c(lpMainTreatmentEffect, lpSubgroupTreatmentEffects)
  
  # Obtain coefficient estimates for SE
  sigmaHat <- t(enrichedSEMatrix) %*% MASS::ginv(enrichedSEMatrix %*% t(enrichedSEMatrix)) %*% c(lpMainSE^2, lpSubgroupSEs^2)
  
  # Transported effects on linear predictor scale
  lpEffect <- t(betaHat) %*% c(1, emTargetProps)
  lpVar <- t(sigmaHat) %*% emTargetVarsAllTerms
  
  # Obtain effects on original scale
  if (link == "log") {
    effect <- exp(lpEffect)
    var <- lpVar * lpEffect
  } else {
    effect <- lpEffect
    var <- lpVar
  }
  
  result <- list(effect = effect,
                 link = link,
                 var = var,
                 effectModifiers = effectModifiers,
                 mainTreatmentEffect = mainTreatmentEffect,
                 mainSE = mainSE,
                 subgroupTreatmentEffects = subgroupTreatmentEffects,
                 subgroupSEs = subgroupSEs,
                 corrStructure = corrStructure,
                 studySampleSize = studySampleSize,
                 aggregateStudyData = aggregateStudyData,
                 targetData = targetData)
  
  class(result) <- "transportInterpolated"
  
  return(result)
}

#' @title Summarize results of transportability analysis using interpolated g-computation
#' 
#' @description
#' Returns summary object which contains the transported effect estimate, organizes subgroup effect estimates into a data frame and produces summary statistics of effect modifiers in the target data, if not already provided.
#'
#' @rdname summary.transportInterpolated
#'
#' @param object Result from \code{transportInterpolated} function
#' @param ... Further arguments from previous function or to pass to next function
#'
#' @return
#' Returns a \code{summary.transportInterpolated} object containing the following components:
#' * \code{effect}, \code{link}, \code{mainTreatmentEffect} and \code{mainSE} as in the \code{transportInterpolated} object
#' * \code{se}: Estimated standard error of the effect estimate
#' * \code{subgroupEffects}: Subgroup effect estimates and their estimated standard errors in the original study, organized into a data frame with specific effect modifier and subgroup information
#' * \code{aggregateStudyData}: Summary statistics of effect modifiers in study data
#' * \code{aggregateTargetData}: Summary statistics of effect modifiers in target data
#' 
#' @export
#' 
#' @md
summary.transportInterpolated <- function (object, ...) {
  transportInterpolatedResult <- object
  
  effect <- transportInterpolatedResult$effect
  effectModifiers <- transportInterpolatedResult$effectModifiers
  # Organize subgroup analyses into table
  subgroupStudyTable <- data.frame(effectModifier = rep(effectModifiers, each = 2),
                                  subgroup = rep(c(1,0), times = length(effectModifiers)),
                                  effect = transportInterpolatedResult$subgroupTreatmentEffects,
                                  se = transportInterpolatedResult$subgroupSEs)
  
  aggregateStudyData <- transportInterpolatedResult$aggregateStudyData
  names(aggregateStudyData) <- effectModifiers
  
  m <- length(effectModifiers)
  targetData <- transportInterpolatedResult$targetData
  
  if (is.data.frame(targetData)) {
    # Calculate aggregate statistics of target data if IPD
    emTargetProps <- apply(targetData[, effectModifiers], 2, function(x) mean(as.numeric(x)))
  } else {
    # Extract aggregate statistics of target data if already aggregate
    emTargetProps <- targetData[effectModifiers]
  }
  
  names(emTargetProps) <- effectModifiers
  
  summaryResult <- list(effect = effect,
                        link = transportInterpolatedResult$link,
                        se = sqrt(transportInterpolatedResult$var),
                        mainTreatmentEffect = transportInterpolatedResult$mainTreatmentEffect,
                        mainSE = transportInterpolatedResult$mainSE,
                        subgroupEffects = subgroupStudyTable,
                        aggregateStudyData = aggregateStudyData,
                        aggregateTargetData = emTargetProps)
  
  class(summaryResult) <- "summary.transportInterpolated"
  
  return(summaryResult)
}

#' @rdname summary.transportInterpolated
#'
#' @param x \code{summary.transportInterpolated} object
#' @param out Output stream
#' @param ... Further arguments from previous function or to pass to next function
#'
#' @export
print.summary.transportInterpolated <- function (x, out = stdout(), ...) {
  summaryResult <- x
  
  write(paste0("Transported ATE: ", summaryResult$effect), out)
  write(paste0("Standard error: ", summaryResult$se), out)
  write(paste0("Link function: ", summaryResult$link), out)
  write(paste0("Source study treatment effect: ", summaryResult$mainTreatmentEffect), out)
  write(paste0("Source study standard error: ", summaryResult$mainSE), out)
  write("Subgroup source treatment effects: ", out)
  print(summaryResult$subgroupEffects, out)
  write("Source data summary: ", out)
  print(summaryResult$aggregateStudyData, out)
  write("Target data summary: ", out)
  print(summaryResult$aggregateTargetData, out)
}

#' @title Visually represent results of transportability analysis using interpolated g-computation
#'
#' @description
#' This function is a wrapper for \code{modelsummary::modelplot} to plot effect estimates in a transportability analysis using interpolated g-computation. Note that the transported variance estimate is used in this function.
#'
#' @param x Result from \code{transportInterpolated} function
#' @param ... Further arguments from previous function or to pass to next function
#'
#' @return
#' A \code{ggplot} object showing the estimates and confidence intervals of the transported effect estimate.
#' 
#' @export
#'
plot.transportInterpolated <- function (x, ...) {
  transportInterpolatedResult <- x
  ti <- data.frame(term = "treatment",
                   estimate = transportInterpolatedResult$effect,
                   conf.low = transportInterpolatedResult$effect - 1.96 * sqrt(transportInterpolatedResult$var),
                   conf.high = transportInterpolatedResult$effect + 1.96 * sqrt(transportInterpolatedResult$var))
  
  gl <- data.frame(n = nrow(transportInterpolatedResult$targetData))
  
  ms <- list(tidy = ti, glance = gl)
  class(ms) <- "modelsummary_list"
  
  resultPlot <- modelsummary::modelplot(ms)
  
  return(resultPlot)
}

#' @title Check validity of interpolated g-computation result object
#' 
#' @description
#' A simple helper function that validates whether the components of the given \code{transportInterpolated} object are of the correct types.
#'
#' @param x Result object from \code{transportInterpolated} function
#'
#' @return A boolean indicating whether all components of \code{transportInterpolated} object have the correct types.
#' 
#' @export
is.transportInterpolated <- function (x) {
  return(is.numeric(x$effect) &
          is.numeric(x$var) &
          x$link %in% c("identity", "log") &
          is.character(x$effectModifiers) &
          is.numeric(x$mainTreatmentEffect) &
          is.numeric(x$mainSE) &
          is.numeric(x$subgroupTreatmentEffects) &
          is.numeric(x$subgroupSEs) &
          is.matrix(x$corrStructure) &
          is.numeric(x$studySampleSize) &
          is.numeric(x$aggregateStudyData) &
          (is.numeric(x$targetData) | is.data.frame(x$targetData)))
}
