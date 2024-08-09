#' @title
#' Transportability analysis using interpolated g-computation
#'
#' @description
#' Estimates the average treatment effect (ATE) using interpolated g-computation in a generalizability or transportability analysis. In particular, the estimators should be unbiased for the ATE in the superpopulation or the target population, respectively.
#'
#' @param link Defaults to \code{"identity"}, which corresponds to absolute treatment effects for continuous responses. The \code{"log"} option accommodates relative treatment effects such as relative risk, odds ratio and hazard ratio.
#' @param effectModifiers Vector of strings indicating effect modifiers to adjust for
#' @param mainTreatmentEffect Estimate of ATE in original study
#' @param mainSE Standard error of estimator of ATE in original study
#' @param subgroupTreatmentEffects Estimates of subgroup ATEs in original study. Please provide subgroup ATEs in the order of effect modifiers listed in \code{effectModifiers}, and provide the ATE of the subgroup whose proportion is provided in \code{targetData} first in each pair
#' @param subgroupSEs 
#' @param corrStructure 
#' @param studySampleSize 
#' @param targetData 
#'
#' @return
#' @export
#'
#' @examples
transportInterpolated <- function (link = c("identity", "log"),
                           effectModifiers,
                           mainTreatmentEffect,
                           mainSE,
                           subgroupTreatmentEffects,
                           subgroupSEs,
                           corrStructure = NULL,
                           studySampleSize,
                           targetData) {
  if (length(subgroupTreatmentEffects) != length(effectModifiers) * 2 |
      length(subgroupSEs) != length(effectModifiers) * 2)
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
  }
  
  m <- length(effectModifiers)
  
  if (nrow(targetData) > 1) {
    # Calculate aggregate statistics of target data if IPD
    if (!all(effectModifiers %in% names(targetData)))
      stop("Please provide effect modifier names as in target data.")
    emTargetProps <- apply(targetData[[effectModifiers]], 2, mean)
    corrStructure <- cor(targetData[[effectModifiers]])
  } else {
    # Extract aggregate statistics of target data if already aggregate
    if (!all(effectModifiers %in% names(targetData)))
      stop("Please provide effect modifier names as in target data.")
    emTargetProps <- targetData[[effectModifiers]]
    if (is.null(corrStructure)) {
      warning("Correlation structure not provided while aggregate target data is provided. Defaulting to independent correlation structure.")
      corrStructure <- diag(m)
    }
  }
  
  emTargetVars <- (emTargetProps * (1 - emTargetProps)) / studySampleSize
  
  # Construct enriched data matrix for treatment effect
  enrichedTEMatrix <- matrix(double((2*m + 1) * (m + 1)), ncol = m + 1)
  enrichedTEMatrix[, 1] <- 1 # Intercept column
  enrichedTEMatrix[1, (2:(m+1))] <- emTargetProps # Overall TE row
  for (j in 1:m) enrichedTEMatrix[1 + j, 2 * j] <- 1 # Filling in 1s at corresponding rows of subgroup treatment effect
  for (i in 1:m) { # BLUP
    enrichedTEMatrix[2*i, - (i + 1)] <- corrStructure[i, -i] * sqrt(emTargetVars[i]) / sqrt(emTargetVars[-i]) * (1 - emTargetProps[-i]) + emTargetProps[-i]
    enrichedTEMatrix[2*i + 1, - (i + 1)] <- corrStructure[i, -i] * sqrt(emTargetVars[i]) / sqrt(emTargetVars[-i]) * (0 - emTargetProps[-i]) + emTargetProps[-i]
  }
  
  # Construct enriched data matrix for SE
  enrichedSEMatrix <- matrix(double((2*m + 1) * ((m^2 + 3*m + 2) / 2)), nrow = 2*m + 1)
  enrichedSEMatrix[, m + 1] <- enrichedTEMatrix^2 # Square columns
  enrichedSEMatrix[, (m + 2):(2*m + 1)] <- 2 * enrichedTEMatrix[, -1] # Double columns
  for (j in 1:(m-1)) {
    for (i in ((j+1):m)) { # Entrywise products of columns
      enrichedSEMatrix[, 2*m + 1 + (m * (m - 1))/2 - ((m - j) * (m - j + 1))/2 + (i - j)] <- 2 * enrichedTEMatrix[, j] * enrichedTEMatrix[, i]
    }
  }
  
  # Construct corresponding EM configuration for target data
  emTargetVarsAllTerms <- double((m^2 + 3*m + 2) / 2)
  emTargetVarsAllTerms[, m + 1] <- c(1, emTargetProps)^2
  emTargetVarsAllTerms[, (m + 2):(2*m + 1)] <- 2 * emTargetProps
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
  lpEffect <- betaHat %*% c(1, emTargetProps)
  lpVar <- sigmaHat %*% emTargetVarsAllTerms
  
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
                 targetData = targetData)
  
  class(result) <- "transportInterpolated"
  
  return(result)
}

#' Title
#'
#' @param object 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
summary.transportInterpolated <- function (object, ...) {
  transportInterpolatedResult <- object
  
  effect <- transportInterpolatedResult$effect
  effectModifiers <- transportInterpolatedResult$effectModifiers
  # Organize subgroup analyses into table
  subgroupStudyTable <- data.frame(effectModifier = rep(effectModifiers, each = 2),
                                  subgroup = rep(c(1,0), times = length(effectModifiers)),
                                  effect = transportInterpolatedResult$subgroupTreatmentEffects,
                                  se = transportInterpolatedResult$subgroupSEs)
  
  m <- length(effectModifiers)
  targetData <- transportInterpolatedResult$targetData
  
  if (nrow(targetData) > 1) {
    # Calculate aggregate statistics of target data if IPD
    emTargetProps <- apply(targetData[[effectModifiers]], 2, mean)
  } else {
    # Extract aggregate statistics of target data if already aggregate
    emTargetProps <- targetData[[effectModifiers]]
  }
  
  names(emTargetProps) <- effectModifiers
  
  summaryResult <- list(effect = effect,
                        link = transportInterpolatedResult$link,
                        se = sqrt(transportInterpolatedResult$var),
                        sourceEffect = transportInterpolatedResult$mainTreatmentEffect,
                        sourceSE = transportInterpolatedResult$mainSE,
                        subgroupEffects = subgroupStudyTable,
                        aggregateTargetData = emTargetProps)
}

#' Title
#'
#' @param x 
#' @param out 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
print.summary.transportInterpolated <- function (x, out, ...) {
  summaryResult <- x
  
  write(paste0("Transported ATE: ", summaryResult$effect), out)
  write(paste0("Standard error: ", summaryResult$se), out)
  write(paste0("Link function: ", summaryResult$log), out)
  write(paste0("Source study treatment effect: ", summaryResult$sourceEffect), out)
  write(paste0("Source study standard error: ", summaryResult$sourceSE), out)
  write("Subgroup source treatment effects: ", out)
  print(summaryResult$subgroupEffects, out)
  write("target data summay: ", out)
  print(summaryResult$aggregateTargetData, out)
}

#' Title
#'
#' @param x 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
plot.transportInterpolated <- function(x, ...) {
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