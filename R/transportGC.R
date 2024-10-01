#' @title Transportability analysis using g-computation
#' 
#' @description
#' Estimates the coefficients of a marginal structural model (MSM) using g-computation in a generalizability or transportability analysis. In particular, the estimators should be unbiased for the coefficients in the superpopulation or the target population, respectively.
#' 
#' @param effectType Type of effect desired for the ATE: \code{"meanDiff"} for mean difference, \code{"rr"} for relative risk, \code{"or"} for odds ratio, and \code{"hr"} for hazard ratio.
#' @param preparedModel A \code{transportGCPreparedModel} object. This is obtained by using the \code{transportGCPreparedModel} function to fit an outcome model using the study data.
#' @param targetData A target dataset.
#' @param bootstrapNum Number of bootstrap datasets to simulate to obtain robust variance estimators
#' 
#' @details
#' The expected workflow is as follows:
#' 
#' \enumerate{
#'  \item{A researcher who wants to perform a generalizability/transportability analysis collects data from the target population.}
#'  \item{They then request the owner of the study data from which they wish to generalize/transport to provide an outcome model fitted using the study data.}
#'  \item{The owner of the study data runs the \code{transportGCPreparedModel} function on the study data to obtain a \code{transportGCPreparedModel} object which contains the fitted outcome model.}
#'  \item{The owner of the study data provides the \code{transportGCPreparedModel} object to the researcher, perhaps in a \code{.rds} file.}
#'  \item{The researcher uses this function and the provided \code{transportGCPreparedModel} object to perform the analysis using g-computation.}
#' }
#' 
#' Since model-fitting objects in \code{R} often contain the data used to fit the model, the \code{transportGCPreparedModel} function wipes this data in the model-fitting object and keeps additional information about the name of the response variable, the name of the treatment variable and the levels of treatment. This is to comply with government regulations regarding access and integration of data sources from different countries.
#' 
#' The MSM-fitting functions do not provide correct standard errors as-is. Bootstrap is used to calculate robust variance estimates of the MSM coefficient estimators. Note that these standard errors are only valid conditional on the observed study data because it is not possible to resample the study data when access to it is restricted.
#'
#' @return
#' A \code{transportGC} object with the following components:
#' * \code{effect}: Calculated ATE
#' * \code{effectType}: Type of effect calculated
#' * \code{var}: Estimated variance of ATE estimator, calculated using bootstrap
#' * \code{preparedModel}: The \code{transportGCPreparedModel} object used to estimate the ATE
#' * \code{bootstrapNum}: Integer indicating number of bootstrap datasets simulated to calculate robust variance estimators.
#' 
#' @export
#' 
#' @md
transportGC <- function (effectType = c("meanDiff", "rr", "or", "hr"),
                         preparedModel,
                         targetData, bootstrapNum = 500) {
  transportGCResult <- transportGCFit(effectType,
                                      preparedModel,
                                      targetData)
  
  if (!preparedModel$wipe) {
    studyData <- preparedModel$studyData
    treatmentLevels <- preparedModel$treatmentLevels
    treatmentGroupData <- list()
    treatmentVector <- as.character(studyData[[preparedModel$treatment]])
    for (level in treatmentLevels) {
      treatmentGroupData[[level]] <- studyData[treatmentVector == level, ]
    }
    nLevels <- sapply(treatmentGroupData, nrow)
    names(nLevels) <- treatmentLevels
  }
  
  nTarget <- nrow(targetData)
  
  bootstrapEstimates <- sapply(1:bootstrapNum,
                                 function (x) {
                                   if (preparedModel$wipe) {
                                    targetBoot <- targetData[sample.int(n = nTarget, replace = T), ]
                                   
                                    resultBoot <- transportGCFit(effectType,
                                                                preparedModel,
                                                                targetData = targetBoot)
                                   } else {
                                      treatmentGroupBoot <- list()
                                      for (level in treatmentLevels) {
                                        nSample <- nLevels[level]
                                        treatmentGroupBoot[[level]] <- treatmentGroupData[[level]][sample.int(nSample, replace = T), ]
                                      }
                                      for (level in treatmentLevels) {
                                        if (!exists("studyBoot")) studyBoot <- treatmentGroupBoot[[level]]
                                        else studyBoot <- rbind(studyBoot, treatmentGroupBoot[[level]])
                                      }
                                     
                                      targetBoot <- targetData[sample.int(n = nTarget, replace = T), ]
                                     
                                      preparedBoot <- suppressWarnings(transportGCPreparedModel(preparedModel$formula,
                                                                               response = preparedModel$response,
                                                                               treatment = preparedModel$treatment,
                                                                               treatmentLevels = preparedModel$treatmentLevels,
                                                                               family = preparedModel$family,
                                                                               method = preparedModel$method,
                                                                               studyData = studyBoot,
                                                                               wipe = F))
                                      resultBoot <- transportGCFit(effectType,
                                                                  preparedBoot,
                                                                  targetData = targetBoot)
                                   }
                                   
                                   return(resultBoot$effect)
                                 })
  
  # varMatrix <- stats::var(bootstrapEstimates)
  # colnames(varMatrix) <- rownames(varMatrix) <- c(names(transportGCResult$msm$coefficients), names(transportGCResult$msm$zeta))
  varMatrix <- stats::var(bootstrapEstimates)
  # if (inherits(preparedModel$outcomeModel, "polr")) {
  #   names(transportGCResult$effect) <- colnames(varMatrix) <- rownames(varMatrix) <- preparedModel$responseLevels
  # }
  transportGCResult$var <- varMatrix
  
  transportGCResult$bootstrapNum <- bootstrapNum
  
  return(transportGCResult)
}

transportGCFit <- function (effectType = c("meanDiff", "rr", "or", "hr"),
                         preparedModel,
                         targetData) {
  # Error checking of consistency from preparedModel object
  if (!is.transportGCPreparedModel(preparedModel)) stop("Please provide a transportGCPreparedModel object.")
  
  # Extract names of response and treatment variables as well as outcome model object
  response <- preparedModel$response
  treatment <- preparedModel$treatment
  outcomeModel <- preparedModel$outcomeModel
  treatmentLevels <- preparedModel$treatmentLevels
  responseLevels <- preparedModel$responseLevels
  
  # Create data frames which are the target dataset, but with the treatment column all equal to 0 or 1
  targetDataReferenceList <- list()
  for (level in treatmentLevels) {
    targetDataReferenceList[[level]] <- targetData
    targetDataReferenceList[[level]][[treatment]] <- level
  }
  
  # Combine the above data frames together. Counterfactuals are predicted with this data frame
  targetDataCounterfactualFrame <- Reduce(function(d1, d2) rbind(d1, d2), targetDataReferenceList)
  targetDataCounterfactualFrame[[treatment]] <- as.factor(targetDataCounterfactualFrame[[treatment]])
  
  # Predict counterfactual outcomes
  # TODO: adapt to coxph since it doesn't give survival time predictions - see Lee (2024) (transporting survival of hiv...) and Chen & Tsiatis (2001)
  if (!inherits(outcomeModel, "coxph") & !inherits(outcomeModel, "polr")) {
    targetDataCounterfactualFrame[[response]] <- stats::predict(outcomeModel,
                                                       newdata = targetDataCounterfactualFrame,
                                                       type = "response")
  } else if (inherits(outcomeModel, "polr")) {
    # predict.polr with type = "probs" returns a matrix. It's done separately to avoid dealing with column names. We're assuming the user has already set the order of levels of the response to what they want.
    # predictorOnlyFormula <- paste0("~ ", as.character(preparedModel$formula) |>
    #                                  strsplit(split = "~") |>
    #                                  unlist() |>
    #                                  utils::tail(n = 1)) |> stats::as.formula()
    # 
    # counterfactualModelFrame <- stats::model.matrix(predictorOnlyFormula, targetDataCounterfactualFrame)[,-1]
    # 
    # counterfactualModelFrame <- counterfactualModelFrame[, names(outcomeModel$coefficients)]
    
    targetDataCounterfactualFrame[[response]] <- stats::predict(outcomeModel,
                                                                newdata = targetDataCounterfactualFrame,
                                                                type = "class")
  } else if (inherits(outcomeModel, "coxph")) {
    nCounterfactual <- nrow(targetDataCounterfactualFrame)
    targetDataCounterfactualFrame[[response]] <- double(nCounterfactual)
    
    # Get max survival time for survmean
    maxTime <- survival:::survmean(survival::survfit(preparedModel$outcomeModel, newdata = targetDataCounterfactualFrame), rmean = "common")$end.time
    
    # For each observation in counterfactual frame, calculate the fitted survival curve
    for (i in 1:nCounterfactual) {
      counterfactualSurvCurve <- survival::survfit(preparedModel$outcomeModel, newdata = targetDataCounterfactualFrame[i, , drop = F])
      targetDataCounterfactualFrame[[response]][i] <- survival:::survmean(counterfactualSurvCurve, rmean = maxTime)$matrix["rmean"]
    }
  }
  
  effectType <- match.arg(effectType, c("meanDiff", "rr", "or", "hr"))
  
  # Calculate ATE based on desired effect type; except for polr, in which distributional causal effects are calculated
  if (inherits(outcomeModel, "polr")) {
    targetModel <- MASS::polr(stats::as.formula(paste0(response, " ~ ", treatment)), data = targetDataCounterfactualFrame, method = preparedModel$outcomeModel$method)
    if (effectType == "meanDiff") effect <- targetModel$coefficients[1]
    else if (effectType == "or") effect <- exp(targetModel$coefficients[1])
  } else if (effectType != "hr") {
    treatmentMean <- mean(targetDataCounterfactualFrame[[response]][targetDataCounterfactualFrame[[treatment]] == treatmentLevels[2]])
    controlMean <- mean(targetDataCounterfactualFrame[[response]][targetDataCounterfactualFrame[[treatment]] == treatmentLevels[1]])
    if (effectType == "meanDiff") effect <- treatmentMean - controlMean
    else if (effectType == "rr") effect <- treatmentMean / controlMean
    else if (effectType == "or") effect <- (treatmentMean / (1 - treatmentMean)) / (controlMean / (1 - controlMean))
  } else {
    effectFormula <- stats::as.formula(paste0("survival::Surv(", response, ", rep(1, nrow(targetDataCounterfactualFrame))) ~ ", treatment))
    effect <- exp(survival::coxph(effectFormula, data = targetDataCounterfactualFrame)$coefficients[1])
  }
  
  transportGCResult <- list(effect = effect,
                            effectType = effectType,
                            preparedModel = preparedModel,
                            data = targetData)
  
  class(transportGCResult) <- "transportGC"
  
  return(transportGCResult)
  
}

#' @title Summarize results of a fitted MSM object using g-computation
#' 
#' @description
#' Returns summary object which contains summary objects for the MSM and the outcome model, as well as information about response and treatment variables. In the MSM summary object, the correct variance estimators are calculated.
#' 
#' @rdname summary.transportGC
#'
#' @param object Result from \code{transportGC} function
#' @param ... Further arguments from previous function or to pass to next function
#'
#' @return
#' The \code{summary.transportGC} function returns a \code{summary.transportGC} object containing the following components:
#' * \code{msmSummary}: Summary object of MSM with correct variance estimates. 
#' * \code{preparedModelSummary}: Summary object of outcome model, provided only for information. No conclusions should be drawn from the outcome model.
#' * \code{response}: String indicating response variable name.
#' * \code{treatment}: String indicating treatment variable name.
#' * \code{treatmentLevels}: Vector of strings indicating levels of treatment variable
#' 
#' @export
summary.transportGC <- function (object, ...) {
  transportGCResult <- object
  
  preparedModel <- transportGCResult$preparedModel
  
  preparedModelSummary <- summary(preparedModel$outcomeModel)
  response <- preparedModel$response
  if (inherits(preparedModel, "polr")) responseLevels <- preparedModel$responseLevels
  treatment <- preparedModel$treatment
  treatmentLevels <- preparedModel$treatmentLevels
  
  # If model is glm, calculate and replace correct SEs
  
  # msm <- transportGCResult$msm
  # 
  # msmSummary <- summary(msm)
  # 
  # if (inherits(msmSummary, "summary.glm")) {
  #   if (!is.null(msm$var)) msmSummary$cov.scaled <- msm$var
  #   msmSummary$cov.unscaled <- msmSummary$cov.scaled / msmSummary$dispersion
  #   msmSummary$coefficients[, 2] <- sqrt(diag(msmSummary$cov.scaled))
  #   msmSummary$coefficients[, 3] <- msmSummary$coefficients[, 1] / msmSummary$coefficients[, 2]
  #   if (msmSummary$family$family == "gaussian") msmSummary$coefficients[, 4] <- 2 * stats::pt(abs(msmSummary$coefficients[, 3]), msmSummary$df[2], lower.tail = F)
  #   else msmSummary$coefficients[, 4] <- 2 * stats::pnorm(abs(msmSummary$coefficients[, 3]), lower.tail = F)
  # }
  # 
  # # Same for polr
  # 
  # if (inherits(msmSummary, "summary.polr")) {
  #   if (!is.null(msm$var)) msmSummary$coefficients[, 2] <- sqrt(diag(msm$var))
  #   msmSummary$coefficients[, 3] <- msmSummary$coefficients[, 1] / msmSummary$coefficients[, 2]
  # }
  
  summaryTransportGC <- list(effect = transportGCResult$effect,
                    effectType = transportGCResult$effectType,
                    var = transportGCResult$var,
                    preparedModelSummary = preparedModelSummary,
                    response = response,
                    treatment = treatment,
                    treatmentLevels = treatmentLevels)
  
  if (inherits(preparedModel$outcomeModel, "polr")) summaryTransportGC$responseLevels <- responseLevels
  
  class(summaryTransportGC) <- "summary.transportGC"
  
  return(summaryTransportGC)
}

#' @rdname summary.transportGC
#'
#' @param x \code{summary.transportGC} object.
#' @param out Output stream.
#' @param ... Further arguments from previous function or to pass to next function
#'
#' @export
#'
#'
print.summary.transportGC <- function (x, out = stdout(), ...) {
  summaryTransportGC <- x
  
  write(paste0("Response: ", summaryTransportGC$response), out)
  write(paste0("Treatment: ", summaryTransportGC$treatment), out)
  write(paste0("Effect type: ", summaryTransportGC$effectType), out)
  write(paste0("ATE estimate: ", summaryTransportGC$effect), out)
  write(paste0("Standard error: ", sqrt(summaryTransportGC$var)), out)
  
  write("Fitted outcome model:", out)
  print(summaryTransportGC$preparedModelSummary, out)
}

#' @title Visually represent results of transportability analysis using g-computation
#' 
#' @description
#' This function is a wrapper for \code{modelsummary::modelplot} to plot the coefficient estimates in a transportability analysis using g-computation. Note that the correct variance estimates are used in this function.
#' 
#' @param x Result from \code{transportGC} function.
#' @param ... Further arguments from previous function or to pass to next function
#'
#' @return
#' A \code{ggplot} object showing the estimates and confidence intervals of the MSM coefficients.
#' 
#' @export
plot.transportGC <- function (x, ...) {
  transportGCResult <- x
  ti <- data.frame(term = transportGCResult$preparedModel$treatment,
                   estimate = transportGCResult$effect,
                   conf.low = transportGCResult$effect - 1.96 * sqrt(transportGCResult$var),
                   conf.high = transportGCResult$effect + 1.96 * sqrt(transportGCResult$var))
  
  gl <- data.frame(n = nrow(transportGCResult$targetData))
  
  ms <- list(tidy = ti, glance = gl)
  class(ms) <- "modelsummary_list"
  
  resultPlot <- modelsummary::modelplot(ms)
  
  return(resultPlot)
}

#' @title Check validity of g-computation result object
#' 
#' @description
#' A simple helper function that validates whether the components of the given \code{transportGC} object are of the correct types.
#' 
#' @param transportGCResult Result object from \code{transportGC} function
#'
#' @return A boolean indicating whether all components of \code{transportGC} object have the correct types.
#' 
#' @export
#'
is.transportGC <- function (transportGCResult) {
  return(is.numeric(transportGCResult$effect) & is.numeric(transportGCResult$var) &
           is.transportGCPreparedModel(transportGCResult$preparedModel) &
           is.data.frame(transportGCResult$data))
}