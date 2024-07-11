#' @title Transportability analysis using IOPW
#' 
#' @description
#' Estimates the coefficients of a marginal structural model (MSM) using IP weighting in a generalizability or transportability analysis. In particular, the estimators should be unbiased for the coefficients in the superpopulation or the target population, respectively. Currently, this function only supports the case when there are exactly two treatment levels.
#' 
#'
#' @param msmFormula A formula for the MSM to be fitted, which usually includes the outcome, the treatment and any effect modifiers.
#' @param propensityScoreModel Either a formula or a \code{glm} object representing the model for treatment assignment given covariates.
#' @param participationModel Either a formula or a \code{glm} object representing the model for study participation given effect modifiers. If a formula and the "data" argument is a list, then add participation as the left-hand side of the formula.
#' @param propensityWeights Vector of custom weights balancing covariates between treatments. Providing them will override the formula or model provided by \code{propensityScoreModel}. This vector should have as any entries as the sample size of the study data.
#' @param participationWeights Vector of custom weights balancing effect modifiers between study and target populations. Providing them will override the formula or model provided by \code{participationModel}. This vector should have as any entries as the sample size of the study data.
#' @param treatment String indicating name of treatment variable. If \code{NULL}, it will be auto-detected from \code{propensityScoreModel} if provided; otherwise it will remain \code{NULL}. Note that when using custom weights, \code{treatment} should be provided so that \code{summary.transportIP} and \code{plot.transportIP} works.
#' @param participation String indicating name of participation variable. If \code{NULL}, it will be auto-detected from \code{participationModel} if provided; otherwise it will remain \code{NULL}. Note that when using custom weights, \code{participation} should be provided so that \code{summary.transportIP} and \code{plot.transportIP} works.
#' @param response String indicating name of response variable. If \code{NULL}, it will be auto-detected form \code{msmFormula}.
#' @param family Either a \code{family} function as used for \code{glm}, or one of \code{c("polr", "coxph", "survreg")}.
#' @param method Link function used for \code{polr}, one of \code{c("logistic", "probit", "loglog", "cloglog", "cauchit")}.
#' @param data Either a single data frame containing merged study and target datasets, or a list containing the study dataset and the target dataset. Note that if participationModel is a glm object, the datasets would have been merged, so provide the merged dataset containing response, treatment, covariates controlled for in the original study, study participation and effect modifiers if this is the case. Make sure to code treatment and participation as 0-1 or TRUE-FALSE, with 1 and TRUE representing treatment group and study data, respectively.
#' @param transport A boolean indicating whether a generalizability analysis (false) or a transportability analysis (true) is done.
#' @param bootstrapNum Number of bootstrap datasets to simulate to obtain robust variance estimate.
#'
#' @details
#' The function fits models of treatment assignment and study participation in order to calculate the weights used to fit the MSM. For each of these models, if a formula is provided, logistic regression is used by default. If a \code{glm} object is provided, the function extracts the necessary weights from the object. The function does not support other weighting methods, so if they are required, provide custom weights.
#' 
#' The MSM-fitting functions do not provide correct standard errors as-is. Bootstrap is used to calculate robust bootstrap variance estimators of the parameter estimators. The function replaces the variance component in \code{summary.glm}, \code{coxph} and \code{survreg} with the robust variance estimators directly. This does not seem to behave well with \code{predict.glm} yet, but prediction is not of primary interest in a generalizability or transportability analysis.
#'
#' @return
#' A \code{transportIP} object containing the following components:
#' * \code{msm}: Raw model fit object for MSM of class \code{glm}, \code{survreg} and \code{coxph}, with the correct variance estimators appropriately replaced. If of class \code{glm}, it will have an extra \code{var} component containing the correct variance estimates.
#' * \code{propensityScoreModel}: Model of treatment assignment, \code{NULL} if not provided and custom propensity weights are used
#' * \code{participationModel}: Model of study participation, \code{NULL} if not provided and custom propensity weights are used
#' * \code{propensityWeights}: Propensity weights used
#' * \code{participationWeights}: Participation weights used
#' * \code{finalWeights}: Weights used to fit MSM
#' * \code{customPropensity}: Boolean indicating whether custom propensity weights are used
#' * \code{customParticipation}: Boolean indicating whether custom participation weights are used
#' * \code{treatment}: String indicating variable name of treatment
#' * \code{participation}: String indicating variable name of participation
#' * \code{response}: String indicating variable name of response
#' * \code{data}: Data provided in \code{data} argument. Either a list containing study data and target data or a data frame containing both.
#' 
#' @export
#'
#' @md
transportIP <- function (msmFormula,
                         propensityScoreModel = NULL,
                         participationModel = NULL,
                         propensityWeights = NULL,
                         participationWeights = NULL,
                         treatment = NULL,
                         participation = NULL,
                         response = NULL,
                         family = stats::gaussian,
                         method = c("logistic", "probit", "loglog", "cloglog", "cauchit"),
                         data, transport = T, bootstrapNum = 500) {
  transportIPResult <- transportIPFit(msmFormula,
                                      propensityScoreModel,
                                      participationModel,
                                      propensityWeights,
                                      participationWeights,
                                      treatment,
                                      participation,
                                      response,
                                      family, method, data, transport)
  
  # Correct variance estimates by performing bootstrap, resampling study and target data separately
  
  if (!transportIPResult$customPropensity & !transportIPResult$customParticipation) {
    # Extract treatment study data, control study data and target data
    propensityScoreModel <- transportIPResult$propensityScoreModel
    studyData <- propensityScoreModel$data
    treatmentLevels <- levels(as.factor(propensityScoreModel$y))
    treatmentGroupData <- list()
    for (level in treatmentLevels) {
      treatmentGroupData[[level]] <- studyData[as.character(propensityScoreModel$y) == level, ]
    }
    nLevels <- sapply(treatmentGroupData, nrow)
    names(nLevels) <- treatmentLevels
    
    participationModel <- transportIPResult$participationModel
    targetData <- participationModel$data[participationModel$y == 0 | participationModel$y == F, ]
    if (transportIPResult$response %in% names(targetData)) targetData[[transportIPResult$response]] <- NULL
    nTarget <- nrow(targetData)
    
    bootstrapEstimates <- t(sapply(1:bootstrapNum,
                            function (x) {
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
                              
                              resultBoot <- transportIPFit(msmFormula,
                                                        propensityScoreModel,
                                                        participationModel,
                                                        propensityWeights,
                                                        participationWeights,
                                                        treatment,
                                                        participation,
                                                        response,
                                                        family, data = list(studyBoot, targetBoot), transport)
                              
                              
                              # Add on intercept estimates from polr case. This is okay because concatenating with a NULL does nothing
                              return(c(resultBoot$msm$coefficients, resultBoot$msm$zeta))
                            }))
    
    # Still okay outside of polr cases because concatenating with a NULL does nothing.
    varMatrix <- stats::var(bootstrapEstimates)
    colnames(varMatrix) <- rownames(varMatrix) <- c(names(transportIPResult$msm$coefficients), names(transportIPResult$msm$zeta))
    transportIPResult$msm$var <- varMatrix
  } else {
    warning("Custom weights are being used. Variance estimates may be biased.")
  }
  
  transportIPResult$bootstrapNum = bootstrapNum
  
  return(transportIPResult)
}


transportIPFit <- function(msmFormula,
                        propensityScoreModel = NULL,
                        participationModel = NULL,
                        propensityWeights = NULL,
                        participationWeights = NULL,
                        treatment = NULL,
                        participation = NULL,
                        response = NULL,
                        family = stats::gaussian,
                        method = c("logistic", "probit", "loglog", "cloglog", "cauchit"),
                        data, transport = T) {
  
  # Auto-detect response, treatment and participation if provided
  if (is.null(response)) response <- all.vars(msmFormula)[1]

  if (is.null(treatment) & !is.null(propensityScoreModel)) {
    if (is.glm(propensityScoreModel)) treatment <- all.vars(propensityScoreModel$formula)[1]
    else treatment <- all.vars(propensityScoreModel)[1]
  }

  if (is.null(participation) & !is.null(participationModel)) {
    if (is.glm(participationModel)) participation <- all.vars(participationModel$formula)[1]
    else participation <- all.vars(participationModel)[1]
  }
  
  # Extract study and target data if provided separately using the fact that only study data would have responses
  studyData <- NULL
  targetData <- NULL
  treatmentIndex <- NULL
  participationIndex <- NULL
  
  if (!is.data.frame(data)) {
    if (response %in% names(data[[1]])) {
      studyData <- data[[1]]
      targetData <- data[[2]]
    }
    else{
      studyData <- data[[2]]
      targetData <- data[[1]]
    }
  } else {
    if (!is.null(treatment)) treatmentIndex <- which(names(data) == treatment)
    if (!is.null(participation)) participationIndex <- which(names(data) == participation)
  }
  

  # If formula is provided for treatment and participation models, fit models ourselves
  if (inherits(propensityScoreModel, "formula")) {
    # Fit model ourselves then reassign it to propensityScoreModel
    if (!is.data.frame(data)) propensityScoreModel <- stats::glm(propensityScoreModel, data = studyData, family = stats::binomial())
    else propensityScoreModel <- stats::glm(propensityScoreModel, data = data[data[,participationIndex] == 1 | data[,participationIndex] == T,], family = stats::binomial())
  }
  
  if (inherits(participationModel, "formula")) {
    # Fit model ourselves then reassign it to participationModel
    if (!is.data.frame(data)) {
      # Handle case when study and target data are not yet merged
      effectModifiers <- all.vars(participationModel)[-1]
      studyParticipationData <- studyData[, names(studyData) %in% c(effectModifiers, "participation"), drop = F]
      targetParticipationData <- targetData[, names(targetData) %in% c(effectModifiers, "participation"), drop = F]
      if (!(participation %in% names(studyParticipationData))) {
        studyParticipationData$participation <- 1
        participation <- "participation"
      }
      if (!(participation %in% names(targetParticipationData))) {
        targetParticipationData$participation <- 0
        participation <- "participation"
      }
      participationData <- rbind(studyParticipationData, targetParticipationData)
      participationIndex <- which(names(participationData) == participation)
      participationData$participation <- as.factor(participationData$participation)
      participationModel <- stats::glm(participationModel, data = participationData, family = stats::binomial())
    }
    else participationModel <- stats::glm(participationModel, data = data, family = stats::binomial())
  }
  
  # If not using custom weights, both propensity and participation models need to be ready at this point
  stopifnot(is.glm(propensityScoreModel) | is.null(propensityScoreModel), is.glm(participationModel) | is.null(participationModel))
  
  # Track custom weights
  
  customPropensity <- F
  if (!is.null(propensityWeights)) {
    warning("Custom propensity weights are being used. Please ensure that these weights are meaningful.")
    customPropensity <- T
  }
  
  customParticipation <- F
  if (!is.null(participationWeights)){
    warning("Custom participation weights are being used. Please ensure that these weights are meaningful.")
    customParticipation <- T
  }
  
  if (!is.null(propensityScoreModel) & !is.null(propensityWeights)) warning("Both propensity model and custom weights are provided, using custom weights.")
  
  if (!is.null(participationModel) & !is.null(participationWeights)) warning("Both participation model and custom weights are provided, using custom weights.")
  
  # Calculate actual weights to be used
  
  if (is.null(propensityWeights)) propensityWeights <- obtainWeights(propensityScoreModel, type = "probability")
  if (is.null(participationWeights)) participationWeights <- obtainWeights(participationModel, type = ifelse(transport, "odds", "probability"))
  
  # Makeshift solution to account for generalizability analysis
  
  if (!transport & length(participationWeights) > length(propensityWeights))
    participationWeights <- participationWeights[participationModel$data[, participationIndex] == 1 |
                                                   participationModel$data[, participationIndex] == T]
  
  finalWeights <-  propensityWeights * participationWeights
              
  if (!customPropensity) if (!(treatment %in% all.vars(msmFormula)[-1])) stop("Treatment is not included in MSM.")
  
  # Fit MSM
  
  # Extract study data because data frames don't behave well with ifelse
  if (!is.data.frame(data)) toAnalyze <- studyData
  else toAnalyze <- data[data[,participationIndex] == 1 | data[,participationIndex] == T,]
  
  # The model fitting functions require weights to be part of the data frame
  toAnalyze$finalWeights <- finalWeights
  
  if (is.character(family)) {
    if (family == "coxph") {
    model <- survival::coxph(msmFormula, data = toAnalyze, weight = finalWeights)
    } else if (family == "survreg") {
    model <- survival::survreg(msmFormula, data = toAnalyze, weight = finalWeights)
    } else if (family == "polr") {
    model <- MASS::polr(msmFormula, data = toAnalyze, weights = finalWeights, method = method)
    }
  } else {
    model <- stats::glm(msmFormula, family = family, data = toAnalyze, weight = finalWeights)
  }
  
  # NOTE: this model object by itself has wrong SEs. The right ones are calculated by sandwich::vcovBS. Should we handle this here or in summary.transportIP?
  # Answer: for glm objects, it should be handled in summary.transportIP. For survreg and coxph, it should be handled in this function
  # Addendum: sandwich::vcovBS doesn't seem to be giving correct variance estimates. Bootstrap is implemented in the wrapper function.
  
  transportIPResult <- list(msm = model,
                            propensityScoreModel = propensityScoreModel,
                            participationModel = participationModel,
                            propensityWeights = propensityWeights,
                            participationWeights = participationWeights,
                            finalWeights = finalWeights,
                            customPropensity = customPropensity,
                            customParticipation = customParticipation,
                            treatment = treatment,
                            participation = participation,
                            response = response,
                            data = data,
                            transport = transport)
  
  class(transportIPResult) <- "transportIP"
  
  return(transportIPResult)
}

# Helper function that extracts weights from models
obtainWeights <- function(model, type = c("probability", "odds")) {
  type <- match.arg(type, c("probability", "odds"))
  
  if (type == "probability") {
    return(ifelse(model$y == T | model$y == 1, 1 / model$fitted.values, 1 / (1 - model$fitted.values)))
  } else {
    participationProb <- stats::predict(model, newdata = model$data[model$y == 1 | model$y == T,], type = "response")
    return((1-participationProb) / participationProb)
  }
}

# Helper function that detects glms
is.glm <- function(x) inherits(x, "glm")

#' @title Summarize results of a fitted MSM using the IOPW approach
#' 
#' @description
#' Returns summary object which contains a summary of the fitted MSM as well as pre- and post-weighting standardized mean differences (SMDs).
#' 
#' @rdname summary.transportIP
#'
#' @param object Result from \code{transportIP} function
#' @param covariates Vector of strings indicating names of covariates in propensity model
#' @param effectModifiers Vector of strings indicating names of effect modifiers in participation model
#' @param ... Further arguments from previous function or to pass to next function
#'
#' @return
#' The \code{summary.transportIP} function returns a \code{summary.transportIP} object containing the following components:
#' * \code{propensitySMD}: Table of unweighted and weighted SMDs of covariates between treatment groups. Only propensity weights are used.
#' * \code{participationSMD}: Table of unweighted and weighted SMDs of effect modifiers between study data and target data. Only participation weights are used.
#' * \code{msmSummary}: Summary object of model object for MSM. The correct variance estimators are set here for \code{glm}, whereas they are set in \code{transportIP} for \code{survreg} and \code{coxph}.
#'
#' @export
summary.transportIP <- function(object, covariates = NULL, effectModifiers = NULL, ...) {
  transportIPResult <- object
  
  treatment <- transportIPResult$treatment
  participation <- transportIPResult$participation
  
  # Calculate SMDs for covariates
  
  # Extract study data and weights
  if (!is.null(transportIPResult$propensityScoreModel)) {
    propensityScoreModel <- transportIPResult$propensityScoreModel
    propensityFormula <- propensityScoreModel$formula
    studyData <- propensityScoreModel$data
    if (is.null(covariates)) covariates <- all.vars(propensityFormula)[-1]
  } else {
    data <- transportIPResult$data
    if (is.data.frame(data)) {
      participationIndex <- which(names(data) == participation)
      studyData <- data[data[, participationIndex] == 1 |
                        data[, participationIndex] == T, ]
    } else {
      if (transportIPResult$response %in% names(data[[1]])) {
        studyData <- data[[1]]
        targetData <- data[[2]]
      }
      else{
        studyData <- data[[2]]
        targetData <- data[[1]]
      }
    }
  }
  propensityWeights <- transportIPResult$propensityWeights
  treatmentIndex <- which(names(studyData) == treatment)
  
  # Calculate SMDs for covariates (should optimize)
  prePropensitySMD <- sapply(covariates, function (covariate) as.double(smd::smd(x = studyData[, which(names(studyData) == covariate)],
                                                                       g = studyData[, treatmentIndex])$estimate))
  prePropensityBalance <- data.frame(variable = covariates, smd = abs(prePropensitySMD), method = rep("Observed", length(covariates)))
  postPropensitySMD <- sapply(covariates, function (covariate) as.double(smd::smd(x = studyData[, which(names(studyData) == covariate)],
                                                                    g = studyData[, treatmentIndex],
                                                                    w = propensityWeights)$estimate))
  postPropensityBalance <- data.frame(variable = covariates, smd = abs(postPropensitySMD), method = rep("Weighted", length(covariates)))
  
  propensityBalance <- rbind(prePropensityBalance, postPropensityBalance)
  
  # Calculate SMDs for effect modifiers
  
  if (!is.null(transportIPResult$participationModel)) {
    participationModel <- transportIPResult$participationModel
    participationFormula <- participationModel$formula
    participationData <- participationModel$data
    participationIndex <- which(names(participationData) == participation)
    if (is.null(effectModifiers)) effectModifiers <- all.vars(participationFormula)[-1]
  } else {
    if (is.data.frame(data)) {
      participationData <- data
    } else {
      # This only happens when user provides custom participation weights and separate study and target data
      studyDataEM <- studyData[, names(studyData) %in% effectModifiers, drop = F]
      studyDataEM$participation <- 1
      targetDataEM <- targetData[, names(targetData) %in% effectModifiers, drop = F]
      targetDataEM$participation <- 0
      participationData <- rbind(studyDataEM, targetDataEM)
      participation <- "participation"
      participationIndex <- which(names(participationData) == participation)
      participationData$participation <- as.factor(participationData$participation)
    }
  }
  # Only addresses transportability - need to implement generalizability later
  if (transportIPResult$transport) {
    participationWeights <- rep(1, nrow(participationData))
    participationWeights[participationData[, participationIndex] == 1 | participationData[, participationIndex] == T] <- transportIPResult$participationWeights
  } else if (!is.null(participationModel)) {
    participationProb <- stats::predict(participationModel, newdata = participationData, type = "response")
    ifelse(participationData[, participationIndex] == 1 | participationData[, participationIndex] == T, 1 / participationProb, 1 / (1 - participationProb))
  } else {
    warning("Custom participation weights used for generalizability analysis. Defaulting to transportability analysis because I don't know what weights you use for target data.")
    participationWeights <- rep(1, nrow(participationData))
    participationWeights[participationData[, participationIndex] == 1 | participationData[, participationIndex] == T] <- transportIPResult$participationWeights
  }
  
  # Should also optimize
  preParticipationSMD <- sapply(effectModifiers, function (effectModifier) as.double(smd::smd(x = participationData[, which(names(participationData) == effectModifier)],
                                                                                 g = participationData[, participationIndex])$estimate))
  preParticipationBalance <- data.frame(variable = effectModifiers, smd = abs(preParticipationSMD), method = rep("Observed", length(effectModifiers)))
  postParticipationSMD <- sapply(effectModifiers, function (effectModifier) as.double(smd::smd(x = participationData[, which(names(participationData) == effectModifier)],
                                                                                              g = participationData[, participationIndex],
                                                                                              w = participationWeights)$estimate))
  postParticipationBalance <- data.frame(variable = effectModifiers, smd = abs(postParticipationSMD), method = rep("Weighted", length(effectModifiers)))
  participationBalance <- rbind(preParticipationBalance, postParticipationBalance)
  
  # If model is glm, calculate and replace correct SEs
  
  msm <- transportIPResult$msm
  
  msmSummary <- summary(msm)
  
  if (inherits(msmSummary, "summary.glm")) {
    if (!is.null(msm$var)) msmSummary$cov.scaled <- msm$var
    msmSummary$cov.unscaled <- msmSummary$cov.scaled / msmSummary$dispersion
    msmSummary$coefficients[, 2] <- sqrt(diag(msmSummary$cov.scaled))
    msmSummary$coefficients[, 3] <- msmSummary$coefficients[, 1] / msmSummary$coefficients[, 2]
    if (msmSummary$family$family == "gaussian") msmSummary$coefficients[, 4] <- 2 * stats::pt(abs(msmSummary$coefficients[, 3]), msmSummary$df[2], lower.tail = F)
    else msmSummary$coefficients[, 4] <- 2 * stats::pnorm(abs(msmSummary$coefficients[, 3]), lower.tail = F)
  }
  
  # Same for polr
  
  if (inherits(msmSummary, "summary.polr")) {
    if (!is.null(msm$var)) msmSummary$coefficients[, 2] <- sqrt(diag(msm$var))
    msmSummary$coefficients[, 3] <- msmSummary$coefficients[, 1] / msmSummary$coefficients[, 2]
  }
  
  summaryTransportIP <- list(propensitySMD = propensityBalance,
                             participationSMD = participationBalance,
                             msmSummary = msmSummary)
  
  class(summaryTransportIP) <- "summary.transportIP"
  
  return(summaryTransportIP)
}

#' @rdname summary.transportIP
#'
#' @param x \code{summary.transportIP} object.
#' @param out Output stream.
#' @param ... Further arguments from previous function or to pass to next function
#'
#' @export
#'
#' 
print.summary.transportIP <- function(x, out = stdout(), ...) {
  summaryTransportIP <- x
  
  write("Absolute SMDs of covariates between treatments before and after weighting:", out)
  print(summaryTransportIP$propensitySMD, out)
  
  write("Absolute SMDs of effect modifiers between study and target populations before and after weighting:", out)
  print(summaryTransportIP$participationSMD, out)
  
  write("MSM results:", out)
  print(summaryTransportIP$msmSummary, out)
}

#' @title Plot graphs relevant to transportability analysis using IOPW
#'
#' @description
#' Plot graphs for assessment of covariate balance and results in a IOPW analysis. This function currently supports mirrored histograms, SMD plots and model coefficient plots.
#'
#' @param x Result from \code{transportIP} function
#' @param type One of \code{"propensityHist", "propensitySMD", "participationHist", "participationSMD", "msm"}. \code{Hist} produces mirrored histograms of estimated probability of treatment between treatment groups (for \code{propensity}), or of estimated probability of participation between study and target data (for \code{participation}). \code{SMD} produces SMD plots of covariates between treatment groups (for \code{propensity}) or effect modifiers between study and target data (for \code{participation}). \code{msm} produces plots showing confidence intervals for the model coefficients, which should have the correct standard errors.
#' @param bins Number of bins for propensity score/participation probability histograms. This is only used for \code{Hist}.
#' @param covariates Vector of strings indicating names of covariates in propensity model
#' @param effectModifiers Vector of strings indicating names of effect modifiers in participation model
#' @param ... Further arguments from previous function or to pass to next function
#'
#' @return
#' A \code{ggplot} object which contains the desired plot.
#'
#' @export
#'
#' @importFrom rlang .data
#' @importFrom stats density
plot.transportIP <- function(x, type = "propensityHist", bins = 50, covariates = NULL, effectModifiers = NULL, ...) {
  transportIPResult <- x
  summaryTransportIP <- summary(transportIPResult, covariates = covariates, effectModifiers = effectModifiers)
  resultPlot <- NULL
  
  # Match argument
  type <- match.arg(type, c("propensityHist", "propensitySMD", "participationHist", "participationSMD", "msm"))
  
  if (type == "propensityHist") {
    # Mirrored histogram of propensity scores
    if (!transportIPResult$customPropensity) {
      propensityModel <- transportIPResult$propensityScoreModel
      studyData <- propensityModel$data
      studyData$propensityScore <- propensityModel$fitted.values
      treatmentVar <- transportIPResult$treatment
      treatmentGroups <- levels(factor(studyData[[treatmentVar]]))
      treatmentGroupIdx <- list()
      for (i in 1:length(treatmentGroups)) {
        treatmentGroupIdx[[i]] <- which(studyData[[treatmentVar]] == treatmentGroups[i])
      }
      resultPlot <- ggplot2::ggplot(data = studyData, mapping = ggplot2::aes(.data$propensityScore)) +
        halfmoon::geom_mirror_histogram(ggplot2::aes(group = .data[[!!treatmentVar]], fill = .data[[!!treatmentVar]]), bins = bins) +
        ggplot2::stat_density(data = studyData[treatmentGroupIdx[[1]],], ggplot2::aes(x = .data$propensityScore, y = -ggplot2::after_stat(density)), color = "black", alpha = 0) + 
        ggplot2::stat_density(data = studyData[treatmentGroupIdx[[2]],], ggplot2::aes(x = .data$propensityScore, y = ggplot2::after_stat(density)), color = "black", alpha = 0) +
        ggplot2::ylab("Count")
    } else {
      stop("Custom propensity weights were used. Please plot your previously estimated propensity scores using the halfmoon package, if desired.")
    }
  } else if (type == "propensitySMD") {
    # SMD plot of covariates
      propensitySMD <- summaryTransportIP$propensitySMD
      resultPlot <- ggplot2::ggplot(propensitySMD, ggplot2::aes(x = .data$variable, y = .data$smd, group = .data$method, color = .data$method)) + 
        ggplot2::geom_line() +
        ggplot2::geom_point() +
        ggplot2::geom_hline(yintercept = 0.1) +
        ggplot2::coord_flip() +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.key = ggplot2::element_blank())
    } else if (type == "participationHist") {
    # Mirrored histogram of probability of participation
      if (!transportIPResult$customParticipation) {
        participationModel <- transportIPResult$participationModel
        allData <- participationModel$data
        allData$participationScore <- participationModel$fitted.values
        participationVar <- transportIPResult$participation
        participationGroups <- levels(factor(allData[[participationVar]]))
        participationGroupIdx <- list()
        for (i in 1:length(participationGroups)) {
          participationGroupIdx[[i]] <- which(allData[[participationVar]] == participationGroups[i])
        }
        resultPlot <- ggplot2::ggplot(allData, ggplot2::aes(.data$participationScore)) +
          halfmoon::geom_mirror_histogram(ggplot2::aes(group = .data[[!!participationVar]], fill = .data[[!!participationVar]]), bins = bins) +
          ggplot2::stat_density(data = allData[participationGroupIdx[[1]],], ggplot2::aes(x = .data$participationScore, y = -ggplot2::after_stat(density)), color = "black", alpha = 0) + 
          ggplot2::stat_density(data = allData[participationGroupIdx[[2]],], ggplot2::aes(x = .data$participationScore, y = ggplot2::after_stat(density)), color = "black", alpha = 0) +
          ggplot2::ylab("Count")
      } else {
        stop("Custom participation weights were used. Please plot your previously estimated participation scores using the halfmoon package, if desired.")
      }
    } else if (type == "participationSMD") {
    # SMD plots of effect modifiers
      participationSMD <- summaryTransportIP$participationSMD
      resultPlot <- ggplot2::ggplot(participationSMD, ggplot2::aes(x = .data$variable, y = .data$smd, group = .data$method, color = .data$method)) + 
        ggplot2::geom_line() +
        ggplot2::geom_point() +
        ggplot2::geom_hline(yintercept = 0.1) +
        ggplot2::coord_flip() +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.key = ggplot2::element_blank())
    } else if (type == "msm") {
      # Coefficient plots
      resultPlot <- modelsummary::modelplot(transportIPResult$msm,
                                            vcov = list(transportIPResult$msm$var),
                                            coef_omit = "(Intercept)")
    }
  
  return(resultPlot)
}

#' @title Check validity of IOPW result object
#' 
#' @description
#' A simple helper function that validates whether the components of the given \code{transportIP} object are of the correct types.
#' 
#' @param transportIPResult Result object from \code{transportIP} function
#'
#' @return A boolean indicating whether all components of \code{transportIP} object have the correct types.
#' 
#' @export
#'
is.transportIP <- function (transportIPResult) {
  return((inherits(transportIPResult$msm, "glm") | inherits(transportIPResult$msm, "coxph") | inherits(transportIPResult$msm, "survreg")) &
           (inherits(transportIPResult$propensityScoreModel, "glm") | is.null(transportIPResult$propensityScoreModel)) & 
           (inherits(transportIPResult$participationModel, "glm") | is.null(transportIPResult$participationModel)) &
           inherits(transportIPResult$propensityWeights, "numeric") & 
           inherits(transportIPResult$participationWeights, "numeric") &
           inherits(transportIPResult$finalWeights, "numeric") &
           inherits(transportIPResult$customPropensity, "logical") &
           inherits(transportIPResult$customParticipation, "logical") &
           (inherits(transportIPResult$treatment, "character") | is.null(transportIPResult$treatment)) &
           (inherits(transportIPResult$participation, "character") | is.null(transportIPResult$participation)) &
           (inherits(transportIPResult$data, "data.frame") | length(transportIPResult$data) == 2) &
           inherits(transportIPResult, "transportIP"))
}