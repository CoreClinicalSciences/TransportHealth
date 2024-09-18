#' @title Prepare an outcome model object for \code{transportGC}
#'
#' @description
#' An outcome model needs to be fitted using the study data to perform transportability analysis using g-computation. However, \code{glm} and functions in \code{survival} typically contain the data used to fit the model in their respective result objects. This function removes most components containing the data from these result objects to comply with data sharing regulations. The party with sole access to the study data may use only this function and provide the results object (possibly in a .rds file) to others who request it.
#' 
#' @param outcomeModel Either a formula or a \code{glm}, \code{survreg} or \code{coxph} object representing the outcome model.
#' @param response String indicating name of response variable. If \code{NULL}, it will be auto-detected from \code{outcomeModel}.
#' @param responseLevels For ordinal responses, vector of strings indicating levels of response variable in the study data. If \code{NULL} and \code{polr} is used, it will be auto-detected using \code{response} and \code{studyData}.
#' @param treatment String indicating name of treatment variable. This argument is required.
#' @param treatmentLevels Vector of strings indicating levels of treatment variable in the study data. If \code{NULL}, it will be auto-detected using \code{treatment} and \code{studyData}.
#' @param family Either a family function as used for \code{glm}, or one of \code{c("coxph", "survreg")}. Only required if \code{outcomeModel} is a formula.
#' @param method Link function used for \code{polr}, one of \code{c("logistic", "probit", "loglog", "cloglog", "cauchit")}. Only required if \code{outcomeModel} is a formula and \code{polr} is used.
#' @param studyData Data frame of the study data.
#' @param wipe Logical indicating whether original study data should be wiped from outcome model-fitting object.
#'
#' @return
#' 
#' A \code{transportGCPreparedModel} object containing the following components:
#' * \code{outcomeModel}: The fitted outcome model with all components containing the data used to fit the model removed
#' * \code{response}: String indicating name of response variable
#' * \code{treatment}: String indicating name of treatment variable
#' * \code{treatmentLevels}: Vector of strings indicating levels of treatment variable
#' * \code{family}: The \code{family} argument provided
#' * \code{wipe}: The \code{wipe} argument provided
#' 
#' @export
#' 
#' @md
transportGCPreparedModel <- function(outcomeModel,
                             response = NULL,
                             responseLevels = NULL,
                             treatment,
                             treatmentLevels = NULL,
                             family = stats::gaussian,
                             method = c("logistic", "probit", "loglog", "cloglog", "cauchit"),
                             studyData = NULL,
                             wipe = T,
                             formula = NULL) {
  
  # Extract response variable
  
  if (inherits(outcomeModel, "formula")) {
    formula <- outcomeModel
  }
  
  if (inherits(outcomeModel, "glm")) {
    if (is.null(formula)) formula <- outcomeModel$formula
    if (is.null(studyData)) studyData <- outcomeModel$data
  }
  
  if (is.null(response)) {
    if (inherits(outcomeModel, "formula")) response <- all.vars(outcomeModel)[1]
    else response <- all.vars(formula)[1]
  }
  
  # Extract treatment variable information. There is no way of detecting treatment from the outcome model formula because it often has covariates in it as well.
  if (is.null(treatment)) stop("Treatment variable name is not specified.")

  if (!is.null(studyData)) treatmentLevels <- levels(studyData[[treatment]])
  else treatmentLevels <- levels(outcomeModel$data[[treatment]])
  
  # If not already a model object, fit the outcome model ourselves. Recall that only the party that has access to the study data should run this function and provide its output to the other party.
  if (inherits(outcomeModel, "formula")) {
    if (is.character(family)) {
      if (family == "coxph") {
        outcomeModel <- survival::coxph(outcomeModel, data = studyData, model = !wipe)
      } else if (family == "survreg") {
        outcomeModel <- survival::survreg(outcomeModel, data = studyData, model = !wipe)
      } else if (family == "polr") {
        outcomeModel <- MASS::polr(outcomeModel, data = studyData, method = method, model = !wipe)
      }
    } else {
        outcomeModel <- stats::glm(outcomeModel, family = family, data = studyData, model = !wipe)
    }
  }
  
  # Extract response variable information when polr is used
  if (inherits(outcomeModel, "polr")) {
    if (!is.null(studyData)) responseLevels <- levels(studyData[[response]])
    else responseLevels <- levels(outcomeModel$data[[response]])
  }
  
  # Erase study data from outcome model object
  if (wipe) {
    if (inherits(outcomeModel, "glm")) {
      outcomeModel$residuals <- outcomeModel$fitted.values <- outcomeModel$linear.predictors <-
        outcomeModel$weights <- outcomeModel$y <- outcomeModel$x <- outcomeModel$model <-
        outcomeModel$data <- outcomeModel$offset <- outcomeModel$xlevels <- NULL
    } else if (inherits(outcomeModel, "survreg")) {
      outcomeModel$linear.predictors <- outcomeModel$means <- outcomeModel$y <- outcomeModel$x <-
        outcomeModel$model <- NULL
    # } else if (inherits(outcomeModel, "coxph")) {
    #   outcomeModel$linear.predictors <- outcomeModel$residuals <- outcomeModel$n <- outcomeModel$nevent <-
    #     outcomeModel$concordance <- outcomeModel$means <-
    #     outcomeModel$model <- NULL
    } else if (inherits(outcomeModel, "polr")) {
      outcomeModel$fitted.values <- outcomeModel$n <- outcomeModel$nobs <- outcomeModel$lp <-
        outcomeModel$model <- NULL
    }
  }
  
  preparedModel <- list(outcomeModel = outcomeModel,
                        response = response,
                        responseLevels = responseLevels,
                        treatment = treatment,
                        treatmentLevels = treatmentLevels,
                        family = family,
                        wipe = wipe,
                        formula = formula)
  
  if (!wipe) preparedModel$studyData <- studyData
  
  class(preparedModel) <- "transportGCPreparedModel"
  
  return(preparedModel)
}

#' @title Check validity of prepared model object for g-computation
#' 
#' @description
#' A simple helper function that validates whether the components of the given \code{transportGCPreparedModel} object are of the correct types.
#' 
#' @param preparedModel Output object from \code{transportGCPreparedModel} function
#'
#' @return
#' A boolean indicating whether all components of \code{transportGCPreparedModel} object have the correct types.
#' 
#' @export
is.transportGCPreparedModel <- function (preparedModel) {
  return((inherits(preparedModel$outcomeModel, "glm") | inherits(preparedModel$outcomeModel, "survreg") | inherits(preparedModel$outcomeModel, "coxph") | inherits(preparedModel$outcomeModel, "polr")) &
           is.character(preparedModel$response) & is.character(preparedModel$treatment) & is.character(preparedModel$treatmentLevels) & is.logical(preparedModel$wipe))
}