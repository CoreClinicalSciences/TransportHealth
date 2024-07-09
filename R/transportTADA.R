#' TransportTADA
#' 
#' @description
#' Estimates the coefficients of a marginal structural model (MSM) using Target Aggregate Data Adjustment (TADA) Transportability Analysis in a transportability analysis.
#'
#' @param msmFormula A formula for the MSM to be fitted, which usually includes the outcome, the treatment and any effect modifiers.
#' @param propensityScoreModel Either a formula or a \code{glm} object representing the model for treatment assignment given covariates.
#' @param matchingCovariates A vector of user-specified covariates which are used for the data matching to obtain participation weights when no custom weights provided.
#' @param propensityWeights Vector of custom weights balancing covariates between treatments. Providing them will override the formula or model provided by \code{propensityScoreModel}. This vector should have as any entries as the sample size of the study data.
#' @param participationWeights Vector of custom weights balancing effect modifiers between study and target populations. Providing them will override the formula or model provided by \code{participationModel}. This vector should have as any entries as the sample size of the study data.
#' @param treatment String indicating name of treatment variable. If \code{NULL}, it will be auto-detected from \code{propensityScoreModel} if provided; otherwise it will remain \code{NULL}. Note that when using custom weights, \code{treatment} should be provided so that \code{summary.transportTADA} and \code{plot.transportTADA} works.
#' @param response String indicating name of response variable. If \code{NULL}, it will be auto-detected form \code{msmFormula}.
#' 
#' @param family Either a \code{family} function as used for \code{glm}, or one of \code{c("coxph", "survreg")}.
#' 
#' @param studyData The individual patient-level data of study population
#' @param aggregateTargetData The aggregate-level data of target population
#'
#' @return
#' A \code{transportTADA} object containing the following components:
#' * \code{msm}: Raw model fit object for MSM of class \code{glm}, \code{survreg} and \code{coxph}, with the correct variance estimators appropriately replaced. If of class \code{glm}, it will have an extra \code{var} component containing the correct variance estimates.
#' * \code{propensityScoreModel}: Model of treatment assignment, \code{NULL} if not provided and custom propensity weights are used.
#' * \code{propensityWeights}: Propensity weights used. When not \code{NULL}, it will be used as custom inputs from user. When \code{NULL}, it will be obtained by logistic regression.
#' * \code{participationWeights}: Participation Weights used. When not \code{NULL}, it will be used as custom inputs from user. When \code{NULL}, it will be obtained by Method of Moments -based method.
#' * \code{finalWeights}: Weights used to fit MSM
#' * \code{customPropensity}: Boolean indicating whether custom propensity weights are used
#' * \code{customParticipation}: Boolean indicating whether custom participation weights are used
#' * \code{treatment}: String indicating variable name of treatment
#' * \code{response}: String indicating variable name of response
#' * \code{studyData}: The original studyData
#' * \code{aggregateTargetData}: The original aggregateTargetData
#' * \code{centeredStudyData}: The centered data for obtaining participation weights with MoM by MAIC::estimate_weights function
#' 
#' @export
#'
#' @examples
#' 
#' @md

transportTADA() <- function(msmFormula, 
                            
                            propensityScoreModel = NULL, 
                            matchingCovariates = NULL, # User-specified matching covariates inputs
                            
                            propensityWeights = NULL, # vector of (custom) propensity weights
                            participationWeights = NULL, # vector of (custom) participation weights
                            
                            treatment = NULL, # string, name of treatment
                            response = NULL, # string, name of response
                            
                            family = stats::gaussian, # stats::gaussian, # any available family for glm such as "gaussian", OR, "coxph" / "survreg"
                            
                            studyData, # data of study population (studyData): N rows, data.frame with responses and variables
                            aggregateTargetData  # data of target population (aggregateTargetData): 1 row, data.frame with only aggregate variables

                            ) {
  
  if (is.null(response)) response <- all.vars(msmFormula)[1] 
  
  if (is.null(treatment) & !is.null(propensityScoreModel)) { 
    if (is.glm(propensityScoreModel)) treatment <- all.vars(propensityScoreModel$formula)[1] 
    else treatment <- all.vars(propensityScoreModel)[1] 
  }
  
  # Check input formats
  stopifnot(is.data.frame(studyData), is.data.frame(aggregateTargetData)) # input format checking
  stopifnot(response %in% names(studyData)) # there must be "response" in studyData dataset
  stopifnot(nrow(aggregateTargetData) == 1) # aggregateTargetData can only be 1 row data.frame
  
  # Determine matching covariates
  if (!is.null(matchingCovariates)) {
    # Check if the user-specified matching covariates exist in both studyData and aggregateTargetData
    validCov <- intersect(matchingCovariates, names(studyData))
    validCov <- intersect(validCov, names(aggregateTargetData))
    
    # If there are any covariates in the user input that don't match, give a warning and remove them
    if (length(validCov) != length(matchingCovariates)) {
      removedCov <- setdiff(matchingCovariates, validCov)
      cat("The following user-specfied matching covariates were not found in neither studyData or aggregateTargetData and have been removed:", 
          toString(removedCov), "\n")
      matchingCovariates <- validCov
    }
    
  } else {
    # If no user input for matching covariates, use all possible matching covariates from both datasets
    matchingCovariates <- intersect(names(studyData), names(aggregateTargetData)) # initial covariates could be used for matching
    
    # **** Quang - this behavior is error-prone. Let's discuss
    # Richard - the point here is some possible variables existing in both datasets, but not meaningful at all for the analysis (like Index) 
    
  }
  
  # If all user-specified matching covariates are invalid, raise an error and stop execution
  if (length(matchingCovariates) == 0) { # Here should be OR validCov???
    stop("All user-specified matching covariates are not found in neither studyData or aggregateTargetData, execution stopped!")
  }
  
  #  If formula is provided for treatment and participation models, fit models ourselves
  if (inherits(propensityScoreModel, "formula")) { # if it is in the format of formula 
    propensityScoreModel <- stats::glm(propensityScoreModel, 
                                       data = studyData, 
                                       family = stats::binomial()) # logistic regression
  }
  
  stopifnot(is.glm(propensityScoreModel) | is.null(propensityScoreModel))
  
  
  
  # -------- Propensity Weights -----------
  # Track Custom Propensity Weights
  customPropensity <- F
  if (!is.null(propensityWeights)) {
    warning("Custom propensity weights are being used. Please ensure that these weights are meaningful.")
    customPropensity <- T
  }
  if (!is.null(propensityScoreModel) & !is.null(propensityWeights)) warning("Both propensity model and custom weights are provided, using custom weights.")
  
  # obtain propensity weights when no custom weights
  if (is.null(propensityWeights)) 
    propensityWeights <- obtainPropensityWeights(propensityScoreModel, type = "probability")
  
  # -------- Participation Weights -----------
  # Track Custom Participation Weights
  customParticipation <- F
  if (!is.null(participationWeights)) {
    warning("Custom participation weights are being used. Please ensure that these weights are meaningful.")
    customParticipation <- T
  }
  
  # obtain participation weights when no custom weights
  
  # Richard - Update with a separate "centeredStudyData"
  if (is.null(participationWeights)) {
    studyVars <- names(studyData)
    centeredStudyData <- studyData  
    
    for (var in studyVars) {
      meanVar <- var
      # **** Quang - this code makes an assumption that the SDs are named Variable_sd. Let's discuss this
      # Richard - Here we at least need to have a fixed format to ensure the sd inputs for user, it could be a must setting, just mention at the function read.me
      sdVar <- paste(var, "sd", sep = "_")
      
      if (meanVar %in% names(aggregateTargetData)) {
        centeredStudyData[[var]] <- studyData[[var]] - aggregateTargetData[[meanVar]] 
      }
      
      # *** continuous variables only: matching both the mean and the sd together: is.numeric(studyData[[var]]) 
      if (meanVar %in% names(aggregateTargetData) && sdVar %in% names(aggregateTargetData)) {
        squaredVarName <- paste(var, "squared_centered", sep = "_")
        centeredStudyData[[squaredVarName]] <- studyData[[var]]^2 - (aggregateTargetData[[meanVar]]^2 + aggreagteTargetData[[sdVar]]^2)
      }
    }
    
    interventionData <- centeredStudyData
    matchCovCentered <- c(names(interventionData), squaredVarName)
    
    participationWeightCalculation <- MAIC::estimate_weights(intervention_data = interventionData,
                                                             matching_vars = matchCovCentered,
                                                             method = "BFGS")
    
    participationWeights <- participationWeightCalculation$weights
  }

  # final weights: propensityWeights accounts for the confounding and participationWeights accounts for the EM distribution difference between populations
  finalWeights <-  propensityWeights * participationWeights
  
  # ------------------------------------------------------------------------------------------------------------- #
  
  if (!customPropensity) 
    if (!(treatment %in% all.vars(msmFormula)[-1])) stop("Treatment is not included in MSM!")
  
  # Fit MSM
  # Variance of coeffs are corrected in this method for coxph() and survreg() with sandwich::vcovBS(model)
  # Variance of coeffs are corrected in summary.transportTADA for glm()
  
  toAnalyze <- studyData # keep it as studyData
  
  # The model fitting functions require weights to be part of the data frame
  toAnalyze$finalWeights <- finalWeights
  
  if (is.character(family) & !(family %in% c("coxph", "survreg"))) {
    stop("Please check the family input and set it as one of the following in character: gaussian, coxph, survreg.")
  }
  
  if (is.character(family)) {
    
    family <- match.arg(family, c("coxph", "survreg"))
    
    if (family == "coxph") {
      model <- survival::coxph(msmFormula, 
                               data = toAnalyze, 
                               weight = finalWeights)
      model$var <- sandwich::vcovBS(model)
    } 
    else if (family == "survreg") {
      model <- survival::survreg(msmFormula, 
                                 data = toAnalyze, 
                                 weight = finalWeights)
      model$var <- sandwich::vcovBS(model)
    }
  } else {
      model <- stats::glm(msmFormula, 
                          family = family, 
                          data = toAnalyze, 
                          weight = finalWeights)
    }
  
  transportTADAResult <- list(msm = model,
                              propensityScoreModel = propensityScoreModel,

                              propensityWeights = propensityWeights,
                              participationWeights = participationWeights,
                              
                              finalWeights = finalWeights,
                              
                              customPropensity = customPropensity,
                              customParticipation = customParticipation,
                              
                              treatment = treatment,
                              response = response,
                              
                              studyData = studyData,
                              aggregateTargetData = aggregateTargetData,
                              
                              centeredStudyData = centeredStudyData
                              
                              )
  
  class(transportTADAResult) <- "transportTADA"
  
  return(transportTADAResult)
  
  # ------------- The End of the Function ------------- #
}

# Helper function that extracts weights from models
# *** In TADA ***: this function is only used for the propensity weights (inverse-prob for studyData)
obtainPropensityWeights <- function(model, type = c("probability")) {
  type <- match.arg(type, c("probability"))
  
  if (type == "probability") {
    return(ifelse(model$y == T | model$y == 1, # here, model$y == 1 equals to A = 1
                  1 / model$fitted.values,# 1 / Pr(A = 1 | L)^hat
                  1 / (1 - model$fitted.values))) # 1 / Pr(A = 0 | L)^hat
  } 
}

# Helper function that detects glms
is.glm <- function(x) {inherits(x, "glm")}


# Summary 

# Plot

# is.transportTADA






































# Please remove this line. I'll go over test data generation with you once we have more progress on this module.
# Put it in the testthat

# ***Test dataset generator: TADA***
# 
# # Generate test data for TransportTADA
# generateTestData_TADA <- function() {
#   expit <- function(x) 1/(1+exp(-x))
# 
#   # Generate study data: studyData
#   nStudy <- 1000
#   
#   sexStudy <- rbinom(nStudy, 1, 0.5) # Male is 1, so female is baseline
#   stressStudy <- rbinom(nStudy, 1, 0.4) # Stressed is 1
#   med2Study <- rbinom(nStudy, 1, 0.1) # 1 means taking other med
#   percentBodyFatStudy <- rnorm(nStudy, 28 - 13 * sexStudy, 2)
#   med1Study <- rbinom(nStudy, 1, expit(0.2 * sexStudy - 0.02 * percentBodyFatStudy + 0.1 * stressStudy))
#   sysBloodPressureStudy <- rnorm(nStudy, 100 + 5 * sexStudy + 0.5 * percentBodyFatStudy + 5 * stressStudy -
#                                    5 * med1Study + med1Study * (-5 * med2Study + 7 * stressStudy))
#   
#   # Put all variables together
#   studyData <- data.frame(sysBloodPressure = sysBloodPressureStudy, 
#                           med1 = as.factor(med1Study), 
#                           sex = as.factor(sexStudy), 
#                           stress = as.factor(stressStudy), 
#                           med2 = as.factor(med2Study), 
#                           percentBodyFat = percentBodyFatStudy)
#   
#   
#   
#   # Generate target data: aggregateTargetData
#   nTarget <- 3000
#   
#   sexTarget <- rbinom(nTarget, 1, 0.3) # Male is 1, so female is baseline
#   stressTarget <- rbinom(nTarget, 1, 0.7) # Stressed is 1
#   med2Target <- rbinom(nTarget, 1, 0.3) # 1 means taking other med
#   percentBodyFatTarget <- rnorm(nTarget, 26 - 12 * sexTarget, 2)
#   
#   # Put all variables together
#   targetData <- data.frame(sex = mean(sexTarget), 
#                            stress = mean(stressTarget),
#                            med2 = mean(med2Target), 
#                            percentBodyFat = mean(percentBodyFatTarget), # mean
#                            percentBodyFat_sd = sd(percentBodyFatTarget)) #sd: applying by matching X^2 in studyData with (sd^2-mean^2) in aggregateTargetData
#   
#   return(list(studyData = studyData, 
#               aggregateTargetData = targetData))
# }
# 
