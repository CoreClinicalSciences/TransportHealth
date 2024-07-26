#' TransportTADA
#' 
#' @description
#' Estimates the coefficients of a marginal structural model (MSM) using Target Aggregate Data Adjustment (TADA) Transportability Analysis in a transportability analysis.
#'
#' @param msmFormula A formula for the MSM to be fitted, which usually includes the outcome, the treatment and any effect modifiers.
#' @param propensityScoreModel Either a formula or a \code{glm} object representing the model for treatment assignment given covariates.
#' @param matchingCovariates A vector of user-specified covariates which are used for the data matching to obtain participation weights when no custom weights provided.
#' @param propensityWeights Vector of custom weights balancing covariates between treatments. Providing them will override the formula or model provided by \code{propensityScoreModel}. This vector should have as any entries as the sample size of the study data.
#' @param participationWeights Vector of custom weights balancing effect modifiers between study and target populations. This vector should have as any entries as the sample size of the study data.
#' @param treatment String indicating name of treatment variable. If \code{NULL}, it will be auto-detected from \code{propensityScoreModel} if provided; otherwise it will remain \code{NULL}. Note that when using custom weights, \code{treatment} should be provided so that \code{summary.transportTADA} and \code{plot.transportTADA} works.
#' @param response String indicating name of response variable. If \code{NULL}, it will be auto-detected form \code{msmFormula}.
#' 
#' @param family Either a \code{family} function as used for \code{glm} such as stats::gaussian(), or one of \code{c("coxph", "survreg")}.
#' 
#' @param studyData The individual patient-level data (IPD) of study population. Ensure that: 1. All binary variables to be used in the matching should be coded 1 and 0;
#' @param aggregateTargetData The aggregate-level data (AgD) of target population. Ensure that: 1. Name columns of mean of continuous variables or proportion of binary variable baselines exactly the same as the column names in the study (IPD) data; 2. Only continuous variables are allowed to consider matching standard deviation (SD) and name the SD column as "variable_sd" in the aggregateTargetData.
#'
#' @return
#' A \code{transportTADA} object containing the following components:
#' * \code{msm}: Raw model fit object for MSM of class \code{glm}, \code{survreg} and \code{coxph}, with the correct variance estimators appropriately replaced. If of class \code{glm}, it will have an extra \code{var} component containing the correct variance estimates.
#' * \code{propensityScoreModel}: Model of treatment assignment, \code{NULL} if not provided and custom propensity weights are used.
#' * \code{propensityWeights}: Propensity weights used. When not \code{NULL}, it will be used as custom inputs from user. When \code{NULL}, it will be obtained by logistic regression.
#' * \code{participationWeights}: Participation Weights used. When not \code{NULL}, it will be used as custom inputs from user. When \code{NULL}, it will be obtained by Method of Moments -based method.
#' * \code{finalWeights}: Weights used to fit MSM.
#' * \code{customPropensity}: Boolean indicating whether custom propensity weights are used.
#' * \code{customParticipation}: Boolean indicating whether custom participation weights are used.
#' * \code{treatment}: String indicating variable name of treatment.
#' * \code{response}: String indicating variable name of response.
#' * \code{studyData}: The original individual patient-level data (IPD) of study population.
#' * \code{aggregateTargetData}: The aggregate-level data (AgD) of target population.
#' * \code{centeredStudyData}: The centered study data for obtaining participation weights with MoM by MAIC::estimate_weights function.
#' 
#' @export
#' 
#' @md

transportTADA <- function(msmFormula, 
                            
                            propensityScoreModel = NULL, 
                            matchingCovariates = NULL, # User-specified matching covariates inputs
                            
                            propensityWeights = NULL, # vector of (custom) propensity weights
                            participationWeights = NULL, # vector of (custom) participation weights
                            
                            treatment = NULL, # string, name of treatment
                            response = NULL, # string, name of response
                            
                            family = stats::gaussian, # any available family for glm such as "gaussian", OR, "coxph" / "survreg"
                            
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
  
  
  # Notice: we need to consider three types of AgD: MEAN & SD for continuous; PROP for binary & ordinal
  # For the user inputs: we need the studyData and aggregateTargetData with all-upper-case column names (SEX, SEX_PROP, AGE_MEAN, AGE_SD)
  # For the binary: the reference level is necessary; for the ordinal: the reference level might be unnecessary
  # All we need to do is to check the validity of user inputs: 
  
  # Determine matching covariates
  # Note that we don't consider the case that user skips the inputs of matching covariates (Reasons at the end of this section)
  if (!is.null(matchingCovariates)) {
    # Check if the user-specified matching covariates exist in both studyData and aggregateTargetData
    validCov <- intersect(matchingCovariates, names(studyData))
    
    validCov <- grepl(paste("^", validCov, sep = "", collapse = "|"), names(aggregateTargetData))
    # for the paste(): if the validCov being c(age, sex), then it gives regular expression ^age|^sex
    # for the grepl(): matches regex pattern within the names of columns in AgD 
    
    # If there are any covariates in the user input that don't match, give a warning and remove them
    if (length(validCov) != length(matchingCovariates)) {
      removedCov <- setdiff(matchingCovariates, validCov)
      cat("The following user-specfied matching covariates were not found in neither studyData or aggregateTargetData and have been removed:",
          toString(removedCov), "\n")
      
      matchingCovariates <- validCov

    # Notice for user to double check the validity of matching variates ready-to-use
    warning(cat("The following covariates are being used for matching:", toString(matchingCovariates),
                ". Please ensure that these covariates are meaningful", "\n"))
    }
  } 
    # # If no user input for matching covariates, use all possible matching covariates from both datasets

    # # **** Quang - this behavior is error-prone. Let's discuss
    # # Richard - the point here is some possible variables existing in both datasets, but not meaningful at all for the analysis (like Index) 
    # # Richard - updated with user notice with the selected ready-to-use covariates, since it's not easy to comprehensively consider input formats and we let it to users. 
  
  

  
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
  
  # obtain participation weights when no custom weights (pending in individual testing scripts)
  if (is.null(participationWeights)) {
    # ----- Participation Weight Estimation ----- #
    
    dummizeIPD <- dummizeIPD(studyData)
    processedAgD <- processedAgD(aggregateTargetData)
    centerIPD <- centerIPD(IPD = dummizeIPD, AgD = processedAgD)
    
    centeredColnames <- grep("_CENTERED", colnames(centerIPD), value = TRUE) # add suffix for centered variables
    
    participationWeightCalculation <- estimateWeights(data = centerIPD,
                                                      centeredColnames = centeredColnames)
    
    participationWeights <- participationWeightCalculation$weights
  }

  
  
  
  # final weights: propensityWeights accounts for the confounding and participationWeights accounts for the EM distribution difference between populations
  finalWeights <-  propensityWeights * participationWeights
  
  # ------------------------------------------------------------------------------------------------------------- #
  
  if (!customPropensity) 
    if (!(treatment %in% all.vars(msmFormula)[-1])) stop("Treatment is not included in MSM!")
  
  # Fit MSM
  # Variance of coeffs are corrected in this method for coxph() and survreg()
  # Variance of coeffs are corrected in summary.transportTADA for glm()
  
  # Use toAnalyze to keep the original studyData safe and clean
  toAnalyze <- studyData
  
  # The model fitting functions require weights to be part of the data frame
  toAnalyze$finalWeights <- finalWeights
  
  if (is.character(family)) {
    
    if (!(family %in% c("coxph", "survreg"))) {
      stop("Please check the family input and set it as one of the following in character: coxph, survreg.")
      }
    
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
                              
                              centeredStudyData = centerIPD
                              
                              )
  
  class(transportTADAResult) <- "transportTADA"
  
  return(transportTADAResult)
  
}
# ------------- The End of the Function ------------- #


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
