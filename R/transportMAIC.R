# Transport MAIC

# Version Date: June 27th 2024
# Richard Yan @Core Clinical Sciences

transportMAIC() <- function(msmFormula, 
                            
                            propensityScoreModel = NULL, 
                            match_cov = NULL, # user-specified matching covariates inputs
                            
                            propensityWeights = NULL, # vector of weights
                            participationWeights = NULL, # vector of weights
                            
                            treatment = NULL, # string, name of treatment
                            response = NULL, # string, name of response
                            
                            family = stats::gaussian, # any available family for glm such as "gaussian", OR, "coxph" / "survreg"
                            
                            IPD, # data of study population (IPD): N rows, data.frame with responses and variables
                            AgD  # data of target population (AgD): 1 row, data.frame with only aggregate variables

                            ) {
  
  if (is.null(response)) response <- all.vars(msmFormula)[1] 
  
  if (is.null(treatment) & !is.null(propensityScoreModel)) { 
    if (is.glm(propensityScoreModel)) treatment <- all.vars(propensityScoreModel$formula)[1] 
    else treatment <- all.vars(propensityScoreModel)[1] 
  }

  studyData <- NULL
  targetData <- NULL
  
  stopifnot(is.data.frame(IPD), is.data.frame(AgD)) # input format checking
  stopifnot(!(response %in% names(IPD))) # there must be "response" in IPD dataset
  stopifnot(!(nrow(AgD) == 1)) # AgD can only be 1 row data.frame
  
  # Check input formats
  stopifnot(is.data.frame(IPD), is.data.frame(AgD))
  stopifnot('response' %in% names(IPD))
  stopifnot(nrow(AgD) == 1)
  
  # Determine matching covariates
  if (!is.null(match_cov)) {
    # Check if the user-specified matching covariates exist in both IPD and AgD
    valid_cov <- intersect(match_cov, names(IPD))
    valid_cov <- intersect(valid_cov, names(AgD))
    
    # If there are any covariates in the user input that don't match, give a warning and remove them
    if (length(valid_cov) != length(match_cov)) {
      removed_cov <- setdiff(match_cov, valid_cov)
      cat("The following user-specfied matching covariates were not found in neither IPD or AgD and have been removed:", toString(removed_cov), "\n")
      match_cov <- valid_cov
    }
    
  } else {
    # If no user input for matching covariates, use all possible matching covariates from both datasets
    match_cov <- intersect(names(IPD), names(AgD)) # initial covariates could be used for matching
  }
  
  # If all user-specified matching covariates are invalid, raise an error and stop execution
  if (length(match_cov) == 0) {
    stop("All user-specified matching covariates are not found in neither IPD or AgD, execution stopped!")
  }
  
  studyData <- IPD[ ,names(IPD) %in% c("response", match_cov)] 
  targetData <- AgD[ ,names(AgD) %in% match_cov] 
  
  
  #  If formula is provided for treatment and participation models, fit models ourselves
  if (inherits(propensityScoreModel, "formula")) { # if it is in the format of formula 
    propensityScoreModel <- stats::glm(propensityScoreModel, 
                                       data = studyData, 
                                       family = stats::binomial()) # logistic regression
  }
  
  stopifnot(is.glm(propensityScoreModel) | is.null(propensityScoreModel))
  
  # Track custom weights
  customPropensity <- F
  if (!is.null(propensityWeights)) {
    warning("Custom propensity weights are being used. Please ensure that these weights are meaningful.")
    customPropensity <- T
  }

  if (!is.null(propensityScoreModel) & !is.null(propensityWeights)) warning("Both propensity model and custom weights are provided, using custom weights.")
  
  if (is.null(propensityWeights)) 
    propensityWeights <- obtainPropensityWeights(propensityScoreModel, type = "probability")
    
  
  if (is.null(participationWeights)) # MoM but not the same function, so not the same here (might be helpful to write another helper to do so)
   {    
     study_vars <- names(studyData)
      for (var in study_vars) {
        mean_var <- var  
        sd_var <- paste(var, "sd", sep = "_")  
      
        if (mean_var %in% names(targetData)) {
          studyData[var] <- studyData[var] - targetData[[mean_var]]
        }
      
        if (mean_var %in% names(targetData) && sd_var %in% names(targetData)) {
          squared_var_name <- paste(var, "squared_centered", sep = "_")
          studyData[[squared_var_name]] <- studyData[[var]]^2 - (targetData[[mean_var]]^2 + targetData[[sd_var]]^2)
        }
      }
    
      intervention_data <- studyData
      match_cov_centered <- c(names(intervention_data), squared_centered_vars)
    
      participationWeightCalculation <- MAIC::estimate_weights(intervention_data = intervention_data,
                                                               matching_vars = match_cov_centered,
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
  # Variance of coeffs are corrected in summary.transportMAIC for glm()
  
  toAnalyze <- studyData
  
  # The model fitting functions require weights to be part of the data frame
  toAnalyze$finalWeights <- finalWeights # add new column
  
  # *** NOTE: a reminder that we just need to use the finalWeights to match the studyData(IPD) only. It is the transportability all about ***
  
  if (is.character(family)) {
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
  } 
  else {  ## *** debug *** when [family] has an input with character but not any correct type of distribution we need a warning to notify the user (## for example when user actually input cox instead of coxph it would be led to glm but not the thing user actually wants [family = family might also lead problem])
      model <- stats::glm(msmFormula, 
                          family = family, 
                          data = toAnalyze, 
                          weight = finalWeights)
      # NOTE: this model object by itself has wrong SEs. The right ones are calculated by sandwich::vcovBS. Should we handle this here or in summary.transportIP?
      # Answer: for glm objects, it should be handled in summary.transportMAIC. For survreg and coxph, it should be handled in this function
      # it just follows the glm() and glm.summary() structure where they get the outputs
    }

  
  transportMAICResult <- list(msm = model,
                              propensityScoreModel = propensityScoreModel,
                              #participationModel = participationModel,
                              
                              propensityWeights = propensityWeights,
                              participationWeights = participationWeights,
                              
                              finalWeights = finalWeights,
                              
                              customPropensity = customPropensity,
                              #customParticipation = customParticipation,
                              
                              treatment = treatment,
                              #participation = participation,
                              response = response,
                              
                              IPD = IPD,
                              AgD = AgD
                              
                              )
  
  class(transportMAICResult) <- "transportMAIC"
  
  return(transportMAICResult)
  
  # ------------- The End of the Function ------------- #
}



# Helper function that extracts weights from models
# *** In MAIC ***: this function is only used for the propensity weights (inverse-prob for IPD)
obtainPropensityWeights <- function(model, type = c("probability") {
  type <- match.arg(type, c("probability")
  
  if (type == "probability") {
    return(ifelse(model$y == T | model$y == 1, # here, model$y == 1 equals to A = 1
                  1 / model$fitted.values,# 1 / Pr(A = 1 | L)^hat
                  1 / (1 - model$fitted.values))) # 1 / Pr(A = 0 | L)^hat
  } 
}

# Helper function that detects glms
is.glm <- function(x) inherits(x, "glm")


# ***Test dataset generator: MAIC***

# Generate test data for TransportMAIC
generateTestData_MAIC <- function() {
  expit <- function(x) 1/(1+exp(-x))

  # Generate study data: IPD
  nStudy <- 1000
  
  sexStudy <- rbinom(nStudy, 1, 0.5) # Male is 1, so female is baseline
  stressStudy <- rbinom(nStudy, 1, 0.4) # Stressed is 1
  med2Study <- rbinom(nStudy, 1, 0.1) # 1 means taking other med
  percentBodyFatStudy <- rnorm(nStudy, 28 - 13 * sexStudy, 2)
  med1Study <- rbinom(nStudy, 1, expit(0.2 * sexStudy - 0.02 * percentBodyFatStudy + 0.1 * stressStudy))
  sysBloodPressureStudy <- rnorm(nStudy, 100 + 5 * sexStudy + 0.5 * percentBodyFatStudy + 5 * stressStudy -
                                   5 * med1Study + med1Study * (-5 * med2Study + 7 * stressStudy))
  
  # Put all variables together
  studyData <- data.frame(sysBloodPressure = sysBloodPressureStudy, 
                          med1 = as.factor(med1Study), 
                          sex = as.factor(sexStudy), 
                          stress = as.factor(stressStudy), 
                          med2 = as.factor(med2Study), 
                          percentBodyFat = percentBodyFatStudy)
  
  
  
  # Generate target data: AgD
  nTarget <- 3000
  
  sexTarget <- rbinom(nTarget, 1, 0.3) # Male is 1, so female is baseline
  stressTarget <- rbinom(nTarget, 1, 0.7) # Stressed is 1
  med2Target <- rbinom(nTarget, 1, 0.3) # 1 means taking other med
  percentBodyFatTarget <- rnorm(nTarget, 26 - 12 * sexTarget, 2)
  
  # Put all variables together
  targetData <- data.frame(sex = mean(sexTarget), 
                           stress = mean(stressTarget),
                           med2 = mean(med2Target), 
                           percentBodyFat = mean(percentBodyFatTarget), # mean
                           percentBodyFat_sd = sd(percentBodyFatTarget)) #sd: applying by matching X^2 in IPD with (sd^2-mean^2) in AgD
  
  return(list(IPD = studyData, 
              AgD = targetData))
}

