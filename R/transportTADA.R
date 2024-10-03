#' @title  Transportability analysis using TADA
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
#' @param family Either a \code{family} function as used for \code{glm} such as stats::gaussian(), or one of \code{c("coxph", "survreg", "polr")}.
#' @param method Link function used for \code{polr}, one of \code{c("logistic", "probit", "loglog", "cloglog", "cauchit")}.
#' @param exOpt A list with components \code{propensity}, \code{participation} and \code{final}. Each component specifies whether weights should be trimmed or truncated. Use the functions \code{trim} and \code{trunc} to specify trimming/truncation. Note that only truncation is supported for final weights.
#' @param studyData The individual participant data (IPD) of study population.
#' @param aggregateTargetData The aggregate-level data (AgD) of target population. Ensure that: 1. Name columns of mean of continuous variables or proportion of binary variable baselines exactly the same as the column names in the study (IPD) data; 2. Only continuous variables are allowed to consider matching standard deviation (SD) and name the SD column as "variable_SD" in the aggregateTargetData.
#'
#' @details
#' The function fits models of treatment assignment and study participation in order to calculate the weights used to fit the MSM. For each of these models, if a formula is provided, logistic regression is used by default. If a \code{glm} object is provided, the function extracts the necessary weights from the object. The function does not support other weighting methods, so if they are required, provide custom weights.
#' 
#' The MSM-fitting functions do not provide correct standard errors as-is. Bootstrap is used to calculate robust bootstrap variance estimators of the parameter estimators. The function replaces the variance component in \code{summary.glm}, \code{coxph} and \code{survreg} with the robust variance estimators directly. This does not seem to behave well with \code{predict.glm} yet, but prediction is not of primary interest in a transportability analysis.
#' 
#' Ensure the binary variables are labelled as 0-1 format
#' 
#' @return
#' A \code{transportTADA} object containing the following components:
#' * \code{msm}: Raw model fit object for MSM of class \code{glm}, \code{survreg} and \code{coxph}, with the correct variance estimators appropriately replaced. If of class \code{glm}, it will have an extra \code{var} component containing the correct variance estimates.
#' * \code{propensityScoreModel}: Model of treatment assignment, \code{NULL} if not provided and custom propensity weights are used.
#' * \code{propensityWeights}: Propensity weights used. When not \code{NULL}, it will be used as custom inputs from user. When \code{NULL}, it will be obtained by logistic regression.
#' * \code{participationWeights}: Participation Weights used. When not \code{NULL}, it will be used as custom inputs from user. When \code{NULL}, it will be obtained by Method of Moments-based method.
#' * \code{participationWeightCalculation}: The object contains participation weights calculation results.
#' * \code{participationWeightSummary}: The summary results related to participation weights such as effect sample sizes and so on as the information for the visualization legend.
#' * \code{finalWeights}: Weights used to fit MSM.
#' * \code{customPropensity}: Boolean indicating whether custom propensity weights are used.
#' * \code{customParticipation}: Boolean indicating whether custom participation weights are used.
#' * \code{treatment}: String indicating the variable name of treatment.
#' * \code{response}: String indicating the variable name of response.
#' * \code{studyData}: The original individual participant data (IPD) of study population.
#' * \code{aggregateTargetData}: The aggregate-level data (AgD) of target population.
#' * \code{processedIPD}: The individual participant data (IPD) after analytically pre-processing in format and calibration. All the ordinal variables have been dummized. \code{NULL} when \code{customParticipation} is \code{TRUE}.
#' * \code{processedAgD}: The aggregate-level data (AgD) after analytically pre-processing in format and calibration. \code{NULL} when \code{customParticipation} is \code{TRUE}.
#' * \code{matchingCovariates}: A vector of user-specified covariates which are used for the data matching to obtain participation weights when no custom weights provided.
#' * \code{centeredStudyData}: The data frame with both the processed study data and centered study data. 
#' * \code{exOpt}: Provided \code{exOpt} argument.
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
                          
                          method = c("logistic", "probit", "loglog", "cloglog", "cauchit"),
                          exOpt = list(propensity = NULL,
                                       participation = NULL,
                                       final = NULL),
                          
                          studyData, # data of study population (studyData): N rows, data.frame with responses and variables
                          aggregateTargetData  # data of target population (aggregateTargetData): 1 row, data.frame with only aggregate variables
                          
                          ) {
  # TODO: bugfix for 1 EM
  
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
    
    validCov_TF <- grepl(paste("^", validCov, sep = "", collapse = "|"), names(aggregateTargetData))
    # for the paste(): if the validCov being c(age, sex), then it gives regular expression ^age|^sex
    # for the grepl(): matches regex pattern within the names of columns in AgD and gives TRUE FALSE
    
    matchedAgDNames <- names(aggregateTargetData)[validCov_TF]
    
    prefixes <- unique(sub("_.*", "", matchedAgDNames))
    
    validCov <- validCov[validCov %in% prefixes]
    
    # If there are any covariates in the user input that don't match, give a warning and remove them
    if (length(validCov) != length(matchingCovariates)) {
      removedCov <- setdiff(matchingCovariates, validCov)
      cat("The following user-specfied matching covariates were not found in neither studyData or aggregateTargetData and have been removed:",
          toString(removedCov), "\n")
    }

    matchingCovariates <- validCov
    
    # Notice for user to double check the validity of matching covariates ready-to-use
    warning(cat("The following covariates are being used for matching:", toString(matchingCovariates),
                ". Please ensure that these covariates are meaningful", "\n"))
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
  if (is.null(propensityWeights)) {
    propensityWeights <- obtainWeights(propensityScoreModel, type = "probability")
    
    if (!is.null(exOpt$propensity)) {
      propensityOpt <- exOpt$propensity
      if (is.trunc(propensityOpt)) propensityWeights <- truncWeights(propensityWeights, propensityOpt)
      else if (is.trim(propensityOpt)) {
        retainedObsIndPropensity <- trimInd(propensityScoreModel, propensityOpt)
        propensityScoreModel <- stats::update(propensityScoreModel,
                                              formula. = propensityScoreModel$formula,
                                              data = cbind(propensityScoreModel$data,
                                                           data.frame(retainedObsIndPropensity = retainedObsIndPropensity)),
                                              weights = retainedObsIndPropensity)
        propensityWeights <- obtainWeights(propensityScoreModel, type = "probability") * retainedObsIndPropensity
      }
    }
  }
  
  
  # -------- Participation Weights -----------
  
  # detect exist ordinal variables: factor, more than 2 levels
  matchingCovariatesSubset <- studyData[matchingCovariates]
  
  # format checking before dummizing the ordinal variables, and ensure 0-1 numeric binary matching covariates
  binaryMatchingVariables <- c()
  
  matchingCovariatesSubset[] <- lapply(names(matchingCovariatesSubset), function(varName) {
    x <- matchingCovariatesSubset[[varName]]
    uniqueVals <- unique(x)
    
    if (length(uniqueVals) == 2) {
      # For matching variables like Gender: "Male" / "Female"
      # For matching variables like Gender but labelled as Male = "3" and Female = "5" (the label has to be 0-1) 
      if ((is.numeric(x) && !all(uniqueVals %in% c(0, 1))) || !is.numeric(x)) {
        binaryMatchingVariables <- c(binaryMatchingVariables, varName)
      }
    }
    
    # For ordinal variables like Season: "Spring" / "Summer" / "Fall" / "Winter"
    if (is.character(x) && length(uniqueVals) > 2) { x <- factor(x) }
    
    x
  })
  
  # Warning user to manually update inputs
  if (length(binaryMatchingVariables) > 0) {
    stop(sprintf("Please manually verify and possibly convert to 0-1 numeric scale the following binary matching covariates before further analysis: (%s)", 
                 paste(binaryMatchingVariables, collapse = ", ")))
  }
  
  
  ordinalCovariates <- names(matchingCovariatesSubset)[sapply(matchingCovariatesSubset, function(x) is.factor(x) && nlevels(x) > 2)]
  
  
  if(length(ordinalCovariates) > 0){
    dummizeCols = ordinalCovariates
    dummizeIPD <- dummizeIPD(studyData, 
                             dummizeCols = dummizeCols)
  } else {
    dummizeIPD <- studyData
  }
  
  processedAgD <- processedAgD(aggregateTargetData)
  centeredIPD <- centerIPD(IPD = dummizeIPD, AgD = processedAgD)
  
  centeredColnames <- grep("_CENTERED", colnames(centeredIPD), value = TRUE) # add suffix for centered variables
  

  # Track Custom Participation Weights
  customParticipation <- F
  if (!is.null(participationWeights)) {
    warning("Custom participation weights are being used. Please ensure that these weights are meaningful.")
    customParticipation <- T
    

    # ------- participationWeightCalculation based on custom participation weights

    wt <- participationWeights
    wt_rs <- (wt / sum(wt)) * sum(stats::complete.cases(matchingCovariatesSubset)) # avoid NAs in matchingCovariatesSubset
    data <- data.frame(weights = wt, 
                       scaledWeights = wt_rs)

    participationWeightCalculation <- list(
      data = data,
      ESS = sum(wt)^2 / sum(wt^2) # effective sample size
    )

    # ------- participationWeightSummary based on custom participation weights

    # participationWeightSummary <- calculateWeightsLegend(weightedData = participationWeightCalculation)
    summaryESS <- round(participationWeightCalculation$ESS, 2)
    ESS_reduction <- round((1 - summaryESS/length(wt)) * 100, 2)
    wt_median <- round(stats::median(wt), 4)
    wt_scaled_median <- round(stats::median(wt_rs), 4)
    
    participationWeightSummary <- list(
      ESS = summaryESS,
      ESS_reduction = ESS_reduction,
      wt_median = wt_median,
      wt_scaled_median = wt_scaled_median,
      nr_na = 0 # unable to get based on custom participation weights
    )

  }
  # obtain participation weights when no custom weights (pending in individual testing scripts)

  # obtain participation weights when no custom weights
  else {
    # ----- Participation Weight Estimation ----- #
    
    # # detect exist ordinal variables: factor, more than 2 levels
    # matchingCovariatesSubset <- studyData[matchingCovariates]
    # 
    # # format checking before dummizing the ordinal variables, and ensure 0-1 numeric binary matching covariates
    # binaryMatchingVariables <- c()
    # 
    # matchingCovariatesSubset[] <- lapply(names(matchingCovariatesSubset), function(varName) {
    #   x <- matchingCovariatesSubset[[varName]]
    #   uniqueVals <- unique(x)
    # 
    #   if (length(uniqueVals) == 2) {
    #     # For matching variables like Gender: "Male" / "Female"
    #     # For matching variables like Gender but labelled as Male = "3" and Female = "5" (the label has to be 0-1) 
    #     if ((is.numeric(x) && !all(uniqueVals %in% c(0, 1))) || !is.numeric(x)) {
    #       binaryMatchingVariables <- c(binaryMatchingVariables, varName)
    #     }
    #   }
    #   
    #   # For ordinal variables like Season: "Spring" / "Summer" / "Fall" / "Winter"
    #   if (is.character(x) && length(uniqueVals) > 2) { x <- factor(x) }
    #   
    #   x
    # })
    # 
    # # Warning user to manually update inputs
    # if (length(binaryMatchingVariables) > 0) {
    #   stop(sprintf("Please manually verify and possibly convert to 0-1 numeric scale the following binary matching covariates before further analysis: (%s)", 
    #                 paste(binaryMatchingVariables, collapse = ", ")))
    # }
    # 
    # 
    # ordinalCovariates <- names(matchingCovariatesSubset)[sapply(matchingCovariatesSubset, function(x) is.factor(x) && nlevels(x) > 2)]
    # 
    # 
    # if(length(ordinalCovariates) > 0){
    #   dummizeCols = ordinalCovariates
    #   dummizeIPD <- dummizeIPD(studyData, 
    #                            dummizeCols = dummizeCols)
    # }
    # 
    # 
    # 
    # processedAgD <- processedAgD(aggregateTargetData)
    # centeredIPD <- centerIPD(IPD = dummizeIPD, AgD = processedAgD)
    # 
    # centeredColnames <- grep("_CENTERED", colnames(centeredIPD), value = TRUE) # add suffix for centered variables
    
  
  # obtain participation weights when no custom weights (pending in individual testing scripts)
    # ----- Participation Weight Estimation ----- #
    
    # dummizeIPD <- dummizeIPD(studyData, ordinalCovariates)
    # processedAgD <- processedAgD(aggregateTargetData)
    # centerIPD <- centerIPD(IPD = dummizeIPD, AgD = processedAgD)
    # 
    # centeredColnames <- grep("_CENTERED", colnames(centerIPD), value = TRUE) # add suffix for centered variables
    
    participationWeightCalculation <- estimateWeights(data = centeredIPD,
                                                      centeredColnames = centeredColnames)
    
    participationWeights <- unlist(participationWeightCalculation$data["weights"])
    
    participationWeightSummary <- calculateWeightsLegend(weightedData = participationWeightCalculation)
    
  }
  
  if (!is.null(exOpt$participation)) {
    if (is.trunc(exOpt$participation)) participationWeights <- truncWeights(participationWeights, exOpt$participation)
    else if (is.trim(exOpt$participation)) warning("trim argument provided for participation weights. Note that only truncating participation weights is supported.")
  }

  # final weights: propensityWeights accounts for the confounding and participationWeights accounts for the EM distribution difference between populations
  finalWeights <-  propensityWeights * participationWeights
  
  if (!is.null(exOpt$final)) {
    if (is.trunc(exOpt$final)) finalWeights <- truncWeights(finalWeights, exOpt$final)
    else if (is.trim(exOpt$final)) warning("trim argument provided for final weights. Note that only truncating final weights is supported.")
  }
  
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

    family <- match.arg(family, c("coxph", "survreg", "polr"))
    
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
    } else if (family == "polr") {
      method <- match.arg(method, c("logistic", "probit", "loglog", "cloglog", "cauchit"))
      model <- MASS::polr(msmFormula, data = toAnalyze, weights = finalWeights, method = method)
      model$var <- sandwich::vcovBS(model)
    }
  } else {

      model <- stats::glm(msmFormula, 
                          family = family, 
                          data = toAnalyze, 
                          weight = finalWeights)
      model$var <- sandwich::vcovBS(model)
    }
  
  transportTADAResult <- list(msm = model, # glm / survreg / coxph / 
                              propensityScoreModel = propensityScoreModel, # glm / NULL

                              propensityWeights = propensityWeights,
                              participationWeights = participationWeights,
                              
                              participationWeightCalculation = participationWeightCalculation, # list, for plot
                              participationWeightSummary = participationWeightSummary, # list, for plot
                              
                              finalWeights = finalWeights,
                              
                              customPropensity = customPropensity, # T/F
                              customParticipation = customParticipation, #T/F
                              
                              treatment = treatment, # chr
                              response = response, # chr
                              
                              exOpt = exOpt,
                              
                              studyData = studyData, # data.frame
                              aggregateTargetData = aggregateTargetData, # data.frame
                              processedIPD = centeredIPD[,grep("_CENTERED$", names(centeredIPD), invert = TRUE, value = TRUE)],
                              processedAgD = processedAgD, # data.frame
                              
                              matchingCovariates = matchingCovariates, # chr
                              
                              centeredStudyData = centeredIPD # data.frame

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

# ----- Helper Functions for Participation Weights ----- #

# User input AgD preprocessing
processedAgD <- function(aggregateTargetData) {
  
  rawAgD <- as.data.frame(aggregateTargetData)
  
  # make all column names to be capital letters to avoid different style (Sex/sex to SEX)
  names(rawAgD) <- toupper(names(rawAgD))
  
  # must have info about total patients N
  if (!"N" %in% names(rawAgD)) {
    stop("Error: Column 'N' for total patient number is missing in the dataset.")
  }
  
  # must be without NA
  if (any(is.na(rawAgD))) {
    stop("Error: Ensure the input dataset has no NA before processing.")
  }
  
  # define column name patterns[-]
  legalSuffix <- c("MEAN", "MEDIAN", "SD", "COUNT", "PROP")
  
  candidateAgD <- names(rawAgD) 
  
  selectCondition1 <- grepl("_", candidateAgD, fixed = TRUE) # variables whose name with "_" would be labelled TRUE for selection
  selectCondition2 <- sapply(candidateAgD, function(candidateAgD) {
    tmpName <- unlist(strsplit(candidateAgD, split = "_"))
    tmpName[length(tmpName)] # this deployment is robust to the cases that there are multiple _ in the column name
  })
  selectCondition2 <- (selectCondition2 %in% legalSuffix) # the last split section of variables in other_colnames by "_"
  selectedCandidates <- candidateAgD[selectCondition1 & selectCondition2] # until this row: get all the column whose column names with "_" and legal suffix: c("MEAN", "MEDIAN", "SD", "COUNT", "PROP")
  
  useAgD <- rawAgD[, c("N", selectedCandidates), drop = FALSE] # selected AgD columns with suffix c("MEAN", "MEDIAN", "SD", "COUNT", "PROP")
  
  # warning about ignored columns with illegal naming
  if (!all(candidateAgD %in% selectedCandidates)) {
    warning(paste0(
      "Following columns are ignored since it does not follow the naming conventions:",
      paste(setdiff(candidateAgD, selectedCandidates), collapse = ",")
    ))
  }
  
  # calculate percentage columns
  countCondition <- grepl("_COUNT$", names(useAgD)) # columns in useAgD with names ending with "_COUNT"
  
  if (any(countCondition)) {
    for (i in which(countCondition)) {
      tmpProp <- useAgD[[i]] / useAgD$N # calculate the proportion with the count and N
      # in case some count are not specified, but proportion are specified, copy over those proportions
      # this also means, in case count is specified, proportion is ignored even it is specified
      propName <- gsub("_COUNT$", "_PROP", names(useAgD)[i])
      if (propName %in% names(useAgD)) {
        tmpProp[is.na(tmpProp)] <- useAgD[is.na(tmpProp), propName]
        names(useAgD)[names(useAgD) == propName] <- paste0(propName, "_redundant")
      }
      useAgD[[i]] <- tmpProp
    }
    names(useAgD) <- gsub("_COUNT$", "_PROP", names(useAgD)) 
  }
  
  useAgD <- useAgD[, !grepl("_redundant$", names(useAgD))] # until this row: always [replace] the proportion based on the new calculated value and remove the existed PROP data. The finalized AgD dataset just has proportion instead of count
  
  # output
  processedAgD <- useAgD
}

# User input IPD preprocessing
# note: dummizeIPD only takes care of ordinal cases. For binary, users have to apply 0-1 label by themselves to avoid the discussion about reference level
dummizeIPD <- function(rawIPD, dummizeCols, refLevel = NULL) { # refLevel: a vector of reference level of EACH ordinal variable in the same order of dummizeCols provided
  for (i in seq_along(dummizeCols)) {
    
    # deal with the ordinal variables needed to be dummized one by one
    rawCols <- rawIPD[[dummizeCols[i]]]
    rawLevels <- stats::na.omit(unique(rawCols))
    
    # dummizeCols is a vector with reference levels in order 
    if (is.null(refLevel) || length(refLevel) < i) { # length(refLevel) < i indicates that user doesn't provide the reference level of the i-th ordinal variables 
      rawCols <- factor(rawCols) # without provided reference level, factorize as default reference level
    } else {
      rawCols <- factor(as.character(rawCols), 
                        levels = c(refLevel[i], setdiff(rawLevels, refLevel[i])))
    }
    
    newCols <- sapply(levels(rawCols)[-1], function(j) {
      as.numeric(rawCols == j) # give 1 or 0 result based on T or F for each level
    })
    
    newCols <- as.data.frame(newCols)
    names(newCols) <- toupper(paste(dummizeCols[i], levels(rawCols)[-1], sep = "_")) # naming as VAR_LEVEL like "SEASON_SPRING" "SEASON_SUMMER"
    
    rawIPD <- cbind(rawIPD, newCols) # combine the new data into the raw IPD data to keep raw data clean
    
    rawIPD <- rawIPD[, !names(rawIPD) %in% dummizeCols, drop = FALSE] # remove the original ordinal variable columns
    
  }
  
  rawIPD # note: the dummized outputs are like (0,1) (1,0) and (0,0)
}

# Centerize the IPD
centerIPD <- function(IPD = dummizeIPD, AgD = processedAgD) {

  # need to remove the original ordinal variable column(s) as converting all variables from factor to numeric
  names(IPD) <- toupper(names(IPD))
  centeredIPD <- IPD
  
  legalSuffix <- c("MEAN", "MEDIAN", "SD", "PROP")
  suffixPat <- paste(paste0("_", legalSuffix, "$"), collapse = "|") # regex pattern: get the suffix out 
  
  useAgD <- AgD[, names(AgD) != "N", drop = FALSE] # no need the N anymore here
  
  # HERE FOR THE continuous variables with _mean and _sd will occur a multiplied results: PERCENTBODYFAT * 2
  paramID <- gsub(suffixPat, "", names(useAgD)) # a vector of variable names, select variables acccording to the regex pattern: SEX_MEAN to SEX
  # paramID <- setdiff(paramID, "N") # here we need to skip "N" since it is not a variable to be discussed
  
  for (j in seq_len(ncol(useAgD))) { 
    if (is.na(useAgD[[j]])) next # skip when NA
    
    # if the AgD of reference level is provided, due to the dummizing rule, skip this column for matching
    if (!( gsub(suffixPat, "", names(useAgD)[j]) %in% names(IPD) ))  next
    
    ipdParam <- paramID[j]  
    
    # factor: as.numeric() 0 1 to 1 2? [first to char then to numeric]
    
    if (grepl("_PROP$", names(useAgD)[j])) {
      centeredIPD[[paste0(ipdParam, "_", "CENTERED")]] <- as.numeric( as.character(IPD[[ipdParam]]) ) - useAgD[[j]]
    } else if (grepl("_MEAN$", names(useAgD)[j])) {
      centeredIPD[[paste0(ipdParam, "_", "CENTERED")]] <- as.numeric( IPD[[ipdParam]] ) - useAgD[[j]]
    } else if (grepl("_MEDIAN$", names(useAgD)[j])) {
      centeredIPD[[paste0(ipdParam, "_MEDIAN_", "CENTERED")]] <- ( as.numeric(IPD[[ipdParam]]) > useAgD[[j]]) - 0.5
    } else if (grepl("_SD$", names(useAgD)[j])) {
      centeredIPD[[paste0(ipdParam, "_SQUARED_", "CENTERED")]] <- ( as.numeric(IPD[[ipdParam]]) ^2) - (useAgD[[j]] ^2 + (useAgD[[paste0(ipdParam, "_MEAN")]] ^2))
    }
  }
  
  return(centeredIPD)
}

# achieve general weights estimation (applied in estimateWeights() )
optimiseWeights <- function(matrix,
                            par = rep(0, ncol(matrix)), # initial parameters
                            method = "BFGS",
                            maxit = 300, # maximum iteration times
                            trace = 0, # trace info outputs
                            ...) {
  if (!all( is.numeric(par) || is.finite(par), length(par) == ncol(matrix) ) ) {
    stop("Error: par must be a numeric vector with finite values of length equal to the number of columns in 'matrix'.")
  }
  
  optResults <- stats::optim(
    par = par,
    fn = function(alpha, X) sum(exp(X %*% alpha)), # function to be minimized
    gr = function(alpha, X) colSums(sweep(X, 1, exp(X %*% alpha), "*")), # gradient
    X = matrix,
    method = method,
    control = list(maxit = maxit, # maximum iteration times
                   trace = trace, ...) # trace info outputs
  )
  if (optResults$convergence != 0) { # convergence checking
    warning("optim() did not converge. ", optResults$message, "\nSee ?optim for more information on convergence code: ", optResults$convergence)
  }
  
  list(
    opt = optResults,
    alpha = optResults$par, # optimal parameters
    wt = exp(matrix %*% optResults$par) # weights
  )
}

# estimate matching weights for population matching
estimateWeights <- function(data,
                            centeredColnames = NULL,
                            startValue = 0,
                            method = "BFGS",
                            nBootIteration = NULL,
                            setSeedBoot = 20240725,
                            ...) {
  # pre check
  checkA <- is.data.frame(data)
  if (!checkA) {
    stop("Error: please ensure the 'data' input is a data.frame during participation weights estimation. hint: estimateWeights()")
  }
  
  checkB <- (!is.null(centeredColnames))
  if (checkB && is.numeric(centeredColnames)) { # when user inputs are column index
    
    checkB1 <- any(centeredColnames < 1 | centeredColnames > ncol(data))
    if (checkB1) { stop("Error: specified `centeredColnames` are out of bound during participation weights estimation. hint: estimateWeights()")}
  } 
  
  else if (checkB && is.character(centeredColnames)) { # when user inputs are column names
    
    checkB2 <- !all(centeredColnames %in% names(data))
    if (checkB2) { stop("Error: one or more specified `centeredColnames` are not found in 'data' input during participation weights estimation. hint: estimateWeights()")}
  } 
  
  else { stop("Error: 'centeredColnames' should be either a numeric or character vector. during participation weights estimation. hint: estimateWeights()")}
  
  # check the format of data itself. Needs to be all numeric for calculation
  checkC <- sapply(centeredColnames, function(ii) { !is.numeric(data[[ii]])})
  if (any(checkC)) {
    stop(paste0(
      "Error: following columns of 'data' are not numeric for the participation weights calculation:",
      paste(which(checkC), collapse = ",")
    ))
  }
  
  
  # prepare data for optimization
  if (is.null(centeredColnames)) 
    centeredColnames <- seq_len(ncol(data)) # an integer sequence with length of ncol(data)
  
  EM <- data[, centeredColnames, drop = FALSE] # effect modifiers
  naEM <- apply(EM, 1, function(EM) any(is.na(EM)))
  nMissing <- sum(naEM)
  rowsWithMissing <- which(naEM)
  EM <- as.matrix(EM[!naEM, , drop = FALSE])
  
  
  # estimate weights
  opt1 <- optimiseWeights(matrix = EM, 
                          par = rep(startValue, ncol(EM)), 
                          method = method, ...)
  
  alpha <- opt1$alpha
  wt <- opt1$wt
  wt_rs <- (wt / sum(wt)) * nrow(EM)
  
  
  # bootstrapping
  # outboot <- if (is.null(nBootIteration)) {
  #   bootSeed <- NULL
  #   bootStrata <- NULL
  #   NULL
  # } else {
  # Make sure to leave '.Random.seed' as-is on exit
  # oldSeed <- globalenv()$.Random.seed
  # on.exit(suspendInterrupts(set_random_seed(oldSeed)))
  # set.seed(setSeedBoot)
  # 
  # rowid_in_data <- which(!ind)
  # arms <- factor(data$ARM[rowid_in_data])
  # boot_statistic <- function(d, w) optimise_weights(d[w, ], 
  #                                                   par = alpha, 
  #                                                   method = method, ...)$wt[, 1]
  # boot_out <- boot(EM,
  #                  statistic = boot_statistic,
  #                  R = nBootIteration,
  #                  strata = arms)
  
  # boot_array <- array(dim = list(nrow(EM), 2, nBootIteration))
  # dimnames(boot_array) <- list(sampled_patient = NULL, 
  #                              c("rowid", "weight"), 
  #                              bootstrap_iteration = NULL)
  # boot_array[, 1, ] <- t(boot.array(boot_out, TRUE))
  # boot_array[, 2, ] <- t(boot_out$t)
  # bootSeed <- boot_out$seed
  # bootStrata <- boot_out$strata
  # boot_array
  #}
  
  # append weights to data
  data$weights <- NA
  data$weights[!naEM] <- wt
  
  data$scaledWeights <- NA
  data$scaledWeights[!naEM] <- wt_rs
  
  if (is.numeric(centeredColnames)) centeredColnames <- names(data)[centeredColnames]
  
  # Output
  outputs <- list(
    data = data,
    centeredColnames = centeredColnames,
    nMissing = nMissing,
    ESS = sum(wt)^2 / sum(wt^2), # effective sample size
    opt = opt1$opt,
    # boot = outboot,
    # bootSeed = bootSeed,
    # bootStrata = bootStrata,
    rowsWithMissing = rowsWithMissing
  )
  
  class(outputs) <- c("estimateWeights", "list")
  outputs
}

# legend for participation plot legend
calculateWeightsLegend <- function(weightedData) {
  
  ESS <- weightedData$ESS
  wt <- weightedData$data$weights
  wt_scaled <- weightedData$data$scaledWeights
  
  # calculate sample size and exclude NA from wt
  nr_na <- sum(is.na(wt))
  n <- length(wt) - nr_na
  wt <- stats::na.omit(wt)
  wt_scaled <- stats::na.omit(wt_scaled)
  
  # calculate ESS reduction and median weights
  ESS_reduction <- (1 - (ESS / n)) * 100
  wt_median <- stats::median(wt)
  wt_scaled_median <- stats::median(wt_scaled)
  
  list(
    ESS = round(ESS, 2),
    ESS_reduction = round(ESS_reduction, 2),
    wt_median = round(wt_median, 4),
    wt_scaled_median = round(wt_scaled_median, 4),
    nr_na = nr_na
  )
}

# function to calculate weighted SD based on weighted.mean() base function
weighted.sd <- function(x, w) {sqrt(sum(w * (x - stats::weighted.mean(x, w)) ^2, na.rm = TRUE) / sum(w)) }

#' @title Summarize results of a fitted MSM using the TADA approach
#' 
#' @description
#' Returns summary object which contains a summary of the fitted MSM, pre- and post-weighting standardized mean differences (SMDs) and pre-post matching table including aggregate level summary of effect modifiers of interest.
#' 
#' @rdname summary.transportTADA
#'
#' @param object Result from \code{transportTADA} function
#' @param covariates Vector of strings indicating names of covariates in propensity model
#' @param effectModifiers Vector of strings indicating names of effect modifiers in participation model. When it is \code{NULL}, the \code{prePostTable} shows all the results
#' @param ... Further arguments from previous function or to pass to next function
#'
#' @return
#' The \code{summary.transportTADA} function returns a \code{summary.transportTADA} object containing the following components:
#' * \code{propensitySMD}: Table of unweighted and weighted SMDs of covariates between treatment groups. Only propensity weights are used.
#' * \code{prePostTable}: Table of unweighted and weighted aggregate level summary as well as pre- and post- difference values of effect modifiers of interest between study data and target data.
#' * \code{msmSummary}: Summary object of model object for MSM. The correct variance estimators are set here for \code{glm}, whereas they are set in \code{transportTADA} for \code{survreg} and \code{coxph}.
#'
#' @export
summary.transportTADA <- function(object, 
                                  covariates = NULL, 
                                  effectModifiers = NULL, ...) { 
  transportTADAResult <- object
  
  treatment <- transportTADAResult$treatment
  matchingCovariates <- toupper(transportTADAResult$matchingCovariates)
  
  originalIPD <- transportTADAResult$studyData
  processedIPD <- transportTADAResult$processedIPD
  processedAgD <- transportTADAResult$processedAgD
  finalWeights <- transportTADAResult$finalWeights
  
  # customParticipation <- transportTADAResult$customParticipation # custom participation weights indicator: T/F
  
  # Calculate SMDs for covariates
  # Extract study data and weights
  if (!is.null(transportTADAResult$propensityScoreModel)) {
    
    propensityScoreModel <- transportTADAResult$propensityScoreModel
    propensityFormula <- propensityScoreModel$formula
    
    studyData <- propensityScoreModel$data
    
    if (is.null(covariates)) covariates <- all.vars(propensityFormula)[-1]
    
  } 
  else {
    
    studyData <- transportTADAResult$studyData
    aggregateTargetData <- transportTADAResult$aggregateTargetData
    
  }
  propensityWeights <- transportTADAResult$propensityWeights
  
  treatmentIndex <- which(names(studyData) == treatment)
  
  # Calculate SMDs for covariates (should optimize)
  if (length(covariates != 0)) {
    prePropensitySMD <- sapply(covariates, function (covariate) as.double(smd::smd(x = studyData[, which(names(studyData) == covariate)],
                                                                                   g = studyData[, treatmentIndex])$estimate))
    
    prePropensityBalance <- data.frame(variable = covariates, 
                                       smd = abs(prePropensitySMD), 
                                       method = rep("Observed", length(covariates)))
    
    postPropensitySMD <- sapply(covariates, function (covariate) as.double(smd::smd(x = studyData[, which(names(studyData) == covariate)],
                                                                                    g = studyData[, treatmentIndex],
                                                                                    w = propensityWeights)$estimate))
    
    postPropensityBalance <- data.frame(variable = covariates, 
                                        smd = abs(postPropensitySMD), 
                                        method = rep("Weighted", length(covariates)))
    
    propensityBalance <- rbind(prePropensityBalance, postPropensityBalance)
  } else propensityBalance <- NULL

  # -------- TADA: pre-post weighting result table ----------
  
  # effectModifiers -> matchingCovariates
  # just use the effectModifiers to control to the certain rows which of user interests (shown as a subset)
  # when effectModifiers is NULL, directly use matchingCovariates as all
  
  if (is.null(effectModifiers)) effectModifiers <- matchingCovariates # all upper cases already
  effectModifiers <- toupper(effectModifiers) # all upper cases
  
  # get the prefixes in rownames of processedAgD and processedIPD to match the effect modifiers as user inputs
  extractPrefix <- function(name) {sapply(strsplit(name, "_"), `[`, 1)}
  
  if(!is.null(processedAgD) & !is.null(processedIPD)){ 
    
    prefixesAgD <- extractPrefix(names(processedAgD))
    prefixesIPD <- extractPrefix(names(processedIPD))
    commonPrefixes <- intersect(prefixesAgD, prefixesIPD)
    
    commonPrefixes <- intersect(commonPrefixes, effectModifiers)
    
    prefixesPattern <- paste(commonPrefixes, collapse = "|")
    
    NAgD <- processedAgD$N
    NIPD <- nrow(processedIPD)
    
    userAgD <- processedAgD[, sapply(names(processedAgD), function(name) grepl(prefixesPattern, name)), drop = F]
    userIPD <- processedIPD[, sapply(names(processedIPD), function(name) grepl(prefixesPattern, name)), drop = F]
    
  }
  
  
  
  # ---- pre weighting IPD summary
  userIPDModified <- userIPD
  userIPDModified[] <- lapply(userIPD, function(x) if(is.factor(x)) as.numeric(as.character(x)) else x) # note: here the factors in userIPD could only be the binary or ordinal variables after dummized, thus are 0-1 scale
  userIPDSummaryList <- list()
  for (name in names(userIPDModified)) {
    if (all(userIPDModified[[name]] %in% c(0, 1))) {
      userIPDSummaryList[[paste(name, "PROP", sep = "_")]] <- mean(userIPDModified[[name]])
    } else if (is.numeric(userIPDModified[[name]])) {
      userIPDSummaryList[[paste(name, "MEAN", sep = "_")]] <- mean(userIPDModified[[name]])
      userIPDSummaryList[[paste(name, "SD", sep = "_")]] <- stats::sd(userIPDModified[[name]])
      userIPDSummaryList[[paste(name, "MEDIUM", sep = "_")]] <- stats::median(userIPDModified[[name]])
    }
  }
  
  # formatting
  userIPDSummaryDataframe <- as.data.frame(userIPDSummaryList)
  transposedUserIPDSummaryDataframe <- t(userIPDSummaryDataframe)
  transposedUserIPDSummaryDataframe <- as.data.frame(transposedUserIPDSummaryDataframe)
  transposedUserIPDSummaryDataframe$EMs <- rownames(transposedUserIPDSummaryDataframe)
  rownames(transposedUserIPDSummaryDataframe) <- NULL
  transposedUserIPDSummaryDataframe <- transposedUserIPDSummaryDataframe[c("EMs", names(transposedUserIPDSummaryDataframe)[1:(ncol(transposedUserIPDSummaryDataframe)-1)])]
  colnames(transposedUserIPDSummaryDataframe) <- c("EMs", "Results")
  
  summaryIPD <- transposedUserIPDSummaryDataframe
  
  # ---- AgD summary
  # formatting
  transposedUserAgDSummaryDataframe <- t(userAgD)
  transposedUserAgDSummaryDataframe <- as.data.frame(transposedUserAgDSummaryDataframe)
  transposedUserAgDSummaryDataframe$EMs <- rownames(transposedUserAgDSummaryDataframe)
  rownames(transposedUserAgDSummaryDataframe) <- NULL
  transposedUserAgDSummaryDataframe <- transposedUserAgDSummaryDataframe[c("EMs", names(transposedUserAgDSummaryDataframe)[1:(ncol(transposedUserAgDSummaryDataframe)-1)])]
  colnames(transposedUserAgDSummaryDataframe) <- c("EMs", "Results")
  
  summaryAgD <- transposedUserAgDSummaryDataframe
  
  # ---- post weighting IPD summary
  weights <- finalWeights
  
  # calculate weighted median based on Hmisc package
  # library(Hmisc)
  
  postWeightingIPDSummaryList <- list()
  
  for (name in names(userIPDModified)) {
    originalEMs <- userIPDModified[[name]]
    if (all(originalEMs %in% c(0, 1))) {
      postWeightingIPDSummaryList[[paste(name, "PROP", sep = "_")]] <- stats::weighted.mean(originalEMs, weights)
    } else if (is.numeric(originalEMs)) {
      postWeightingIPDSummaryList[[paste(name, "MEAN", sep = "_")]] <- stats::weighted.mean(originalEMs, weights)
      postWeightingIPDSummaryList[[paste(name, "SD", sep = "_")]] <- weighted.sd(originalEMs, weights)
      postWeightingIPDSummaryList[[paste(name, "MEDIUM", sep = "_")]] <- Hmisc::wtd.quantile(originalEMs, weights, probs = 0.5)
    }
  }
  
  # formatting
  postWeightingIPDSummaryDataframe <- as.data.frame(postWeightingIPDSummaryList)
  transposedPostWeightingIPDSummaryDataframe <- t(postWeightingIPDSummaryDataframe)
  transposedPostWeightingIPDSummaryDataframe <- as.data.frame(transposedPostWeightingIPDSummaryDataframe, stringsAsFactors = FALSE)
  transposedPostWeightingIPDSummaryDataframe$EMs <- rownames(transposedPostWeightingIPDSummaryDataframe)
  rownames(transposedPostWeightingIPDSummaryDataframe) <- NULL
  names(transposedPostWeightingIPDSummaryDataframe)[1] <- "Results"
  transposedPostWeightingIPDSummaryDataframe <- transposedPostWeightingIPDSummaryDataframe[c("EMs", names(transposedPostWeightingIPDSummaryDataframe)[1:(ncol(transposedPostWeightingIPDSummaryDataframe)-1)])]
  colnames(transposedPostWeightingIPDSummaryDataframe) <- c("EMs", "Results")
  
  summaryIPDWeighted <- transposedPostWeightingIPDSummaryDataframe
  
  # merge tables together
  merge <- merge(summaryIPD, summaryAgD, by = "EMs", all = FALSE)
  finalMerge <- merge(merge, summaryIPDWeighted, by = "EMs", all = FALSE)
  colnames(finalMerge) <- c("EMs", "IPD (pre-weighting)", "AgD", "IPD (post-weighting)")
  
  finalMerge[] <- lapply(finalMerge, function(x) if(is.numeric(x)) round(x, 3) else x)
  prePostTable <- finalMerge
  
  
  colnames(prePostTable) <- c("Effect Modifiers", 
                              paste0("Study ", "(N = ", NIPD, ")"), 
                              paste0("Target ", "(N = ", NAgD, ")"),
                              "Study (Post-Weighting)")
  
  # SMD-like summary data
  prePostTable$"Pre Weighting Difference" <- prePostTable[,3] - prePostTable[,2]
  prePostTable$"Post Weighting Difference" <- prePostTable[,3] - prePostTable$`Study (Post-Weighting)`
 
  # Final Table
  prePostTable <- prePostTable[,c("Effect Modifiers", 
                                  paste0("Study ", "(N = ", NIPD, ")"), 
                                  paste0("Target ", "(N = ", NAgD, ")"),
                                  "Study (Post-Weighting)",
                                  "Pre Weighting Difference",
                                  "Post Weighting Difference")]
  
  # --------------------------------------------------------- 

  # If model is glm, calculate and replace correct SEs
  
  msm <- transportTADAResult$msm
  
  msmSummary <- suppressWarnings(summary(msm))
  
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
  
  
  
  
  summaryTransportTADA <- list(propensitySMD = propensityBalance,
                               
                               prePostTable = prePostTable,
                               
                               msmSummary = msmSummary)
  
  class(summaryTransportTADA) <- "summary.transportTADA"
  
  return(summaryTransportTADA)
}


#' @rdname summary.transportTADA
#'
#' @param x \code{summary.transportTADA} object.
#' @param out Output stream.
#' @param ... Further arguments from previous function or to pass to next function
#'
#' @export
#'
print.summary.transportTADA <- function(x, out = stdout(), ...) {
  summaryTransportTADA <- x
  
  write("Absolute SMDs of covariates between treatments before and after weighting:", out)
  print(summaryTransportTADA$propensitySMD, out)
  
  write("Aggregate level data summary of effect modifiers of interests before and after weighting:", out)
  print(summaryTransportTADA$prePostTable, out)
  
  write("MSM results:", out)
  print(summaryTransportTADA$msmSummary, out)
}

#' @title Plot graphs relevant to transportability analysis using TADA
#'
#' @description 
#' Plot graphs for assessment of covariate balance and results in a TADA analysis. This function currently supports mirrored histograms, SMD plots and model coefficient plots.
#' 
#' @param x Result from \code{transportTADA} function
#' @param type One of \code{"propensityHist", "propensitySMD", "participationHist", "msm"}. \code{Hist} produces mirrored histograms of estimated probability of treatment between treatment groups (for \code{propensity}), or of estimated participation weights between study and target data (for \code{participation}). \code{SMD} produces SMD plots of covariates between treatment groups (for \code{propensity}). \code{msm} produces plots showing confidence intervals for the model coefficients, which should have the correct standard errors.
#' @param bins Number of bins for propensity score/participation weight histograms. This is only used for \code{Hist}.
#' @param maxWeight Maximum participation weight value, i.g., range, that user wants to show in the \code{participationHist}.
#' @param covariates Vector of strings indicating names of covariates in propensity model
#' @param effectModifiers Vector of strings indicating names of effect modifiers in participation model
#' @param ... Further arguments from previous function or to pass to next function
#'
#' @return A \code{ggplot} object which contains the desired plot.
#' 
#' @export
#'
#' @importFrom rlang .data
#' @importFrom stats density
plot.transportTADA <- function(x, type = "propensityHist", bins = 50, maxWeight = NULL, covariates = NULL, effectModifiers = NULL, ...) {
  transportTADAResult <- x
  summaryTransportTADA <- summary(transportTADAResult, covariates = covariates, effectModifiers = effectModifiers)
  resultPlot <- NULL
  
  # Match argument
  type <- match.arg(type, c("propensityHist", "propensitySMD", 
                            "participationHist", 
                            "msm"))
  
  if (type == "propensityHist") {
    # Mirrored histogram of propensity scores
    if (!transportTADAResult$customPropensity) {
      propensityModel <- transportTADAResult$propensityScoreModel
      studyData <- propensityModel$data
      studyData$propensityScore <- propensityModel$fitted.values
      treatmentVar <- transportTADAResult$treatment
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
    } 
    else 
      { stop("Custom propensity weights were used. Please plot your previously estimated propensity scores using the halfmoon package, if desired.")}
  } 
  
  else if (type == "propensitySMD") {
    # SMD plot of covariates
    propensitySMD <- summaryTransportTADA$propensitySMD
    if (is.null(propensitySMD)) stop("No covariates to plot, please provide some.")
    resultPlot <- ggplot2::ggplot(propensitySMD, ggplot2::aes(x = .data$variable, y = .data$smd, group = .data$method, color = .data$method)) + 
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::geom_hline(yintercept = 0.1) +
      ggplot2::coord_flip() +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.key = ggplot2::element_blank())
  } 

  else if (type == "participationHist") {
    
    weightedData <- transportTADAResult$participationWeightCalculation
    weightsStat <- transportTADAResult$participationWeightSummary
    
    main_title <- c("Participation Weights", "Scaled Participation Weights")
    bin_col <- "#2FB9AB"
    vline_col <- "red"
    bins <- bins

    wt_data0 <- weightedData$data[, c("weights", "scaledWeights")]
    colnames(wt_data0) <- main_title
    wt_data <- utils::stack(wt_data0)
    wt_data$median <- ifelse(wt_data$ind == main_title[1],
                             weightsStat$wt_median, 
                             weightsStat$wt_scaled_median
    )
    
    # Legend Data
    lab <- with(weightsStat, {
      lab <- c(paste0("Median = ", wt_median), paste0("Median = ", wt_scaled_median))
      lab <- paste0(lab, "\nESS = ", ESS, "\nReduction% = ", ESS_reduction)
      if (nr_na > 0) lab <- paste0(lab, "\n#Missing Weights = ", nr_na)
      lab
    })
    legend_data <- data.frame(ind = main_title, lab = lab)
    
    resultPlot <- ggplot2::ggplot(wt_data) +
      ggplot2::geom_histogram(ggplot2::aes(x = .data[["values"]]), bins = bins, color = bin_col, fill = bin_col) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = .data[["median"]]),
                          color = vline_col,
                          linetype = "dashed") +
      ggplot2::theme_bw() +
      ggplot2::facet_wrap(~ind, ncol = 1) +
      ggplot2::geom_text(data = legend_data,
        ggplot2::aes(label = .data[["lab"]]), x = Inf, y = Inf, hjust = 1, vjust = 1, size = 3) +
      ggplot2::theme(
        axis.title = ggplot2::element_text(size = 12),
        axis.text = ggplot2::element_text(size = 12)) +
      ggplot2::ylab("Frequency") +
      ggplot2::xlab("Weight")
    
    if (!is.null(maxWeight)) {resultPlot <- resultPlot + ggplot2::xlim(c(0, maxWeight))}
  } 
  
  else if (type == "msm") {
    # Coefficient plots
    resultPlot <- modelsummary::modelplot(transportTADAResult$msm,
                                          vcov = list(transportTADAResult$msm$var),
                                          coef_omit = "(Intercept)")
  }
  
  return(resultPlot)
}



#' @title Check validity of TADA result object
#'
#' @description 
#' A simple helper function that validates whether the components of the given \code{transportTADA} object are of the correct types.
#'
#' @param transportTADAResult Result object from \code{transportTADA} function
#'
#' @return A boolean indicating whether all components of \code{transportTADA} object have the correct types.
#' 
#' @export
#'
is.transportTADA <- function (transportTADAResult) {
  return((inherits(transportTADAResult$msm, "glm") | inherits(transportTADAResult$msm, "coxph") | inherits(transportTADAResult$msm, "survreg")) &
           
         (inherits(transportTADAResult$propensityScoreModel, "glm") | is.null(transportTADAResult$propensityScoreModel)) & 
           
         inherits(transportTADAResult$propensityWeights, "numeric") & 
         inherits(transportTADAResult$participationWeights, "numeric") &

         inherits(transportTADAResult$participationWeightSummary, "list")  & 
         inherits(transportTADAResult$participationWeightCalculation, "list")  & 
           
         inherits(transportTADAResult$finalWeights, "numeric") &
           
         inherits(transportTADAResult$customPropensity, "logical") &
         inherits(transportTADAResult$customParticipation, "logical") &
           
         (inherits(transportTADAResult$treatment, "character") | is.null(transportTADAResult$treatment)) &
         (inherits(transportTADAResult$response, "character") | is.null(transportTADAResult$response)) &
           
         inherits(transportTADAResult$studyData, "data.frame")  &
         inherits(transportTADAResult$aggregateTargetData, "data.frame")  &
         (inherits(transportTADAResult$processedIPD, "data.frame") | is.null(transportTADAResult$processedIPD)) &
         (inherits(transportTADAResult$processedAgD, "data.frame") | is.null(transportTADAResult$processedAgD)) &
           
         inherits(transportTADAResult$matchingCovariates, "character")  &

         (inherits(transportTADAResult$centeredStudyData, "data.frame") | is.null(transportTADAResult$centeredStudyData)) &
         
         inherits(transportTADAResult, "transportTADA"))
}
