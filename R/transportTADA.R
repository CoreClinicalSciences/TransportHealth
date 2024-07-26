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
    
    # ----- Helper Functions ----- #
    
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
        tmp[length(tmpName)] # this deployment is robust to the cases that there are multiple _ in the column name
      })
      selectCondition2 <- (selectCondition2 %in% legalSuffix) # the last split section of variables in other_colnames by "_"
      selectedCandidates <- candidateAgD[selectCondition1 & selectCondition2] # until this row: get all the column whose column names with "_" and legal suffix: c("MEAN", "MEDIAN", "SD", "COUNT", "PROP")
      
      useAgD <- rawAgD[, c(selectedCandidates), drop = FALSE] # selected AgD columns with suffix c("MEAN", "MEDIAN", "SD", "COUNT", "PROP")
      
      # warning about ignored columns with illegal naming
      if (!all(candidateAgD %in% selectedCandidates)) {
        warning(paste0(
          "following columns are ignored since it does not follow the naming conventions:",
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
      useAgD
    }
    
    # User input IPD preprocessing
    dummizeIPD <- function(rawIPD, dummizeCols, refLevel = NULL) {
      for (i in seq_along(dummizeCols)) {
        
        rawCols <- rawIPD[[dummizeCols[i]]]
        rawLevels <- na.omit(unique(rawCols))
        
        # dummizeCols is a vector with reference levels in order 
        if (is.null(refLevel) || length(refLevel) < i) {
          rawCols <- factor(rawCols)
        } else {
          rawCols <- factor(as.character(rawCols), 
                            levels = c(refLevel[i], setdiff(rawLevels, refLevel[i])))
        }
        
        newCols <- sapply(levels(rawCols)[-1], function(j) {
          as.numeric(rawCols == j)
        })
        
        newCols <- as.data.frame(newCols)
        names(newCols) <- toupper(paste(dummizeCols[i], levels(rawCols)[-1], sep = "_")) # naming as VAR_LEVEL like "SEASON_SPRING" "SEASON_SUMMER"
        
        rawIPD <- cbind(rawIPD, newCols) # combine the new data into the raw AgD data to keep raw data clean
      }
      
      rawIPD
    }
    
    # Centerize the IPD
    centerIPD <- function(IPD, AgD) {
      centeredIPD <- IPD
      
      legalSuffix <- c("MEAN", "MEDIAN", "SD", "PROP")
      suffixPat <- paste(paste0("_", legalSuffix, "$"), collapse = "|") # regex pattern: get the suffix out 
      
      useAgD <- AgD
      paramID <- gsub(suffixPat, "", names(useAgD)) # a vector of variable names, select variables acccording to the regex pattern: SEX_MEAN to SEX
      
      for (j in seq_len(ncol(useAgD))) { 
        if (is.na(useAgD[[j]])) next # skip when NA
        
        ipdParam <- paramID[j]
        
        if (grepl("_MEAN$|_PROP$", names(useAgD)[j])) {
          centeredIPD[[paste0(ipdParam, "_", "CENTERED")]] <- IPD[[ipdParam]] - useAgD[[j]]
        } else if (grepl("_MEDIAN$", names(useAgD)[j])) {
          centeredIPD[[paste0(ipdParam, "_MEDIAN_", "CENTERED")]] <- (IPD[[ipdParam]] > useAgD[[j]]) - 0.5
        } else if (grepl("_SD$", names(useAgD)[j])) {
          centeredIPD[[paste0(ipdParam, "_SQUARED_", "CENTERED")]] <- (IPD[[ipdParam]] ^2) - (useAgD[[j]] ^2 + (useAgD[[paste0(ipdParam, "_MEAN")]] ^2))
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
      
      optResults <- optim(
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
        stop("Error: please ensure the 'data' input is a data.frame.")
      }
      
      checkB <- (!is.null(centeredColnames))
      if (checkB && is.numeric(centeredColnames)) { # when user inputs are column index
        
        checkB1 <- any(centeredColnames < 1 | centeredColnames > ncol(data))
        if (checkB1) { stop("Error: specified `centeredColnames` are out of bound.")}
      } 
      
      else if (checkB && is.character(centeredColnames)) { # when user inputs are column names
        
        checkB2 <- !all(centeredColnames %in% names(data))
        if (checkB2) { stop("Error: one or more specified `centeredColnames` are not found in 'data'.")}
      } 
      
      else { stop("Error: 'centeredColnames' should be either a numeric or character vector.")}
      
      # check the format of data itself. Needs to be all numeric for calculation
      checkC <- sapply(centeredColnames, function(ii) { !is.numeric(data[[ii]])})
      if (any(checkC)) {
        stop(paste0(
          "Error: following columns of 'data' are not numeric for the calculation:",
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
        rowsWithMissing = rowsWithMissing
      )
      
      class(outputs) <- c("estimateWeights", "list")
      outputs
    }
    
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
  # Variance of coeffs are corrected in this method for coxph() and survreg() with sandwich::vcovBS(model)
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
    # if (!is.glm.family(family)) {
    #   stop("Please check the 'family' input. Ensure it is either a character specifying 'coxph' or 'survreg', or a valid glm family object like stats::gaussian().")
    # }
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

# Helper function that check the validity of user input family for glm()
# is.glm.family <- function(family) {
#   required_components <- c("family", "link", "linkinv", "linkfun", "variance")
#   all(required_components %in% names(family))
# }


# Dummy Variables
# Since the dplyr and tidyr updates frquently, try to use some alternative instead
# |> instead of %>% for the R package
# processedData <- studyData
# 
# for (var in studyVars) {
#   is_ordinal <- length(unique(studyData[[var]])) > 2 & length(unique(studyData[[var]])) <= 8
#   
#   if (is_ordinal) {
#     tempData <- studyData %>%
#       select(ID, all_of(var)) %>%
#       mutate(!!var := factor(!!sym(var))) %>%
#       pivot_wider(names_from = !!sym(var), values_from = !!sym(var),
#                   values_fn = list(~ if_else(!is.na(.), 1, 0)),
#                   values_fill = list(~ 0)) %>%
#       rename_with(~paste0(var, "_", .), -ID) 
#     
#     processedData <- processedData %>%
#       left_join(tempData, by = "ID")
#   }
# }
# 
# ordinalVars <- studyData |>
#   dplyr::summarise(across(everything(), ~length(unique(.)) > 2 & length(unique(.)) <= 8)) |>
#   dplyr::select_if(is.logical) |>
#   names()
# 
# # When there is at least one (non-binary) ordinal variable: convert it to dummy variable and name them by: variable_category
# if (length(ordinalVars) > 0) {
#   ordinalVarsDummies <- studyData |>
#     dplyr::mutate(across(all_of(ordinalVars), ~factor(.))) |>
#     tidyr::pivot_wider(names_from = all_of(ordinalVars), values_from = all_of(ordinalVars),
#                 values_fn = list(~ if_else(!is.na(.), 1, 0)),
#                 values_fill = list(~ 0))
#} 





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
