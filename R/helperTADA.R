# ----- Helper Functions for TADA ----- #

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
    rawLevels <- stats::na.omit(unique(rawCols))
    
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