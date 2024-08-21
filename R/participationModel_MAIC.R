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

test_data <- generateTestData_MAIC()
IPD <- test_data$IPD
AgD <- test_data$AgD

participationModel <- participation ~ stress + med2

if (inherits(participationModel, "formula")) {
  
  participationModel <- glm(participationModel, family = binomial(), data = IPD)
} else if (inherits(participationModel, "glm")) {
  
  participationModel <- participationModel
} else {
  stop("participationModel must be either a formula or a fitted glm object!")
}

# Define the function to compute weights based on the model parameters beta and participationModel
obtainParticipationWeights <- function(beta, participationModel) {
  # Define the variables that modify the effect
  effectModifiers <- all.vars(participationModel)[-1]
  participation <- "participation"  # Specify the name for the participation indicator
  
  # Check and ensure the IPD dataset contains the participation indicator
  if (!(participation %in% names(IPD))) {
    IPD$participation <- 1  # Assume participation if missing
  }
  # Check and ensure the AgD dataset contains the participation indicator
  if (!(participation %in% names(AgD))) {
    AgD$participation <- 0  # Assume non-participation if missing
  }
  
  # Extract relevant data from IPD and AgD
  IPDParticipationData <- IPD[, names(IPD) %in% c(effectModifiers, participation)]
  participationData$participation <- as.factor(participationData$participation)  # Convert participation to a factor type
  
  # Update and fit the model, construct the design matrix X
  participationModel <- update(participationModel, data = participationData)
  X <- model.matrix(participationModel)
  weights <- exp(X %*% beta)  # Calculate weights using the formula
  
  # Compute target means and the weighted squared differences as the optimization objective
  target_means <- colMeans(X[-1, , drop = FALSE])  # Calculate the means for the target data
  objective_value <- sum((X - target_means)^2 * weights)  # Calculate the weighted squared differences
  
  # Return weights instead of the objective function value
  return(list(weights = weights, 
              objective_value = objective_value))
}


# Initialize beta parameters
initial_beta <- rep(0, length(coef(participationModel)))

# Use the optimization function to estimate the optimal beta parameters
result <- optim(initial_beta, 
                obtainParticipationWeights, 
                participationModel = participationModel, 
                method = "BFGS")

# Output the estimated beta values and the calculated weights
print("Estimated beta values:")
print(result$par)
print("Calculated weights for each individual:")
print(result$value$weights)
