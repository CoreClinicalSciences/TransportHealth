set.seed(20240501)

expit <- function(x) 1/(1+exp(-x))

generateTestData <- function() {
  # Generate study data
  nStudy <- 1000
  sexStudy <- rbinom(nStudy, 1, 0.5) # Male is 1, so female is baseline
  stressStudy <- rbinom(nStudy, 1, 0.4) # Stressed is 1
  med2Study <- rbinom(nStudy, 1, 0.1) # 1 means taking other med
  percentBodyFatStudy <- rnorm(nStudy, 28 - 13 * sexStudy, 2)
  med1Study <- rbinom(nStudy, 1, expit(0.2 * sexStudy - 0.02 * percentBodyFatStudy + 0.1 * stressStudy))
  sysBloodPressureStudy <- rnorm(nStudy, 100 + 5 * sexStudy + 0.5 * percentBodyFatStudy + 5 * stressStudy -
                                   5 * med1Study + med1Study * (-5 * med2Study + 7 * stressStudy))
  
  # Put all variables together
  studyData <- data.frame(sysBloodPressure = sysBloodPressureStudy, med1 = as.factor(med1Study), sex = as.factor(sexStudy), stress = as.factor(stressStudy), med2 = as.factor(med2Study), percentBodyFat = percentBodyFatStudy)
  
  # Generate target data
  nTarget <- 1500
  sexTarget <- rbinom(nTarget, 1, 0.3) # Male is 1, so female is baseline
  stressTarget <- rbinom(nTarget, 1, 0.7) # Stressed is 1
  med2Target <- rbinom(nTarget, 1, 0.3) # 1 means taking other med
  percentBodyFatTarget <- rnorm(nTarget, 26 - 12 * sexTarget, 2)
  
  # Put all variables together
  targetData <- data.frame(sex = as.factor(sexTarget), stress = as.factor(stressTarget), med2 = as.factor(med2Target), percentBodyFat = percentBodyFatTarget)
  
  return(list(studyData = studyData, targetData = targetData))
}

testData <- generateTestData()
studyData <- testData$studyData
targetData <- testData$targetData

usethis::use_data(studyData, overwrite = TRUE)
usethis::use_data(targetData, overwrite = TRUE)