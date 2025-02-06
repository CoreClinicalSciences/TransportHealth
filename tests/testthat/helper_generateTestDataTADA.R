expit <- function(x) 1/(1+exp(-x))

# ---------- Generate test data for TADA for TransportHealth ----------
generateTestDataTADA <- function() {
  
  # ---------- Generate study data in IPD ----------
  nStudy <- 1000
  
  sexStudy <- rbinom(nStudy, 1, 0.5) # Male is 1, so female is baseline
  stressStudy <- rbinom(nStudy, 1, 0.4) # Stressed is 1
  med2Study <- rbinom(nStudy, 1, 0.1) # 1 means taking other med
  percentBodyFatStudy <- rnorm(nStudy, 28 - 13 * sexStudy, 2)
  toxicGradeStudy <- sample(c("Low", "Medium", "High"), size = nStudy, replace = TRUE, prob = c(0.7, 0.15, 0.15))
  
  # treatment
  med1Study <- rbinom(nStudy, 1, expit(0.2 * sexStudy - 0.02 * percentBodyFatStudy + 0.1 * stressStudy))
  # response
  sysBloodPressureStudy <- rnorm(nStudy, 100 + 5 * sexStudy + 0.5 * percentBodyFatStudy + 5 * stressStudy -
                                   5 * med1Study + med1Study * (-5 * med2Study + 7 * stressStudy))
  
  htStudy <- sysBloodPressureStudy >= 114
  htStageStudy <- integer(nStudy)
  for (i in 1:nStudy) {
    if (sysBloodPressureStudy[i] <= 107) htStageStudy[i] <- 1
    else if (sysBloodPressureStudy[i] <= 114) htStageStudy[i] <- 2
    else htStageStudy[i] <- 3
  }
  htStageStudy <- as.factor(htStageStudy)
  osStudy <- rexp(nStudy, 0.01 * exp(5 * sexStudy + 0.5 * percentBodyFatStudy + 5 * stressStudy -
                                       5 * med1Study + med1Study * (-5 * med2Study + 7 * stressStudy)))
  
  
  # Put all variables together
  studyData <- data.frame( sysBloodPressure = sysBloodPressureStudy, # response
                           med1 = as.factor(med1Study), # treatment
                           
                           sex = as.factor(sexStudy), 
                           stress = as.factor(stressStudy), 
                           med2 = as.factor(med2Study), 
                           toxicGrade = as.factor(toxicGradeStudy),
                           percentBodyFat = percentBodyFatStudy,
                           ht = htStudy, htStage = htStageStudy, os = osStudy)
  
  
  
  # ---------- Generate target data in AgD ----------
  nTarget <- 1500
  
  sexTarget <- rbinom(nTarget, 1, 0.3) # Male is 1, so female is baseline
  stressTarget <- rbinom(nTarget, 1, 0.7) # Stressed is 1
  med2Target <- rbinom(nTarget, 1, 0.3) # 1 means taking other med
  percentBodyFatTarget <- rnorm(nTarget, 26 - 12 * sexTarget, 2)
  toxicGradeTarget <- sample(c("Low", "Medium", "High"), size = nTarget, replace = TRUE, prob = c(0.45, 0.35, 0.20))

  # Put all variables together
  targetData <- data.frame(N = nTarget,
                           
                           # binary
                           sex_COUNT = sum(sexTarget), # count of Male
                           sex_PROP = 0.93, # an incorrect value to test the replacement function, the new prop should be calculated based on 
                           
                           stress_PROP = sum(stressTarget)/nTarget, # proportion of stress
                           
                           med2_PROP = sum(med2Target)/nTarget, 
                           
                           # ordinal with three subclass
                           toxicGrade_LOW_PROP = as.numeric(table(toxicGradeTarget)["Low"]) / nTarget,
                           toxicGrade_Medium_COUNT = as.numeric(table(toxicGradeTarget)["Medium"]),
                           # We don't need to prop for reference level candidate
                           # like: for male and female we just need to provide one of two and it's informative
                           # OR: user could provide all AgD for all levels and we could detect the ref level by the first element of factor IPD
                           toxicGrade_high_PROP = as.numeric(table(toxicGradeTarget)["High"]) / nTarget,
                           
                           # continuous
                           percentBodyFat_MEAN = mean(percentBodyFatTarget),
                           percentBodyFat_SD = sd(percentBodyFatTarget),
                           percentBodyFat_MEDIAN = median(percentBodyFatTarget)
                           )
  

  return(list(studyData = studyData, 
              aggregateTargetData = targetData))
}



