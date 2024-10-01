generateTestDataInterpolated <- function() {
  testData <- generateTestData()
  studyData <- testData$studyData
  targetData <- testData$targetData
  studyData$percentBodyFatDicho <- as.numeric(studyData$percentBodyFat >= 17)
  targetData$percentBodyFatDicho <- as.numeric(targetData$percentBodyFat >= 17)
  
  nStudy <- nrow(studyData)
  
  mainTreatmentEffect <- with(studyData, mean(sysBloodPressure[med1 == 1]) - mean(sysBloodPressure[med1 == 0]))
  mainSE <- with(studyData, sqrt((var(sysBloodPressure[med1 == 1]) * (length(which(med1 == 1)) - 1) +
                            var(sysBloodPressure[med1 == 0]) * (length(which(med1 == 0)) - 1)) / (nStudy - 2)) * 
                            sqrt((1/length(which(med1 == 1))) + (1/length(which(med1 == 0)))))
  
  
  effectModifiers <- c("med2", "percentBodyFatDicho")
  
  aggregateStudyData <- apply(studyData[, effectModifiers], 2, function(x) {x |> as.numeric() |> mean()})
  
  subgroupLevels <- rep(c(1,0), times = length(effectModifiers))
  
  subgroupTreatmentEffects <- sapply(1:length(subgroupLevels), function(i) {
                              with(studyData[studyData[[effectModifiers[ceiling(i/2)]]] == subgroupLevels[i],],
                                   mean(sysBloodPressure[med1 == 1]) - mean(sysBloodPressure[med1 == 0])) 
                              })
  
  subgroupSEs <- sapply(1:length(subgroupLevels), function(i) {
    with(studyData[studyData[[effectModifiers[ceiling(i/2)]]] == subgroupLevels[i],],
         sqrt((var(sysBloodPressure[med1 == 1]) * (length(which(med1 == 1)) - 1) +
              var(sysBloodPressure[med1 == 0]) * (length(which(med1 == 0)) - 1)) / (nStudy - 2)) * 
           sqrt((1/length(which(med1 == 1))) + (1/length(which(med1 == 0))))
    )
  })
  
  targetData$n <- nrow(targetData)
  
  aggregateTargetData <- apply(targetData, 2, function(x) mean(as.numeric(x)))
  names(aggregateTargetData) <- names(targetData)
  
  return(list(mainTreatmentEffect = mainTreatmentEffect,
              mainSE = mainSE,
              effectModifiers = effectModifiers,
              subgroupTreatmentEffects = subgroupTreatmentEffects,
              subgroupSEs = subgroupSEs,
              nStudy = nStudy,
              aggregateStudyData = aggregateStudyData,
              targetData = targetData,
              aggregateTargetData = aggregateTargetData))
}