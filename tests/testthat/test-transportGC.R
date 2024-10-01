test_that("Outcome model provided as formula", {
  set.seed(20240918)
  testData <- generateTestData()
  expect_no_error(preparedModel <- transportGCPreparedModel(outcomeModel = sysBloodPressure ~ med1 + sex + stress + percentBodyFat + med2 + med1:stress + med1:med2,
                                            treatment = "med1",
                                            family = stats::gaussian,
                                            studyData = testData$studyData))
  
  expect_true(is.transportGCPreparedModel(preparedModel))
  
  expect_no_error(transportGCResult <- transportGC("meanDiff",
                                   preparedModel,
                                   testData$targetData))
  
  expect_true(inherits(transportGCResult, "transportGC"))
  
  expect_no_error(preparedModel <- transportGCPreparedModel(outcomeModel = sysBloodPressure ~ med1 + sex + stress + percentBodyFat + med2 + med1:stress + med1:med2,
                                                            treatment = "med1",
                                                            family = stats::gaussian,
                                                            studyData = testData$studyData,
                                                            wipe = F))
  
  expect_true(is.transportGCPreparedModel(preparedModel))
  
  expect_no_error(transportGCResult <- transportGC("meanDiff",
                                                   preparedModel,
                                                   testData$targetData))
  
  expect_true(is.transportGC(transportGCResult))
  
  expect_no_error(transportGCSummary <- summary(transportGCResult))
  
  expect_true(is.numeric(transportGCSummary$effect))
  expect_true(is.numeric(transportGCSummary$var))
  expect_true(inherits(transportGCSummary$preparedModelSummary, c("summary.glm", "summary.survreg", "summary.coxph", "summary.polr")))
  expect_true(is.character(transportGCSummary$response))
  expect_true(is.character(transportGCSummary$treatment))
  expect_true(is.character(transportGCSummary$treatmentLevels))
  
  expect_no_error(transportGCPlot <- plot(transportGCResult))
  
  expect_true(ggplot2::is.ggplot(transportGCPlot))
})

test_that("Outcome model provided as glm", {
  testData <- generateTestData()
  
  outcomeModel <- glm(sysBloodPressure ~ med1 + sex + stress + percentBodyFat + med2 + med1:stress + med1:med2,
                      data = testData$studyData,
                      family = gaussian)
  
  expect_no_error(preparedModel <- transportGCPreparedModel(outcomeModel = outcomeModel,
                                                            treatment = "med1"))
  
  expect_true(is.transportGCPreparedModel(preparedModel))
  
  expect_no_error(transportGCResult <- transportGC("meanDiff",
                                                   preparedModel,
                                                   testData$targetData))
  
  expect_true(inherits(transportGCResult, "transportGC"))
  
  expect_no_error(preparedModel <- transportGCPreparedModel(outcomeModel = outcomeModel,
                                                            treatment = "med1",
                                                            wipe = F))
  
  expect_true(is.transportGCPreparedModel(preparedModel))
  
  transportGCResult <- transportGC("meanDiff",
                                                   preparedModel,
                                                   testData$targetData)
  
  expect_true(is.transportGC(transportGCResult))
  
  expect_no_error(transportGCSummary <- summary(transportGCResult))
  
  expect_true(is.numeric(transportGCSummary$effect))
  expect_true(is.numeric(transportGCSummary$var))
  expect_true(inherits(transportGCSummary$preparedModelSummary, c("summary.glm", "summary.survreg", "summary.coxph", "summary.polr")))
  expect_true(is.character(transportGCSummary$response))
  expect_true(is.character(transportGCSummary$treatment))
  expect_true(is.character(transportGCSummary$treatmentLevels))
  
  expect_no_error(transportGCPlot <- plot(transportGCResult))
  
  expect_true(ggplot2::is.ggplot(transportGCPlot))
})

test_that("Other outcomes", {
  set.seed(20240918)
  testData <- generateTestData()
  studyData <- testData$studyData
  targetData <- testData$targetData
  preparedModel <- transportGCPreparedModel(outcomeModel = ht ~ med1 + sex + stress + percentBodyFat + med2 + med1:stress + med1:med2,
                                                            treatment = "med1",
                                                            family = binomial,
                                                            studyData = studyData,
                                            wipe = F)
  
  expect_true(is.transportGCPreparedModel(preparedModel))
  
  expect_no_error(transportGCResult <- transportGC("or",
                                                   preparedModel,
                                                   targetData))
  
  expect_true(inherits(transportGCResult, "transportGC"))
  
  expect_no_error(transportGCResult <- transportGC("rr",
                                                   preparedModel,
                                                   targetData))
  
  expect_true(inherits(transportGCResult, "transportGC"))
  
  preparedModel <- transportGCPreparedModel(outcomeModel = htStage ~ med1 + sex + stress + percentBodyFat + med2 + med1:stress + med1:med2,
                                                            treatment = "med1",
                                                            family = "polr",
                                                            studyData = studyData,
                                            wipe = F)
  
  expect_true(is.transportGCPreparedModel(preparedModel))
  
  expect_no_error(transportGCResult <- transportGC("or",
                                                   preparedModel,
                                                   targetData))
  
  expect_true(inherits(transportGCResult, "transportGC"))
  
  # preparedModel <- transportGCPreparedModel(outcomeModel = survival::Surv(os, rep(1, nrow(testData$studyData))) ~ med1 + sex + stress + percentBodyFat + med2 + med1:stress + med1:med2,
  #                                                           treatment = "med1",
  #                                                           family = "coxph",
  #                                                           studyData = studyData)
  # 
  # expect_true(is.transportGCPreparedModel(preparedModel))
  # 
  # transportGCResult <- transportGC("hr",
  #                                                  preparedModel,
  #                                                  targetData, bootstrapNum = 1)
  # 
  # expect_true(inherits(transportGCResult, "transportGC"))
  
  preparedModel <- transportGCPreparedModel(outcomeModel = survival::Surv(os, rep(1, nrow(testData$studyData))) ~ med1 + sex + stress + percentBodyFat + med2 + med1:stress + med1:med2,
                                                            treatment = "med1",
                                                            family = "survreg",
                                                            studyData = studyData,
                                            wipe = F)
  
  expect_true(is.transportGCPreparedModel(preparedModel))
  
  expect_no_error(transportGCResult <- transportGC("hr",
                                                   preparedModel,
                                                   targetData))
  
  expect_true(inherits(transportGCResult, "transportGC"))
})
