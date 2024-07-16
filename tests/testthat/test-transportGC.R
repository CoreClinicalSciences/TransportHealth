test_that("Outcome model provided as formula", {
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
