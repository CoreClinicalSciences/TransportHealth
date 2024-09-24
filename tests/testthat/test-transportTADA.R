test_that("Scenario 1: without custom weights, formula provided for propensityScoreModel", {
  set.seed(20240730)
  data <- generateTestDataTADA()
  
  expect_no_error(testResult <- suppressWarnings(transportTADA(msmFormula = sysBloodPressure ~ med1,
                                                               propensityScoreModel = med1 ~ sex + percentBodyFat + stress,
                                                
                                                               matchingCovariates = c("sex", "stress", "med2", "toxicGrade", "percentBodyFat"),
                                                
                                                               propensityWeights = NULL, 
                                                               participationWeights = NULL, 
                                                
                                                               treatment = NULL, # med1Study
                                                               response = NULL, # sysBloodPressureStudy
                                                
                                                               family = gaussian,
                                                
                                                               studyData = data$studyData,
                                                               aggregateTargetData = data$aggregateTargetData)))
  
  expect_true(is.transportTADA(testResult))

  expect_no_error(testSummary <- summary(testResult))
  
  expect_true(inherits(testSummary, "summary.transportTADA")) 
  expect_true(inherits(testSummary$propensitySMD, "data.frame"))
  expect_true(inherits(testSummary$prePostTable, "data.frame"))
  expect_true(inherits(testSummary$msmSummary, "summary.glm"))
  
  expect_no_error(testPlot <- plot(testResult, type = "propensityHist"))
  expect_true(ggplot2::is.ggplot(testPlot))
  
  expect_no_error(testPlot <- plot(testResult, type = "propensitySMD"))
  expect_true(ggplot2::is.ggplot(testPlot))
  
  expect_no_error(testPlot <- plot(testResult, maxWeight = 3, type = "participationHist"))
  expect_true(ggplot2::is.ggplot(testPlot))
  
  expect_no_error(testPlot <- plot(testResult, type = "msm"))
  expect_true(ggplot2::is.ggplot(testPlot))
  
  testPlot + ggplot2::scale_fill_manual(values = c("#073660", "#2FB9AB")) + ggplot2::theme(text = ggplot2::element_text(size = 15))
  
})

test_that("Scenario 2: without custom weights, glm provided for propensityScoreModel", {
  set.seed(20240730)
  data <- generateTestDataTADA()
  
  propensityScoreModel <- glm(med1 ~ sex + percentBodyFat + stress, data = data$studyData, family = binomial())
  
  expect_no_error(testResult <- suppressWarnings(transportTADA(msmFormula = sysBloodPressure ~ med1,
                                                               propensityScoreModel = propensityScoreModel,
                                                               
                                                               matchingCovariates = c("sex", "stress", "med2", "toxicGrade", "percentBodyFat"),
                                                               
                                                               propensityWeights = NULL, 
                                                               participationWeights = NULL, 
                                                               
                                                               treatment = NULL, # med1Study
                                                               response = NULL, # sysBloodPressureStudy
                                                               
                                                               family = gaussian,
                                                               
                                                               studyData = data$studyData,
                                                               aggregateTargetData = data$aggregateTargetData)))
  
  expect_true(is.transportTADA(testResult))
  
  expect_no_error(testSummary <- summary(testResult))
  
  expect_true(inherits(testSummary, "summary.transportTADA")) 
  expect_true(inherits(testSummary$propensitySMD, "data.frame"))
  expect_true(inherits(testSummary$prePostTable, "data.frame"))
  expect_true(inherits(testSummary$msmSummary, "summary.glm"))
  
  expect_no_error(testPlot <- plot(testResult, type = "propensityHist"))
  expect_true(ggplot2::is.ggplot(testPlot))
  
  expect_no_error(testPlot <- plot(testResult, type = "propensitySMD"))
  expect_true(ggplot2::is.ggplot(testPlot))
  
  expect_no_error(testPlot <- plot(testResult, maxWeight = 3, type = "participationHist"))
  expect_true(ggplot2::is.ggplot(testPlot))
  
  expect_no_error(testPlot <- plot(testResult, type = "msm"))
  expect_true(ggplot2::is.ggplot(testPlot))
  
  testPlot + ggplot2::scale_fill_manual(values = c("#073660", "#2FB9AB")) + ggplot2::theme(text = ggplot2::element_text(size = 15))
  
})

test_that("Scenario 3: customize both propensity weights and participation weights", {
  set.seed(20240730)
  data <- generateTestDataTADA()
  
  expect_no_error(testResult <- suppressWarnings(transportTADA(msmFormula = sysBloodPressure ~ med1,
                                                               propensityScoreModel = med1 ~ sex + percentBodyFat + stress,
                                                               
                                                               matchingCovariates = c("sex", "stress", "med2", "toxicGrade", "percentBodyFat"),
                                                               
                                                               propensityWeights = rep(1, nrow(data$studyData)), 
                                                               participationWeights = rep(1, nrow(data$studyData)), 
                                                               
                                                               treatment = "med1", # med1Study
                                                               response = NULL, # sysBloodPressureStudy
                                                               
                                                               family = gaussian,
                                                               
                                                               studyData = data$studyData,
                                                               aggregateTargetData = data$aggregateTargetData))) #|> expect_warning() |> expect_warning()
  
  expect_true(is.transportTADA(testResult))
  
  expect_no_error(testSummary <- summary(testResult, covariates = c("sex", "percentBodyFat", "stress"), effectModifiers = c("stress", "med2")))
  
  expect_true(inherits(testSummary, "summary.transportTADA")) 
  expect_true(inherits(testSummary$propensitySMD, "data.frame"))
  expect_true(inherits(testSummary$prePostTable, "data.frame"))
  expect_true(inherits(testSummary$msmSummary, "summary.glm"))
  
  expect_error(testPlot <- plot(testResult, type = "propensityHist")) # expect error when custom propensity weights 
  #expect_true(ggplot2::is.ggplot(testPlot))
  
  expect_no_error(testPlot <- plot(testResult, type = "propensitySMD"))
  expect_true(ggplot2::is.ggplot(testPlot))
  
  expect_no_error(testPlot <- suppressWarnings(plot(testResult, maxWeight = 3, type = "participationHist")))
  expect_true(ggplot2::is.ggplot(testPlot))
  
  expect_no_error(testPlot <- plot(testResult, type = "msm"))
  expect_true(ggplot2::is.ggplot(testPlot))
  
  testPlot + ggplot2::scale_fill_manual(values = c("#073660", "#2FB9AB")) + ggplot2::theme(text = ggplot2::element_text(size = 15))
  
})


# rlang::last_trace()
# dummizeIPD(dummizeCols): auto detecting based on is.factor and the total number of subtypes. Only need to dummize when subtypes are greater than 2


