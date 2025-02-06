test_that("Scenario 1: aggregate study data, IPD target data", {
  set.seed(20240813)
  testDataInterpolated <- generateTestDataInterpolated()
  
  expect_no_warning(testResult <- with(testDataInterpolated,
                                   transportInterpolated(link = "identity",
                                                        effectModifiers = effectModifiers,
                                                        mainTreatmentEffect = mainTreatmentEffect,
                                                        mainSE = mainSE,
                                                        subgroupTreatmentEffects = subgroupTreatmentEffects,
                                                        subgroupSEs = subgroupSEs,
                                                        studySampleSize = nStudy,
                                                        aggregateStudyData = aggregateStudyData,
                                                        targetData = targetData)))
  
  expect_true(is.transportInterpolated(testResult))
  
  expect_no_warning(testSummary <- summary(testResult))
  
  expect_true(is.data.frame(testSummary$subgroupEffects))
  expect_false(is.data.frame(testSummary$aggregateStudyData))
  expect_false(is.data.frame(testSummary$aggregateTargetData))
  
  expect_no_warning(print(testSummary))
  
  expect_no_warning(testPlot <- plot(testResult))
  expect_true(ggplot2::is.ggplot(testPlot))
})

test_that("Scenario 2: aggregate study data, aggregate target data", {
  set.seed(20240813)
  testDataInterpolated <- generateTestDataInterpolated()
  
  expect_no_warning(testResult <- with(testDataInterpolated,
                                       transportInterpolated(link = "identity",
                                                             effectModifiers = effectModifiers,
                                                             mainTreatmentEffect = mainTreatmentEffect,
                                                             mainSE = mainSE,
                                                             subgroupTreatmentEffects = subgroupTreatmentEffects,
                                                             subgroupSEs = subgroupSEs,
                                                             corrStructure = diag(length(effectModifiers)),
                                                             studySampleSize = nStudy,
                                                             aggregateStudyData = aggregateStudyData,
                                                             targetData = aggregateTargetData)))
  
  expect_warning(testResult2 <- with(testDataInterpolated,
                                       transportInterpolated(link = "identity",
                                                             effectModifiers = effectModifiers,
                                                             mainTreatmentEffect = mainTreatmentEffect,
                                                             mainSE = mainSE,
                                                             subgroupTreatmentEffects = subgroupTreatmentEffects,
                                                             subgroupSEs = subgroupSEs,
                                                             studySampleSize = nStudy,
                                                             aggregateStudyData = aggregateStudyData,
                                                             targetData = aggregateTargetData)))
  
  expect_true(is.transportInterpolated(testResult))
  
  expect_no_warning(testSummary <- summary(testResult))
  
  expect_true(is.data.frame(testSummary$subgroupEffects))
  expect_false(is.data.frame(testSummary$aggregateStudyData))
  expect_false(is.data.frame(testSummary$aggregateTargetData))
  
  expect_no_warning(print(testSummary))
  
  expect_no_warning(testPlot <- plot(testResult))
  expect_true(ggplot2::is.ggplot(testPlot))
})
