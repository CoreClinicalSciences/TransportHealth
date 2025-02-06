#' Example datasets for transportability analysis
#' 
#' A pair of example datasets that represent a use case for transportability analysis.
#' 
#' @format Two data frames representing data from an RCT on the effects of an antihypertensive medication on blood pressure (\code{studyData}; 1000 rows) and data from an observational study (\code{targetData}; 1500 rows). The variables are:
#' \describe{
#'  \item{sysBloodPressure}{Systolic blood pressure of patient}
#'  \item{med1}{0-1 indicator for uptake of antihypertensive medication vs. placebo}
#'  \item{sex}{Sex of patient}
#'  \item{percentBodyFat}{Percentage body fat of patient}
#'  \item{stress}{0-1 indicator for whether patient is generally highly stressed}
#'  \item{med2}{0-1 indicator for uptake of another medication that potentially has interaction with the antihypertensive medication}
#' }
#' 
#' @name mockData
"studyData"

#' @rdname mockData
"targetData"