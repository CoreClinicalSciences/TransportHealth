# TransportHealth

**TransportHealth** is an R package, which was developed by [Core Clinical Sciences](https://www.coreclinicalsciences.com/), that provides functions to conduct transportability and generalizability analyses. 

Transportability and generalizability analyses are types of causal inference methods that allow us to quantitatively assess external validity of evidence from randomized clinical trials (RCTs) and other diverse research studies. We assess how applicable the findings from original studies are to settings that the original studies were not conducted in or to the broader population that was not included in the original studies. In transportability analyses, we aim to transport the findings from study sample to an external target population by adjusting for the different distribution of effect modifiers between the study and the target sample. Generalizability analyses are similar but concerns the case where the study sample is a subset of the target population. 

**The scope of the current version of our package includes the following analyses.**

- Inverse probability (IP) weighting for mergeable individual patient-level datasets (IPDs) of original and target studies

- G-computation for unmergable IPDs of original and target studies

**For the future scope, we are currently developing the following methods.**

- Target Aggregate Data Adjustment (TADA) method that can transport findings from the original IPD study to aggregate (summary-level) data of a target study

- Methods that can transport aggregate data from the original study to IPD from a target study

- Methods that can transport aggregate data from the original study to aggregate-level data from a target study data


## Getting Started

If you are just getting started with **TransportHealth**, we recommend starting with the tutorial vignettes. These are examples provided as part of our package documentation. Future papers, as they are published as preprints or after being peer reviewed, will be added here, as they become available. 

## Installation

Install the latest development version from [GitHub](CoreClinicalSciences/TransportHealth)

``` r
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("CoreClinicalSciences/TransportHealth")
```

## Citing TransportHealth

To cite `TransportHealth`, please see [here](https://coreclinicalsciences.github.io/TransportHealth/authors.html#citation)

