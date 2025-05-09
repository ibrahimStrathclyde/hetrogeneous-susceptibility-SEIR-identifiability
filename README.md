# SEIR Model Parameter Identifiability Analysis

This repository contains the R code and analysis used in the contribution to "Some statistical aspects of the Covid-19 Pandemic" (Wood et al., 2025) by M. Gabriela M. Gomes, Ibrahim Mohammed, and Chris Robertson.

## Overview

The main focus of this repository is to demonstrate that:

1. Parameter identifiability issues in SEIR models are not specific to models with heterogeneity parameters.
2. Using data from multiple concurrent epidemics can improve parameter identifiability in both homogeneous and heterogeneous models.

We use maximum likelihood estimation (MLE) to fit various SEIR models to simulated epidemiological data, assessing parameter correlations and the accuracy of parameter recovery.

## Contents

This repository is organized as follows:
seir-parameter-identifiability/
├── README.md                             # Repository documentation
├── LICENSE                               # MIT License
├── .gitignore                            # Git ignore rules
├── install_packages.R                    # Script to install required packages
├── R/                                    # Core functions
│   ├── MaxLik_fit_functions.R            # Functions for model fitting
│   ├── utility_functions.R               # Utility functions
│   └── visualization_functions.R         # Plotting functions
├── scripts/                              # Analysis scripts
│   ├── 1_single_epidemic_homogeneous.R   # Homogeneous model, single epidemic
│   ├── 2_single_epidemic_heterogeneous.R # Heterogeneous model, single epidemic
│   ├── 3_dual_epidemic_homogeneous.R     # Homogeneous model, dual epidemics
│   └── 4_dual_epidemic_heterogeneous.R   # Heterogeneous model, dual epidemics
├── results/                              # Output data (will be generated)
│   └── .gitkeep
├── figures/                              # Output figures (will be generated)
│   └── .gitkeep
└── paper/                                # Publication files
└── RSS_Discussion_Wood.pdf           # Discussion paper


The scripts folder contains the main analysis files that should be run in numerical order. The R folder contains supporting functions that are used by these scripts. The results and figures folders will be populated when the scripts are executed.

## Models

### Homogeneous SEIR Model
The standard SEIR model without individual variation, where NPIs (non-pharmaceutical interventions) are modeled through a time-varying transmission parameter.

### Heterogeneous SEIR Model
An SEIR model that incorporates individual variation in susceptibility, modeled using a gamma distribution. This model better captures the decrease in susceptibility of the remaining population as the epidemic progresses.

In our implementation, individual variation in susceptibility is parameterized by:
- **v**: Coefficient of variation of the susceptibility distribution
- **λ**: Defined as λ = 1 + v², which appears in our model equations and determines the rate at which herd immunity is approached

The model dynamics are given by:

dS/dt = -c(t)R₀γIS^λ
dE/dt = c(t)R₀γIS^λ - δE
dI/dt = δE - γI
where S, E, and I are susceptible, exposed, and infectious population fractions, δ is the progression rate from exposed to infectious, γ is the recovery rate, R₀ is the basic reproduction number, and c(t) represents time-varying intervention strength.

## Simulation Parameters

In our analysis, we simulate epidemics with specific initial conditions to generate synthetic data for model fitting. The key parameters include:

### Initial Conditions Relationship

For all simulations, we maintain the relationship between initial exposed (E₀) and infectious (I₀) individuals as:I₀ = E₀/2.5
This ratio was derived mathematically by solving the equations of E and I early in the initial growth of the epidemic.

### Multiple Epidemic Scenarios

For the dual epidemic analyses, we simulate pairs of epidemics with different initial conditions to represent outbreaks peaking at different times. The specific cases used in our paper figures are:

- **Case 1**: Base scenario with large difference in initial conditions
  - E₀(epidemic 1) = 1040, I₀(epidemic 1) = 416
  - E₀(epidemic 2) = 60, I₀(epidemic 2) = 24

- **Case 2**: Medium difference in initial conditions
  - E₀(epidemic 1) = 520, I₀(epidemic 1) = 208
  - E₀(epidemic 2) = 60, I₀(epidemic 2) = 24

- **Case 3**: Small difference in initial conditions
  - E₀(epidemic 1) = 130, I₀(epidemic 1) = 52
  - E₀(epidemic 2) = 60, I₀(epidemic 2) = 24

In all cases, the second epidemic has fixed initial conditions, while the first epidemic's size varies. This allows us to investigate how the ratio of initial prevalence between epidemics affects parameter identifiability.

## Requirements

To run this code, you'll need:

- R version 4.0.0 or higher
- The following R packages:
  - tidyverse
  - deSolve
  - GGally
  - gridExtra
  - grid

You can install these packages with:

```r
install.packages(c("tidyverse", "deSolve", "GGally", "gridExtra", "grid"))


Usage

1. Clone this repository:
git clone https://github.com/yourusername/seir-parameter-identifiability.git

2. Run the analysis scripts in numerical order:

Rscript scripts/1_single_epidemic_homogeneous.R
Rscript scripts/2_single_epidemic_heterogeneous.R
Rscript scripts/3_dual_epidemic_homogeneous.R

3. Alternatively, open the scripts in RStudio and run them interactively.

Findings

1. Parameter correlations are a feature of both homogeneous and heterogeneous SEIR models, not specifically introduced by heterogeneity parameters.
2. Fitting models to concurrent epidemics with different initial conditions significantly reduces parameter correlations, particularly between:

Heterogeneity parameter (v) and intervention strength (c_value2)
Basic reproduction number (R0) and intervention timing (t0)


3. This approach enables more reliable parameter estimation, addressing identifiability concerns that have limited the use of heterogeneous models in epidemiological policy research.


Citation
If you use this code or findings in your research, please cite:
Gomes, M. G. M., Mohammed, I., & Robertson, C. (2025). M. Gabriela M. Gomes, Ibrahim Mohammed and Chris Robertson's contribution to the Discussion of 'Some statistical aspects of the Covid-19 Pandemic' by Wood et al. Journal of the Royal Statistical Society Series A: Statistics in Society, 1-3.


License
This project is licensed under the MIT License - see the LICENSE file for details.
Contributors

M. Gabriela M. Gomes (University of Strathclyde, NOVA School of Science and Technology)
Ibrahim Mohammed (University of Strathclyde, Abubakar Tafawa Balewa University, Bauchi, Nigeria)
Chris Robertson (University of Strathclyde, Public Health Scotland)