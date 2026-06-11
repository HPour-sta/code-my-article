# Code for "Hybrid Analytical-Numerical Algorithm for Stochastic Frontier Models" article.

**Authors:** Hadis Pouresmaili, Reza Pourmousa

## Repository Structure

code-my-article/

├── code/

│ ├── 00_hybrid_algorithm.R

│ ├── 01_simulation_NK2.R

│ ├── 02_simulation_extreme.R

│ ├── 03_empirical_application.R

│ └── 04_sensitivity_seed.R

## Requirements

- R >= 4.2.0
- Packages: `maxLik`, `truncnorm`, `pracma`, `randtoolbox`

## Installation

 ```r
 install.packages(c("maxLik", "truncnorm", "pracma", "randtoolbox"))
 ```
 
 ## Run the Code
 
 ```r
 source("code/00_hybrid_algorithm.R")
 source("code/01_simulation_NK2.R")
 source("code/02_simulation_extreme.R")
 source("code/03_empirical_application.R")
 source("code/04_sensitivity_seed.R")
 ```
 
 ## Contact

 Hadis.Pouresmaili@math.uk.ac.ir
