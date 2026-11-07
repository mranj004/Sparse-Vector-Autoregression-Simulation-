\# Sparse Vector Autoregression Simulation



\## Project Goal

This GitHub repository aims to see how different innovation covariance structures affect the estimation accuracy of sparse Vector Autoregressive (VAR) models.



In this specific repository, we generate data using a known, stable, sparse autoregressive coefficient matrix A in which follows a VAR(1) model



$$

Y\_t = AY\_{t-1) + \\epsilon\_t \\quad \\epsilon\_t \\sim N(0, \\Sigma\_{\\epsilon}

$$



Inspired by Basu and Michailidis (2015), this repository uses **Block-I, Block-II, and Toeplitz** for the innovation covariance structures and correlation parameters $\\rho = \\{0.5, 0.7, 0.9\\}$



For each combination of covariance structure and correlation strength, observations are generated. These observations are then used to estimate the coefficient matrix by being stored in $\\hat{A}$ with row-wise LASSO regression.



The main goal is to see how well $\\hat{A} \\approx A$ under different innovation covariance structures and correlation strengths, while holding the true coefficient matrix fixed.



\## Methodology



Each row of $A$ is estimated using the LASSO estimator



$$
\\hat{a\_i} = \\arg\\min\_{a \\in \\mathbb{R}^d} \\left{\\frac{1}{2T}\\sum\_{t=2}^T(y\_{ti} - a^{\\top{y\_{t-1}}}^2 + \\lambda||a||\_1 \\right} 

$$



Just like in the paper, the relative Frobenius estimation error is used, which is defined by 

$$
\\frac{||\\hat{A} - A||\_F}{||A||\_F}

$$



MSE, precision, recall, and F1 score is also calculated to assess coefficient support recovery.







\## Repository Structure



| File | Description |

|---|---|

| `00\_requirements.R` | Loads the required R packages |

| `01\_generate\_A.R` | Generates a stable, sparse coefficient matrix \\(A\\) |

| `02\_generate\_errors.R` | Generates Gaussian innovations with Block-I, Block-II, or Toeplitz covariance |

| `03\_generate\_VAR.R` | Generates VAR(1) observations using \\(A\\) and the innovations |

| `04\_estimate\_A.R` | Estimates \\(\\hat A\\) using row-wise LASSO |

| `05\_evaluate\_A.R` | Evaluates coefficient and support recovery |

| `06\_run\_simulations.R` | Runs repeated experiments across covariance structures and correlation strengths |



\## Reproducing the Simulations



Kindly install the required R packages and run:



```r

source("00\_requirements.R")

source("06\_run\_simulations.R")



sim <- runSim(

&#x20; T = 300,

&#x20; d = 10,

&#x20; spar = 0.8,

&#x20; rho\_A = 0.8,

&#x20; lambda = 0.05,

&#x20; n\_rep = 30,

&#x20; rhoVal = c(0.5, 0.7, 0.9),

&#x20; seedA = 293,

&#x20; seedSim = 123

)



print(sim$summary)

```



The simulation generates results for nine experimental conditions, corresponding to three covariance structures and three correlation strengths.



Results are saved in the `results/` directory.



\## Reference

Basu, S., \& Michailidis, G. (2015). Regularized estimation in sparse high-dimensional time series models. \*The Annals of Statistics, 43\*(4), 1535–1567. https://doi.org/10.1214/15-AOS1315

