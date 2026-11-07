# Sparse Vector Autoregression Simulation

## Project Goal

This GitHub repository aims to see how different Gaussian innovation covariance structures affect the estimation accuracy of sparse Vector Autoregressive (VAR) models. We compare two LASSO-based estimation procedures, which is a fixed regularization parameter $\lambda=0.05$ from \texttt{glmnet} and a data-adaptive regularization parameter selected using \texttt{sparsevar}

For the simulation, we generate data using a known, stable, sparse autoregressive coefficient matrix $A$ in which follows a VAR(1) model

$$
Y\_t = AY\_{t-1) + \\epsilon\_t \\quad \\epsilon\_t \\sim N(0, \\Sigma\_{\\epsilon}
$$

Inspired by Basu and Michailidis (2015), we consider **Block-I, Block-II, and Toeplitz** for the innovation covariance structures with the correlation parameters 
$$
\rho = \\{0.5, 0.7, 0.9\\}
$$

For each combination of covariance structure and correlation strength, we generate observations from the same true coefficient matrix $A$ and these observations are then used to estimate two coefficient matrices
- $\hat{A}_{\mathrm{fixed}}$: A row-wise LASSO using $\texttt{glmnet}$ with $\lambda = 0.05$
- $\hat{A}_{\mathrm{adaptive}}$: LASSO using $\textttt{sparsevar}$ with $\lambda$ selected via time-slice validation

The main goal is to see how each estimated matrix approximates the true matrix $A$ and whether the data-adaptive estimation proceduce improves coefficient recovery under different innovation covariance structures

## Methodology
Each row of $A$ is estimated using the LASSO estimator

### 1. Data Generation
We first generate a sparse coefficient matrix $A$ with a spectral radius less than one to esnure the stability of the VAR(1) process

The matrix $A$ remains fixed throughout the experiment. But the innovation covariance matrix $\Sigma_{\epsilon}$ changes according to the selected covariance structure and correlation parameter. 

For each condition, we simulate the observations $Y$ using the VAR(1) data generating process

### 2. Coefficient Estimation
We compare two LASSO-based estimation procedures.

**Fixed-$\lambda$ LASSO ($\texttt{glmnet}$)**
$$
\\hat{a\_i} = \\arg\\min\_{a \\in \\mathbb{R}^d} \\left{\\frac{1}{2T}\\sum\_{t=2}^T(y\_{ti} - a^{\\top{y\_{t-1}}}^2 + \\lambda||a||\_1 \\right} 
$$
where $\lambda = 0.05$ is held fixed across the experimental conditions.

**Adaptive-$\lambda$ LASSO ($\texttt{sparsevar}$)**
We also estiamte the VAR(1) coefficient matrix using the $\texttt{sparsevar}$ package with an L1 penalty. Unlike the fixed estimator, this approach selects its regularization parameter from the data via time-slice validation.

Both of these procedures use the same simulated observation within each experimental condition
 
### 3. Evaluation Metrics
Just like in the paper, the relative Frobenius estimation error is used, which is defined by 

$$
\frac{|\\hat{A} - A|\_F}{|A|\_F}
$$

Smaller values imply that the estimated coefficient matrix is closer to the true matrix

MSE, precision, recall, and F1 score is also calculated to assess coefficient support recovery.

## Repository Structure
| File | Description |

|---|---|

| `00\_requirements.R` | Loads the required R packages |

| `01\_generate\_A.R` | Generates a stable, sparse coefficient matrix \\(A\\) |

| `02\_generate\_errors.R` | Generates Gaussian innovations with Block-I, Block-II, or Toeplitz covariance |

| `03\_generate\_VAR.R` | Generates VAR(1) observations using \\(A\\) and the innovations |

| `04\_estimate\_A.R` | Estimates \\(\\hat A\\) using row-wise LASSO |

| `05\_evaluate\_A.R` | Evaluates coefficient and support recovery |

| `06\_run\_simulation.R` | Runs repeated experiments across covariance structures and correlation strengths |



## Reproducing the Simulations

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
The experiment consists of nine covariance conditions (three covariance strucures and three correlation strengths), evaluated using two estimation procedures

Results are saved in the `results/` directory.

## Reference
Basu, S., \& Michailidis, G. (2015). Regularized estimation in sparse high-dimensional time series models. \*The Annals of Statistics, 43\*(4), 1535–1567. https://doi.org/10.1214/15-AOS1315

