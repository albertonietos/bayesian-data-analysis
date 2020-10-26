## ----setup, include=FALSE--------------------------------------------------------------------------------------------------------------------------
# This chunk just sets echo = TRUE as default (i.e. print all code)
knitr::opts_chunk$set(echo = TRUE)


## ---- warning=FALSE--------------------------------------------------------------------------------------------------------------------------------
library(aaltobda)
library(rstan)
library(bayesplot)

# rstan options
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Import the data
data("bioassay")


## --------------------------------------------------------------------------------------------------------------------------------------------------
code_bioassay <- "bioassay_model.stan"
writeLines(readLines(code_bioassay))


## --------------------------------------------------------------------------------------------------------------------------------------------------
bioassay_data <- list(N = length(bioassay$x),
                      x = bioassay$x,
                      n = bioassay$n,
                      y = bioassay$y,
                      mu = c(0, 10),
                      Sigma = matrix(data = c(4, 10, 10, 100), nrow = 2, ncol = 2))


## --------------------------------------------------------------------------------------------------------------------------------------------------
max_iters <- 2000
n_chains <- 5
fit <- stan(file = 'bioassay_model.stan', data = bioassay_data, chains = n_chains, 
            iter = max_iters, warmup = floor(max_iters/2))


## --------------------------------------------------------------------------------------------------------------------------------------------------
monitor(fit)


## ---- fig.width=5, fig.height=5--------------------------------------------------------------------------------------------------------------------
p <- mcmc_scatter(fit, pars = c('theta[1]', 'theta[2]'), alpha = 0.2)
(p + labs(title = "Scatter plot of alpha and beta draws",
         x = expression(alpha),
         y = expression(beta)) + stat_density_2d(color = "black", size = .5) 
  + xlim(-4, 10) + ylim(-10, 40))

