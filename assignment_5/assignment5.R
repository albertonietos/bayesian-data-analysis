## ----setup, include=FALSE-------------------------------------------------------------------------------------------------------------------------------------
# This chunk just sets echo = TRUE as default (i.e. print all code)
knitr::opts_chunk$set(echo = TRUE)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
library(aaltobda)
library(rstan)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
data("bioassay")

density_ratio <-
  function(alpha_propose,
           alpha_previous,
           beta_propose,
           beta_previous,
           x,
           y,
           n)
  {
    # multivariate normal prior distribution parameters
    mu_prior <- c(0, 10)
    sigma_prior <-
      matrix(data = c(4, 10, 10, 100),
             nrow = 2,
             ncol = 2)
    
    # Compute log-posterior of proposal distribution
    log_likelihood_proposal <-
      bioassaylp(alpha_propose, beta_propose, x, y, n)
    log_prior_proposal <- dmvnorm(c(alpha_propose, beta_propose),
                                  mu_prior, sigma_prior, log = TRUE)
    log_posterior_proposal <-
      log_likelihood_proposal + log_prior_proposal
    
    # Compute log-posterior of previous distribution
    log_likelihood_previous <-
      bioassaylp(alpha_previous, beta_previous, x, y, n)
    log_prior_previous <- dmvnorm(c(alpha_previous, beta_previous),
                                  mu_prior,
                                  sigma_prior,
                                  log = TRUE)
    log_posterior_previous <-
      log_likelihood_previous + log_prior_previous
    
    # Compute the ratio of densities
    return (exp(log_posterior_proposal - log_posterior_previous))
    
  }


density_ratio(
  alpha_propose = 1.89,
  alpha_previous = 0.374,
  beta_propose = 24.76,
  beta_previous = 20.04,
  x = bioassay$x,
  y = bioassay$y,
  n = bioassay$n
)
density_ratio(
  alpha_propose = 0.374,
  alpha_previous = 1.89,
  beta_propose = 20.04,
  beta_previous = 24.76,
  x = bioassay$x,
  y = bioassay$y,
  n = bioassay$n
)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
Metropolis_bioassay <- function(sigma_alpha, sigma_beta, max_iters) {
  # First, we sample a (normal) proposal distribution as a starting point
  mu <- c(0, 10)
  sigma <- matrix(data = c(4, 10, 10, 100),
                  nrow = 2,
                  ncol = 2)
  theta_0 <- rmvnorm(n = 1, mu, sigma)
  
  # Set the starting point for alpha and beta
  alpha <- numeric(length = max_iters)
  beta <- numeric(length = max_iters)
  alpha[1] <- theta_0[1, 1]
  beta[1] <- theta_0[1, 2]
  
  # Run n times
  for (i in 2:max_iters) {
    alpha_previous <- alpha[[i - 1]]
    beta_previous <- beta[[i - 1]]
    
    # Sample proposal parameters from the proposal distribution
    alpha_propose <-
      rnorm(1, mean = alpha_previous, sd = sigma_alpha)
    beta_propose <- rnorm(1, mean = beta_previous, sd = sigma_beta)
    
    # Compute the density ratio
    r <- density_ratio(
      alpha_propose,
      alpha_previous,
      beta_propose,
      beta_previous,
      bioassay$x,
      bioassay$y,
      bioassay$n
    )
    
    # Set the parameters for this iteration
    if (runif(1) <= min(r, 1)) {
      # Update parameters with proposed parameters
      alpha[i] <- alpha_propose
      beta[i] <- beta_propose
    } else {
      # Update parameters with previous parameters
      alpha[i] <- alpha_previous
      beta[i] <- beta_previous
    }
  }
  
  # return list of alpha and beta for all iterations
  list(alpha = alpha, beta = beta)
}

# Metropolis_bioassay(sigma_alpha = 1, sigma_beta = 3, max_iters = 3000)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
max_iters <- 3000
N_chains <- 10
warm_up_length <- floor(max_iters / 2)

alpha_chains <-
  matrix(nrow = N_chains, ncol = max_iters - warm_up_length)
beta_chains <-
  matrix(nrow = N_chains, ncol = max_iters - warm_up_length)
for (i in 1:N_chains) {
  chains <- Metropolis_bioassay(sigma_alpha = 1,
                                sigma_beta = 5,
                                max_iters = max_iters)
  alpha_chains[i,] <- chains$alpha[(warm_up_length + 1):max_iters]
  beta_chains[i,] <- chains$beta[(warm_up_length + 1):max_iters]
}


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(
  alpha_chains[1,],
  type = "l",
  col = 1,
  xlab = "Chain length",
  ylab = "alpha"
)
for (i in 2:N_chains) {
  lines(
    alpha_chains[i,],
    type = "l",
    col = i,
    xlab = "Chain length",
    ylab = "alpha"
  )
}
title("Convergence of chains for alpha")


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(
  beta_chains[1,],
  type = "l",
  col = 1,
  xlab = "Chain length",
  ylab = "beta"
)
for (i in 2:N_chains) {
  lines(
    beta_chains[i,],
    type = "l",
    col = i,
    xlab = "Chain length",
    ylab = "beta"
  )
}
title("Convergence of chains for beta")


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
Rhat(alpha_chains)
Rhat(beta_chains)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
 
plot(alpha_chains[1,], beta_chains[1,], col = 1, xlab = "alpha", ylab = "beta")
for (i in 2:N_chains) {
  points(alpha_chains[i,], beta_chains[i,], col = i, xlab = "alpha", ylab = "beta")
}
title("Draws of alpha and beta")

