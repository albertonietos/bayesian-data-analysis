## ----setup, include=FALSE-----------------------------------------------------------------------------------------
# This chunk just sets echo = TRUE as default (i.e. print all code)
knitr::opts_chunk$set(echo = TRUE)


## -----------------------------------------------------------------------------------------------------------------
library(aaltobda)
library(ggplot2)
library(grid)
library(gridExtra)


## -----------------------------------------------------------------------------------------------------------------
S <- 4000  # number of draws
data("bioassay_posterior")
tail(bioassay_posterior)


## -----------------------------------------------------------------------------------------------------------------
# Means
mu_alpha_post <- mean(bioassay_posterior$alpha)
mu_beta_post <- mean(bioassay_posterior$beta)
# MCSEs
mcse_mu_alpha_post <- sqrt(var(bioassay_posterior$alpha) / S)
mcse_mu_beta_post <- sqrt(var(bioassay_posterior$beta) / S)


## -----------------------------------------------------------------------------------------------------------------
# Quantiles
quantile_alpha_post <-
  quantile(x = bioassay_posterior$alpha, probs = c(0.05, 0.95))
quantile_beta_post <-
  quantile(x = bioassay_posterior$beta, probs = c(0.05, 0.95))
# MCSEs of quantiles of alpha
mcse_quantile_alpha_post_lower <-
  mcse_quantile(draws = bioassay_posterior$alpha, prob = 0.05)
mcse_quantile_alpha_post_upper <-
  mcse_quantile(draws = bioassay_posterior$alpha, prob = 0.95)
# MCSEs of quantiles of beta
mcse_quantile_beta_post_lower <-
  mcse_quantile(draws = bioassay_posterior$beta, prob = 0.05)
mcse_quantile_beta_post_upper <-
  mcse_quantile(draws = bioassay_posterior$beta, prob = 0.95)


## -----------------------------------------------------------------------------------------------------------------
data("bioassay")
log_importance_weights <- function(alpha, beta) {
  return (bioassaylp(alpha, beta, bioassay$x, bioassay$y, bioassay$n))
}

alpha <- c(1.896, -3.6, 0.374, 0.964, -3.123, -1.581)
beta <- c(24.76, 20.04, 6.15, 18.65, 8.16, 17.4)
round(log_importance_weights(alpha, beta), 2)


## -----------------------------------------------------------------------------------------------------------------
normalized_importance_weights <- function(alpha, beta) {
  # Compute the log ratios
  log_imp_ratios <- log_importance_weights(alpha, beta)
  # Exponentiate the log ratios
  imp_ratios <- exp(log_imp_ratios)
  # Normalize the ratios (sum of all is 1)
  imp_ratios_normlz <- imp_ratios / sum(imp_ratios)
}
round(normalized_importance_weights(alpha = alpha, beta = beta), 3)


## -----------------------------------------------------------------------------------------------------------------
# Bivariate normal distribution prior
mu_alpha_prior <- 0
sd_alpha_prior <- 2
mu_beta_prior <- 10
sd_beta_prior <- 10
rho <- 0.6

mu <- c(mu_alpha_prior, mu_beta_prior)
sigma <-
  array(c(
    c(sd_alpha_prior ^ 2, rho * sd_alpha_prior * sd_beta_prior),
    c(rho * sd_alpha_prior * sd_beta_prior, sd_beta_prior ^ 2)
  ),
  dim = c(2, 2))

# Number of draws
S <- 4000
# Random draws from bivariate normal distribution
sample <- data.frame(rmvnorm(S, mu, sigma))
colnames(sample) <- c("alpha", "beta")

# Compute the normalized importance ratios
ratios <- normalized_importance_weights(sample$alpha, sample$beta)

# Plot the histogram
hist(ratios, xlab = "Normalized importance weights", main = "Histogram of normalized weights")


## -----------------------------------------------------------------------------------------------------------------
S_eff <- function(alpha, beta) {
  # Compute the normalized importance ratios
  ratios <- normalized_importance_weights(alpha, beta)
  return (1 / sum(ratios ^ 2))
}
round(S_eff(alpha = sample$alpha, beta = sample$beta), 3)


## -----------------------------------------------------------------------------------------------------------------
posterior_mean <- function(alpha, beta) {
  ratios <- normalized_importance_weights(alpha, beta)
  posterior_mean_alpha <- sum(alpha * ratios)
  posterior_mean_beta <- sum(beta * ratios)
  return (c(posterior_mean_alpha, posterior_mean_beta))
}
round(posterior_mean(alpha = alpha, beta = beta),3)

# Number of draws
S <- 4000
# Random draws from bivariate normal distribution
sample <- data.frame(rmvnorm(S, mu, sigma))
colnames(sample) <- c("alpha", "beta")

# Compute the posterior mean
round(posterior_mean(alpha = sample$alpha, beta = sample$beta),3)

# Compute the MCSEs of the posterior mean
posterior_mean_mcse <- function(alpha, beta) {
  S_effect <- S_eff(alpha = alpha, beta = beta)
  ratios <- normalized_importance_weights(alpha = alpha, beta = beta)
  
  # Expectation of the square of theta -> E[theta^2]
  E_sqrd_theta <- c(sum(alpha^2 * ratios), sum(beta^2 * ratios))
  # Expectation of theta squared -> E[theta]^2
  E_theta_sqrd <- posterior_mean(alpha = alpha, beta = beta)^2
  # variance = E[theta^2] - E[theta]^2
  variance <- E_sqrd_theta - E_theta_sqrd
  
  mcse_mu_post <- sqrt(variance / S_effect)
  return(mcse_mu_post)
}

posterior_mean_mcse(alpha = sample$alpha, beta = sample$beta)

