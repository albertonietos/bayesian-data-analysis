## ----setup, include=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------
# This chunk sets echo = TRUE as default, that is print all code.
# knitr::opts_chunk$set can be used to set other notebook generation options, too.
# include=FALSE inside curly brackets makes this block not be included in the pdf.
knitr::opts_chunk$set(echo = TRUE)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(aaltobda)
data("algae")

# test data
algae_test <- c(0, 1, 1, 0, 0, 0)

n <- length(algae)
n
y <- sum(algae)
y


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
alpha <- 2
beta <- 10
beta_point_est <- function(prior_alpha, prior_beta, data) {
  n <- length(data)
  y <- sum(data)
  
  return((prior_alpha+y)/(prior_alpha+prior_beta+n))
}
beta_point_est(prior_alpha = alpha, prior_beta = beta, data = algae)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
beta_interval <- function(prior_alpha, prior_beta, data, prob) {
  n <- length(data)
  y <- sum(data)
  
  # posterior parameters
  posterior_alpha <- prior_alpha + y
  posterior_beta <- prior_beta + n - y
  
  # lower and upper intervals
  lower <- (1-prob)/2.0
  upper <- 1-lower
  
  return (qbeta(c(lower, upper), posterior_alpha, posterior_beta))
}
beta_interval(prior_alpha = 2, prior_beta = 10, data = algae, prob = 0.9)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
beta_low <- function(prior_alpha, prior_beta, data, pi_0) {
  n <- length(data)
  y <- sum(data)
  
  # posterior parameters
  posterior_alpha <- prior_alpha + y
  posterior_beta <- prior_beta + n - y
  
  return(pbeta(pi_0, posterior_alpha, posterior_beta))
}
beta_low(prior_alpha = 2, prior_beta = 10, data = algae, pi_0 = 0.2)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot_posterior <- function(prior_alpha, prior_beta, data) {
  n <- length(data)
  y <- sum(data)
  
  # posterior parameters
  posterior_alpha <- prior_alpha + y
  posterior_beta <- prior_beta + n - y
  
  x <- seq(0, 1, length = 10000)
  density <- dbeta(x, prior_alpha, prior_beta)
  plot(x, density, type = "l", col = "blue", xlab = "p(algae)", ylab = "Density", ylim = c(0, 25), main = paste("Beta(", prior_alpha, ",", prior_beta, ")"))
  
  density <- dbeta(x, posterior_alpha, posterior_beta)
  lines(x, density, type = "l", col = "green")
  legend("topright", legend = c("Prior", "Posterior"),
         col = c("blue", "green"), lty = 1, cex = 0.8)
}
plot_posterior(prior_alpha = 1, prior_beta = 1, data = algae)
plot_posterior(prior_alpha = 2, prior_beta = 10, data = algae)
plot_posterior(prior_alpha = 2, prior_beta = 30, data = algae)
plot_posterior(prior_alpha = 20, prior_beta = 100, data = algae)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
n <- length(algae)
y <- sum(algae)
        
prior_alpha = 1
prior_beta = 1       
x <- seq(0, 0.3, length = 10000)
density <- dbeta(x, prior_alpha + y, prior_beta + n - y)
plot(x, density, type = "l", col = "blue", xlab = "p(algae)", ylab = "Posterior density", ylim = c(0, 25))

prior_alpha = 2
prior_beta = 10
density <- dbeta(x, prior_alpha + y, prior_beta + n - y)
lines(x, density, type = "l", col = "green")

prior_alpha = 2 
prior_beta = 30
density <- dbeta(x, prior_alpha + y, prior_beta + n - y)
lines(x, density, type = "l", col = "red")

prior_alpha = 20
prior_beta = 100
density <- dbeta(x, prior_alpha + y, prior_beta + n - y)
lines(x, density, type = "l", col = "purple")

legend("topright", legend = c("Beta(1,1)", "Beta(2,10)", "Beta(2,30)", "Beta(20,100)"),
       col = c("blue", "green", "red", "purple"), lty = 1, cex = 0.8)

