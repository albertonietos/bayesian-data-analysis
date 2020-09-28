## ----setup, include=FALSE-----------------------------------------------------
# This chunk just sets echo = TRUE as default (i.e. print all code)
knitr::opts_chunk$set(echo = TRUE)


## -----------------------------------------------------------------------------
library(aaltobda)
data("windshieldy1")
head(windshieldy1)


## -----------------------------------------------------------------------------
windshieldy_test <- c(13.357, 14.928, 14.896, 14.820)


## -----------------------------------------------------------------------------
n <- length(windshieldy1)
y <- sum(windshieldy1)
s2 <- var(windshieldy1)
my <- mean(windshieldy1)


## -----------------------------------------------------------------------------
mu_point_est <- function(data) {
  n <- length(data)
  df <- n - 1
  
  mean <- mean(data)
  s2 <- var(data)
  scale <- s2 / n

  point_est <- mean(rtnew(100000, df, mean = mean, scale = scale))
  return(point_est)
}
mu_point_est(data = windshieldy1)


## -----------------------------------------------------------------------------
mu_interval <- function(data, prob) {
  n <- length(data)
  df <- n - 1

  mean <- mean(data)
  s2 <- var(data)
  scale <- s2 / n

  # lower and upper intervals
  lower <- (1-prob)/2.0
  upper <- 1-lower

  sample<-rt(100000, df)*sqrt(scale)+mean

  return (quantile(sample, c(lower, upper)))
  }
  mu_interval(data = windshieldy1, prob = 0.95)


## ---- fig.height=3, fig.width=4-----------------------------------------------
plot_density <- function(data){
  n <- length(data)
  df <- n - 1

  mean <- mean(data)
  s2 <- var(data)
  scale <- s2 / n
  
  y <- seq(from = 13, to = 16, length.out = 1000)
  density <- dtnew(y, df, mean = mean, scale = scale)
  plot(y, density, type = "l")
  title("Probability density function of y")
}
plot_density(windshieldy1)


## -----------------------------------------------------------------------------
mu_pred_point_est <- function(data) {
  n <- length(data)
  df <- n - 1
  
  mean <- mean(data)
  s2 <- var(data)
  scale <- sqrt(1+(1/n))*sqrt(s2)

  point_est <- mean(rtnew(100000, df, mean = mean, scale = scale))
  return(point_est)
}
mu_pred_point_est(data = windshieldy1)


## -----------------------------------------------------------------------------
mu_pred_interval <- function(data, prob) {
  n <- length(data)
  df <- n - 1

  mean <- mean(data)
  s2 <- var(data)
  scale <- sqrt(1+(1/n))*sqrt(s2)

  # lower and upper intervals
  lower <- (1-prob)/2.0
  upper <- 1-lower

  sample <- rtnew(100000, df, mean = mean, scale = scale)

  return (quantile(sample, c(lower, upper)))
  }
  mu_pred_interval(data = windshieldy1, prob = 0.95)


## ---- fig.height=3, fig.width=4-----------------------------------------------
plot_pred_density <- function(data){
  n <- length(data)
  df <- n - 1

  mean <- mean(data)
  s2 <- var(data)
  scale <- sqrt(1+(1/n))*sqrt(s2)
  
  y <- seq(from = 8, to = 22, length.out = 1000)
  density <- dtnew(y, df, mean = mean, scale = scale)
  plot(y, density, type = "l")
  title("Probability density function of y")
}
plot_pred_density(windshieldy1)


## -----------------------------------------------------------------------------
p0 <- rbeta(100000, 40, 636)
p1 <- rbeta(100000, 23, 659)

posterior_odds_ratio_point_est <- function (p0, p1) {
  num <- p1 / (1 - p1)
  den <- p0 / (1 - p0)
  odds_ratio <- num / den
  
  return (mean (odds_ratio))
}

posterior_odds_ratio_point_est(p0 = p0, p1 = p1)


## -----------------------------------------------------------------------------
posterior_odds_ratio_interval <- function (p0, p1, prob) {
  num <- p1 / (1 - p1)
  den <- p0 / (1 - p0)
  odds_ratio <- num / den
  
  # lower and upper intervals
  lower <- (1-prob)/2.0
  upper <- 1-lower

  return (quantile(odds_ratio, c(lower, upper)))
}
posterior_odds_ratio_interval(p0 = p0, p1 = p1, prob = 0.95)


## ---- fig.height=3, fig.width=4-----------------------------------------------
num <- p1 / (1 - p1)
den <- p0 / (1 - p0)
odds_ratio <- num / den
hist(odds_ratio)


## ---- fig.height=3, fig.width=4-----------------------------------------------
n_0 <- 674
y_0 <- 39

n_1 <- 680
y_1 <- 22

prior_alpha = 1
prior_beta = 1       
x <- seq(0, 0.3, length = 10000)
density_0 <- dbeta(x, prior_alpha + y_0, prior_beta + n_0 - y_0)
density_1 <- dbeta(x, prior_alpha + y_1, prior_beta + n_1 - y_1)
plot(x, density_0, type = "l", col = "blue", xlab = "probability", ylab = "Posterior density", 
     ylim = c(0, 60), xlim = c(0, 0.12))


lines(x, density_1, type = "l", col = "green")
legend("topright", legend = c("Control", "Treatment"),
       col = c("blue", "green"), lty = 1, cex = 0.8)


## ---- fig.height=5, fig.width=7-----------------------------------------------
prior_alpha = 1
prior_beta = 1       
x <- seq(0, 0.3, length = 10000)
density_0 <- dbeta(x, prior_alpha + y_0, prior_beta + n_0 - y_0)
density_1 <- dbeta(x, prior_alpha + y_1, prior_beta + n_1 - y_1)
plot(x, density_0, type = "l", col = "blue", xlab = "p", ylab = "Posterior density", 
     ylim = c(0, 60), xlim = c(0, 0.1))
lines(x, density_1, type = "l", lty = 2, col = "blue")


prior_alpha = 2
prior_beta = 10
density_0 <- dbeta(x, prior_alpha + y_0, prior_beta + n_0 - y_0)
density_1 <- dbeta(x, prior_alpha + y_1, prior_beta + n_1 - y_1)
lines(x, density_0, type = "l", lty = 2, col = "green")
lines(x, density_1, type = "l", lty = 2, col = "green")



legend("topright", legend = c("p0 w/ Beta(1,1)", "p1 w/ Beta(1,1)", "p0 w/ Beta(2,10)", 
                              "p1 w/ Beta(2,10)"),
       col = c("blue", "blue","green", "green"), lty = c(1,2, 1, 2), cex = 0.8)
title("Posterior distributions for p0 and p1 for different choice of priors")


## -----------------------------------------------------------------------------
mu <- function (data) {
  n <- length(data)
  df <- n - 1
  
  mean <- mean(data)
  s2 <- var(data)
  scale <- s2 / n

  return (rtnew(100000, df, mean = mean, scale = scale))
}
mu_1 <- mu(windshieldy1)
mu_2 <- mu(windshieldy2)

mu_d_est <- mean(mu_1-mu_2)
mu_d_est


## -----------------------------------------------------------------------------
prob <- 0.95

# lower and upper intervals
lower <- (1-prob)/2.0
upper <- 1-lower

(quantile(mu_1-mu_2, c(lower, upper)))



## ---- fig.height=3, fig.width=4-----------------------------------------------
hist(mu_1-mu_2)

