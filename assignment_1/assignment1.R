## ----setup, include=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------
# This chunk just sets echo = TRUE as default (i.e. print all code)
knitr::opts_chunk$set(echo = TRUE)


## ---- fig.height=4, fig.width=4-------------------------------------------------------------------------------------------------------------------------------------------
mu <- 0.2
var <- 0.01
alpha <- mu * (((mu*(1-mu))/var) - 1)
beta <- (alpha*(1-mu))/mu
x <- seq(0, 1, length=100)
y <- dbeta(x, alpha, beta)
plot(x, y)


## ---- fig.height=4, fig.width=4-------------------------------------------------------------------------------------------------------------------------------------------
n <- 1000
sample <- rbeta(n, alpha, beta)
hist(sample)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mu_bar <- mean(sample)
var_bar <- var(sample)
print(paste("Sample mean:", format(mu_bar)))
print(paste("Sample variance:", format(var_bar)))


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
lower <- quantile(sample, probs = 0.025)
upper <- quantile(sample, probs = 0.975)
print(paste("The central 95% probability interval of the distribution from the drawn samples"))
print(paste("is [", format(lower), ",", format(upper), "]"))


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
boxes <- matrix(c(2,4,1,5,1,3), ncol = 2,
                dimnames = list(c("A", "B", "C"), c("red", "white")))
boxes
p_red <- function(boxes) {
  p_A <- 0.4
  p_B <- 0.1
  p_C <- 0.5
  p <- p_A*boxes[1,1]/sum(boxes[1,]) + p_B*boxes[2,1]/sum(boxes[2,]) + p_C*boxes[3,1]/sum(boxes[3,])
  return(p)
}
p_red(boxes = boxes)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
p_box <- function(boxes) {
  p_A <- 0.4
  p_B <- 0.1
  p_C <- 0.5
  p_A.red <- p_A*boxes[1,1]/sum(boxes[1,])/p_red(boxes)
  p_B.red <- p_B*boxes[2,1]/sum(boxes[2,])/p_red(boxes)
  p_C.red <- p_C*boxes[3,1]/sum(boxes[3,])/p_red(boxes)
  return(c(p_A.red, p_B.red, p_C.red))
}
p_box(boxes = boxes)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
p_identical_twin <- function(fraternal_prob, identical_prob) {
  p_it.and.tb <- identical_prob * 0.5
  p_ft.and.tb <- fraternal_prob * 0.5 * 0.5
  return(p_it.and.tb/(p_it.and.tb + p_ft.and.tb))
}
p_identical_twin(fraternal_prob = 1/150, identical_prob = 1/400)

