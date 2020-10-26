// The input data is a vector 'x' of length 'N'.
data {
  int<lower=0> N;     // total number of rows
  vector[N] x;        // dose level
  int<lower=0> n[N];  // number of animals
  int y[N];  // number of deaths (positive outcomes)
  
  // parameters of prior distribution
  vector[2] mu;
  matrix[2, 2] Sigma;
  
}

// The parameters accepted by the model. Our model
// accepts the 2D parameter 'theta'.
parameters {
  vector[2] theta;
}

transformed parameters {
  vector[N] rel;
  rel = theta[1] + theta[2]*x;  // dose-response relation
}

// The model to be estimated.
model {
  theta ~ multi_normal(mu, Sigma);  // prior
  y ~ binomial_logit(n, rel);       // likelihood
}

