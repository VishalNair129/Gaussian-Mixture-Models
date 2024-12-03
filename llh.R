llh <- function(params, data, covars) {
  loglik = 0
  for (i in 1:length(data)) {
    W = exp(params[1] + params[2] * covars[i, 1]) / (1 + exp(params[1] + params[2] * covars[i, 1]))
    q = exp(params[11]) / (1 + exp(params[11]))
    mu = params[3] + params[4] * covars[i, 1] + params[5] * covars[i, 2] + params[6] * covars[i, 3] + 
      params[7] * covars[i, 4] + params[8] * covars[i, 5] + params[9] * covars[i, 6]
    sigma = params[10]
    
    # Use W consistently in the log-likelihood calculation
    loglik = loglik + log(W * dnorm(data[i], mean = mu, sd = sigma) + (1 - W) * dnorm(data[i], mean = mu / q, sd = sigma / q))
  }
  return((-1) * loglik) 
}

