#simulate from wrapped cauchy distribution

rwcauchy <- function(n, mu, rho){
  
  if (rho == 0) 
    result <- runif(n, 0, 2 * pi)
  else if (rho == 1) 
    result <- rep(mu, n)
  else {
    scale <- -log(rho)
    result <- rcauchy(n, mu, scale)%%(2 * pi)
  }
  return(result)
  
}

#simulate from abe-ley distribution
rweiSSVM <- function(n, beta = 1, alpha = 1, mu = pi, kappa = 1, lambda = 0){
  
  theta <- c()
  x <- c()
  
  #simulate theta
  for(i in 1:n){
    
    u <- runif(1,0,1)
    
    theta.1 <- rwcauchy(1, mu, tanh(kappa/2))
    
    if(u < (1 + lambda*sin(theta.1 - mu))/2){
      theta[i] <- theta.1 
    }else{
      theta[i] <- -theta.1
    }
    
    shape <- beta * (1-tanh(kappa)*cos(theta-mu))^(1/alpha)
    #compute scale parameter for different parametrization,
    #shape = (scale_parWeibull in R)^(-alpha) = 1/(scale_parWeibull in R)^alpha dus 
    #scale_parWeibull in R = kth root(1/shape) dus
    #scale_parWeibull in R = (1/b)^(1/alpha)
    scale <- (1/shape)^(1/alpha) 
    
    x[i] <- rweibull(1, shape = alpha, scale = scale)
    
  }
  
  out <- cbind(theta %% (2*pi), x)
  
  colnames(out) <- c("theta", "x")
  
  return(as.data.frame(out))
  
}
