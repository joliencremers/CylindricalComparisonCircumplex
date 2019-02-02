# Posterior Sampling CL-PN
#
# This script contains functions to sample from the posterior of the CL-PN model for cylindrical data
#





##### define helper functions #####


## Functions to sample from posterior distributions ##

# Sample Gamma and Sigma, the regression coefficients and residual variance of 
# the regression on the linear outcome y.

# input

# N = sample size
# X = design matrix
# y = linear outcome
# sig = residual variance (starting value / previous iteration)
# p.gamma.mu = prior mean vector for the regression coefficients
# p.gamma.prec = prior precision for the regression coefficients
# p.sig.shape = prior shape parameter for the residual variance
# p.sig.scale = prior scale parameter for the residual variance

post.gamma.sigma <- function(N, X, y, sigma, p.gamma.mu, p.gamma.prec, p.sigma.shape, p.sigma.scale){
  
  prec.post <- t(X)%*%X + p.gamma.prec
  var.post <- solve(prec.post)
  mu.post <- solve(prec.post) %*% ((p.gamma.prec%*%p.gamma.mu) + (t(X)%*%y))
  
  gamma <- t(mvrnorm(1, mu.post, sigma*var.post))
  
  shape.post <- N/2 + p.sigma.shape
  scale.post <- 0.5*t(y-X%*%t(gamma))%*%(y-X%*%t(gamma)) + p.sigma.scale
  
  sigma <- 1/rgamma(1, shape.post , scale.post)

  
  return(list("gamma" = gamma, "sigma" = sigma) )
  
}

# Sample Beta, the regression coefficients of one component (I, II) of the 
# regression on the circular outcome theta.

# input

# Z = design matrix
# Y = (cos(theta), sin(theta))*r (augmented circular outcome)
# p.beta.mu = prior mean for the regression coefficients
# p.beta.prec = prior precision for the regression coefficients


post.beta <- function(Z, Y, p.beta.mu, p.beta.prec){
  
  mu.post <- solve(t(Z)%*%Z + p.beta.prec) %*% (p.beta.prec%*%p.beta.mu + t(Z)%*%Y)
  var.post <- solve(t(Z)%*%Z + p.beta.prec)
  
  beta <- mvrnorm(1, mu.post, var.post)
  
}

# Sampling r, the latent lengths

# input

# theta = circular outcome
# beta.I = regression coefficients for the first component (starting value / previous iteration)
# beta.II = regression coefficients for the second component (starting value / previous iteration)
# ZI = design matrix for the prediction of the cosine component of the circular outcome
# ZII = design matrix for the prediction of the sine component of the circular outcome
# N = sample size
# r = r value (starting value / previous iteration)

post.r <- function(theta, beta.I, beta.II, ZI, ZII, N, r){
  
  mub1 <- ZI%*%beta.I
  mub2 <- ZII%*%beta.II
  
  Dbd <- cos(theta)*mub1 + sin(theta)*mub2
  
  for(i in 1:N){
    
    m <- runif(1, 0, 1)*exp(-0.5*(r[i]-Dbd[i])^2)
    u <- runif(1, 0, 1)
    r1 <- Dbd[i] + max(-Dbd[i], -sqrt(-2*log(m)))
    r2 <- Dbd[i] + sqrt(-2*log(m))
    r[i] <- sqrt((r2^2 - r1^2)*u + r1^2)
    
  }
  
  return(r)
  
}


## Functions to compute linear and circular likelihoods ##


# loglikelihood linear outcome
# input
# see descriptions above

llik.lin <- function(y, X, sigma, gamma){
  
  log(1) - log(sqrt(2*pi*sigma^2)) + sum((y-X%*%t(gamma))^2/(2*sigma^2))
  
}

# loglikelihood circular outcome
# input
# see descriptions above

llik.circ <- function(theta, ZI, ZII, beta.I, beta.II, r){
  
  mub1 <- ZI%*%beta.I
  mub2 <- ZII%*%beta.II
  
  mu <- cbind(mub1, mub2)
  ll <- length(theta) * (log(1) - log(2*pi))
  
  for(i in 1: length(theta)){
    
    Dbd <- cos(theta[i])*mub1[i] + sin(theta[i])*mub2[i]
    
    ll <- ll + (-t(mu[i,])%*%mu[i,])/2 + 0.5*(-r[i]^2 + 2*r[i]*Dbd)
    
  }
 
  return(ll)
  
}





##### Define the wrapper function for sampling from the CL-GPN model #####

# input

# theta = circular outcome
# y = linear outcome
# X = design matrix for the linear outcome
# ZI = design matrix for the first component of the circular outcome
# ZII = design matrix for the second component of the circular outcome
# theta.hold = circular outcome in holdout set
# y.hold = linear outcome in holdout set
# X.hold = design matrix for the linear outcome in holdout set
# ZI = design matrix for the first component of the circular outcome in holdout set
# ZII = design matrix for the second component of the circular outcome in holdout set
# its = amount of iterations

# output

# Gamma = posterior samples for Gamma (regression parameters linear outcome)
# BI = posterior samples for BI (regression parameters first component circular outcome)
# BII = posterior samples for BI (regression parameters second component circular outcome)
# Sigma = posterior samples Sigma (error variance linear outcome) 
# ll.circ = circular loglikelihood
# ll.lin = linear loglikelihood
# theta_pred = posterior predictive values for theta (using actual covariate values)
# y_pred = posterior predictive values for y (using actual covariate values)
# theta_pred.hold = posterior predictive values for theta in the holdout set (using actual covariate values)
# y_pred.hold = posterior predictive values for y in the holdout set (using actual covariate values)

CLPN <- function(theta, y, X, ZI, ZII, its, theta.hold, y.hold, X.hold, ZI.hold, ZII.hold){
  
  #Determine the sample size
  N <- nrow(X)
  N.hold <- length(theta.hold)
  
  #Define prior parameters
  p.gamma.mu <- rep(0, 4)
  p.gamma.prec <- diag(4)*0.0001
  
  p.beta.mu <- rep(0, 2)
  p.beta.prec <- diag(2)*0.0001
  
  p.sigma.scale <- 0.001
  p.sigma.shape <- 0.001
  
  #Set starting values
  Beta.I <- c(0,0)
  Beta.II <- c(0,0)
  r <- rep(1, N)
  r.hold <- rep(1, N.hold)
  #r <- y
  sigma <- c(1)
  gamma <- c(0,0,0,0)
  
  #Compute the "latent" outcome (circular)
  Y <- cbind(cos(theta)*r, sin(theta)*r)
  
  #Create empty results vectors
  Gamma <- matrix(NA, its, 4)
  BI <- matrix(NA, its, 2)
  BII <- matrix(NA, its, 2)
  Sigma <- rep(NA, its)
  theta_pred <- matrix(NA, N, its)
  y_pred <- matrix(NA, N, its)
  theta_pred.hold <- matrix(NA, N.hold, its)
  y_pred.hold <- matrix(NA, N.hold, its)
  
  ll.lin <- rep(NA, its)
  ll.circ <- rep(NA, its)

  #Start the MCMC chain
  for(i in 1:its){
    
    #Sample posterior parameters
    gs <- post.gamma.sigma(N, X, y, sigma, p.gamma.mu, p.gamma.prec, p.sigma.shape, p.sigma.scale)
    
    gamma <- gs$gamma
    sigma <- gs$sigma
    
    Beta.I <- post.beta(ZI, Y[,1], p.beta.mu, p.beta.prec)
    Beta.II <- post.beta(ZII, Y[,2], p.beta.mu, p.beta.prec)
    
    #Sample r
    r <- post.r(theta, Beta.I, Beta.II, ZI, ZII, N, r)
    r.hold <- post.r(theta.hold, Beta.I, Beta.II, ZI.hold, ZII.hold, N.hold, r.hold)

    #compute predicted values
    for(j in 1:N){
      Y_pred <- mvrnorm(1, c(ZI[j,]%*%Beta.I, ZII[j,]%*%Beta.II), diag(2))
      theta_pred[j,i] <- atan2(Y_pred[2], Y_pred[1])
      y_pred[j,i] <- rnorm(1, X[j,]%*%t(gamma), sqrt(sigma))
      
    }

    #compute predicted values holdout set
    for(j in 1:N.hold){
      Y_pred <- mvrnorm(1, c(ZI.hold[j,]%*%Beta.I, ZII.hold[j,]%*%Beta.II), diag(2))
      theta_pred.hold[j,i] <- atan2(Y_pred[2], Y_pred[1])
      y_pred.hold[j,i] <- rnorm(1, X.hold[j,]%*%t(gamma), sqrt(sigma))
      
    }
    
    
    #Compute the "latent" outcome (circular)
    Y <- cbind(cos(theta), sin(theta))*r
    
    #Update the design matrix for the linear outcome
    X <- cbind(X[,1], Y, X[,4])
    X.hold <- cbind(X.hold[,1], cos(theta.hold)*r.hold, sin(theta.hold)*r.hold, X.hold[,4])
    
    #Save results in matrices
    Gamma[i,] <- gamma
    BI[i,] <- Beta.I
    BII[i,] <- Beta.II
    Sigma[i] <- sqrt(sigma)
    
    ll.lin[i] <- llik.lin(y.hold, X.hold, sqrt(sigma), gamma)
    ll.circ[i] <- llik.circ(theta.hold, ZI.hold, ZII.hold, Beta.I, Beta.II, r.hold)
    
    

    
  }
  
  return(list("Gamma"= Gamma, "BI" = BI, "BII" = BII, "Sigma" = Sigma,
              "ll.circ" = ll.circ, "ll.lin" = ll.lin,
              "theta_pred" = theta_pred, "y_pred" = y_pred,
              "theta_pred.hold" = theta_pred.hold, "y_pred.hold" = y_pred.hold))
  
}

