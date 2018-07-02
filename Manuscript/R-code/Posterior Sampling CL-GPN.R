#Functions for posterior distributions

#Sampling Gamma and Sigma, the regression coefficients and residual variance of 
#the regression on the linear outcome y.

#N = sample size
#X = design matrix
#y = linear outcome
#sig = residual variance (starting value or value from previous iteration)
#p.gamma.mu = prior mean vector for the regression coefficients
#p.gamma.prec = prior precision for the regression coefficients
#p.sig.shape = prior shape parameter for the residual variance
#p.sig.scale = prior scale parameter for the residual variance

post.gamma.sig <- function(N, X, y, sig, p.gamma.mu, p.gamma.prec, p.sig.shape, p.sig.scale){
  
  prec.post <- t(X)%*%X + p.gamma.prec
  var.post <- solve(prec.post)
  mu.post <- solve(prec.post) %*% ((p.gamma.prec%*%p.gamma.mu) + (t(X)%*%y))
  
  gamma <- t(mvrnorm(1, mu.post, sig*var.post))
  
  shape.post <- N/2 + p.sig.shape
  scale.post <- 0.5*t(y-X%*%t(gamma))%*%(y-X%*%t(gamma)) + p.sig.scale
  
  sig <- 1/rgamma(1, shape.post , scale.post)
  
  return(list("gamma" = gamma, "sig" = sig) )
  
}

#Sampling rho, element of sigma, the variance-covariance matrix of the circular 
#outcome

post.rho <- function(Y, Z, Beta.I, Beta.II, tau, p.rho.mu, p.rho.var){
  
  num <- tau^(-1) * t(Y[,1] - Z%*%Beta.I)%*%(Y[,2] - Z%*%Beta.II) + p.rho.mu*p.rho.var^(-1)
  den <- tau^(-1) * t(Y[,2] - Z%*%Beta.II)%*%(Y[,2] - Z%*%Beta.II) + p.rho.var^(-1)
  
  mu.rho.post <- num/den
  var.rho.post <- 1/(tau^(-1)*(t(Y[,2] - Z%*%Beta.II)%*%(Y[,2] - Z%*%Beta.II)) + p.rho.var^(-1))
  
  rho <- rnorm(1, mu.rho.post, var.rho.post)
  
}

#Sampling tau, element of sigma, the variance-covariance matrix of the circular 
#outcome

post.tau <- function(N, Y, Z, Beta.I, Beta.II, rho, p.tau.shape, p.tau.scale){
  
tau.shape.post <- N/2 + p.tau.shape
comp <- Y[,1]-((Z%*%Beta.I) + rho*(Y[,2] - Z%*%Beta.II))
tau.scale.lik <- 0.5*sum(t(comp)%*%comp)
tau.scale.post <- tau.scale.lik + p.tau.scale
  
tau <- 1/rgamma(1, tau.shape.post, tau.scale.post)
  
}

#Sampling Beta, the regression coefficients of one component (I, II) of the 
#regression on the circular outcome theta.

#p = the amount of predictors in Z
#Z = design matrix
#Y = (cos(theta), sin(theta))*r
#beta = the regression coefficients (previous iteration/staring values)
#p.beta.var = prior variance for the regression coefficients, mean is set to 0
#sigma = the variance-covariance matrix (previous iteration/starting values)


post.beta <- function(p, Z, Y, beta, p.beta.var, sigma){
  
  Zj <- Z[,p]

  betak <- as.matrix(beta[,-p])
  
  var.post <- solve(sum(Zj^2)*solve(sigma) + p.beta.var^(-1)*diag(2))
  
  Zbetak <- t(rbind(betak[1,]%*%Z[,-p], betak[2,]%*%Z[,-p]))
  
  mu.post <- var.post%*%solve(sigma)%*%colSums(-Zj*Zbetak + Zj*Y)

  beta <- mvrnorm(1, mu.post, var.post)
  
}

#Sampling r, the latent lengths

#theta = circular outcome
#beta.I = regression coefficients for the first component
#beta.II = regression coefficients for the second component
#ZI = design matrix for the first component
#ZII = design matrix for the second component
#N = sample size
#r = r value (starting value or from the previous iteration)
#sigma = the variance-covariance matrix of the circular outcome

post.r <- function(theta, beta.I, beta.II, Z, N, r, sigma){
  
  mub1 <- Z%*%beta.I
  mub2 <- Z%*%beta.II
  
  mu <- cbind(mub1, mub2)
  
  u <- cbind(cos(theta), sin(theta))

  for(i in 1:N){
    
    A <- t(u[i,])%*%solve(sigma)%*%(u[i,])
    
    B <- t(u[i,])%*%solve(sigma)%*%(mu[i,])
    
    m <- runif(1, 0, exp(-0.5*A*(r[i]-B/A)^2))
    unif <- runif(1, 0, 1)
    r1 <- B/A + max(-B/A, -sqrt((-2*log(m))/A))
    r2 <- B/A + sqrt((-2*log(m))/A)
    r[i] <- sqrt((r2^2 - r1^2)*unif + r1^2)
    
  }
  
  return(r)
  
}

llik.lin <- function(y, X, sig, gamma){
  
  log(1) - log(sqrt(2*pi*sig^2)) + sum((y-X%*%t(gamma))^2/(2*sig^2))
  
}

llik.circ <- function(theta, Z, beta.I, beta.II, r, tau, Sigma){

  mub1 <- Z%*%beta.I
  mub2 <- Z%*%beta.II
  
  mu <- cbind(mub1, mub2)
  u <- cbind(cos(theta), sin(theta))
  
  ll <- 0
  
  for(i in 1: length(theta)){
    
    Dbd <- cbind(r[i]*cos(theta[i]) - mub1[i], r[i]*sin(theta[i]) - mub2[i])
    
    ll <- ll + log(r[i]) + log(1) - log(2*pi+sqrt(tau)) - (Dbd%*%solve(Sigma)%*%t(Dbd))/(2*tau)

  }
  
  return(ll)
  
}

#MCMC sampler

#theta = circular outcome
#y = linear outcome
#X = design matrix for the linear outcome
#Z = design matrix for the circular outcome
#its = amount of iterations
#p = amount of predictors in Z

CLGPN <- function(theta, y, X, Z, its, p, theta.hold, y.hold, X.hold, Z.hold){
  
  #Determine the sample size
  N <- nrow(X)
  N.hold <- length(theta.hold)
  
  #Define prior parameters
  p.gamma.mu <- rep(0, 4)
  p.gamma.prec <- diag(4)*0.0001
  
  p.beta.var <- 10^5
  
  p.sig.scale <- 0.001
  p.sig.shape <- 0.001
  
  p.tau.scale <- 0.01
  p.tau.shape <- 0.01
  
  p.rho.mu <- 0
  p.rho.var <- 10000
  
  #Set starting values
  Beta.I <- c(0,0)
  Beta.II <- c(0,0)
  r <- rep(1, N)
  r.hold <- rep(1, N.hold)
  #r <- y
  sig <- c(1)
  gamma <- c(0,0,0,0)
  tau <- 1
  rho <- 0
  sigma <- matrix(c(tau + rho^2, rho, rho, 1), 2, 2)
  
  #Compute the "latent" outcome (circular)
  Y <- cbind(cos(theta)*r, sin(theta)*r)
  
  #Create empty results vectors
  Gamma <- matrix(NA, its, 4)
  BI <- matrix(NA, its, 2)
  BII <- matrix(NA, its, 2)
  Sig <- rep(NA, its)
  Sigma <- array(NA, dim = c(2,2,its))
  theta_pred <- matrix(NA, N, its)
  y_pred <- matrix(NA, N, its)
  theta_pred.hold <- matrix(NA, N.hold, its)
  y_pred.hold <- matrix(NA, N.hold, its)
  theta_pred_min <- rep(NA, its)
  theta_pred_max <- rep(NA, its)
  theta_pred_mean <- rep(NA, its)

  ll.lin <- rep(NA, its)
  ll.circ <- rep(NA, its)
  
  #Start the MCMC chain
  for(i in 1:its){

    #Sample posterior parameters
    gs <- post.gamma.sig(N, X, y, sig, p.gamma.mu, p.gamma.prec, p.sig.shape, p.sig.scale)
    
    gamma <- gs$gamma
    sig <- gs$sig
    
    beta <- rbind(Beta.I, Beta.II)
    
    for(j in 1:p){
      
      beta[,j] <- post.beta(j, Z, Y, beta, p.beta.var, sigma)
      
    }
    
    Beta.I <- beta[1,]
    Beta.II <- beta[2,]
    
    rho <-   post.rho(Y, Z, Beta.I, Beta.II, tau, p.rho.mu, p.rho.var) 
    tau <- post.tau(N, Y, Z, Beta.I, Beta.II, rho, p.tau.shape, p.tau.scale)
    
    sigma <- matrix(c(tau + rho^2, rho, rho, 1), 2, 2)
    
    #Sample r (only in case of non-circumplex data)  
    r <- post.r(theta, Beta.I, Beta.II, Z, N, r, sigma)
    r.hold <- post.r(theta.hold, Beta.I, Beta.II, Z.hold, N.hold, r.hold, sigma)

    
    #compute predicted values regression
    Y_pred <- mvrnorm(1, c(c(1, min(Z[,2]))%*%Beta.I,  c(1, min(Z[,2]))%*%Beta.II), sigma)
    theta_pred_min[i] <- atan2(Y_pred[2], Y_pred[1])

    Y_pred <- mvrnorm(1, c(c(1, max(Z[,2]))%*%Beta.I,  c(1, max(Z[,2]))%*%Beta.II), sigma)
    theta_pred_max[i] <- atan2(Y_pred[2], Y_pred[1])

    Y_pred <- mvrnorm(1, c(c(1, mean(Z[,2]))%*%Beta.I,  c(1, mean(Z[,2]))%*%Beta.II), sigma)
    theta_pred_mean[i] <- atan2(Y_pred[2], Y_pred[1])

    #compute predicted values
    for(j in 1:N){
      Y_pred <- mvrnorm(1, c(Z[j,]%*%Beta.I, Z[j,]%*%Beta.II), sigma)
      theta_pred[j,i] <- atan2(Y_pred[2], Y_pred[1])
      y_pred[j,i] <- rnorm(1, X[j,]%*%t(gamma), sqrt(sig))
      
    }
    
    #compute predicted values holdout data
    for(j in 1:N.hold){
      Y_pred <- mvrnorm(1, c(Z.hold[j,]%*%Beta.I, Z.hold[j,]%*%Beta.II), sigma)
      theta_pred.hold[j,i] <- atan2(Y_pred[2], Y_pred[1])
      y_pred.hold[j,i] <- rnorm(1, X.hold[j,]%*%t(gamma), sqrt(sig))
      
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
    Sig[i] <- sqrt(sig)
    Sigma[,,i] <- sigma
    
    ll.lin[i] <- llik.lin(y.hold, X.hold, sqrt(sig), gamma)
    #ll.circ[i] <- llik.circ(theta.hold, Z.hold, Beta.I, Beta.II, y.hold, tau, sigma)
    ll.circ[i] <- llik.circ(theta.hold, Z.hold, Beta.I, Beta.II, r.hold, tau, sigma)

  }
  
  return(list("Gamma"= Gamma, "BI" = BI, "BII" = BII, "Sig" = Sig, "Sigma" = Sigma,
              "ll.circ" = ll.circ, "ll.lin" = ll.lin,
              "theta_pred" = theta_pred, "y_pred" = y_pred,
              "theta_pred.hold" = theta_pred.hold, "y_pred.hold" = y_pred.hold,
              "theta_pred.min" = theta_pred_min,  "theta_pred.max" = theta_pred_max,
              "theta_pred.mean" = theta_pred_mean))
  
}

