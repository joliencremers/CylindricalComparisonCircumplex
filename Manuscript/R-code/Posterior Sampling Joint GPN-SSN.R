require(MCMCpack)
require(truncnorm)
#require(LaplacesDemon)
require(bayesm)
require(mvtnorm)

#Constraint matrix function (for 1 circular and 1 linear variable)

cmat <- function(var){diag(3)*(c(sqrt(var), sqrt(var), 1))}

#sample from matrix normal (adapted from LaplacesDemon, function is.symmetric.matrix does not work properly; floating point problems...)

rmatrixnorm <- function (M, U, V) 
{
  if (missing(M)) 
    stop("Matrix M is missing.")
  if (!is.matrix(M)) 
    M <- matrix(M)
  if (missing(U)) 
    stop("Matrix U is missing.")
  if (!is.matrix(U)) 
    U <- matrix(U)
    if (missing(V)) 
    stop("Matrix V is missing.")
  if (!is.matrix(V)) 
    V <- matrix(V)
  if (nrow(M) != nrow(U)) 
    stop("Dimensions of M and U are incorrect.")
  if (ncol(M) != ncol(V)) 
    stop("Dimensions of M and V are incorrect.")
  n <- nrow(U)
  k <- ncol(V)
  Z <- matrix(rnorm(n * k), n, k)
  X <- M + t(chol(U)) %*% Z %*% chol(V)
  return(X)
}

#Functions for posterior distributions

#Sample Beta, a matrix of regression coefficients for the circular and linear 
#variables

#Y = the circular and linear outcomes
#X = a design matrix
#Beta = Beta matrix (previous iteration/starting value)
#Kappa = Kappa matrix (previous iteration/starting value)
#Sigma  = Sigma matrix (previous iteration/starting value)
#p.Beta = 
#p.Kappa = 

post.Beta <- function(Y, X, Beta, Sigma, p.Beta, p.Kappa){

Kappa.post <-  t(X)%*%X + p.Kappa 
Beta.post <- solve(Kappa.post) %*% (t(X)%*%Y + p.Kappa%*%p.Beta) 

Beta <- rmatrixnorm(Beta.post, solve(Kappa.post), Sigma)

return(Beta)

}


post.Sigma <- function(Y, X, N, Beta, Sigma, p.Psi, p.nu, p.Beta, p.Kappa){
  
Psi.post <- p.Psi + t(Y - X%*%Beta)%*%(Y - X%*%Beta) + t(Beta - p.Beta)%*%p.Kappa%*%(Beta - p.Beta)
nu.post <- p.nu + N  

Sigma <- rwishart(nu.post, solve(Psi.post))$IW

return(Sigma)

}


post.lambda <- function(rawdata, X, D, Sigma, Beta, p.Omega, p.gamma, m, c, N){
  
  Sig.y <- Sigma[(c+1):m,(c+1):m]
  Sig.wy <- Sigma[1:c,(c+1):m]
  Sig.w <- Sigma[1:c, 1:c]
  mu.y <- X%*%Beta[,(c+1):m]
  mu.w <- X%*%Beta[,1:c]
  y <- as.matrix(rawdata[,(c+1):m])
  
  
  Sigma.ycondw <- Sig.y-t(Sig.wy)%*%solve(Sig.w)%*%Sig.wy
  mu.ycondw <- t(mu.y) + Sig.wy%*%solve(Sig.w)%*%t(rawdata[,1:c] - mu.w)
  
  
  Omega.post <- matrix(0, m-c, m-c)
  gamma.post <- matrix(0, m-c, m-c)
  
  for(t in 1:N){
    
    Omega.post <- Omega.post + 
                  diag(m-c)*D[t,]%*%solve(Sigma.ycondw)%*% diag(m-c)*D[t,] + 
                  solve(p.Omega)
    
    gamma.post <- gamma.post + 
                  diag(m-c)*D[t,]%*%solve(Sigma.ycondw)%*%(y[t,]-mu.ycondw[,t]) + 
                  solve(p.Omega)%*%p.gamma
    
  }
  
  Omega.post <- solve(Omega.post)
  gamma.post <- Omega.post %*% gamma.post
  
  lambda <- rnorm(1, gamma.post, sqrt(Omega.post))
  
  return(lambda)
  
}



post.d <- function(X, rawdata, Beta, D, lambda, Sigma, m, c, N){
  
  Sig.y <- Sigma[(c+1):m,(c+1):m]
  Sig.wy <- Sigma[1:c,(c+1):m]
  Sig.w <- Sigma[1:c, 1:c]
  mu.y <- X%*%Beta[,(c+1):m]
  mu.w <- X%*%Beta[,1:c]
  y <- as.matrix(rawdata[,(c+1):m])
  
  Sigma.ycondw <- Sig.y-t(Sig.wy)%*%solve(Sig.w)%*%Sig.wy
  mu.ycondw <- t(mu.y) + Sig.wy%*%solve(Sig.w)%*%t(rawdata[,1:c] - mu.w)
  
  var.d <- solve(t(diag(m-c)*lambda)%*%solve(Sigma.ycondw)%*%(diag(m-c)*lambda) + diag(m-c))
  
  for(t in 1:N){
    
    mean.d <- var.d%*%t(diag(m-c)*lambda)%*%solve(Sigma.ycondw)%*%(y[t,]-mu.ycondw[,t])
    D[t,] <- rtruncnorm(1, a = 0, mean = mean.d, sd = sqrt(var.d))
    
  }
  
  return(D)
  
}


post.r <- function(X, rawdata, r, Sigma, Beta, m, p, N){
  
  
  for(i in seq(1, 2*p, 2)){
    
    u <- rawdata[,i:(i+1)]
    Sig.wi <- matrix(Sigma[i:(i+1),i:(i+1)], 2, 2)
    Sig.nwi <- matrix(Sigma[-(i:(i+1)),-(i:(i+1))], m-2, m-2)
    Sig.nwiwi <- matrix(Sigma[-(i:(i+1)), i:(i+1)], m-2, 2) #same as Sig.nwiwi (symmetrix matrix)
    mu.nwi <- X%*%Beta[,-(i:(i+1))]
    mu.wi <- X%*%Beta[,i:(i+1)]
    nwi <- matrix(rawdata[,-(i:(i+1))], N, m-2)
    

      for(t in 1:N){
        
        Sigma.wicond <- Sig.wi - t(Sig.nwiwi)%*%solve(Sig.nwi)%*%Sig.nwiwi
        mu.wicond <- mu.wi[t,] + t(Sig.nwiwi)%*%solve(Sig.nwi)%*%(nwi[t,] - mu.nwi[t,])
          
        A <- t(u[t,])%*%solve(Sigma.wicond)%*%u[t,]
        B <- t(u[t,])%*%solve(Sigma.wicond)%*%mu.wicond
        
        v <- runif(1, 0, exp(-0.5*A*(r[t,i]-B/A)^2))
        unif <- runif(1, 0, 1)
        r1 <- B/A + max(-B/A, -sqrt((-2*log(v))/A))
        r2 <- B/A + sqrt((-2*log(v))/A)
        r[t,i] <- sqrt((r2^2 - r1^2)*unif + r1^2)
        
      }
    
    
  }
  
  return(r)
  
}


llik.lin <- function(X, y, rawdata, lambda, D, Sigma, Beta, N){

  ll <- 0
  

  Sig.w <- matrix(Sigma[1:2,1:2], 2, 2)
  Sig.yw <- Sigma[1:2,3]
  Sig.y <- Sigma[3,3]
  
  mu.y <- as.numeric(X%*%Beta[,3])
  mu.w <- X%*%Beta[,1:2]
  
  
  for(t in 1:N){
    
    mean <- mu.y[t] + lambda*D[t, ] + t(Sig.yw)%*%solve(Sig.w)%*%(rawdata[t, 1:2] - mu.w[t,])
    sd <- sqrt(Sig.y + t(Sig.yw)%*%solve(Sig.w)%*%Sig.yw)
    
    ll <- ll + log(dnorm(y[t], mean, sd))
  }
  
  return(ll)
  
}

llik.circ <- function(X, theta, y, rawdata, lambda, D, Sigma, Beta, N){
  
  ll <- 0
  

    Sig.w <- matrix(Sigma[1:2,1:2], 2, 2)
    Sig.yw <- Sigma[1:2,3]
    Sig.y <- Sigma[3,3]
    
    mu.y <- as.numeric(X%*%Beta[,3])
    mu.w <- X%*%Beta[,1:2]
    
    for(t in 1:N){
      
      Sigma.wicond <- Sig.w + Sig.yw%*%solve(Sig.y)%*%t(Sig.yw)
      mu.wicond <- mu.w[t,] + as.numeric(Sig.yw%*%solve(Sig.y)%*%(y[t] - mu.y[t] - lambda*D[t,]))
      
      ll <- ll + log(dmvnorm(rawdata[t,1:2], mu.wicond, Sigma.wicond))

      
    }
  
  return(ll)
  
}


#p = amount of circular variables
#q amount of linear variables

JGPNSSN <- function(theta, y, X, its, p, q, theta.hold, y.hold, X.hold){
  
  N <- length(theta)
  N.hold <- length(theta.hold)
  
  c <- 2*p
  m <- c+q
  k <- ncol(X)
  
  #priors
  p.Kappa <- diag(k)*0.0001
  p.Beta <- matrix(0, k, m)
  
  p.nu <- 1
  p.Psi <- diag(m)*0.0001
  
  p.Omega <- diag(m-c)*10000
  p.gamma  <- 0
  
  #starting values:
  D <- as.matrix(rep(1, (m-c)*N), N, m-c)
  D.hold <- as.matrix(rep(1, (m-c)*N.hold), N.hold, m-c)
  r <- matrix(rep(1, p*N), N, p)
  r.hold <- matrix(rep(1, p*N.hold), N.hold, p)
  Sigma <- diag(m)
  #Kappa <- diag(k)
  Beta <- matrix(0, k, m)
  lambda <- 0
  
  rawdata <- cbind(cos(theta)*r,sin(theta)*r,y)
  Y <- (rawdata) - cbind(matrix(0, N, 2*p), D%*%(diag(m-c)*lambda))
  
  rawdata.hold <- cbind(cos(theta.hold)*r.hold, sin(theta.hold)*r.hold, y.hold)
  Y.hold <- (rawdata.hold) - cbind(matrix(0, N.hold, 2*p), D.hold%*%(diag(m-c)*lambda))
  
  lambda.res <- matrix(NA, its, m-c)
  Beta.res <- array(NA, dim = c(k, m, its))
  Betacon.res <- array(NA, dim = c(k, m, its))
  Sigma.res <- array(NA, dim = c(m, m, its))
  Sigmacon.res <- array(NA, dim = c(m, m, its))
  theta_pred <- matrix(NA, N, its)
  y_pred <- matrix(NA, N, its)
  theta_pred.hold <- matrix(NA, N.hold, its)
  y_pred.hold <- matrix(NA, N.hold, its)
  theta_pred_min <- rep(NA, its)
  theta_pred_max <- rep(NA, its)
  theta_pred_mean <- rep(NA, its)
  
  ll.lin <- rep(NA, its)
  ll.circ <- rep(NA, its)
  
  for(it in 1:its){
    
    
    Beta <- post.Beta(Y, X, Beta, Sigma, p.Beta, p.Kappa)
    Sigma <- post.Sigma(Y, X, N, Beta, Sigma, p.Psi, p.nu, p.Beta, p.Kappa)
    lambda <- post.lambda(rawdata, X, D, Sigma, Beta, p.Omega, p.gamma, m, c, N)
    D <- post.d(X, rawdata, Beta, D, lambda, Sigma, m, c, N)
    D.hold <- post.d(X.hold, rawdata.hold, Beta, D.hold, lambda, Sigma, m, c, N.hold)
    #r <- post.r(X, rawdata, r, Sigma, Beta, m, p, N)
    r <- y #hij lijkt niet te werken omdat (lin var) y == r.... op deze manier werkt het wel
    r.hold <- y.hold
    
    #save results in matrices
    
    lambda.res[it,] <- lambda
    Beta.res[,,it] <- Beta
    Sigma.res[,,it] <- Sigma
    
    conmat <- solve(cmat(Sigma[2,2]))
    Sigmacon.res[,,it] <- conmat%*%Sigma%*%conmat 
    Ypred <- X%*%Beta
    Ypredc <- conmat%*%t(Ypred)
    Betacon.res[,,it] <- solve(t(X)%*%X)%*%t(X)%*%t(Ypredc)

    
    #compute predicted values regression
    Y_pred <- mvrnorm(1, c(1, min(X[,2]))%*%Betacon.res[,,it]+c(0,0,lambda*D[1,]), Sigmacon.res[,,it])
    theta_pred_min[it] <- atan2(Y_pred[2], Y_pred[1])
    
    Y_pred <- mvrnorm(1, c(1, mean(X[,2]))%*%Betacon.res[,,it]+c(0,0,lambda*D[1,]), Sigmacon.res[,,it])
    theta_pred_mean[it] <- atan2(Y_pred[2], Y_pred[1])
    
    Y_pred <- mvrnorm(1, c(1, max(X[,2]))%*%Betacon.res[,,it]+c(0,0,lambda*D[1,]), Sigmacon.res[,,it])
    theta_pred_max[it] <- atan2(Y_pred[2], Y_pred[1])
    
    #compute predicted values
    for(j in 1:N){
      Y_pred <- mvrnorm(1, X[j,]%*%Betacon.res[,,it]+c(0,0,lambda*D[j,]), Sigmacon.res[,,it])
      theta_pred[j,it] <- atan2(Y_pred[2], Y_pred[1])
      y_pred[j,it] <- Y_pred[3]
      
    }
    #compute predicted values for the holdout data
    for(j in 1:N.hold){
      Y_pred <- mvrnorm(1, X.hold[j,]%*%Betacon.res[,,it]+c(0,0,lambda*D.hold[j,]), Sigmacon.res[,,it])
      theta_pred.hold[j,it] <- atan2(Y_pred[2], Y_pred[1])
      y_pred.hold[j,it] <- Y_pred[3]
      
    }
    
    #recompute augmented data
    rawdata <- cbind(cos(theta)*r,sin(theta)*r,y)
    Y <- (rawdata) - cbind(matrix(0, N, 2*p), D%*%(diag(m-c)*lambda))
    
    rawdata.hold <- cbind(cos(theta.hold)*r.hold,sin(theta.hold)*r.hold,y.hold)
    Y.hold <- (rawdata.hold) - cbind(matrix(0, N.hold, 2*p), D.hold%*%(diag(m-c)*lambda))
    
    ll.circ[it] <- llik.circ(X.hold, theta.hold, y.hold, rawdata.hold, lambda, D.hold, Sigmacon.res[,,it], Betacon.res[,,it], N.hold)
    ll.lin[it] <- llik.lin(X.hold, y.hold, rawdata.hold, lambda, D.hold, Sigmacon.res[,,it], Betacon.res[,,it], N.hold)

      }
  
  return(list("lambda"= lambda.res, "Beta" = Beta.res, "Sigma" = Sigma.res,
              "Sigmacon" = Sigmacon.res, 'Betacon' = Betacon.res,
              "ll.circ" = ll.circ, "ll.lin" = ll.lin,
              "theta_pred" = theta_pred, "y_pred" = y_pred,
              "theta_pred.hold" = theta_pred.hold, "y_pred.hold" = y_pred.hold,
              "theta_pred.min" = theta_pred_min,  "theta_pred.max" = theta_pred_max,
              "theta_pred.mean" = theta_pred_mean))
  
}


