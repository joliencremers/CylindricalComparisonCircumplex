# Abe-Ley optimization
#
# This script contains helper functions for optimization of the Abe-Ley model for cylindrical data
#


## define helper functions ##

# input function for optimization

# input
# par = vector of starting values for the parameters
# data = dataframe containing the circular (theta) and linear outcome  (y)

# output = loglikelihood

func.regII <- function(par, data){
  
  beta <- par[1:2]
  gamma <- par[3:4]
  kappa <- par[5]
  alpha <- par[6]
  lambda <- par[7]
  
  y <- data$y
  theta <- data$theta
  x <- as.matrix(data[, c("ax", "bx")])
  z <- as.matrix(data[, c("az", "bz")])
  n <- length(data$y)
  
  p1 <- n*(log(alpha) - log(2*pi*cosh(kappa)))
  p2 <- alpha * sum(log(exp(x%*%gamma)))
  p3 <- sum(log(1 + lambda*sin(theta - beta[1] + 2*atan(z[,2]*beta[2]))))
  p4 <- (alpha - 1) * sum(log(y))
  p5 <- sum((exp(x%*%gamma)*y)^alpha*(1 - tanh(kappa)*cos(theta - beta[1] +  2*atan(z[,2]*beta[2]))))
  
  
  p1+p2+p3+p4-p5
  
  
}

## function to compute the linear loglikelihood ##

# input
# par = vector of parameter values
# data = dataframe containing the circular (theta) and linear outcome  (y)

func.reg.cond.lin <- function(par, data){
  
  
  beta <- par[1:2]
  gamma <- par[3:4]
  kappa <- par[5]
  alpha <- par[6]
  lambda <- par[7]
  
  y <- data$y
  theta <- data$theta
  x <- as.matrix(data[, c("ax", "bx")])
  z <- as.matrix(data[, c("az", "bz")])
  n <- length(data$y)
  
  shape <- exp(x%*%gamma)*(1-tanh(kappa)*cos(theta - (beta[1] +  2*atan(z[,2]*beta[2]))))^(1/alpha)
  
  log(alpha) + sum(log(shape^alpha)) + sum(log(y^(alpha - 1))) - sum((shape*y)^(alpha))

}

## function to compute the circular loglikelihood ##

# input
# par = vector of parameter values
# data = dataframe containing the circular (theta) and linear outcome  (y)

func.reg.cond.circ <- function(par, data){
  
  
  beta <- par[1:2]
  gamma <- par[3:4]
  kappa <- par[5]
  alpha <- par[6]
  lambda <- par[7]
  
  y <- data$y
  theta <- data$theta
  x <- as.matrix(data[, c("ax", "bx")])
  z <- as.matrix(data[, c("az", "bz")])
  n <- length(data$y)
  
  conc <- y^{alpha}*exp(x%*%gamma)^{alpha}*tanh(kappa)
  mu <-  beta[1] +  2*atan(z[,2]*beta[2])
  
  log(1) - sum(log(2*pi*besselI(conc, 0))) + sum(log(1 + lambda * sin(theta - mu))) + 
  sum((exp(x%*%gamma)*y)^(alpha)*tanh(kappa)*cos(theta-mu))
  
}