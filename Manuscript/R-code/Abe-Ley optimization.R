func.regII <- function(par, data){
  
  eta <- par[1:2]
  nu <- par[3:4]
  kappa <- par[5]
  alpha <- par[6]
  lambda <- par[7]
  
  y <- data$y
  theta <- data$theta
  x <- as.matrix(data[, c("ax", "bx")])
  z <- as.matrix(data[, c("az", "bz")])
  n <- length(data$y)
  
  p1 <- n*(log(alpha) - log(2*pi*cosh(kappa)))
  p2 <- alpha * sum(log(exp(x%*%nu)))
  p3 <- sum(log(1 + lambda*sin(theta - eta[1] + 2*atan(z[,2]*eta[2]))))
  p4 <- (alpha - 1) * sum(log(y))
  p5 <- sum((exp(x%*%nu)*y)^alpha*(1 - tanh(kappa)*cos(theta - eta[1] +  2*atan(z[,2]*eta[2]))))
  
  
  p1+p2+p3+p4-p5
  
  
}

func.reg.cond.lin <- function(par, data){
  
  
  eta <- par[1:2]
  nu <- par[3:4]
  kappa <- par[5]
  alpha <- par[6]
  lambda <- par[7]
  
  y <- data$y
  theta <- data$theta
  x <- as.matrix(data[, c("ax", "bx")])
  z <- as.matrix(data[, c("az", "bz")])
  n <- length(data$y)
  
  shape <- exp(x%*%nu)*(1-tanh(kappa)*cos(theta - (eta[1] +  2*atan(z[,2]*eta[2]))))^(1/alpha)
  
  log(alpha) + sum(log(shape^alpha)) + sum(log(y^(alpha - 1))) - sum((shape*y)^(alpha))

}

func.reg.cond.circ <- function(par, data){
  
  
  eta <- par[1:2]
  nu <- par[3:4]
  kappa <- par[5]
  alpha <- par[6]
  lambda <- par[7]
  
  y <- data$y
  theta <- data$theta
  x <- as.matrix(data[, c("ax", "bx")])
  z <- as.matrix(data[, c("az", "bz")])
  n <- length(data$y)
  
  conc <- y^{alpha}*exp(x%*%nu)^{alpha}*tanh(kappa)
  mu <-  eta[1] +  2*atan(z[,2]*eta[2])
  
  log(1) - sum(log(2*pi*besselI(conc, 0))) + sum(log(1 + lambda * sin(theta - mu))) + 
  sum((exp(x%*%nu)*y)^(alpha)*tanh(kappa)*cos(theta-mu))
  
}