kappa.to.mu <- function (kappa, sigma){
  (sigma /sqrt(2)) * (1/kappa - kappa)
}

mu.to.kappa <- function(mu, sigma){
  (sqrt(2*sigma^2 + mu^2) - mu) / (sqrt(2) * sigma)
}

xi.eta <- function(dat, theta, m, Sigma){
  Sigma.inv <- solve(Sigma)
  qy <- mahalanobis(dat, center=theta, cov=Sigma.inv, inverted=TRUE)
  d <- NCOL(dat)
  gamma <- sqrt(2 + t(m) %*% Sigma.inv %*% m)
  delta <- sqrt(qy)
  lambda <- 1 - d/2
  e <- log(delta/gamma)
  f <- log(besselK(x= gamma*delta, nu= lambda-1))
  h <- log(besselK(x= gamma*delta, nu= lambda))
  g <- log(besselK(x= gamma*delta, nu= lambda+1))
  xi <- exp(-e + f - h)
  eta <- exp(e + g - h)
  return(list(xi = xi, eta = eta))
}

Sigma.update <- function(dat, theta, m, xi, eta){
  d <- NCOL(dat)
  n <- NROW(dat)
  Sigma <- matrix(rep(0, d^2), nrow= d)
  for (i in 1:n){
    Sigma <- Sigma + (1/n) * xi[i] * (dat[i,] - theta) %*%
      t(dat[i,] - theta) - (1/n) * eta[i] * m %*% t(m)
  }
  Sigma
}


## theta.update <- function(dat, m, xi, eta){
##     n <- NROW(dat)
##     (apply(xi * dat, 2, sum) - n*m) / sum(xi)
## }

theta.update <- function(dat, xi, eta){
  n <- NROW(dat)
  num <- apply(xi * dat, 2, sum) * sum(eta) - n * apply(dat, 2, sum)
  denom <- sum(xi) * sum(eta) - n^2
  num / denom
}

m.update <- function(dat, theta, xi, eta){
  n <- NROW(dat)
  (apply(dat, 2, sum) - n*theta) / sum(eta)
}


AML.EM <- function(dat, tol=1e-14){
  dat <- as.matrix(dat)
  n <- NROW(dat)
  d <- NCOL(dat)
  N <- 1000
  theta <- matrix(rep(0, d*N), nrow = N)
  m <- matrix(rep(0, d*N), nrow = N)
  for(i in 1:N){
    if (i == 1){
      theta[1,] <- apply(dat,2,median)
      m[1,] <- rep(0,d)
      Sigma <- Sigma.1(dat, theta[1,], m[1,])
    }else{
      xi.cur <- xi.eta(dat, theta[i-1,], m[i-1,], Sigma)$xi
      eta.cur <- xi.eta(dat, theta[i-1,], m[i-1,], Sigma)$eta
      theta[i,] <- theta.update(dat, xi.cur, eta.cur)
      m[i,] <- m.update(dat, theta[i,], xi.cur, eta.cur)
      if (is.na(m[i,1]) || is.na(theta[i,1]) == TRUE){
        break
        i = i-1
      }
      test1 <- t(theta[i,] - theta[i-1,]) %*% (theta[i,] - theta[i-1,])
      test2 <- t(m[i,] - m[i-1,]) %*% (m[i,] - m[i-1,])
      if ((test1  < tol) && (test2 < tol))
        break
      Sigma <- Sigma.update(dat, theta[i,], m[i,], xi.cur, eta.cur)
    }
  }
  return(list(theta = ifelse(is.na(theta[i,]) == TRUE, theta[i-1,], theta[i,]),
              Sigma = Sigma,
              m = ifelse(is.na(m[i,] == TRUE), m[i-1,], m[i,]),
              iterations = i,
              allTheta = theta[1:i,],
              allm = m[1:i,]))
}
