kappa_to_mu <- function (kappa, sigma){
  (sigma /sqrt(2)) * (1/kappa - kappa)
}

mu_to_kappa <- function(mu, sigma){
  (sqrt(2 * sigma^2 + mu^2) - mu) / (sqrt(2) * sigma)
}

xi_eta <- function(dat, theta, m, Sigma){
  Sigma_inv <- solve(Sigma)
  qy <- stats::mahalanobis(dat, center = theta, cov = Sigma_inv, inverted = TRUE)
  d <- NCOL(dat)
  gamma <- c(sqrt(2 + t(m) %*% Sigma_inv %*% m))
  delta <- c(sqrt(qy))
  lambda <- 1 - d/2
  e <- log(delta / gamma)
  f <- log(besselK(x = gamma * delta, nu = lambda - 1))
  h <- log(besselK(x = gamma * delta, nu = lambda))
  g <- log(besselK(x = gamma * delta, nu = lambda + 1))
  xi <- exp(-e + f - h)
  eta <- exp(e + g - h)
  return(list(xi = xi, eta = eta))
}

Sigma_1 <- function(dat, theta, m){
  d <- NCOL(dat)
  n <- NROW(dat)
  Sigma <- matrix(rep(0, d^2), nrow = d)
  for (i in 1:n){
    Sigma <- Sigma + (1 / n) * (dat[i,] - theta - m) %*% t(dat[i,] - theta - m)
  }
  Sigma
}


Sigma_update <- function(dat, theta, m, xi, eta){
  d <- NCOL(dat)
  n <- NROW(dat)
  Sigma <- matrix(rep(0, d^2), nrow = d)
  for (i in 1:n){
    Sigma <- Sigma + (1 / n) * xi[i] * (dat[i,] - theta) %*%
      t(dat[i,] - theta) - (1/n) * eta[i] * m %*% t(m)
  }
  Sigma
}


theta_update <- function(dat, xi, eta){
  n <- NROW(dat)
  num <- apply(xi * dat, 2, sum) * sum(eta) - n * apply(dat, 2, sum)
  denom <- sum(xi) * sum(eta) - n^2
  num / denom
}

m_update <- function(dat, theta, xi, eta){
  n <- NROW(dat)
  (apply(dat, 2, sum) - n * theta) / sum(eta)
}


AML_EM <- function(dat, tol = 1e-5){
  dat <- as.matrix(dat)
  n <- NROW(dat)
  d <- NCOL(dat)
  N <- 1000
  theta <- matrix(rep(0, d * N), nrow = N)
  m <- matrix(rep(0, d * N), nrow = N)
  for (i in 1:N) {
    if (i == 1) {
      theta[1, ] <- apply(dat, 2, median)
      m[1, ] <- rep(0, d)
      Sigma <- Sigma_1(dat, theta[1, ], m[1, ])
    } else {
      xi_cur <- xi_eta(dat, theta[i-1, ], m[i-1, ], Sigma)$xi
      eta_cur <- xi_eta(dat, theta[i-1, ], m[i-1, ], Sigma)$eta
      theta[i, ] <- theta_update(dat, xi_cur, eta_cur)
      m[i, ] <- m_update(dat, theta[i, ], xi_cur, eta_cur)
      if (is.na(m[i, 1]) || is.na(theta[i, 1])){
        break
        i = i - 1
      }
      test1 <- t(theta[i, ] - theta[i-1, ]) %*% (theta[i, ] - theta[i-1, ])
      test2 <- t(m[i, ] - m[i-1,]) %*% (m[i, ] - m[i-1, ])
      if (test2 < tol)
        break
      Sigma <- Sigma_update(dat, theta[i, ], m[i, ], xi_cur, eta_cur)
    }
  }
  return(
    list(theta = ifelse(is.na(theta[i, ]), theta[i-1, ], theta[i, ]),
         Sigma = Sigma,
         m = ifelse(is.na(m[i, ]), m[i-1, ], m[i, ]))
  )
}
