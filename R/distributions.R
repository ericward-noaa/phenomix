# library(stats)
# library(gnorm) # https://github.com/maryclare/gnorm/blob/master/R/gnorm.R
#
# qthill
# ddnorm
# ddt
# gnorm::dgnorm()
# ddgnorm
# qgnorm
# qdnorm

#' Calculate the quantile of a Student-t distribution
#'
#' @param p The quantile of the distribution
#' @param mean The mean of the distribution
#' @param sd The std deviation of the distribution
#' @param df The df of the distribution
#' @importFrom stats qt
#' @export
qthill <- function(p, mean, sd, df) {
  qt(p, df) * std_dev + mean
}

#' Calculate the density of a double normal distribution
#'
#' @param x The value to evaluate the density at
#' @param mu The mean of the distribution
#' @param sigma1 The std deviation of the left half of the distribution
#' @param sigma2 The std deviation of the right half of the distribution
#' @importFrom stats dnorm
#' @export
ddnorm <- function(x, mu, sigma1, sigma2) {
  z <- log(2) - log(sigma1 + sigma2)
  if (x < mu) {
    z <- z + dnorm((x - mu) / sigma1, mean = 0, sd = 1, log = TRUE)
  } else {
    z <- z + dnorm((x - mu) / sigma2, mean = 0, sd = 1, log = TRUE)
  }
  return(z)
}

#' Calculate the density of a double Student-t distribution
#'
#' @param x The value to evaluate the density at
#' @param mu The mean of the distribution
#' @param sigma1 The std deviation of the left half of the distribution
#' @param sigma2 The std deviation of the right half of the distribution
#' @param tdf_1 The df deviation of the left half of the distribution
#' @param tdf_2 The df deviation of the right half of the distribution
#' @importFrom stats dt
#' @export
ddt <- function(x, mu, sigma1, sigma2, tdf_1, tdf_2) {
  z <- log(2) - log(sigma1 + sigma2)
  if (x < mu) {
    # Scaling the x value and adjusting for the degrees of freedom
    scaled_x <- (x - mu) / sigma1
    z <- z + dt(scaled_x, df = tdf_1, log = TRUE) - log(sigma1)
  } else {
    # Scaling the x value and adjusting for the degrees of freedom
    scaled_x <- (x - mu) / sigma2
    z <- z + dt(scaled_x, df = tdf_2, log = TRUE) - log(sigma2)
  }
  return(z)
}

#' Calculate the density of a double general normal distribution
#'
#' @param x The value to evaluate the density at
#' @param mu The mean of the distribution
#' @param alpha1 The alpha parameter of the left half of the distribution
#' @param alpha2 The alpha parameter of the right half of the distribution
#' @param beta1 The beta parameter of the left half of the distribution
#' @param beta2 The beta parameter of the right half of the distribution
#' @param sigma1 The sigma of the left half of the distribution
#' @param sigma2 The sigma of the right half of the distribution
#' @importFrom gnorm dgnorm
#' @export
ddgnorm <- function(x, mu, alpha1, alpha2, beta1, beta2, sigma1, sigma2) {
  z <- log(2) - log(sigma1 + sigma2)
  if (x < mu) {
    z <- z + dgnorm(x, mu, alpha1, beta1)
  } else {
    z <- z + dgnorm(x, mu, alpha2, beta2)
  }
  return(z)
}

#' Calculate the quantile of a double normal distribution
#'
#' @param p The value to evaluate the quantile function at
#' @param mu The mean of the distribution
#' @param sigma1 The sigma of the left half of the distribution
#' @param sigma2 The sigma of the right half of the distribution
#' @importFrom stats qnorm
#' @export
qdnorm <- function(p, mu, sigma1, sigma2) {
  r <- sigma1 / (sigma1 + sigma2)
  if (p < r) {
    z <- mu + sigma1 * qnorm(0.5 * p * (sigma1 + sigma2) / sigma1, 0, 1)
  } else {
    z <- mu + sigma2 * qnorm(0.5 * ((sigma1 + sigma2) * (1 + p) - 2 * sigma1) / sigma2, 0, 1)
  }
  return(z)
}

#' Calculate the quantile of a double Student-t distribution
#'
#' @param p The quantile to evaluate the function at
#' @param mu The mean of the distribution
#' @param sigma1 The std deviation of the left half of the distribution
#' @param sigma2 The std deviation of the right half of the distribution
#' @param tdf_1 The df deviation of the left half of the distribution
#' @param tdf_2 The df deviation of the right half of the distribution
#' @export
qdt <- function(p, mu, sigma1, sigma2, tdf_1, tdf_2) {
  r <- sigma1 / (sigma1 + sigma2)
  if (p < r) {
    z <- mu + sigma1 * qthill(0.5 * p * (sigma1 + sigma2) / sigma1, tdf_1, 0, 1)
  } else {
    z <- mu + sigma2 * qthill(0.5 * ((sigma1 + sigma2) * (1 + p) - 2 * sigma1) / sigma2, tdf_2, 0, 1)
  }
  return(z)
}

#' Calculate the quantile of a double general normal distribution
#'
#' @param p The value to evaluate the quantile function at
#' @param mu The mean of the distribution
#' @param sigma1 The sigma of the left half of the distribution
#' @param sigma2 The sigma of the right half of the distribution
#' @param beta_ratio_1 The beta ratio the left half of the distribution
#' @param beta_ratio_2 The beta ratio of the right half of the distribution
#' @param beta_1 The beta parameter of the left half of the distribution
#' @param beta_2 The beta parameter of the right half of the distribution
#' @importFrom gnorm qgnorm
#' @export
qdgnorm <- function(p, mu, sigma1, sigma2, beta_ratio_1, beta_ratio_2, beta_1, beta_2) {
  r <- sigma1 / (sigma1 + sigma2)
  if (p < r) {
    z <- mu + sigma1 * qgnorm(0.5 * p * (sigma1 + sigma2) / sigma1, mu, sigma1 * beta_ratio_1, beta_1)
  } else {
    z <- mu + sigma2 * qgnorm(0.5 * ((sigma1 + sigma2) * (1 + p) - 2 * sigma1) / sigma2, mu, sigma2 * beta_ratio_2, beta_2)
  }
  return(z)
}
