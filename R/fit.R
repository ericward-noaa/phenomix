#' Fitting function to be called by user
#'
#' This function creates a list of parameters, sets up TMB object and attempts to
#' do fitting / estimation
#'
#' @param data_list A list of data, as output from create_data
#' @param silent Boolean passed to TMB::MakeADFun, whether to be verbose or not (defaults to FALSE)
#' @param inits Optional named list of parameters for starting values, defaults to NULL
#' @param control Optional control list for stats::nlminb. For arguments see ?nlminb. Defaults to eval.max=2000, iter.max=1000, rel.tol=1e-10. For final model runs, the rel.tol should be even smaller
#' @param limits Whether to include limits for stats::nlminb. Can be a list of (lower, upper), or TRUE to use suggested hardcoded limits. Defaults to NULL,
#' where no limits used.
#' @param fit_model Whether to fit the model. If not, returns a list including the data, parameters, and initial values. Defaults to TRUE
#' @importFrom stats runif rnorm
#' @importFrom methods is
#' @import RTMB
#' @import gnorm
#' @export
#' @examples
#' data(fishdist)
#'
#' # example of fitting fixed effects, no trends, no random effects
# set.seed(2)
# datalist <- create_data(fishdist[which(fishdist$year > 1970), ],
#   asymmetric_model = FALSE,
#   est_mu_re = FALSE, est_sigma_re = FALSE
# )
# fit <- fit(datalist, silent=TRUE)
#' #
#' # # example of model with random effects in means only, and symmetric distribution
# set.seed(1)
# datalist <- create_data(fishdist[which(fishdist$year > 1970),], asymmetric_model = FALSE,
#                      est_sigma_re = FALSE)
# fit <- fit(datalist)
#' # # example of model with random effects in variances
# set.seed(1)
# datalist <- create_data(fishdist[which(fishdist$year > 1970),], asymmetric_model = TRUE,
#                          est_mu_re = TRUE)
# fit <- fit(datalist)
#' #
#' # # example of model with poisson response
# set.seed(1)
# datalist <- create_data(fishdist[which(fishdist$year > 1970),], asymmetric_model = FALSE,
#                         est_sigma_trend=FALSE, est_mu_trend=FALSE, est_mu_re = TRUE,
#                         family="poisson")
#' # fit <- fit(datalist)
fit <- function(data_list,
                silent = FALSE,
                inits = NULL,
                control = list(eval.max = 2000, iter.max = 1000, rel.tol = 1e-10),
                limits = NULL, # can also be a list, or TRUE
                fit_model = TRUE) {

  # create list of parameter starting values -- used in both
  # asymmetric and symmetric models
  parameters <- list(
    theta = rnorm(n = data_list$nLevels, log(mean(data_list$y[which(!is.na(data_list$y))]))),
    b_mu = rep(0, ncol(data_list$mu_mat)),#runif(ncol(data_list$mu_mat),0,10), #rep(0, ncol(data_list$mu_mat)),
    log_sigma_mu_devs = 0,
    mu_devs = rep(0, data_list$nLevels),
    b_sig1 = rep(1, ncol(data_list$sig_mat)),
    b_sig2 = rep(1, ncol(data_list$sig_mat)),
    log_sigma1_sd = 0,
    sigma1_devs = rep(0, data_list$nLevels),
    log_sigma2_sd = 0,
    sigma2_devs = rep(0, data_list$nLevels),
    log_obs_sigma = 0.0,
    log_tdf_1 = 0,
    log_tdf_2 = 0,
    log_beta_1 = 0.1,
    log_beta_2 = 0.1
  )
  parameters$b_mu[1] <- mean(data_list$x[which(!is.na(data_list$y))])
  #CV <- 0.1
  #parameters$b_sig1[1] <- 1# CV * parameters$b_mu[1]
  #parameters$b_sig2[1] <- parameters$b_sig1[1]
  if (data_list$family == 1) {
    # for gaussian, don't log-transform theta
    parameters$theta <- rep(0, data_list$nLevels) # rnorm(n = data_list$nLevels, mean(data_list$y))
  }

  # If inits is included, use that instead of parameters
  # if (!is.null(inits)) parameters <- inits

  # Mapping off params as needed:
  tmb_map <- list()
  if (data_list$family %in% c(2, 4)) {
    # don't include obs_sigma for poisson or binomial
    tmb_map <- c(tmb_map, list(log_obs_sigma = as.factor(NA)))
  }

  if (data_list$asymmetric == 0) {
    # map off pars not needed
    tmb_map <- c(tmb_map, list(
      b_sig2 = rep(as.factor(NA), ncol(data_list$sig_mat)),
      log_sigma2_sd = as.factor(NA),
      sigma2_devs = rep(as.factor(NA), data_list$nLevels)
    ))
  }

  if (data_list$tail_model == 0) {
    # then fit gaussian model and map off both parameters
    tmb_map <- c(tmb_map, list(
      log_tdf_1 = as.factor(NA),
      log_tdf_2 = as.factor(NA),
      log_beta_1 = as.factor(NA),
      log_beta_2 = as.factor(NA)
    ))
  }
  if (data_list$tail_model == 1) {
    # fit the t-model
    if (data_list$asymmetric == 0) {
      # then map off the tdf2, because model is symmetric
      tmb_map <- c(tmb_map, list(
        log_tdf_2 = as.factor(NA),
        log_beta_1 = as.factor(NA),
        log_beta_2 = as.factor(NA)
      ))
    } else {
      tmb_map <- c(tmb_map, list(
        log_beta_1 = as.factor(NA),
        log_beta_2 = as.factor(NA)
      ))
    }
  }
  if (data_list$tail_model == 2) {
    # then fit gaussian model and map off both parameters
    if (data_list$asymmetric == 0) {
      # then map off the beta2, because model is symmetric
      tmb_map <- c(tmb_map, list(
        log_beta_2 = as.factor(NA),
        log_tdf_1 = as.factor(NA),
        log_tdf_2 = as.factor(NA)
      ))
    } else {
      tmb_map <- c(tmb_map, list(
        log_tdf_1 = as.factor(NA),
        log_tdf_2 = as.factor(NA)
      ))
    }
  }
  if (data_list$asymmetric == 1) {
    if (data_list$share_shape == 1) {
      if (data_list$tail_model == 1) {
        # map off 2nd nu parameter
        tmb_map <- c(tmb_map, list(
          log_tdf_2 = as.factor(NA)
        ))
      }
      if (data_list$tail_model == 2) {
        tmb_map <- c(tmb_map, list(
          log_beta_2 = as.factor(NA)
        ))
      }
    }
  }

  if (data_list$nLevels == 1) {
    data_list$est_mu_re <- 0
    data_list$est_sigma_re <- 0
  }

  # if don't estimate mean random effects map them off
  if (data_list$est_mu_re == 0) {
    tmb_map <- c(tmb_map, list(
      log_sigma_mu_devs = as.factor(NA),
      mu_devs = rep(as.factor(NA), data_list$nLevels)
    ))
  }

  if (data_list$est_sigma_re == 0) {
    tmb_map <- c(tmb_map, list(
      log_sigma1_sd = as.factor(NA),
      log_sigma2_sd = as.factor(NA),
      sigma1_devs = rep(as.factor(NA), data_list$nLevels),
      sigma2_devs = rep(as.factor(NA), data_list$nLevels)
    ))
  } else {
    if (data_list$asymmetric == 0) {
      tmb_map <- c(tmb_map, list(
        log_sigma2_sd = as.factor(NA),
        sigma2_devs = rep(as.factor(NA), data_list$nLevels)
      ))
    }
  }

  random <- ""
  if (data_list$est_mu_re == 1) {
    random <- c(random, "mu_devs")
  }
  if (data_list$est_sigma_re == 1) {
    random <- c(random, "sigma1_devs")
    if (data_list$asymmetric == TRUE) random <- c(random, "sigma2_devs")
  }
  random <- random[-1]
  if (length(random) == 0) random <- NULL

  # perhaps it'd be better to split this function out?
  f <- function(parms) {

    RTMB::getAll(data_list, parms, warn=FALSE)
    nll<-0
    y <- RTMB::OBS(y) # optional, https://cran.r-project.org/web/packages/RTMB/vignettes/RTMB-introduction.html
    # derived parameters
    obs_sigma <- RTMB:::getValues(exp(log_obs_sigma))
    tdf_1 <- exp(log_tdf_1) + 2
    tdf_2 <- exp(log_tdf_2) + 2
    beta_1 <- exp(log_beta_1)
    beta_2 <- exp(log_beta_2)
    if(share_shape==1) {
      tdf_2 <- tdf_1
      beta_2 <- beta_1
    }

    if(use_t_prior & tail_model == 1) {
      #https:#statmodeling.stat.columbia.edu/2015/05/17/do-we-have-any-recommendations-for-priors-for-student_ts-degrees-of-freedom-parameter/
        nll <- nll + dgamma(RTMB:::getValues(tdf_1), nu_prior[1], nu_prior[2], log=TRUE)
        if(asymmetric == 1) nll <- nll + dgamma(RTMB:::getValues(tdf_2), nu_prior[1], nu_prior[2], log=TRUE)
    }
    if(use_beta_prior & tail_model == 2) {
      # value ~ 2-3 is normal
      nll <- nll + dgamma(RTMB:::getValues(beta_1), beta_prior[1], beta_prior[2], log=TRUE)
      if(asymmetric == 1) nll <- nll + dgamma(RTMB:::getValues(beta_2), beta_prior[1], beta_prior[2], log=TRUE)
    }

    sigma1 <- sigma2 <- mu <- alpha1 <- alpha2 <- lower25 <- upper75 <- range <- rep(NA, nLevels)
    n <- length(y)

    # calculations for beta for gnorm dist if implemented
    beta_ratio <- rep(NA,2)
    if(tail_model == 2) {
      beta_ratio[1] <- sqrt(exp(lgamma(1.0/RTMB:::getValues(beta_1))) / exp(lgamma(3.0/RTMB:::getValues(beta_1))))
      if(asymmetric == 1) {
        beta_ratio[2] <- sqrt(exp(lgamma(1.0/RTMB:::getValues(beta_2))) / exp(lgamma(3.0/RTMB:::getValues(beta_2))))
      }
    }

    # random effects components
    for(i in 1:nLevels) {
      # random effects contributions of mean and sigma1
      if(est_mu_re==1) {
        nll <- nll + dnorm(RTMB:::getValues(mu_devs[i]), 0, exp(RTMB:::getValues(log_sigma_mu_devs)), log=TRUE)
      }
      if(est_sigma_re==1) {
        nll <- dnorm(RTMB:::getValues(sigma1_devs[i]), 0, exp(RTMB:::getValues(log_sigma1_sd)), log=TRUE)
        if(asymmetric == 1) {
          nll <- dnorm(RTMB:::getValues(sigma2_devs[i]), 0, exp(RTMB:::getValues(log_sigma2_sd)), log=TRUE)
        }
      }
    }

    # fixed effects for mu and both sigmas
    mu <- mu_mat * RTMB:::getValues(b_mu)
    sigma1 <- sig_mat * RTMB:::getValues(b_sig1)
    if(asymmetric == 1) sigma2 <- sig_mat * RTMB:::getValues(b_sig2)

    for(i in 1:nLevels) {

      if(est_sigma_re == 1) {
        sigma1[i] <- sigma1[i] + RTMB:::getValues(sigma1_devs[i])
        if(asymmetric == 1) sigma2[i] <- sigma2[i] + RTMB:::getValues(sigma2_devs[i])
      }

      # calculate alphas if the gnorm model is used
      if(tail_model == 2) {
        alpha1[i] <- sigma1[years[i]]*beta_ratio[1]
        if(asymmetric==1) {
          alpha2[i] <- sigma2[years[i]]*beta_ratio[2]
        }
      }

      # trend in in normal space, e.g. not log-linear
      if(est_mu_re == 1) {
        mu[i] <- mu[i] + RTMB:::getValues(mu_devs[i])
      }

      # this is all for calculating quantiles of normal distributions
      if(tail_model==0) {
        if(asymmetric == 0) {
          lower25[i] <- qnorm(0.25, mu[i], sigma1[i])
          upper75[i] <- qnorm(0.75, mu[i], sigma1[i])
        } else {
          lower25[i] <- qdnorm(0.25, mu[i], sigma1[i], sigma2[i])
          upper75[i] <- qdnorm(0.75, mu[i], sigma1[i], sigma2[i])
        }
      }
      # this is all for calculating quantiles of Student-t distributions
      if(tail_model==1) {
        if(asymmetric == 0) {
          lower25[i] <- qthill(0.25, RTMB:::getValues(tdf_1), mu[i], sigma1[i])
          upper75[i] <- qthill(0.75, RTMB:::getValues(tdf_1), mu[i], sigma1[i])
        } else {
          lower25[i] <- qdt(0.25, mu[i], sigma1[i], sigma2[i], RTMB:::getValues(tdf_1), RTMB:::getValues(tdf_2))
          upper75[i] <- qdt(0.75, mu[i], sigma1[i], sigma2[i], RTMB:::getValues(tdf_1), RTMB:::getValues(tdf_2))
        }
      }

      # this is all for calculating quantiles of Student-t distributions
      if(tail_model==2) {
        if(asymmetric == 0) {
          lower25[i] <- qgnorm(0.25, mu[i], sigma1[i]*beta_ratio[1], RTMB:::getValues(beta_1))
          upper75[i] <- qgnorm(0.75, mu[i], sigma1[i]*beta_ratio[1], RTMB:::getValues(beta_1))
        } else {
          lower25[i] <- qdgnorm(0.25, mu[i], sigma1[i], sigma2[i], beta_ratio[1], beta_ratio[2], RTMB:::getValues(beta_1), RTMB:::getValues(beta_2))
          upper75[i] <- qdgnorm(0.75, mu[i], sigma1[i], sigma2[i], beta_ratio[1], beta_ratio[2], RTMB:::getValues(beta_1), RTMB:::getValues(beta_2))
        }
      }
      range[i] = upper75[i] - lower25[i]
    }

    # this is for the predictions
    log_dens <- pred <- rep(NA, n)
    for(i in 1:n) {
      if(asymmetric == 1) {
        # model is asymmetric, left side smaller / right side bigger
        if(tail_model==0) {
          log_dens[i] <- ddnorm(x[i], mu[years[i]], sigma1[years[i]], sigma2[years[i]])
        }
        if(tail_model==1) {
          log_dens[i] <- ddt(x[i], mu[years[i]], sigma1[years[i]], sigma2[years[i]], RTMB:::getValues(tdf_1), RTMB:::getValues(tdf_2))
        }
        if(tail_model==2) {
          log_dens[i] <- ddgnorm(x[i], mu[years[i]], alpha1[years[i]], alpha2[years[i]], RTMB:::getValues(beta_1), RTMB:::getValues(beta_2), sigma1[years[i]], sigma2[years[i]])
        }
        pred[i] <- log_dens[i] + RTMB:::getValues(theta[years[i]])

      } else {
        if(tail_model==0) {
          # model is symmetric around mu, gaussian tails
          log_dens[i] <- dnorm(x[i], mu[years[i]], sigma1[years[i]], log=TRUE)
        } else {
          if(tail_model==1) {
            # model is symmetric around mu, student-t tails
            log_dens[i] <- dt((x[i] - mu[years[i]]) / sigma1[years[i]], RTMB:::getValues(tdf_1), log=TRUE) - log(sigma1[years[i]])
          } else {
            # gnorm, copied from maryclare/gnorm
            # alpha = sqrt( var * gamma(1/beta) / gamma(3/beta) ), alpha = sigma(1)*beta_ratio(1)
            log_dens[i] <- dgnorm(x[i], mu[years[i]], alpha1[years[i]], RTMB:::getValues(beta_1))
          }
        }
        pred[i] <- log_dens[i] + RTMB:::getValues(theta[years[i]])
      }
    }

    # this is for the cumulative annual predictions
    year_log_tot <- rep(0, nLevels)
    year_tot <- rep(0, nLevels)
    for(i in 1:nLevels) {

      for(t in 1:366) {

        if(asymmetric == 1) {
          # model is asymmetric, left side smaller / right side bigger
          if(tail_model==0) {
            dens <- ddnorm(t, mu[i], sigma1[i], sigma2[i])
          }
          if(tail_model==1) {
            dens <- ddt(t, mu[i], sigma1[i], sigma2[i], RTMB:::getValues(tdf_1), RTMB:::getValues(tdf_2))
          }
          if(tail_model==2) {
            dens <- ddgnorm(t, mu[i], alpha1[i], alpha2[i], RTMB:::getValues(beta_1), RTMB:::getValues(beta_2), sigma1[i], sigma2[i])
          }
          year_log_tot[i] <- year_log_tot[i] + dens + RTMB:::getValues(theta[i])
          year_tot[i] <- year_tot[i] + exp(dens + RTMB:::getValues(theta[i]))

        } else {
          if(tail_model==0) {
            # model is symmetric around mu, gaussian tails
            dens <- dnorm(t, mu[i], sigma1[i], log=TRUE)
          } else {
            if(tail_model==1) {
              # model is symmetric around mu, student-t tails
              dens <- dt((t - mu[i]) / sigma1[i], RTMB:::getValues(tdf_1), log=TRUE) - log(sigma1[i])
            } else {
              # gnorm, copied from maryclare/gnorm
              # alpha = sqrt( var * gamma(1/beta) / gamma(3/beta) ), alpha = sigma(1)*beta_ratio(1)
              dens <- dgnorm(t, mu[i], alpha1[i], RTMB:::getValues(beta_1))
            }
          }
          year_log_tot[i] <- year_log_tot[i] + dens + RTMB:::getValues(theta[i])
          year_tot[i] <- year_tot[i] + exp(dens + RTMB:::getValues(theta[i]))
        }
      }
    }

    # this is the likelihood
    s1 <- 0
    s2 <- 0

    if(family==1) {
      # gaussian, both data and predictions in log space
      nll <- nll + sum(dnorm(y, pred, obs_sigma, log=TRUE))
    }
    if(family==2) {
      for(i in 1:n) {
        if(pred[i] > 20) pred[i] <- 20 # not needed
        nll <- dpois(y[i], exp(pred[i]), log=TRUE)
      }
    }
    if(family==3) {
      for(i in 1:n) {
        s1 = pred[i]
        #s2 = s1 + pow(s1, Type(2))*obs_sigma
        s2 <- 2 * s1 - RTMB:::getValues(log_obs_sigma) # log(var - mu)
        nll <- nll + RTMB::dnbinom_robust(y[i], s1, s2, log=TRUE)
      }
    }
    if(family==4) {
      for(i in 1:n) {
        nll <- nll + RTMB::dbinom_robust(y[i], 1, pred[i], log=TRUE)
      }
    }
    if(family==5) {
      # lognormal, both data and predictions in log space
      nll <- sum(dnorm(log(y), pred, obs_sigma, log=TRUE))
    }
    # ADREPORT section
    RTMB::ADREPORT(theta) # nuisance parameter
    RTMB::ADREPORT(sigma1) # sigma, LHS
    RTMB::ADREPORT(mu) # mean of curves by year
    RTMB::ADREPORT(b_mu) # sigma, LHS
    RTMB::ADREPORT(b_sig1) # mean of curves by year
    RTMB::ADREPORT(year_tot) # mean of curves by year
    RTMB::ADREPORT(year_log_tot) # mean of curves by year

    if(family %in% c(2,4) == FALSE) {
      RTMB::ADREPORT(obs_sigma) # obs sd (or phi, NB)
    }
    RTMB::ADREPORT(pred) # predictions in link space (log)

    RTMB::ADREPORT(lower25) # lower quartile
    RTMB::ADREPORT(upper75) # upper quartile
    RTMB::ADREPORT(range) # diff between upper and lower quartiles
    if(tail_model==1) {
      RTMB::ADREPORT(tdf_1) # tdf for LHS
    }
    if(tail_model==2) {
      RTMB::ADREPORT(beta_1) # tdf for LHS
    }
    if(asymmetric == 1) {
      # these are only reported for asymmetric model
      RTMB::ADREPORT(b_sig2) # mean of curves by year
      if(est_sigma_re==1) {
        RTMB::ADREPORT(sigma2) # same as above, but RHS optionally
      }
      if(tail_model==1) {
        RTMB::ADREPORT(tdf_2)
      }
      if(tail_model==2) {
        RTMB::ADREPORT(beta_2)
      }
    }
    return (-nll)

  }

  obj <- RTMB::MakeADFun(
    func = f,
    #data = data_list,
    parameters = parameters,
    map = tmb_map,
    #DLL = "phenomix",
    random = random,
    silent = silent
  )

  mod_list <- list(
    obj = obj,
    init_vals = obj$par,
    # pars = pars,
    # sdreport = sdreport,
    data_list = data_list
  )

  if (fit_model == TRUE) {
    # Attempt to do estimation
    if (!is.null(inits)) {
      init <- inits
    } else {
      init <- obj$par
    }
    if (is.null(limits)) {
      pars <- stats::nlminb(
        start = init,
        objective = obj$fn,
        gradient = obj$gr,
        control = control
      )
    } else {
      if (is(limits,"list")) {
        lower_limits <- limits$lower
        upper_limits <- limits$upper
      } else {
        lim <- limits(parnames = names(obj$par), max_theta = data_list$max_theta)
        lower_limits <- lim$lower
        upper_limits <- lim$upper
      }

      pars <- stats::nlminb(
        start = init, objective = obj$fn,
        gradient = obj$gr, control = control,
        lower = lower_limits,
        upper = upper_limits
      )
    }

    sdreport <- RTMB::sdreport(obj)
    mod_list <- c(mod_list, list(pars = pars, sdreport = sdreport, tmb_map = tmb_map, tmb_random = random, tmb_parameters = parameters))
  }

  return(structure(mod_list, class = "phenomix"))
}
