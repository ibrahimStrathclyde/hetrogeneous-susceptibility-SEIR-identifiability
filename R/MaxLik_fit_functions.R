# ===============================================================================
# MaxLik_fit_functions.R
#
# This file contains functions for fitting reduced SEIR models using maximum
# likelihood estimation. It includes functions for both homogeneous and
# heterogeneous models, single and dual epidemic fitting.
#
# Authors: M. Gabriela M. Gomes, Ibrahim Mohammed, Chris Robertson
# ===============================================================================

# ===============================================================================
# Single epidemic likelihood functions and parameter fitting
# ===============================================================================

#' Log-likelihood function for Poisson observations with reduced model with NPI
#'
#' @param params Model parameters including R0, v, intervention timing, etc.
#' @param sim.data Observed epidemic data with time and reports columns
#' @param initial_state Initial state vector for the SEIR model
#' @param times Vector of time points (if NULL, uses global 'times' variable)
#' @return Log-likelihood value (positive number)
poisson.loglik.NPI.reduced <- function(params, sim.data, initial_state, times = NULL) {
  # If times is not provided, use the global times variable
  if (is.null(times)) {
    if (!exists("times", envir = .GlobalEnv)) {
      stop("The 'times' variable is not defined in the global environment")
    }
    times <- get("times", envir = .GlobalEnv)
  }
  
  # Integrate the model equations to get the epidemic trajectory
  out <- as.data.frame(ode(
    y = initial_state, 
    times = times, 
    func = Reduced.m_intervene, 
    parms = params
  ))
  
  # Calculate daily incidence from cumulative cases
  Daily_incidence <- c(0, diff(out[, "C"]))
  df <- out %>% mutate(Inc = Daily_incidence)
  
  # Compute model predictions (expected counts)
  # Use a small positive number instead of zero to avoid log(0)
  lambda_ <- ifelse(df[,"Inc"] == 0, 0.0001, df[,"Inc"])
  
  # Validate sim.data
  if (!is.data.frame(sim.data)) {
    sim.data <- as.data.frame(sim.data)
  }
  
  # Check if the "reports" column exists
  if (!"reports" %in% names(sim.data)) {
    stop("The 'reports' column is missing in sim.data")
  }
  
  # Check if the lengths match
  if (nrow(sim.data) != nrow(df)) {
    stop("The lengths of sim.data and df do not match")
  }
  
  # Compute log likelihood using Poisson probability mass function
  loglik <- sum(dpois(
    x = sim.data[,"reports"], 
    lambda = lambda_, 
    log = TRUE
  ))
  
  return(loglik)
}

#' Objective function for optimization (returns negative log-likelihood)
#'
#' @param par Vector of transformed parameters to be estimated
#' @param sim.data Observed epidemic data
#' @return Negative log-likelihood value (for minimization)
f4_optim_reducedm.NPI <- function(par, sim.data) {
  # Transform parameters to their natural scale
  params <- c(
    R0 = exp(par[1]),           # R0 (positive)
    v = exp(par[2]),            # Coefficient of variation (positive)
    t0 = exp(par[3]),           # Intervention start time (positive)
    t1 = t1_spec,               # Ramp-up end time (fixed)
    t2 = t2_spec,               # Intervention end time (fixed)
    t3 = t3_spec,               # Ramp-down end time (fixed)
    c_value1 = c_value1_spec,   # Initial transmission factor (fixed)
    c_value2 = expit(par[4]),   # Intervention strength (0-1)
    c_value3 = c_value3_spec,   # Final transmission factor (fixed)
    rho = rho_spec,             # Relative infectiousness in E compartment (fixed)
    delta = delta_spec,         # Rate of transition from E to I (fixed)
    gamma = gamma_spec,         # Recovery rate (fixed)
    N = N,                      # Population size (fixed)
    tfinal = tfinal_spec        # Final time (fixed)
  )
  
  # Calculate the log-likelihood
  loglik <- poisson.loglik.NPI.reduced(
    params, 
    sim.data = sim.data, 
    initial_state = initial_state
  )
  
  # Return negative log-likelihood for minimization
  return(-loglik)
}

#' Function to fit the reduced model to a single epidemic dataset
#'
#' @param dat Data frame with time and reports columns
#' @return List containing parameter estimates, transformed parameters, and Hessian
fit4_reducedm_loglik.NPI <- function(dat) {
  # Starting values for optimization (on transformed scale)
  # log(2) for R0, log(2) for v, log(12) for t0, logit(0.4) for c_value2
  start_par <- c(log(2), log(2), log(12), logit(0.4))
  
  # Run optimization to find maximum likelihood estimates
  fit1 <- optim(
    par = start_par,
    fn = f4_optim_reducedm.NPI,
    sim.data = dat,
    method = "Nelder-Mead",
    control = list(trace = 0, maxit = 900),
    hessian = TRUE
  )
  
  # Calculate fitted parameters on natural scale
  fittedparams <- c(
    R0 = exp(fit1$par[1]),
    v = exp(fit1$par[2]),
    t0 = exp(fit1$par[3]),
    c_value2 = expit(fit1$par[4]),
    # Calculate AIC: 2*k + 2*minimized_negative_loglik
    # Note: fit1$value is the minimized negative log-likelihood
    AIC = 2 * length(fit1$par) + 2 * fit1$value,
    # Store the minimized negative log-likelihood
    negloglik = fit1$value,
    # Store the maximized log-likelihood
    loglik = -fit1$value,
    # Convergence code (0 indicates successful convergence)
    convergence = fit1$convergence
  )
  
  # Return results
  return(list(
    parms = fittedparams,        # Parameter estimates on natural scale
    trans_parms = fit1$par,      # Parameter estimates on transformed scale
    trans_hessian = fit1$hessian # Hessian matrix (for confidence intervals)
  ))
}

# ===============================================================================
# Dual epidemic likelihood functions and parameter fitting
# ===============================================================================

#' Function to fit the reduced model to two epidemic datasets simultaneously
#'
#' @param par Vector of transformed parameters to be estimated
#' @param sim.data_1 First epidemic dataset
#' @param sim.data_2 Second epidemic dataset
#' @return Combined negative log-likelihood
f4_2epi_optim_reduced.NPI <- function(par, sim.data_1, sim.data_2) {
  # Transform parameters to their natural scale
  params <- c(
    R0 = exp(par[1]),
    v = exp(par[2]),
    t0 = exp(par[3]),
    t1 = t1_spec,
    t2 = t2_spec,
    t3 = t3_spec,
    c_value1 = c_value1_spec,
    c_value2 = expit(par[4]),
    c_value3 = c_value3_spec,
    rho = rho_spec,
    delta = delta_spec,
    gamma = gamma_spec,
    N = N,
    tfinal = tfinal_spec
  )
  
  # Calculate log-likelihood for first epidemic
  loglik1 <- poisson.loglik.NPI.reduced(
    params, 
    sim.data = sim.data_1, 
    initial_state = initial_state_1
  )
  
  # Calculate log-likelihood for second epidemic
  loglik2 <- poisson.loglik.NPI.reduced(
    params, 
    sim.data = sim.data_2, 
    initial_state = initial_state_2
  )
  
  # Sum log-likelihoods (this is valid because the two epidemics are independent)
  loglik = loglik1 + loglik2
  
  # Return negative combined log-likelihood
  return(-loglik)
}

#' Function to fit the reduced model to two epidemic datasets
#'
#' @param dat1 First epidemic dataset (data frame with time and reports columns)
#' @param dat2 Second epidemic dataset (data frame with time and reports columns)
#' @return List containing parameter estimates, transformed parameters, and Hessian
fit4_2epic_reduced.loglik.NPI <- function(dat1, dat2) {
  # Starting values for optimization
  start_par <- c(log(2), log(2), log(12), logit(0.4))
  
  # Run optimization to find maximum likelihood estimates
  fit1 <- optim(
    par = start_par,
    fn = f4_2epi_optim_reduced.NPI,
    sim.data_1 = dat1,
    sim.data_2 = dat2,
    method = "Nelder-Mead",
    control = list(trace = 0, maxit = 400),
    hessian = TRUE
  )
  
  # Calculate fitted parameters on natural scale
  fittedparams <- c(
    R0 = exp(fit1$par[1]),
    v = exp(fit1$par[2]),
    t0 = exp(fit1$par[3]),
    c_value2 = expit(fit1$par[4]),
    # Calculate AIC: 2*k + 2*minimized_negative_loglik
    AIC = 2 * length(fit1$par) + 2 * fit1$value,
    # Store the minimized negative log-likelihood
    negloglik = fit1$value,
    # Store the maximized log-likelihood
    loglik = -fit1$value,
    # Convergence code
    convergence = fit1$convergence
  )
  
  # Return results
  return(list(
    parms = fittedparams,
    trans_parms = fit1$par,
    trans_hessian = fit1$hessian
  ))
}

# ===============================================================================
# Homogeneous model functions (single epidemic)
# ===============================================================================

#' Objective function for homogeneous model with NPI
#'
#' @param par Vector of transformed parameters to be estimated
#' @param sim.data Observed epidemic data
#' @return Negative log-likelihood value
f_optim_reducedm.poisloglikwithNPI <- function(par, sim.data) {
  # Transform parameters to their natural scale
  params <- c(
    R0 = exp(par[1]),
    v = v_spec,  # Fixed at 0 for homogeneous model
    t0 = exp(par[2]), 
    t1 = t1_spec, 
    t2 = t2_spec, 
    t3 = t3_spec,
    c_value1 = c_value1_spec,
    c_value2 = expit(par[3]),
    c_value3 = c_value3_spec,
    rho = rho_spec,
    delta = delta_spec,
    gamma = gamma_spec,
    N = N, 
    tfinal = tfinal_spec
  )
  
  # Calculate the negative log-likelihood
  loglik <- -poisson.loglik.withNPI.reduced(
    params, 
    sim.data = sim.data, 
    initial_state = initial_state
  )
  
  return(loglik)
}

#' Function to fit homogeneous model to a single epidemic with NPI
#'
#' @param dat Data frame with time and reports columns
#' @return List containing parameter estimates, transformed parameters, and Hessian
fit3_hom_1epic_loglikwithNPI <- function(dat) {
  fit <- optim(
    par = c(log(2), log(12), logit(0.2)), 
    fn = f_optim_reducedm.poisloglikwithNPI,
    sim.data = dat,
    method = "Nelder-Mead", 
    control = list(trace = 0, maxit = 300),
    hessian = TRUE
  )
  
  # Calculate fitted parameters on natural scale
  fittedparams <- c(
    R0 = exp(fit$par[1]),
    t0 = exp(fit$par[2]),
    c_value2 = expit(fit$par[3]),
    AIC = 2 * length(fit$par) - 2 * (-fit$value),
    value = fit$value,
    convergence = fit$convergence
  )
  
  return(list(
    parms = fittedparams,
    trans_parms = fit$par,
    trans_hessian = fit$hessian
  ))
}

# ===============================================================================
# Homogeneous model functions (dual epidemic)
# ===============================================================================

#' Objective function for homogeneous model with NPI for two epidemics
#'
#' @param par Vector of transformed parameters to be estimated
#' @param sim.data_1 First epidemic dataset
#' @param sim.data_2 Second epidemic dataset
#' @return Combined negative log-likelihood
f3_optim_reducedm.poisloglikwithNPI <- function(par, sim.data_1, sim.data_2) {
  # Transform parameters to their natural scale
  params <- c(
    R0 = exp(par[1]),
    v = v_spec,  # Fixed at 0 for homogeneous model
    t0 = exp(par[2]), 
    t1 = t1_spec, 
    t2 = t2_spec, 
    t3 = t3_spec,
    c_value1 = c_value1_spec,
    c_value2 = expit(par[3]),
    c_value3 = c_value3_spec,
    rho = rho_spec,
    delta = delta_spec,
    gamma = gamma_spec,
    N = N, 
    tfinal = tfinal_spec
  )
  
  # Calculate log-likelihood for first epidemic
  loglik1 <- -poisson.loglik.withNPI.reduced(
    params, 
    sim.data = sim.data_1,
    initial_state = initial_state_1
  )
  
  # Calculate log-likelihood for second epidemic
  loglik2 <- -poisson.loglik.withNPI.reduced(
    params, 
    sim.data = sim.data_2,
    initial_state = initial_state_2
  )
  
  # Sum negative log-likelihoods
  loglik = loglik1 + loglik2
  
  return(loglik)
}

#' Function to fit homogeneous model to two epidemics with NPI
#'
#' @param dat1 First epidemic dataset
#' @param dat2 Second epidemic dataset
#' @return List containing parameter estimates, transformed parameters, and Hessian
fit3_hom_2epic_loglikwithNPI <- function(dat1, dat2) {
  fit <- optim(
    par = c(log(2), log(12), logit(0.2)), 
    fn = f3_optim_reducedm.poisloglikwithNPI,
    sim.data_1 = dat1,
    sim.data_2 = dat2,
    method = "Nelder-Mead",
    control = list(trace = 0, maxit = 600),
    hessian = TRUE
  )
  
  # Calculate fitted parameters on natural scale
  fittedparams <- c(
    R0 = exp(fit$par[1]),
    t0 = exp(fit$par[2]),
    c_value2 = expit(fit$par[3]),
    AIC = 2 * length(fit$par) - 2 * (-fit$value),
    value = fit$value,
    convergence = fit$convergence
  )
  
  return(list(
    parms = fittedparams,
    trans_parms = fit$par,
    trans_hessian = fit$hessian
  ))
}

# ===============================================================================
# Additional utility functions for likelihood calculations
# ===============================================================================

#' Log-likelihood function for homogeneous model with NPI
#'
#' @param params Model parameters
#' @param sim.data Observed epidemic data
#' @param initial_state Initial state vector for the SEIR model
#' @return Log-likelihood value
poisson.loglik.withNPI.reduced <- function(params, sim.data, initial_state) {
  # Integrate the model equations
  out <- as.data.frame(ode(
    y = initial_state,
    times = times,
    func = Reduced.m_intervene,
    parms = params
  ))
  
  # Calculate daily incidence from cumulative cases
  Daily_incidence <- c(0, diff(out[, "C"]))
  df <- out %>% mutate(Inc = Daily_incidence)
  
  # Ensure values are valid for Poisson likelihood
  lambda_ <- ifelse(df[,"Inc"] == 0, 0.0001, df[,"Inc"])
  
  # Validate sim.data
  if (!is.data.frame(sim.data)) {
    sim.data <- as.data.frame(sim.data)
  }
  
  # Check if the "reports" column exists
  if (!"reports" %in% names(sim.data)) {
    stop("The 'reports' column is missing in sim.data")
  }
  
  # Check if the lengths match
  if (nrow(sim.data) != nrow(df)) {
    stop("The lengths of sim.data and df do not match")
  }
  
  # Compute log likelihood
  loglik <- sum(dpois(
    x = sim.data[,"reports"],
    lambda = lambda_,
    log = TRUE
  ))
  
  return(loglik)
}