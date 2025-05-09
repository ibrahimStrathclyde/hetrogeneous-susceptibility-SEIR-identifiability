# ===============================================================================
# utility_functions.R
# 
# This file contains utility functions used throughout the analysis for:
# - Mathematical transformations (logit, expit)
# - Discretizing gamma distributions for heterogeneous models
# - Initializing model states
# - ODE system definitions for different SEIR models
# 
# Authors: M. Gabriela M. Gomes, Ibrahim Mohammed, Chris Robertson
# ===============================================================================

# ===============================================================================
# Mathematical utility functions
# ===============================================================================

#' Logit transformation
#' 
#' @param p A probability between 0 and 1
#' @return The logit of p: log(p/(1-p))
logit <- function(p) {
  log(p/(1-p))
}

#' Inverse logit (expit) transformation
#' 
#' @param x A real number
#' @return The inverse logit of x: 1/(1+exp(-x))
expit <- function(x) {
  1/(1+exp(-x))
}

# ===============================================================================
# Gamma distribution discretization
# ===============================================================================

#' Discretize a gamma distribution using multiplicative method
#' 
#' This function takes a gamma distribution with shape and rate alpha
#' (mean 1, variance 1/alpha) and discretizes it into n_groups intervals.
#' It uses a multiplicative method to ensure the discretized distribution
#' has exactly the same mean and approximately the same variance as the
#' continuous distribution.
#' 
#' @param n_groups Number of groups for discretization (must be even)


# ===============================================================================
# ODE model functions
# ===============================================================================

#' ODE function for heterogeneous SEIR model with time-varying transmission
#' 
#' @param t Current time point
#' @param y Vector of state variables
#' @param parms List of model parameters
#' @return List containing the derivatives of all state variables
hetsus_model.ct <- function(t, y, parms) {
  x <- parms[grepl("^x", names(parms))]
  parms_2 <- parms[!grepl("^x", names(parms))]
  z_list <- as.list(parms_2)
  z_list$x <- x  # This is necessary as x is a vector within parms
  z_list$t <- t  # x needs to be a named vector within this function 
  z_list$y <- y
  
  with(z_list, {
    S <- y[1:K]
    E <- y[(K+1):(2*K)]
    I <- y[(2*K+1):(3*K)]
    R <- y[(3*K+1):(4*K)]
    C <- y[(4*K+1):(5*K)]
    
    # Define intervention function
    prox <- 1.0  # Default value
    
    if (t <= t0) {
      prox <- c_value1
    } else if (t <= t1) {
      prox <- c_value1 - (c_value1 - c_value2) * (t - t0) / (t1 - t0)
    } else if (t <= t2) {
      prox <- c_value2
    } else if (t <= t3) {
      prox <- c_value3 + (c_value3 - c_value2) * (t - t3) / (t3 - t2)
    } else {
      prox <- c_value3
    }
    
    # Calculate transmission rate with intervention effect
    Beta <- R0*prox/(rho / delta + 1 / gamma)
    
    # SEIR model equations with heterogeneity
    dS <- -Beta*(rho*E+I)*(S/N)^(1+v^2)    
    dE <- Beta*(rho*E+I)*(S/N)^(1+v^2)-delta*E
    dI <- delta*E-gamma*I
    dR <- gamma*I
    dC <- delta*E
    return(list(c(dS, dE, dI, dR, dC)))
  })
}

# ===============================================================================
# Simulation functions
# ===============================================================================

#' Simulate cases using a reduced SEIR model with intervention
#' 
#' This function simulates epidemic data using a reduced SEIR model with
#' intervention (non-pharmaceutical interventions). The model can be homogeneous
#' (v=0) or heterogeneous (v>0).
#' 
#' @param R0 Basic reproduction number
#' @param delta Rate of leaving exposed compartment
#' @param rho Relative infectiousness in E compartment
#' @param gamma Recovery rate
#' @param v Coefficient of variation of susceptibility (0 for homogeneous)
#' @param N Total population size
#' @param E0 Initial number of exposed individuals
#' @param I0 Initial number of infectious individuals
#' @param t0 Time when adaptive behavior begins
#' @param t1 Time when lockdown begins
#' @param t2 Time when lockdown ends
#' @param t3 Time when transmission returns to baseline
#' @param c_value1 Initial transmission factor
#' @param c_value2 Intervention strength
#' @param c_value3 Final transmission factor
#' @param tfinal Final simulation time
#' @return List containing simulated data and plots
simulate_cases_reduced_model <- function(R0, delta, rho, gamma, v, N, E0, I0, 
                                         t0, t1, t2, t3, c_value1, c_value2, c_value3, tfinal) {
  
  params <- list(R0, delta, rho, gamma, v, N, E0, I0, t0, t1, t2, t3, c_value1, c_value2, c_value3, tfinal)
  names(params) <- c("R0", "delta", "rho", "gamma", "v", "N", "E0", "I0", "t0", 
                     "t1", "t2", "t3", "c_value1", "c_value2", "c_value3", "tfinal")
  
  parms <- c(R0, gamma, rho, delta, v, N, t0, t1, t2, t3, c_value1, c_value2, c_value3)
  names(parms) <- c("R0", "gamma", "rho", "delta", "v", "N", "t0", "t1", "t2", "t3", 
                    "c_value1", "c_value2", "c_value3")
  
  # Initial state
  initial_state <- c(S = N-E0-I0, E = E0, I = I0, R = 0, C = 0)
  
  # Time vector for simulation
  times <- seq(0, params$tfinal, by = 1)
  
  # Perform numerical integration
  result <- as.data.frame(ode(y = initial_state, times = times, func = Reduced.m_intervene, parms = parms))
  Daily_incidence <- c(0, diff(result[, "C"]))
  df <- result %>% mutate(Inc = Daily_incidence)
  
  # Generate Poisson-distributed cases
  cases <- rpois(length(df[, "Inc"]), lambda = df$Inc)
  
  # Create data frame of simulated data
  sim_data <- data.frame(time = times, reports = cases)
  
  # Create plots
  simgraph_Daily <- ggplot(sim_data, aes(x = time)) +
    geom_point(aes(y = reports), color = "red") +
    geom_vline(xintercept = t0, linetype = "dashed", color = "black") +
    labs(x = "Time", y = "Population", 
         title = paste0("Simulated data daily cases of covid 19 reduced model, R0= ", 
                        params$R0, ", CV= ", params$v)) +
    theme_minimal() +
    theme(legend.position = "top")
  
  simgraph_Prev <- ggplot(df, aes(x = time)) +
    geom_line(aes(y = I), color = "red") +
    geom_vline(xintercept = t0, linetype = "dashed", color = "black") +
    labs(x = "Time", y = "Population", 
         title = paste0("Simulated data Prevalence cases of covid 19 reduced model, R0= ", 
                        params$R0, ", CV= ", params$v)) +
    theme_minimal() +
    theme(legend.position = "top")
  
  return(list(sim_data = sim_data, simgraph.Prev. = simgraph_Prev, simgraph.Daily = simgraph_Daily))
}


#' ODE function for homogeneous SEIR model with time-varying transmission
#' 
#' @param time Current time point
#' @param y Vector of state variables
#' @param parms List of model parameters
#' @return List containing the derivatives of all state variables
seir.ct <- function(time, y, parms) {
  with(as.list(parms), {
    S <- y[1]
    E <- y[2]
    I <- y[3]
    R <- y[4]
    C <- y[5]
    
    # Define intervention function
    prox <- 1.0  # Default value
    
    if (time <= t0) {
      prox <- c_value1
    } else if (time <= t1) {
      prox <- c_value1 + (c_value2 - c_value1) * (time - t0) / (t1 - t0)
    } else if (time <= t2) {
      prox <- c_value2
    } else if (time <= t3) {
      prox <- c_value3 + (c_value3 - c_value2) * (time - t3) / (t3 - t2)
    } else {
      prox <- c_value3
    }
    
    # SEIR model equations
    dS <- -Beta*prox*(rho*E+I)*(S/N)^(1+v^2)    
    dE <- Beta*prox*(rho*E+I)*(S/N)^(1+v^2)-delta*E
    dI <- delta*E-gamma*I
    dR <- gamma*I
    dC <- delta*E
    return(list(c(dS, dE, dI, dR, dC)))
  })
}

#' ODE function for reduced model with intervention (homogeneous when v=0)
#' 
#' @param t Current time point
#' @param y Vector of state variables
#' @param parms List of model parameters
#' @return List containing the derivatives of all state variables
Reduced.m_intervene <- function(t, y, parms) {
  with(as.list(c(t, y, parms)), {
    S <- y[1]
    E <- y[2]
    I <- y[3]
    R <- y[4]
    C <- y[5]
    
    # Define intervention function
    prox <- 1.0  # Default value
    
    if (t <= t0) {
      prox <- c_value1
    } else if (t <= t1) {
      prox <- c_value1 - (c_value1 - c_value2) * (t - t0) / (t1 - t0)
    } else if (t <= t2) {
      prox <- c_value2
    } else if (t <= t3) {
      prox <- c_value3 + (c_value3 - c_value2) * (t - t3) / (t3 - t2)
    } else {
      prox <- c_value3
    }
    
    # Calculate transmission rate with intervention effect
    
    Beta <- R0*prox/(rho / delta + 1 / gamma)
    
    
    # SEIR model equations
    dS<- -Beta*(rho*E+I)*(S/N)^(1+v^2)    
    dE<- Beta*(rho*E+I)*(S/N)^(1+v^2)-delta*E
    dI<- delta*E-gamma*I
    dR<- gamma*I
    dC<-delta*E
    return(list(c(dS, dE, dI, dR,dC)))
  })
  
}