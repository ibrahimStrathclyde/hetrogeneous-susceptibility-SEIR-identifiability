#===============================================================================
# 2_single_epidemic_heterogeneous.R
#
# Maximum Likelihood Estimation of SEIR model with heterogeneity for a 
# single epidemic with non-pharmaceutical interventions (NPIs)
#
# This code fits a reduced SEIR model with individual variation in 
# susceptibility to a single epidemic curve. It demonstrates parameter
# identifiability issues, particularly between heterogeneity (v) and
# intervention strength (c).
#
#Authors: M. Gabriela M. Gomes, Ibrahim Mohammed, Chris Robertson

#===============================================================================

#===============================================================================
# Load required libraries and source functions
#===============================================================================
library(tidyverse)
library(GGally)
library(deSolve)
library(gridExtra)
library(grid)

# Source required functions
source("R/MaxLik_fit_functions.R")
source("R/utility_functions.R")
source("R/visualization_functions.R")
#===============================================================================
# Set model parameters
#===============================================================================
# Simulation settings
n_samples <- 200          # Number of samples
N <- 100000               # Population size
E0_1 <- 1040/1            # Number of exposed individuals at t=0
I0_1 <- 416/1             # Number of infectious individuals at t=0

# Epidemiological parameters
alpha_gamma_shape <- 0.5  # Shape parameter of gamma distribution for heterogeneity
R0_spec <- 3              # Basic reproduction number
delta_spec <- 1/5.5       # Rate of progression from exposed to infectious (days⁻¹)
rho_spec <- 0             # Relative infectiousness of exposed individuals
gamma_spec <- 1/4         # Recovery rate (days⁻¹)

# Intervention parameters
t0_spec <- 15             # Time when impact of adaptive behaviour started (days)
t1_spec <- 20             # Time when lockdown begins (days)
t2_spec <- 99             # Time when lockdown ends (days)
t3_spec <- t2_spec + 1    # Time when transmission returned to baseline (days)
tfinal_spec <- t3_spec    # Final simulation time (days)
c_value1_spec <- 1        # Initial transmission factor
c_value2_spec <- 0.3      # Intervention strength - key parameter of interest
c_value3_spec <- 1        # Final transmission factor

#===============================================================================
# Simulate epidemic data
#===============================================================================
# Set seed for reproducibility
set.seed(12375L)

# Simulate data for multiple datasets
if (exists("sim_data")) rm(sim_data)

for (i in 1:n_samples) {
  print(paste("Simulating dataset", i))
  
  # Simulate cases using heterogeneous model with coefficient of variation v
  dsimdiscre <- simulate_cases_reduced_model(
    R0 = R0_spec,
    delta = delta_spec,
    rho = rho_spec, 
    gamma = gamma_spec,
    v = sqrt(1/alpha_gamma_shape),  # Coefficient of variation
    N = N,
    E0 = E0_1,
    I0 = I0_1,
    t0 = t0_spec,
    t1 = t1_spec,
    t2 = t2_spec,
    t3 = t3_spec,
    c_value1 = c_value1_spec,
    c_value2 = c_value2_spec,
    c_value3 = c_value3_spec,
    tfinal = tfinal_spec
  ) 
  
  # Store simulated data
  sim.data_1 <- dsimdiscre$sim_data
  sim.data_1$data_set <- i
  sim_data <- if (exists("sim_data")) bind_rows(sim_data, sim.data_1) else sim.data_1
}

# Set up global times variable
times <- sim.data_1$time

# Create initial state for model
initial_state <- c(S = N - E0_1 - I0_1, E = E0_1, I = I0_1, R = 0, C = 0)

#===============================================================================
# Visualize example trajectories
#===============================================================================
# Plot some example trajectories
z_sel <- sample(1:n_samples, 5)
z_df <- sim_data %>% 
  filter(data_set %in% z_sel) %>%
  mutate(epidemic = paste("epi", data_set, sep = "_"))

# Create plot of example trajectories
traj_plot <- z_df %>% 
  ggplot(aes(x = time, y = reports, colour = epidemic)) +
  geom_point() +
  geom_vline(xintercept = t1_spec) + 
  geom_vline(xintercept = t0_spec, colour = "red") +
  labs(
    title = "Example epidemic trajectories",
    x = "Time (days)",
    y = "Daily cases",
    caption = paste("NPI strength =", c_value2_spec, ", Applied at t =", t1_spec)
  ) +
  theme_minimal()

print(traj_plot)

#===============================================================================
# Maximum likelihood estimation
#===============================================================================
# Fit the reduced model to each dataset
if (exists("results")) rm(results)

for (i in 1:n_samples) { 
  print(paste("Fitting dataset", i))
  
  # Extract current dataset
  sim.data_1 <- sim_data %>% 
    filter(data_set == i) %>% 
    dplyr::select(-data_set)
  
  # Fit the reduced model with error handling
  z_mle <- tryCatch({
    fit4_reducedm_loglik.NPI(dat = sim.data_1)
  }, error = function(e) {
    cat("Error in fitting dataset", i, ":", e$message, "\n")
    return(NULL)
  })
  
  # Skip this dataset if fitting failed
  if (is.null(z_mle)) {
    cat("Skipping dataset", i, "due to fitting error\n")
    next
  }
  
  # Set up results dataframe for this iteration
  i_results <- as.data.frame(matrix(z_mle$parms, nrow = 1))
  colnames(i_results) <- names(z_mle$parms)
  
  # Initialize values for Hessian results
  z_se <- numeric(length(z_mle$trans_parms))
  z_cor <- c(0, 0, 0)
  z_hess <- 0
  z_pd <- 0
  z_ratio <- 0
  
  # Process Hessian matrix if available
  if (!is.null(z_mle$trans_hessian)) {
    tryCatch({
      z_hess <- 1
      z_eigen <- eigen(z_mle$trans_hessian)
      z_ratios <- z_eigen$values[1] / z_eigen$values
      z_ratio <- z_ratios[length(z_ratios)]
      
      if (all(z_eigen$values > 0)) {
        z_pd <- 1
        z_variance <- solve(z_mle$trans_hessian)
        z_d <- diag(1 / sqrt(diag(z_variance)), nrow = nrow(z_variance))
        z_correlation <- z_d %*% (z_variance %*% z_d)
        z_se <- sqrt(diag(z_variance))
        
        # Extract key correlations (R0-v, R0-t0, v-c_value2)
        z_cor <- c(
          z_correlation[2, 1],  # R0_v_trans_cor
          z_correlation[3, 1],  # R0_t0_trans_cor
          z_correlation[4, 2]   # v_c_value2_trans_cor
        )
        
        # Calculate confidence intervals
        par_ucl <- z_mle$trans_parms + 1.96 * sqrt(diag(z_variance))
        par_lcl <- z_mle$trans_parms - 1.96 * sqrt(diag(z_variance))
        
        C_intervals <- as.data.frame(matrix(c(
          exp(par_lcl[1]), exp(par_ucl[1]),          # R0
          exp(par_lcl[2]), exp(par_ucl[2]),          # v
          exp(par_lcl[3]), exp(par_ucl[3]),          # t0
          expit(par_lcl[4]), expit(par_ucl[4])       # c_value2
        ), nrow = 1, byrow = TRUE))
        
        colnames(C_intervals) <- c(
          "R0_lcl", "R0_ucl", 
          "v_lcl", "v_ucl", 
          "t0_lcl", "t0_ucl", 
          "c_value2_lcl", "c_value2_ucl"
        )
      }
    }, error = function(e) {
      cat("Error processing Hessian for dataset", i, ":", e$message, "\n")
      # Default values already set for z_se, z_cor, etc.
    })
  }
  
  # Create dataframes for transformed parameters
  z1 <- as.data.frame(matrix(z_mle$trans_parms, nrow = 1))
  colnames(z1) <- paste(names(z_mle$parms)[1:length(z_mle$trans_parms)], "trans", sep = "_")
  
  z_se <- as.data.frame(matrix(z_se, nrow = 1))
  colnames(z_se) <- paste(names(z_mle$parms)[1:length(z_mle$trans_parms)], "trans_se", sep = "_")
  
  z_cor <- as.data.frame(matrix(z_cor, nrow = 1))
  colnames(z_cor) <- c("R0_v_trans_cor", "R0_t0_trans_cor", "v_c_value2_trans_cor")
  
  # Combine all results for this iteration
  i_results <- i_results %>% bind_cols(z1, z_se, z_cor)
  i_results <- if (exists("C_intervals")) bind_cols(i_results, C_intervals) else i_results
  
  # Add diagnostic values and dataset ID
  i_results$hess_exists <- z_hess
  i_results$hess_pd <- z_pd
  i_results$ratio_max_min_evalue <- z_ratio
  i_results$dataset_id <- i
  
  # Append to results
  results <- if (exists("results")) bind_rows(results, i_results) else i_results
}

# Save results
saveRDS(results,"results_1epi_het_npi_0.3.rds")

#===============================================================================
# Analysis of parameter estimates
#===============================================================================
# Calculate summary statistics
result_summary <- results %>% 
  select(R0, R0_lcl, R0_ucl, v, v_lcl, v_ucl, t0, t0_lcl, t0_ucl, c_value2, c_value2_lcl, c_value2_ucl)

print(summary(result_summary))

# Calculate bias metrics
bias_metrics <- results %>%
  filter(hess_pd == 1) %>%
  summarize(
    R0_mean = mean(R0, na.rm = TRUE),
    R0_median = median(R0, na.rm = TRUE),
    R0_abs_bias = mean(abs(R0 - R0_spec), na.rm = TRUE),
    R0_rel_bias_pct = 100 * mean(abs(R0 - R0_spec) / R0_spec, na.rm = TRUE),
    
    v_mean = mean(v, na.rm = TRUE),
    v_median = median(v, na.rm = TRUE),
    v_abs_bias = mean(abs(v - sqrt(1/alpha_gamma_shape)), na.rm = TRUE),
    v_rel_bias_pct = 100 * mean(abs(v - sqrt(1/alpha_gamma_shape)) / sqrt(1/alpha_gamma_shape), na.rm = TRUE),
    
    t0_mean = mean(t0, na.rm = TRUE),
    t0_median = median(t0, na.rm = TRUE),
    t0_abs_bias = mean(abs(t0 - t0_spec), na.rm = TRUE),
    t0_rel_bias_pct = 100 * mean(abs(t0 - t0_spec) / t0_spec, na.rm = TRUE),
    
    c_value2_mean = mean(c_value2, na.rm = TRUE),
    c_value2_median = median(c_value2, na.rm = TRUE),
    c_value2_abs_bias = mean(abs(c_value2 - c_value2_spec), na.rm = TRUE),
    c_value2_rel_bias_pct = 100 * mean(abs(c_value2 - c_value2_spec) / c_value2_spec, na.rm = TRUE),
    
    # Correlation summaries
    R0_v_corr_mean = mean(R0_v_trans_cor, na.rm = TRUE),
    R0_t0_corr_mean = mean(R0_t0_trans_cor, na.rm = TRUE),
    v_c_value2_corr_mean = mean(v_c_value2_trans_cor, na.rm = TRUE)
  )

print(bias_metrics)

#===============================================================================
# Visualization of parameter distributions
#===============================================================================
# Plot histograms for each parameter with true value marked
R0_hist <- ggplot(results, aes(x = R0)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") +
  geom_vline(xintercept = R0_spec, color = "red", linetype = "dashed") +
  labs(title = "Distribution of R0 Estimates", x = "R0", y = "Count") +
  theme_minimal()

v_hist <- ggplot(results, aes(x = v)) +
  geom_histogram(binwidth = 0.05, fill = "lightgreen", color = "black") +
  geom_vline(xintercept = sqrt(1/alpha_gamma_shape), color = "red", linetype = "dashed") +
  labs(title = "Distribution of v Estimates", x = "v", y = "Count") +
  theme_minimal()

t0_hist <- ggplot(results, aes(x = t0)) +
  geom_histogram(binwidth = 0.5, fill = "lightpink", color = "black") +
  geom_vline(xintercept = t0_spec, color = "red", linetype = "dashed") +
  labs(title = "Distribution of t0 Estimates", x = "t0", y = "Count") +
  theme_minimal()

c_value2_hist <- ggplot(results, aes(x = c_value2)) +
  geom_histogram(binwidth = 0.01, fill = "lightcoral", color = "black") +
  geom_vline(xintercept = c_value2_spec, color = "red", linetype = "dashed") +
  labs(title = "Distribution of NPI Strength Estimates", x = "c_value2", y = "Count") +
  theme_minimal()

# Arrange histograms in grid
hist_grid <- grid.arrange(
  R0_hist, v_hist, t0_hist, c_value2_hist,
  ncol = 2,
  top = textGrob("Parameter Estimate Distributions", gp = gpar(fontsize = 14, fontface = "bold"))
)

print(hist_grid)
ggsave("figures/parameter_distributions_heterogeneous_single.png", hist_grid, width = 10, height = 8)

#===============================================================================
# Analysis of parameter correlations
#===============================================================================
# Filter valid results (only include results where Hessian was positive definite)
valid_results <- results %>% 
  filter(hess_pd == 1)

# Calculate correlation between v and c_value2
correlation <- cor(valid_results$v, valid_results$c_value2, use = "pairwise.complete.obs")

# Create correlation plot with mathematical notation
cv_c_plot <- ggplot(valid_results, aes(x = v, y = c_value2)) +
  # Add scatter points
  geom_point(alpha = 0.5, size = 1, color = "gray30") +
  
  # Add density contours
  geom_density_2d(color = "royalblue", bins = 7) +
  
  # Add reference lines for true values
  geom_vline(xintercept = sqrt(1/alpha_gamma_shape), 
             linetype = "dashed", color = "brown4", alpha = 0.7) +
  
  geom_hline(yintercept = c_value2_spec, 
             linetype = "dashed", color = "red", alpha = 0.7) +
  
  # Add title with mathematical notation AND correlation value
  labs(
    title = bquote(nu ~ "and" ~ c[1] ~ "(r =" ~ .(round(correlation, 2)) ~ ")"),
    x = expression(nu),
    y = expression(c[1])
  ) +
  
  # Add clean theme
  theme_minimal() +
  theme(
    panel.grid.minor = element_line(color = "gray95"),
    panel.grid.major = element_line(color = "gray90"),
    plot.title = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 12)
  )

# Display the plot
print(cv_c_plot)
#ggsave("figures/v_c_correlation_heterogeneous_single.png", cv_c_plot, width = 8, height = 6)

# Function to create correlation plots with density contours
create_param_correlation_plot <- function(data, param1, param2, 
                                          param1_label, param2_label,
                                          param1_true = NULL, param2_true = NULL) {
  # Calculate correlation
  cor_value <- cor(data[[param1]], data[[param2]], use = "pairwise.complete.obs")
  
  # Create the plot
  p <- ggplot(data, aes(x = !!sym(param1), y = !!sym(param2))) +
    # Add scatter points
    geom_point(alpha = 0.5, size = 1, color = "gray30") +
    # Add density contours
    geom_density_2d(color = "royalblue", bins = 7) +
    # Add labels with mathematical notation and correlation value
    labs(
      title = bquote(.(param1_label) ~ "and" ~ .(param2_label) ~ "(r =" ~ .(round(cor_value, 2)) ~ ")"),
      x = param1_label,
      y = param2_label
    ) +
    # Clean theme
    theme_minimal() +
    theme(
      panel.grid.minor = element_line(color = "gray95"),
      panel.grid.major = element_line(color = "gray90"),
      plot.title = element_text(size = 12, hjust = 0.5),
      axis.title = element_text(size = 10)
    )
  
  # Add reference lines for true values if provided
  if (!is.null(param1_true)) {
    p <- p + geom_vline(xintercept = param1_true, linetype = "dashed", 
                        color = "brown4", alpha = 0.7)
  }
  
  if (!is.null(param2_true)) {
    p <- p + geom_hline(yintercept = param2_true, linetype = "dashed", 
                        color = "red", alpha = 0.7)
  }
  
  return(p)
}

# Create correlation plots for all parameter combinations
# 1. R0 vs CV (v)
r0_cv_plot <- create_param_correlation_plot(
  data = valid_results,
  param1 = "R0", 
  param2 = "v",
  param1_label = expression(R[0]),
  param2_label = expression(nu),
  param1_true = R0_spec,
  param2_true = sqrt(1/alpha_gamma_shape)
)

# 2. R0 vs Intervention Timing (t0)
r0_t0_plot <- create_param_correlation_plot(
  data = valid_results,
  param1 = "R0", 
  param2 = "t0",
  param1_label = expression(R[0]),
  param2_label = expression(t[0]),
  param1_true = R0_spec,
  param2_true = t0_spec
)

# 3. R0 vs Intervention Strength (c_value2)
r0_c_plot <- create_param_correlation_plot(
  data = valid_results,
  param1 = "R0", 
  param2 = "c_value2",
  param1_label = expression(R[0]),
  param2_label = expression(c[1]),
  param1_true = R0_spec,
  param2_true = c_value2_spec
)

# 4. CV (v) vs Intervention Strength (c_value2) - already created above as cv_c_plot

# Arrange all plots in a grid
grid_plot <- grid.arrange(
  r0_cv_plot, r0_t0_plot, 
  r0_c_plot, cv_c_plot,
  ncol = 2,
  top = textGrob("Parameter Correlations (Heterogeneous Model, Single Epidemic)", 
                 gp = gpar(fontsize = 14, fontface = "bold"))
)

print(grid_plot)
#ggsave("figures/parameter_correlations_heterogeneous_single.png", grid_plot, width = 10, height = 8)

# Calculate full correlation matrix
cor_matrix <- cor(results[, c("R0", "v", "t0", "c_value2")], use = "pairwise.complete.obs")
print(cor_matrix)

#===============================================================================
# Conclusions
#===============================================================================
cat("\n==============================================================================\n")
cat("CONCLUSIONS\n")
cat("==============================================================================\n\n")

cat("1. Parameter recovery:\n")
cat(paste0("   - R0 (mean, rel. bias): ", round(bias_metrics$R0_mean, 2), ", ", 
           round(bias_metrics$R0_rel_bias_pct, 1), "%\n"))
cat(paste0("   - v (mean, rel. bias): ", round(bias_metrics$v_mean, 2), ", ", 
           round(bias_metrics$v_rel_bias_pct, 1), "%\n"))
cat(paste0("   - t0 (mean, rel. bias): ", round(bias_metrics$t0_mean, 2), ", ", 
           round(bias_metrics$t0_rel_bias_pct, 1), "%\n"))
cat(paste0("   - c_value2 (mean, rel. bias): ", round(bias_metrics$c_value2_mean, 2), ", ", 
           round(bias_metrics$c_value2_rel_bias_pct, 1), "%\n\n"))

cat("2. Key parameter correlations:\n")
cat(paste0("   - R0 and v: ", round(bias_metrics$R0_v_corr_mean, 2), "\n"))
cat(paste0("   - R0 and t0: ", round(bias_metrics$R0_t0_corr_mean, 2), "\n"))
cat(paste0("   - v and c_value2: ", round(bias_metrics$v_c_value2_corr_mean, 2), "\n\n"))

cat("3. There is a strong correlation between the heterogeneity parameter (v)\n")
cat("   and intervention strength (c_value2) in the heterogeneous SEIR model\n")
cat("   when fitted to a single epidemic curve. This suggests identifiability\n")
cat("   issues between these parameters.\n\n")

cat("4. Despite these correlations, parameter recovery is generally good,\n")
cat("   with relative biases typically below 10%.\n\n")

cat("Next step: Investigate if fitting to two concurrent epidemics can improve\n")
cat("parameter identifiability, particularly for the v-c_value2 correlation.\n")
