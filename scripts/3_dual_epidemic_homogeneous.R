# ==============================================================================
# 3_dual_epidemic_homogeneous.R
# 
# This script fits a homogeneous SEIR model to two concurrent epidemic datasets
# with non-pharmaceutical interventions (NPIs).
# 
# The main goal is to assess if fitting to two epidemics improves parameter
# identifiability compared to fitting to a single epidemic.
#
# Authors: M. Gabriela M. Gomes, Ibrahim Mohammed, Chris Robertson
# ==============================================================================

# Load necessary libraries
library(tidyverse)
library(GGally)
library(deSolve)
library(gridExtra)
library(grid)

# Source required functions
source("R/MaxLik_fit_functions.R")
source("R/utility_functions.R")
source("R/visualization_functions.R")

# Create results directory if it doesn't exist
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# ==============================================================================
# Set parameters for simulation
# ==============================================================================
n_samples <- 100         # Number of simulated datasets
K <- 20                  # Number of groups for gamma discretization (for discretized model)
N <- 100000              # Population size

# First epidemic parameters
E0_1 <- 1040/1           # Number exposed at t=0 for first epidemic
I0_1 <- 416/1            # Number infected at t=0 for first epidemic

# Second epidemic parameters (smaller initial conditions)
E0_2 <- 60               # Number exposed at t=0 for second epidemic
I0_2 <- 24               # Number infected at t=0 for second epidemic

# Common parameters for both epidemics
v_spec <- 0              # Coefficient of variation (v=0 for homogeneous model)
R0_spec <- 3
delta_spec <- 1/5.5      # Rate of leaving exposed compartment
rho_spec <- 0            # Relative infectiousness in E compartment
gamma_spec <- 1/4        # Recovery rate
t0_spec <- 15            # Time when impact of adaptive behaviour started
t1_spec <- 20            # Time when lockdown begins
t2_spec <- 99            # Time when lockdown ends
t3_spec <- t2_spec + 1   # Time when transmission returned to baseline
tfinal_spec <- t3_spec   # Final simulation time
c_value1_spec <- 1       # Initial transmission factor
c_value2_spec <- 0.3     # Intervention strength
c_value3_spec <- 1       # Final transmission factor

# Set seed for reproducibility
set.seed(12357801L)

# ==============================================================================
# Simulate epidemic data for first epidemic
# ==============================================================================
if (exists("sim_data1")) rm(sim_data1)

for (i in 1:n_samples) {
  print(paste("Simulating first epidemic dataset", i))
  
  # Simulate cases using homogeneous model (v=0)
  dsimdiscre <- simulate_cases_reduced_model(
    R0 = R0_spec,
    delta = delta_spec,
    rho = rho_spec, 
    gamma = gamma_spec,
    v = v_spec,
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
  sim_data1 <- if (exists("sim_data1")) bind_rows(sim_data1, sim.data_1) else sim.data_1
}

# Set up global times variable
times <- sim.data_1$time

# ==============================================================================
# Simulate epidemic data for second epidemic (different initial conditions)
# ==============================================================================
# Set seed for second epidemic simulation
set.seed(98765)

if (exists("sim_data2")) rm(sim_data2)

for (i in 1:n_samples) {
  print(paste("Simulating second epidemic dataset", i))
  
  # Simulate cases using homogeneous model (v=0) but with different initial conditions
  dsimdiscre <- simulate_cases_reduced_model(
    R0 = R0_spec,
    delta = delta_spec,
    rho = rho_spec, 
    gamma = gamma_spec,
    v = v_spec,
    N = N,
    E0 = E0_2,
    I0 = I0_2,
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
  sim.data_2 <- dsimdiscre$sim_data
  sim.data_2$data_set <- i
  sim_data2 <- if (exists("sim_data2")) bind_rows(sim_data2, sim.data_2) else sim.data_2
}

# Create initial states for homogeneous model
initial_state_1 <- c(S = N - E0_1 - I0_1, E = E0_1, I = I0_1, R = 0, C = 0)
initial_state_2 <- c(S = N - E0_2 - I0_2, E = E0_2, I = I0_2, R = 0, C = 0)

# Plot some example trajectories for both epidemics
z_sel <- sample(1:n_samples, 5)

# First epidemic trajectories
z_df1 <- sim_data1 %>% 
  filter(data_set %in% z_sel) %>%
  mutate(epidemic = paste("epi1_", data_set, sep = ""))

# Second epidemic trajectories
z_df2 <- sim_data2 %>% 
  filter(data_set %in% z_sel) %>%
  mutate(epidemic = paste("epi2_", data_set, sep = ""))

# Create plot of example trajectories for first epidemic
traj_plot1 <- z_df1 %>% 
  ggplot(aes(x = time, y = reports, colour = epidemic)) +
  geom_point() +
  geom_vline(xintercept = t1_spec) + 
  geom_vline(xintercept = t0_spec, colour = "red") +
  labs(
    title = "Example first epidemic trajectories (homogeneous)",
    x = "Time (days)",
    y = "Daily cases",
    caption = paste("NPI strength =", c_value2_spec, ", Applied at t =", t0_spec)
  ) +
  theme_minimal()

# Create plot of example trajectories for second epidemic
traj_plot2 <- z_df2 %>% 
  ggplot(aes(x = time, y = reports, colour = epidemic)) +
  geom_point() +
  geom_vline(xintercept = t1_spec) + 
  geom_vline(xintercept = t0_spec, colour = "red") +
  labs(
    title = "Example second epidemic trajectories (homogeneous)",
    x = "Time (days)",
    y = "Daily cases",
    caption = paste("NPI strength =", c_value2_spec, ", Applied at t =", t0_spec)
  ) +
  theme_minimal()

# Display trajectory plots
combined_plot <- grid.arrange(traj_plot1, traj_plot2, ncol = 2)
ggsave("figures/homogeneous_dual_epidemic_trajectories.png", combined_plot, width = 12, height = 6)

# ==============================================================================
# Fit the homogeneous model to both epidemics simultaneously
# ==============================================================================
if (exists("results")) rm(results)

for (i in 1:n_samples) {
  print(paste("Fitting dataset pair", i))
  
  # Extract current datasets
  sim.data_1 <- sim_data1 %>% 
    filter(data_set == i) %>% 
    dplyr::select(-data_set)
  
  sim.data_2 <- sim_data2 %>% 
    filter(data_set == i) %>% 
    dplyr::select(-data_set)
  
  # Fit the homogeneous model with error handling
  z_mle <- tryCatch({
    fit3_hom_2epic_loglikwithNPI(dat1 = sim.data_1, dat2 = sim.data_2)
  }, error = function(e) {
    cat("Error in fitting dataset pair", i, ":", e$message, "\n")
    return(NULL)
  })
  
  # Skip this dataset if fitting failed
  if (is.null(z_mle)) {
    cat("Skipping dataset pair", i, "due to fitting error\n")
    next
  }
  
  # Set up results dataframe for this iteration
  i_results <- as.data.frame(matrix(z_mle$parms, nrow = 1))
  colnames(i_results) <- names(z_mle$parms)
  
  # Initialize values for Hessian results
  z_se <- numeric(length(z_mle$trans_parms))
  z_cor <- c(0, 0)
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
        
        # Extract key correlations (R0-t0, t0-c_value2)
        z_cor <- c(
          z_correlation[2, 1],  # R0_t0_trans_cor
          z_correlation[3, 2]   # t0_c_value2_trans_cor
        )
        
        # Calculate confidence intervals
        par_ucl <- z_mle$trans_parms + 1.96 * sqrt(diag(z_variance))
        par_lcl <- z_mle$trans_parms - 1.96 * sqrt(diag(z_variance))
        
        C_intervals <- as.data.frame(matrix(c(
          expit(par_lcl[1]) * (rho_spec / delta_spec + 1 / gamma_spec), 
          expit(par_ucl[1]) * (rho_spec / delta_spec + 1 / gamma_spec),
          exp(par_lcl[2]), exp(par_ucl[2]),
          expit(par_lcl[3]), expit(par_ucl[3])
        ), nrow = 1, byrow = TRUE))
        
        colnames(C_intervals) <- c(
          "R0_lcl", "R0_ucl", 
          "t0_lcl", "t0_ucl", 
          "c_value2_lcl", "c_value2_ucl"
        )
      }
    }, error = function(e) {
      cat("Error processing Hessian for dataset pair", i, ":", e$message, "\n")
      # Default values already set for z_se, z_cor, etc.
    })
  }
  
  # Create dataframes for transformed parameters
  z1 <- as.data.frame(matrix(z_mle$trans_parms, nrow = 1))
  colnames(z1) <- paste(names(z_mle$parms)[1:length(z_mle$trans_parms)], "trans", sep = "_")
  
  z_se <- as.data.frame(matrix(z_se, nrow = 1))
  colnames(z_se) <- paste(names(z_mle$parms)[1:length(z_mle$trans_parms)], "trans_se", sep = "_")
  
  z_cor <- as.data.frame(matrix(z_cor, nrow = 1))
  colnames(z_cor) <- c("R0_t0_trans_cor", "t0_c_value2_trans_cor")
  
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

# ==============================================================================
# Analyze results
# ==============================================================================

# Calculate summary statistics
result_summary <- results %>% 
  select(R0, R0_lcl, R0_ucl, t0, t0_lcl, t0_ucl, c_value2, c_value2_lcl, c_value2_ucl)

print(summary(result_summary))

# Save results
saveRDS(results, "results/results_2epi_hom_npi_0.3_case1.rds")

# Calculate bias metrics
bias_metrics <- results %>%
  filter(hess_pd == 1) %>%
  summarize(
    R0_mean = mean(R0, na.rm = TRUE),
    R0_median = median(R0, na.rm = TRUE),
    R0_abs_bias = mean(abs(R0 - R0_spec), na.rm = TRUE),
    R0_rel_bias_pct = 100 * mean(abs(R0 - R0_spec) / R0_spec, na.rm = TRUE),
    
    t0_mean = mean(t0, na.rm = TRUE),
    t0_median = median(t0, na.rm = TRUE),
    t0_abs_bias = mean(abs(t0 - t0_spec), na.rm = TRUE),
    t0_rel_bias_pct = 100 * mean(abs(t0 - t0_spec) / t0_spec, na.rm = TRUE),
    
    c_value2_mean = mean(c_value2, na.rm = TRUE),
    c_value2_median = median(c_value2, na.rm = TRUE),
    c_value2_abs_bias = mean(abs(c_value2 - c_value2_spec), na.rm = TRUE),
    c_value2_rel_bias_pct = 100 * mean(abs(c_value2 - c_value2_spec) / c_value2_spec, na.rm = TRUE),
    
    # Correlation summaries
    R0_t0_corr_mean = mean(R0_t0_trans_cor, na.rm = TRUE),
    t0_c_value2_corr_mean = mean(t0_c_value2_trans_cor, na.rm = TRUE)
  )

print(bias_metrics)
write.csv(bias_metrics, "results/bias_metrics_homogeneous_dual_epidemic.csv", row.names = FALSE)

# ==============================================================================
# Create parameter distribution plots
# ==============================================================================

# Plot histograms for each parameter with true value marked
R0_hist <- ggplot(results, aes(x = R0)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") +
  geom_vline(xintercept = R0_spec, color = "red", linetype = "dashed") +
  labs(title = "Distribution of R0 Estimates", x = "R0", y = "Count") +
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
  R0_hist, t0_hist, c_value2_hist,
  ncol = 2,
  top = textGrob("Parameter Estimate Distributions (Homogeneous Model, Two Epidemics)", 
                 gp = gpar(fontsize = 14, fontface = "bold"))
)

print(hist_grid)
ggsave("figures/homogeneous_dual_parameter_distributions.png", hist_grid, width = 10, height = 8)

# ==============================================================================
# Create parameter correlation plots
# ==============================================================================

# Filter valid results (only include results where Hessian was positive definite)
valid_results <- results %>% filter(hess_pd == 1)

# 1. Create correlation plot for R0 and t0
# Calculate correlation
r0_t0_correlation <- cor(valid_results$R0, valid_results$t0, use = "pairwise.complete.obs")

# Create the plot with mathematical notation in title and reference lines
r0_t0_plot <- ggplot(valid_results, aes(x = R0, y = t0)) +
  # Add scatter points
  geom_point(alpha = 0.5, size = 1, color = "gray30") +
  
  # Add density contours
  geom_density_2d(color = "royalblue", bins = 7) +
  
  # Add reference lines for true values
  geom_vline(xintercept = R0_spec, 
             linetype = "dashed", color = "brown4", alpha = 0.7) +
  
  geom_hline(yintercept = t0_spec, 
             linetype = "dashed", color = "red", alpha = 0.7) +
  
  # Add title with mathematical notation AND correlation value
  labs(
    title = bquote(R[0] ~ "and" ~ t[0] ~ "(r =" ~ .(round(r0_t0_correlation, 2)) ~ ")"),
    x = expression(R[0]),
    y = expression(t[0])
  ) +
  
  # Add clean theme
  theme_minimal() +
  theme(
    panel.grid.minor = element_line(color = "gray95"),
    panel.grid.major = element_line(color = "gray90"),
    plot.title = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 12)
  )

# 2. Create correlation plot for t0 and c_value2
# Calculate correlation
t0_c_correlation <- cor(valid_results$t0, valid_results$c_value2, use = "pairwise.complete.obs")

# Create the plot with mathematical notation in title and reference lines
t0_c_plot <- ggplot(valid_results, aes(x = t0, y = c_value2)) +
  # Add scatter points
  geom_point(alpha = 0.5, size = 1, color = "gray30") +
  
  # Add density contours
  geom_density_2d(color = "royalblue", bins = 7) +
  
  # Add reference lines for true values
  geom_vline(xintercept = t0_spec, 
             linetype = "dashed", color = "brown4", alpha = 0.7) +
  
  geom_hline(yintercept = c_value2_spec, 
             linetype = "dashed", color = "red", alpha = 0.7) +
  
  # Add title with mathematical notation AND correlation value
  labs(
    title = bquote(t[0] ~ "and" ~ c[1] ~ "(r =" ~ .(round(t0_c_correlation, 2)) ~ ")"),
    x = expression(t[0]),
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

# 3. Create correlation plot for R0 and c_value2
# Calculate correlation
r0_c_correlation <- cor(valid_results$R0, valid_results$c_value2, use = "pairwise.complete.obs")

# Create the plot with mathematical notation in title and reference lines
r0_c_plot <- ggplot(valid_results, aes(x = R0, y = c_value2)) +
  # Add scatter points
  geom_point(alpha = 0.5, size = 1, color = "gray30") +
  
  # Add density contours
  geom_density_2d(color = "royalblue", bins = 7) +
  
  # Add reference lines for true values
  geom_vline(xintercept = R0_spec, 
             linetype = "dashed", color = "brown4", alpha = 0.7) +
  
  geom_hline(yintercept = c_value2_spec, 
             linetype = "dashed", color = "red", alpha = 0.7) +
  
  # Add title with mathematical notation AND correlation value
  labs(
    title = bquote(R[0] ~ "and" ~ c[1] ~ "(r =" ~ .(round(r0_c_correlation, 2)) ~ ")"),
    x = expression(R[0]),
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

# Arrange the correlation plots
correlation_grid <- grid.arrange(
  r0_t0_plot, t0_c_plot, r0_c_plot,
  ncol = 3,
  top = textGrob("Key Parameter Correlations (Homogeneous Model, Two Epidemics)", 
                 gp = gpar(fontsize = 14, fontface = "bold"))
)

print(correlation_grid)
ggsave("figures/homogeneous_dual_key_correlations.png", correlation_grid, width = 15, height = 5)

# ==============================================================================
# Calculate confidence interval widths
# ==============================================================================
# Calculate confidence interval widths
ci_widths <- results %>%
  filter(hess_pd == 1) %>%
  mutate(
    R0_width = R0_ucl - R0_lcl,
    t0_width = t0_ucl - t0_lcl,
    c_value2_width = c_value2_ucl - c_value2_lcl
  ) %>%
  select(ends_with("_width"))

# Summarize CI widths
ci_summary <- ci_widths %>%
  summarize(
    R0_width_mean = mean(R0_width, na.rm = TRUE),
    t0_width_mean = mean(t0_width, na.rm = TRUE),
    c_value2_width_mean = mean(c_value2_width, na.rm = TRUE),
    
    R0_width_median = median(R0_width, na.rm = TRUE),
    t0_width_median = median(t0_width, na.rm = TRUE),
    c_value2_width_median = median(c_value2_width, na.rm = TRUE)
  )

print(ci_summary)
write.csv(ci_summary, "results/ci_widths_homogeneous_dual_epidemic.csv", row.names = FALSE)

# ==============================================================================
# Box plots for all parameters
# ==============================================================================
combined_long <- results %>%
  select(R0, t0, c_value2) %>%
  pivot_longer(cols = everything(), names_to = "Parameter", values_to = "Value")

# Boxplot for all parameters
boxplot <- ggplot(combined_long, aes(x = Parameter, y = Value, fill = Parameter)) +
  geom_boxplot(alpha = 0.7) +
  geom_hline(data = data.frame(
    Parameter = c("R0", "t0", "c_value2"),
    true_value = c(R0_spec, t0_spec, c_value2_spec)
  ), aes(yintercept = true_value), linetype = "dashed", color = "red") +
  labs(
    title = "Boxplot of Parameter Estimates (Homogeneous Model, Two Epidemics)",
    x = "Parameter", 
    y = "Value"
  ) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3")

print(boxplot)
ggsave("figures/homogeneous_dual_parameter_boxplots.png", boxplot, width = 8, height = 6)

# ==============================================================================
# Load single epidemic results for comparison
# ==============================================================================
if (file.exists("results/results_1epi_hom_npi_0.3_case1.rds")) {
  single_results <- readRDS("results/results_1epi_hom_npi_0.3_case1.rds")
  
  # Calculate correlation statistics for single epidemic
  single_cor_stats <- single_results %>%
    filter(hess_pd == 1) %>%
    summarize(
      R0_t0_corr = mean(R0_t0_trans_cor, na.rm = TRUE),
      t0_c_value2_corr = mean(t0_c_value2_trans_cor, na.rm = TRUE),
      R0_c_value2_corr = mean(R0_c_value2_trans_cor, na.rm = TRUE)
    )
  
  # Calculate correlation statistics for dual epidemics
  dual_cor_stats <- results %>%
    filter(hess_pd == 1) %>%
    summarize(
      R0_t0_corr = mean(R0_t0_trans_cor, na.rm = TRUE),
      t0_c_value2_corr = mean(t0_c_value2_trans_cor, na.rm = TRUE)
    )
  
  # Create a data frame for comparison
  cor_comparison <- data.frame(
    Parameter_Pair = c("R0-t0", "t0-c_value2"),
    Single_Epidemic = c(single_cor_stats$R0_t0_corr, single_cor_stats$t0_c_value2_corr),
    Dual_Epidemics = c(dual_cor_stats$R0_t0_corr, dual_cor_stats$t0_c_value2_corr)
  )
  
  print("Correlation comparison between single and dual epidemic fitting:")
  print(cor_comparison)
  write.csv(cor_comparison, "results/correlation_comparison_homogeneous.csv", row.names = FALSE)
  
  # Create a bar plot to compare correlations
  cor_comparison_long <- cor_comparison %>%
    pivot_longer(cols = c(Single_Epidemic, Dual_Epidemics), 
                 names_to = "Fitting_Type", values_to = "Correlation")
  
  correlation_comparison_plot <- ggplot(cor_comparison_long, 
                                        aes(x = Parameter_Pair, y = abs(Correlation), 
                                            fill = Fitting_Type)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    labs(
      title = "Absolute Parameter Correlations: Single vs. Dual Epidemic Fitting",
      subtitle = "Homogeneous SEIR Model",
      x = "Parameter Pair",
      y = "Absolute Correlation",
      fill = "Fitting Type"
    ) +
    scale_fill_brewer(palette = "Set1") +
    theme_minimal() +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  print(correlation_comparison_plot)
  ggsave("figures/homogeneous_correlation_comparison.png", correlation_comparison_plot, 
         width = 8, height = 6)
}

# ==============================================================================
# Print final conclusions
# ==============================================================================
cat("\n==============================================================================\n")
cat("CONCLUSIONS\n")
cat("==============================================================================\n\n")

cat(paste0("1. Parameter recovery was ", 
           ifelse(mean(c(bias_metrics$R0_rel_bias_pct, 
                         bias_metrics$t0_rel_bias_pct, 
                         bias_metrics$c_value2_rel_bias_pct)) < 10, 
                  "good", "moderate"), 
           " with mean relative bias: \n"))
cat(paste0("   - R0: ", round(bias_metrics$R0_rel_bias_pct, 1), "%\n"))
cat(paste0("   - t0: ", round(bias_metrics$t0_rel_bias_pct, 1), "%\n"))
cat(paste0("   - c_value2: ", round(bias_metrics$c_value2_rel_bias_pct, 1), "%\n\n"))

cat("2. Key parameter correlations for dual epidemic fitting:\n")
cat(paste0("   - R0 and t0: ", round(bias_metrics$R0_t0_corr_mean, 2), "\n"))
cat(paste0("   - t0 and c_value2: ", round(bias_metrics$t0_c_value2_corr_mean, 2), "\n\n"))

if (exists("cor_comparison")) {
  cat("3. Correlation comparison between single and dual epidemic fitting:\n")
  cat(paste0("   - R0-t0 correlation: ", 
             round(cor_comparison$Single_Epidemic[1], 2), " (single) vs ", 
             round(cor_comparison$Dual_Epidemics[1], 2), " (dual)\n"))
  cat(paste0("   - t0-c_value2 correlation: ", 
             round(cor_comparison$Single_Epidemic[2], 2), " (single) vs ", 
             round(cor_comparison$Dual_Epidemics[2], 2), " (dual)\n\n"))
  
  # Calculate reduction in correlation magnitude
  reduction_R0_t0 <- (abs(cor_comparison$Single_Epidemic[1]) - abs(cor_comparison$Dual_Epidemics[1])) / 
    abs(cor_comparison$Single_Epidemic[1]) * 100
  
  reduction_t0_c_value2 <- (abs(cor_comparison$Single_Epidemic[2]) - abs(cor_comparison$Dual_Epidemics[2])) / 
    abs(cor_comparison$Single_Epidemic[2]) * 100
  
  cat(paste0("   - Reduction in |R0-t0| correlation: ", round(reduction_R0_t0, 1), "%\n"))
  cat(paste0("   - Reduction in |t0-c_value2| correlation: ", round(reduction_t0_c_value2, 1), "%\n\n"))
  
  cat("4. Fitting to two concurrent epidemics ", 
      ifelse(mean(c(reduction_R0_t0, reduction_t0_c_value2)) > 20, 
             "significantly reduces", "somewhat reduces"),
      " parameter correlations\n")
  cat("   in the homogeneous model, improving parameter identifiability.\n\n")
} else {
  cat("3. Single epidemic results not available for comparison.\n\n")
}

cat("Next steps: Compare these results with the heterogeneous model fit to two epidemics\n")
cat("to assess if dual fitting has similar benefits for models with individual variation.\n")