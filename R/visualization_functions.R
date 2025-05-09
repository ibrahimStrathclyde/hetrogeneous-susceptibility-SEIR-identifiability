# ===============================================================================
# visualization_functions.R
#
# This file contains functions for visualizing results from SEIR model fitting.
# It includes functions for creating parameter distribution plots, correlation
# plots, and other visualizations to assess parameter identifiability.
#
# Authors: M. Gabriela M. Gomes, Ibrahim Mohammed, Chris Robertson
# ===============================================================================

#' Create parameter correlation plot with density contours
#'
#' @param data Data frame containing parameter estimates
#' @param param1 Name of first parameter column
#' @param param2 Name of second parameter column
#' @param param1_label Expression for first parameter label
#' @param param2_label Expression for second parameter label
#' @param param1_true True value of first parameter (optional)
#' @param param2_true True value of second parameter (optional)
#' @return ggplot object
#'
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
    # Add title with mathematical notation AND correlation value
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

#' Create a parameter correlation matrix plot
#'
#' @param data Data frame containing parameter estimates
#' @param params Vector of parameter names to include in plot
#' @param true_values Named list of true parameter values
#' @param title Plot title
#' @return ggpairs object
#'
create_correlation_matrix <- function(data, params, true_values = NULL, title = "Parameter Correlations") {
  plot <- data %>% 
    select(all_of(params)) %>% 
    ggpairs(
      title = title,
      upper = list(continuous = wrap("cor", size = 5)),
      lower = list(continuous = wrap("points", alpha = 0.3))
    )
  
  return(plot)
}

#' Create parameter distribution histograms with true values marked
#'
#' @param data Data frame containing parameter estimates
#' @param param Parameter column name
#' @param param_name Formatted parameter name for plot title
#' @param true_value True parameter value
#' @param binwidth Width of histogram bins
#' @param fill_color Fill color for histogram
#' @return ggplot object
#'
create_parameter_histogram <- function(data, param, param_name, true_value, 
                                       binwidth = NULL, fill_color = "skyblue") {
  if (is.null(binwidth)) {
    # Auto-calculate binwidth if not provided
    binwidth <- (max(data[[param]], na.rm = TRUE) - min(data[[param]], na.rm = TRUE)) / 20
  }
  
  plot <- ggplot(data, aes(x = !!sym(param))) +
    geom_histogram(binwidth = binwidth, fill = fill_color, color = "black") +
    geom_vline(xintercept = true_value, color = "red", linetype = "dashed") +
    labs(title = paste("Distribution of", param_name, "Estimates"), 
         x = param_name, y = "Count") +
    theme_minimal()
  
  return(plot)
}

#' Create a grid of parameter histograms
#'
#' @param data Data frame containing parameter estimates
#' @param params Named list of parameters to plot (names = column names, values = display names)
#' @param true_values Named list of true parameter values
#' @param binwidths Named list of binwidths for each parameter (optional)
#' @param fill_colors Named list of fill colors for each parameter (optional)
#' @param title Grid title
#' @return grid.arrange object
#'
create_parameter_histogram_grid <- function(data, params, true_values, 
                                            binwidths = NULL, fill_colors = NULL, 
                                            title = "Parameter Estimate Distributions") {
  plots <- list()
  
  for (i in seq_along(params)) {
    param <- names(params)[i]
    param_name <- params[[i]]
    
    # Get binwidth and fill color if provided
    bw <- if (!is.null(binwidths) && param %in% names(binwidths)) binwidths[[param]] else NULL
    fc <- if (!is.null(fill_colors) && param %in% names(fill_colors)) fill_colors[[param]] else "skyblue"
    
    plots[[i]] <- create_parameter_histogram(
      data = data,
      param = param,
      param_name = param_name,
      true_value = true_values[[param]],
      binwidth = bw,
      fill_color = fc
    )
  }
  
  grid_plot <- do.call(grid.arrange, c(
    plots, 
    ncol = 2,
    top = textGrob(title, gp = gpar(fontsize = 14, fontface = "bold"))
  ))
  
  return(grid_plot)
}

#' Create a boxplot for parameter estimates
#'
#' @param data Data frame containing parameter estimates
#' @param params Vector of parameter names to include in boxplot
#' @param param_labels Named list of parameter labels
#' @param true_values Named list of true parameter values
#' @param title Plot title
#' @return ggplot object
#'
create_parameter_boxplot <- function(data, params, param_labels = NULL, 
                                     true_values = NULL, title = "Parameter Estimates") {
  # Convert data to long format
  data_long <- data %>%
    select(all_of(params)) %>%
    pivot_longer(cols = everything(), names_to = "Parameter", values_to = "Value")
  
  # Create parameter labels if not provided
  if (is.null(param_labels)) {
    param_labels <- params
    names(param_labels) <- params
  }
  
  # Create the boxplot
  plot <- ggplot(data_long, aes(x = Parameter, y = Value, fill = Parameter)) +
    geom_boxplot(alpha = 0.7) +
    labs(
      title = title,
      x = "",
      y = "Value"
    ) +
    theme_minimal() +
    scale_fill_brewer(palette = "Set3") +
    theme(legend.position = "none")
  
  # Add reference lines for true values if provided
  if (!is.null(true_values)) {
    true_values_df <- data.frame(
      Parameter = names(true_values),
      true_value = unlist(true_values)
    )
    
    plot <- plot + 
      geom_hline(data = true_values_df, 
                 aes(yintercept = true_value), 
                 linetype = "dashed", color = "red")
  }
  
  return(plot)
}

#' Create a comparison plot of parameter correlations for single vs dual epidemic fitting
#'
#' @param single_data Results from single epidemic fitting
#' @param dual_data Results from dual epidemic fitting
#' @param params Vector of parameter pairs to compare
#' @param title Plot title
#' @param subtitle Plot subtitle
#' @return ggplot object
#'
create_correlation_comparison <- function(single_data, dual_data, params, 
                                          title = "Correlation Comparison", 
                                          subtitle = NULL) {
  # Calculate correlations for single epidemic
  single_cors <- sapply(params, function(pair) {
    p1 <- pair[1]
    p2 <- pair[2]
    cor(single_data[[p1]], single_data[[p2]], use = "pairwise.complete.obs")
  })
  
  # Calculate correlations for dual epidemic
  dual_cors <- sapply(params, function(pair) {
    p1 <- pair[1]
    p2 <- pair[2]
    cor(dual_data[[p1]], dual_data[[p2]], use = "pairwise.complete.obs")
  })
  
  # Create data frame for plotting
  param_names <- sapply(params, function(pair) paste(pair[1], "-", pair[2]))
  comparison_df <- data.frame(
    Parameter_Pair = param_names,
    Single_Epidemic = single_cors,
    Dual_Epidemics = dual_cors
  )
  
  # Reshape to long format
  comparison_long <- comparison_df %>%
    pivot_longer(cols = c(Single_Epidemic, Dual_Epidemics), 
                 names_to = "Fitting_Type", values_to = "Correlation")
  
  # Create bar plot
  plot <- ggplot(comparison_long, 
                 aes(x = Parameter_Pair, y = abs(Correlation), 
                     fill = Fitting_Type)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    labs(
      title = title,
      subtitle = subtitle,
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
  
  return(plot)
}

#' Plot example epidemic trajectories
#'
#' @param data Data frame containing epidemic trajectories
#' @param group_var Name of column used to group different trajectories
#' @param time_var Name of column containing time values
#' @param value_var Name of column containing values to plot
#' @param vline_positions Vector of x positions for vertical lines
#' @param vline_colors Vector of colors for vertical lines
#' @param title Plot title
#' @param x_label X-axis label
#' @param y_label Y-axis label
#' @param caption Plot caption
#' @return ggplot object
#'
plot_trajectories <- function(data, group_var, time_var, value_var,
                              vline_positions = NULL, vline_colors = NULL,
                              title = "Epidemic Trajectories", 
                              x_label = "Time (days)", 
                              y_label = "Daily cases",
                              caption = NULL) {
  
  plot <- ggplot(data, aes_string(x = time_var, y = value_var, colour = group_var)) +
    geom_point() +
    labs(
      title = title,
      x = x_label,
      y = y_label,
      caption = caption
    ) +
    theme_minimal()
  
  # Add vertical lines if provided
  if (!is.null(vline_positions)) {
    if (is.null(vline_colors)) {
      vline_colors <- rep("black", length(vline_positions))
    }
    
    for (i in seq_along(vline_positions)) {
      plot <- plot + 
        geom_vline(xintercept = vline_positions[i], 
                   colour = vline_colors[i], 
                   linetype = "dashed")
    }
  }
  
  return(plot)
}