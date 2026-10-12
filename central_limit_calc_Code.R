library(e1071)

# Set reproducible seed
STUDENT_ID_LAST4 <- 2747 
set.seed(STUDENT_ID_LAST4)

B <- 20000
n_range <- 10:40

# 1. Define ALL distributions
distributions <- list(
  "Standard Normal" = list(
    r_fun = function(n) rnorm(n, mean = 0, sd = 1),
    pop_sd = 1
  ),
  "Six-Sided Die" = list(
    r_fun = function(n) sample(1:6, size = n, replace = TRUE),
    pop_sd = sqrt(35/12)
  ),
  "Exponential (lambda=1)" = list(
    r_fun = function(n) rexp(n, rate = 1),
    pop_sd = 1
  ),
  "Binomial Small (n=5, p=0.10)" = list(
    r_fun = function(n) rbinom(n, size = 5, prob = 0.10),
    pop_sd = sqrt(0.45)
  ),
  "Cauchy (Location=0, Scale=1)" = list(
    r_fun = function(n) rcauchy(n, location = 0, scale = 1),
    pop_sd = NA # Infinite/undefined variance
  )
)

# 2. Simulation Engine (Values n=10 to 30)
analyze_distribution <- function(dist_name, dist_info) {
  cat("Processing:", dist_name, "...\n")
  
  summary_table <- data.frame(
    n = n_range,
    sim_mean = numeric(length(n_range)),
    sim_se = numeric(length(n_range)),
    theo_se = numeric(length(n_range)),
    skewness = numeric(length(n_range)),
    kurtosis = numeric(length(n_range)),
    qq_r2 = numeric(length(n_range)) 
  )
  
  sample_means_list <- list()
  
  for (i in seq_along(n_range)) {
    n <- n_range[i]
    means <- replicate(B, mean(dist_info$r_fun(n)))
    
    sample_means_list[[as.character(n)]] <- means
    
    # Calculate metrics, removing extreme values for safety
    summary_table$sim_mean[i] <- mean(means, na.rm = TRUE)
    summary_table$sim_se[i]   <- sd(means, na.rm = TRUE)
    summary_table$theo_se[i]  <- if (!is.na(dist_info$pop_sd)) dist_info$pop_sd / sqrt(n) else NA
    summary_table$skewness[i] <- skewness(means, na.rm = TRUE)
    summary_table$kurtosis[i] <- kurtosis(means, na.rm = TRUE) 
    
    # Calculate LSRL fit for Q-Q plot to store R^2 quantitatively
    qq_vals <- qqnorm(means, plot.it = FALSE)
    valid_qq <- is.finite(qq_vals$x) & is.finite(qq_vals$y)
    if (sum(valid_qq) > 1) {
      qq_lsrl <- lm(qq_vals$y[valid_qq] ~ qq_vals$x[valid_qq])
      summary_table$qq_r2[i] <- summary(qq_lsrl)$r.squared
    } else {
      summary_table$qq_r2[i] <- NA
    }
  }
  
  return(list(name = dist_name, summary = summary_table, draws = sample_means_list))
}

# Run the simulation for ALL distributions
all_results <- list()
for (d in names(distributions)) {
  all_results[[d]] <- analyze_distribution(d, distributions[[d]])
}

# 3. Visual Functions

# A. Generate Individual Plots for n = 10..30 (Safeguarded against non-finite values)
plot_all_n_diagnostics <- function(res_obj) {
  draws_list <- res_obj$draws
  dist_title <- res_obj$name
  
  for (n_val in names(draws_list)) {
    means <- draws_list[[n_val]]
    means_clean <- means[is.finite(means)] # Remove extreme/Inf values for safe plotting
    
    if (length(means_clean) < 2) next
    
    par(mfrow = c(1, 2), mar = c(5.5, 4.5, 4.0, 2.1), oma = c(0, 0, 2.5, 0))
    
    # Compute safe scale metrics (used for plot AND text labels)
    m_val <- mean(means_clean)
    s_val <- sd(means_clean)
    
    h_temp <- hist(means_clean, plot = FALSE)
    dens_y <- density(means_clean)$y
    norm_y <- if (!is.na(s_val) && s_val > 0) dnorm(m_val, mean = m_val, sd = s_val) else 0
    max_y  <- max(c(h_temp$density, dens_y, norm_y), na.rm = TRUE)
    if (!is.finite(max_y) || max_y <= 0) max_y <- 1
    
    # Left: Histogram with Safe Y-Limits
    hist(means_clean, probability = TRUE, col = "skyblue", border = "white",
         ylim = c(0, max_y * 1.15),
         main = paste("Histogram (n =", n_val, ")"),
         xlab = "Sample Mean", ylab = "Density", cex.main = 1.1)
    
    lines(density(means_clean), col = "darkblue", lwd = 2)
    
    if (!is.na(s_val) && s_val > 0) {
      curve(dnorm(x, mean = m_val, sd = s_val), add = TRUE, col = "red", lwd = 2.5, lty = 2)
    }
    
    legend("topright", legend = c("Sample KDE", "Normal Curve"),
           col = c("darkblue", "red"), lty = c(1, 2), lwd = 2, bty = "n", cex = 0.8)
    
    # --> NEW ADDITION: Add Labeled Simulated Mean & SD at the bottom of the Histogram <--
    mtext(paste("Sim. Mean =", round(m_val, 4), " | Sim. SD =", round(s_val, 4)), 
          side = 1, line = 4.2, cex = 1.0, font = 2, col = "darkgreen")
    
    # Right: Q-Q Plot
    qq_data <- qqnorm(means_clean, main = paste("Q-Q Plot (n =", n_val, ")"),
                      xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", cex.main = 1.1)
    
    valid_qq <- is.finite(qq_data$x) & is.finite(qq_data$y)
    if (sum(valid_qq) > 1) {
      qq_lsrl <- lm(qq_data$y[valid_qq] ~ qq_data$x[valid_qq])
      abline(qq_lsrl, col = "red", lwd = 2)
    }
    
    # Kurtosis label at the bottom of the Q-Q plot
    kurt_val <- kurtosis(means_clean)
    mtext(paste("Excess Kurtosis =", round(kurt_val, 4)), 
          side = 1, line = 4.2, cex = 1.0, font = 2, col = "purple")
    
    # Master Header for both plots
    mtext(paste0(dist_title, " — n = ", n_val), side = 3, outer = TRUE, line = 0.5, cex = 1.2, font = 2)
  }
}

# B. Generate Summary Metric Trends (Safeguarded)
plot_quantitative_summary <- function(res_obj) {
  df <- res_obj$summary
  dist_title <- res_obj$name
  
  par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3.5, 1.5), oma = c(0, 0, 2, 0))
  
  safe_trend_plot <- function(x, y, col, ylab, main_title) {
    valid <- is.finite(y)
    if (any(valid)) {
      plot(x[valid], y[valid], type = "b", pch = 19, col = col,
           xlab = "Sample Size (n)", ylab = ylab, main = main_title)
      abline(h = 0, col = "red", lty = 2)
    } else {
      plot(x, rep(0, length(x)), type = "n", xlab = "Sample Size (n)", ylab = ylab, main = main_title)
      mtext("Undefined / Non-finite", side = 3, line = -2, col = "red", cex = 0.8)
    }
  }
  
  safe_trend_plot(df$n, df$skewness, "darkgreen", "Skewness", "Skewness")
  safe_trend_plot(df$n, df$kurtosis, "purple", "Excess Kurtosis", "Kurtosis")
  
  # Standard Error Plot
  valid_se <- is.finite(df$sim_se)
  if (any(valid_se)) {
    plot(df$n[valid_se], df$sim_se[valid_se], type = "b", pch = 19, col = "blue",
         xlab = "Sample Size (n)", ylab = "Standard Error", main = "SE Convergence")
    if (!all(is.na(df$theo_se))) lines(df$n, df$theo_se, col = "red", lty = 2)
  } else {
    plot(df$n, rep(0, length(df$n)), type = "n", xlab = "Sample Size (n)", ylab = "Standard Error", main = "SE Convergence")
    mtext("Undefined (Infinite Variance)", side = 3, line = -2, col = "red", cex = 0.8)
  }
  
  mtext(paste(dist_title, "Summary Trends (n = 10 to 40)"), side = 3, outer = TRUE, font = 2)
  par(mfrow = c(1, 1))
}

# ==============================================================================
# 4. EXECUTION LOOP & PDF OUTPUT
# ==============================================================================

# Reset active plot devices
graphics.off()

# Save all diagnostic plots and dashboards into a clean PDF file
pdf("CLT_summary_Plots.pdf", width = 11, height = 8.5)

for (d_name in names(all_results)) {
  plot_all_n_diagnostics(all_results[[d_name]])
  plot_quantitative_summary(all_results[[d_name]])
}

dev.off()
cat("Process complete! CLT_summary_Plots.pdf generated successfully.\n")
# ==============================================================================
# 5. CONSOLE OUTPUT: Simulated vs. Theoretical SE Differences
# ==============================================================================

cat("\n=== Difference Between Simulated and Theoretical SE ===\n")

for (d_name in names(all_results)) {
  # Skip Cauchy (since theoretical SD is undefined)
  if (grepl("Cauchy", d_name)) next
  
  # Extract the summary dataframe for this distribution
  df <- all_results[[d_name]]$summary
  
  cat("\nDistribution:", d_name, "\n")
  
  # Calculate difference for each n
  diff_vec <- df$sim_se - df$theo_se
  
  out_df <- data.frame(
    n = df$n,
    Simulated_SE   = round(df$sim_se, 5),
    Theoretical_SE = round(df$theo_se, 5),
    Difference     = round(diff_vec, 5)
  )
  
  # Print the main comparison table
  print(out_df, row.names = FALSE)
  
  # Print the mean difference summary line
  mean_diff <- mean(diff_vec, na.rm = TRUE)
  cat("-------------------------------------------------------\n")
  cat("Mean Difference across n = 10 to 30:", round(mean_diff, 5), "\n")
  cat("-------------------------------------------------------\n")
}
cat("\n=======================================================\n")