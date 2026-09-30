library(e1071)

# Set reproducible seed
STUDENT_ID_LAST4 <- 1234 
set.seed(STUDENT_ID_LAST4)

B <- 5000
n_range <- 10:30

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
    
    summary_table$sim_mean[i] <- mean(means)
    summary_table$sim_se[i]   <- sd(means)
    summary_table$theo_se[i]  <- if (!is.na(dist_info$pop_sd)) dist_info$pop_sd / sqrt(n) else NA
    summary_table$skewness[i] <- skewness(means)
    summary_table$kurtosis[i] <- kurtosis(means) 
    
    # Calculate LSRL fit for Q-Q plot to store R^2 quantitatively
    qq_vals <- qqnorm(means, plot.it = FALSE)
    qq_lsrl <- lm(qq_vals$y ~ qq_vals$x)
    summary_table$qq_r2[i] <- summary(qq_lsrl)$r.squared
  }
  
  return(list(name = dist_name, summary = summary_table, draws = sample_means_list))
}

# Run the simulation for ALL distributions
all_results <- list()
for (d in names(distributions)) {
  all_results[[d]] <- analyze_distribution(d, distributions[[d]])
}

# 3. Visual Functions

# A. Generate Individual Plots (Histogram + Q-Q with LSRL & Kurtosis) for n = 10..30
plot_all_n_diagnostics <- function(res_obj) {
  draws_list <- res_obj$draws
  dist_title <- res_obj$name
  
  for (n_val in names(draws_list)) {
    means <- draws_list[[n_val]]
    
    # Increase bottom margin (first value) to 5.5 to fit the Kurtosis label
    par(mfrow = c(1, 2), mar = c(5.5, 4.5, 4.0, 2.1), oma = c(0, 0, 2.5, 0))
    
    # Left: Histogram
    hist(means, probability = TRUE, col = "skyblue", border = "white",
         main = paste("Histogram (n =", n_val, ")"),
         xlab = "Sample Mean", ylab = "Density", cex.main = 1.1)
    lines(density(means), col = "darkblue", lwd = 2)
    
    # Right: Q-Q Plot
    qq_data <- qqnorm(means, main = paste("Q-Q Plot (n =", n_val, ")"),
                      xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", cex.main = 1.1)
    
    # Calculate and draw explicit LSRL on Q-Q Plot
    qq_lsrl <- lm(qq_data$y ~ qq_data$x)
    abline(qq_lsrl, col = "red", lwd = 2)
    
    # Add Labeled Kurtosis Value at the Bottom of the Q-Q Plot
    kurt_val <- kurtosis(means)
    mtext(paste("Excess Kurtosis =", round(kurt_val, 4)), 
          side = 1, line = 4.2, cex = 1.0, font = 2, col = "purple")
    
    # Outer Title
    mtext(paste0(dist_title, " — n = ", n_val), side = 3, outer = TRUE, line = 0.5, cex = 1.2, font = 2)
  }
}

# B. Generate Summary Metric Trends
plot_quantitative_summary <- function(res_obj) {
  df <- res_obj$summary
  dist_title <- res_obj$name
  par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3.5, 1.5), oma = c(0, 0, 2, 0))
  
  # Skewness Trend
  plot(df$n, df$skewness, type = "b", pch = 19, col = "darkgreen",
       xlab = "Sample Size (n)", ylab = "Skewness", main = "Skewness")
  abline(h = 0, col = "red", lty = 2)
  
  # Kurtosis Trend
  plot(df$n, df$kurtosis, type = "b", pch = 19, col = "purple",
       xlab = "Sample Size (n)", ylab = "Excess Kurtosis", main = "Kurtosis")
  abline(h = 0, col = "red", lty = 2)
  
  # SE Trend
  plot(df$n, df$sim_se, type = "b", pch = 19, col = "blue",
       xlab = "Sample Size (n)", ylab = "Standard Error", main = "SE Convergence")
  if (!all(is.na(df$theo_se))) lines(df$n, df$theo_se, col = "red", lty = 2)
  
  mtext(paste(dist_title, "Summary Trends (n = 10 to 30)"), side = 3, outer = TRUE, font = 2)
  par(mfrow = c(1, 1))
}

# ==============================================================================
# 4. EXECUTION LOOP: Render Output for EVERY Distribution
# ==============================================================================

# This loop ensures every distribution listed gets all individual n=10-30 graphs
# and its final summary metric dashboard printed to your R plot viewer.

for (d_name in names(all_results)) {
  # 1. Output Individual Plots (n=10 to 30) with LSRL & Kurtosis Label
  plot_all_n_diagnostics(all_results[[d_name]])
  
  # 2. Output the Summary Trend Dashboard
  plot_quantitative_summary(all_results[[d_name]])
}
# Open PDF device to record all plots into a single file
pdf("CLT_summary_Plots.pdf", width = 11, height = 8.5)

# Run execution loop to generate and save all plots across all distributions
for (d_name in names(all_results)) {
  # 1. Output Individual Diagnostic Plots (n = 10 to 30) with LSRL & Kurtosis Label
  plot_all_n_diagnostics(all_results[[d_name]])
  
  # 2. Output the Summary Trend Dashboard for the distribution
  plot_quantitative_summary(all_results[[d_name]])
}

# Close and finalize the PDF file
dev.off()