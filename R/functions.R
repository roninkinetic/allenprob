# Function to calculate probability for sample proportions
sample_proportion_prob <- function(population_proportion, sample_size, sample_proportion_target) {
  std_error <- sqrt(population_proportion * (1 - population_proportion) / sample_size)
  z_score <- (sample_proportion_target - population_proportion) / std_error
  pnorm(z_score)
}

# Function to calculate the confidence interval for a proportion
calculate_proportion_confidence_interval <- function(successes, sample_size, confidence_level) {
  # Calculate the sample proportion
  p_hat <- successes / sample_size

  # Calculate the standard error
  standard_error <- sqrt((p_hat * (1 - p_hat)) / sample_size)

  # Calculate the critical z-value for the given confidence level
  alpha <- 1 - confidence_level
  z_value <- qnorm(1 - alpha / 2)  # Two-tailed test

  # Calculate the margin of error
  margin_of_error <- z_value * standard_error

  # Calculate the confidence interval
  lower_bound <- p_hat - margin_of_error
  upper_bound <- p_hat + margin_of_error

  # Return the confidence interval as a list
  return(list(lower_bound = lower_bound, upper_bound = upper_bound))
}

# Function to calculate sample mean probability
calculate_sample_mean_probability <- function(population_mean, population_std_dev, sample_size, sample_mean_target) {
  # Calculate standard error of the mean
  std_error <- population_std_dev / sqrt(sample_size)

  # Calculate the z-score
  z_score <- (sample_mean_target - population_mean) / std_error

  # Calculate and return the probability
  probability <- 1 - pnorm(z_score)
  return(probability)
}


# Expo probability between times
calculate_exponential_probability <- function(lambda, time_start, time_end) {
  # Calculate the cumulative probability up to time_end and time_start
  P_T_leq_end <- pexp(time_end, rate = lambda)
  P_T_leq_start <- pexp(time_start, rate = lambda)

  # Calculate the probability between time_start and time_end
  probability <- P_T_leq_end - P_T_leq_start
  return(probability)
}


# Pnorm probability between two bounds
prob_between <- function(lower_bound, upper_bound, mean, sd) {
  # Calculate the probability of being between the bounds
  pnorm(upper_bound, mean = mean, sd = sd) - pnorm(lower_bound, mean = mean, sd = sd)
}

# Calculate confidence intervals with sample mean, sample sd, sample size, and confidence level
calculate_confidence_interval <- function(sample_mean, sample_sd, sample_size, confidence_level) {
  # Calculate the critical z-value for the specified confidence level
  alpha <- 1 - confidence_level  # Calculate the alpha level
  z_value <- qnorm(1 - alpha / 2)  # Use the upper tail for two-tailed test

  # Calculate the standard error
  standard_error <- sample_sd / sqrt(sample_size)

  # Calculate the margin of error
  margin_of_error <- z_value * standard_error

  # Calculate the confidence interval
  lower_bound <- sample_mean - margin_of_error
  upper_bound <- sample_mean + margin_of_error

  # Return the confidence interval as a list
  return(list(lower_bound = lower_bound, upper_bound = upper_bound))
}

# Function to calculate the probability of sample mean being greater than or equal to a given value
calculate_probability_mean_greater_than <- function(population_mean, population_sd, sample_size, sample_mean) {
  # Calculate the standard error of the mean
  standard_error <- population_sd / sqrt(sample_size)

  # Calculate the z-score for the sample mean
  z_score <- (sample_mean - population_mean) / standard_error

  # Calculate the probability using the cumulative distribution function
  probability <- 1 - pnorm(z_score)  # P(X >= sample_mean)

  return(probability)
}

# Function to calculate the probability of sample mean being less than or equal to a given value
calculate_probability_mean_less_than_or_equal <- function(population_mean, population_sd, sample_size, sample_mean) {
  # Calculate the standard error of the mean
  standard_error <- population_sd / sqrt(sample_size)

  # Calculate the z-score for the sample mean
  z_score <- (sample_mean - population_mean) / standard_error

  # Calculate the probability using the cumulative distribution function
  probability <- pnorm(z_score)  # P(X <= sample_mean)

  return(probability)
}

# Function to calculate the probability of a cost being less than a given value
calculate_probability_cost_less_than <- function(mean_cost, sd_cost, threshold_cost) {
  # Calculate the z-score for the given cost threshold
  z_score <- (threshold_cost - mean_cost) / sd_cost

  # Calculate the probability using the cumulative distribution function
  probability <- pnorm(z_score)  # P(X < threshold_cost)

  return(probability)
}

# Function to calculate the probability of a cost being greater than a given value
calculate_probability_cost_greater_than <- function(mean_cost, sd_cost, threshold_cost) {
  # Calculate the z-score for the given cost threshold
  z_score <- (threshold_cost - mean_cost) / sd_cost

  # Calculate the probability using the cumulative distribution function
  probability <- 1 - pnorm(z_score)  # P(X > threshold_cost)

  return(probability)
}

# Generic function to perform a one-sample t-test and calculate the p-value
calculate_p_value_greater_than <- function(data, population_mean, significance_level) {
  # Perform the one-sample t-test
  t_test_result <- t.test(data, mu = population_mean, alternative = "greater")

  # Extract the p-value
  p_value <- t_test_result$p.value

  # Interpret the results based on the significance level
  interpretation <- ifelse(p_value < significance_level,
                           "The sample mean is significantly greater than the population mean.",
                           "The sample mean is not significantly greater than the population mean.")

  return(list(p_value = p_value, interpretation = interpretation))
}


# Function to perform a two-tailed hypothesis test for a proportion
hypothesis_test_proportion <- function(successes, trials, significance_level) {
  # Calculate the sample proportion
  p_hat <- successes / trials

  # Define the null hypothesis proportion
  p_0 <- 0.5

  # Calculate the standard error
  standard_error <- sqrt((p_0 * (1 - p_0)) / trials)

  # Calculate the z-test statistic
  z <- (p_hat - p_0) / standard_error

  # Calculate the p-value for the two-tailed test
  p_value <- 2 * (1 - pnorm(abs(z)))  # two-tailed

  # Determine the critical z-value for the given significance level
  alpha <- significance_level
  z_critical <- qnorm(1 - alpha / 2)

  # Decision rule
  if (abs(z) > z_critical) {
    result <- "Reject the null hypothesis."
  } else {
    result <- "Fail to reject the null hypothesis."
  }

  # Return a list with results
  return(list(z = z, p_value = p_value, result = result))
}


# Function to calculate the confidence interval for a proportion
calculate_proportion_confidence_interval <- function(successes, sample_size, confidence_level) {
  # Calculate the sample proportion
  p_hat <- successes / sample_size

  # Calculate the standard error
  standard_error <- sqrt((p_hat * (1 - p_hat)) / sample_size)

  # Calculate the critical z-value for the given confidence level
  alpha <- 1 - confidence_level
  z_value <- qnorm(1 - alpha / 2)  # Two-tailed test

  # Calculate the margin of error
  margin_of_error <- z_value * standard_error

  # Calculate the confidence interval
  lower_bound <- p_hat - margin_of_error
  upper_bound <- p_hat + margin_of_error

  # Return the confidence interval as a list
  return(list(lower_bound = lower_bound, upper_bound = upper_bound))
}

# Function to find the value of x such that P(X > x) = probability
find_x_value <- function(mean, sd, probability) {
  # Calculate the z-score for the given probability
  z_value <- qnorm(1 - probability)  # Use 1 - probability for the upper tail

  # Calculate the value of x
  x_value <- mean + z_value * sd

  return(x_value)
}
