#' Calculate Probability for Sample Proportions
#'
#' This function calculates the probability of observing a sample proportion
#' greater than or equal to a specified target sample proportion based on
#' a given population proportion and sample size.
#'
#' The calculation is based on the standard normal distribution, and the
#' steps involved are as follows:
#'
#' 1. Calculate the standard error of the sample proportion using the formula:
#'    SE = sqrt((p * (1 - p)) / n)
#'    where p is the population proportion and n is the sample size.
#'
#' 2. Compute the z-score using the formula:
#'    z = (p_hat - p) / SE
#'    where p_hat is the target sample proportion.
#'
#' 3. Use the z-score to calculate the cumulative probability using the standard normal distribution:
#'    P(Z <= z) = pnorm(z)
#'
#' @param population_proportion The proportion in the population (e.g., 0.5).
#' @param sample_size The size of the sample (an integer).
#' @param sample_proportion_target The target sample proportion.
#'
#' @return The probability as a numeric value representing the cumulative probability of the
#'         sample proportion being less than or equal to the target.
#'
#' @examples
#' sample_proportion_prob(0.5, 100, 0.55)
#'
#' @export
sample_proportion_prob <- function(population_proportion, sample_size, sample_proportion_target) {
  # Calculate the standard error
  std_error <- sqrt(population_proportion * (1 - population_proportion) / sample_size)

  # Calculate the z-score
  z_score <- (sample_proportion_target - population_proportion) / std_error

  # Calculate and return the cumulative probability
  pnorm(z_score)
}


#' Calculate Confidence Interval for a Proportion
#'
#' This function calculates the confidence interval for a population proportion
#' based on the number of successes observed in a sample, the sample size,
#' and the desired confidence level.
#'
#' The calculation involves the following steps:
#'
#' 1. Compute the sample proportion (p_hat) using the formula:
#'    p_hat = successes / sample_size
#'
#' 2. Calculate the standard error (SE) of the sample proportion using the formula:
#'    SE = sqrt((p_hat * (1 - p_hat)) / n)
#'    where n is the sample size.
#'
#' 3. Determine the critical z-value for the specified confidence level:
#'    z = qnorm(1 - (alpha / 2))
#'    where alpha is the significance level (e.g., alpha = 1 - confidence_level).
#'
#' 4. Calculate the margin of error using the formula:
#'    ME = z * SE
#'
#' 5. The confidence interval is then given by:
#'    CI = (p_hat - ME, p_hat + ME)
#'
#' @param successes The number of successes in the sample (an integer).
#' @param sample_size The size of the sample (an integer).
#' @param confidence_level The confidence level (e.g., 0.95 for 95% confidence).
#'
#' @return A list containing the lower and upper bounds of the confidence interval:
#'   - `lower_bound`: The lower bound of the confidence interval.
#'   - `upper_bound`: The upper bound of the confidence interval.
#'
#' @examples
#' calculate_proportion_confidence_interval(50, 100, 0.95)
#'
#' @export
calculate_proportion_confidence_interval <- function(successes, sample_size, confidence_level) {
  # Calculate the sample proportion
  p_hat <- successes / sample_size

  # Calculate the standard error
  standard_error <- sqrt((p_hat * (1 - p_hat)) / sample_size)

  # Calculate the significance level
  alpha <- 1 - confidence_level

  # Calculate the critical z-value for the given confidence level
  z_value <- qnorm(1 - alpha / 2)

  # Calculate the margin of error
  margin_of_error <- z_value * standard_error

  # Return the confidence interval as a list
  list(lower_bound = p_hat - margin_of_error, upper_bound = p_hat + margin_of_error)
}

#' Calculate Sample Mean Probability
#'
#' This function calculates the probability of obtaining a sample mean that is at least
#' as extreme as a specified target value, given a population mean, population standard
#' deviation, and sample size. This is useful for hypothesis testing and determining
#' how likely a sample mean is under the null hypothesis.
#'
#' The calculation involves the following steps:
#'
#' 1. Calculate the standard error (SE) of the sample mean using the formula:
#'    SE = sigma / sqrt(n)
#'    where sigma is the population standard deviation and n is the sample size.
#'
#' 2. Compute the z-score for the target sample mean using the formula:
#'    z = (x_bar - mu) / SE
#'    where x_bar is the target sample mean and mu is the population mean.
#'
#' 3. Calculate the probability of obtaining a sample mean at least as extreme as the target
#'    value using the cumulative distribution function of the standard normal distribution:
#'    P(X >= x_bar) = 1 - P(Z <= z) = 1 - pnorm(z)
#'
#' @param population_mean The mean of the population (numeric).
#' @param population_std_dev The standard deviation of the population (numeric).
#' @param sample_size The size of the sample (integer).
#' @param sample_mean_target The target sample mean (numeric).
#'
#' @return The probability as a numeric value representing the likelihood of obtaining a sample
#'         mean at least as extreme as the target.
#'
#' @examples
#' calculate_sample_mean_probability(100, 15, 30, 105)
#'
#' @export
calculate_sample_mean_probability <- function(population_mean, population_std_dev, sample_size, sample_mean_target) {
  # Calculate the standard error
  std_error <- population_std_dev / sqrt(sample_size)

  # Calculate the z-score
  z_score <- (sample_mean_target - population_mean) / std_error

  # Calculate and return the probability of obtaining a sample mean at least as extreme as the target
  1 - pnorm(z_score)
}

#' Calculate Probability Between Two Times for an Exponential Distribution
#'
#' This function calculates the probability of an event occurring between two specified times
#' in an exponential distribution characterized by a rate parameter. The exponential
#' distribution is commonly used to model the time until an event occurs.
#'
#' The calculation involves the following steps:
#'
#' 1. Compute the cumulative distribution function (CDF) at the end time using:
#'    P(X <= t) = pexp(t, rate = lambda)
#'    This gives the probability that the event occurs by time t.
#'
#' 2. Compute the CDF at the start time using:
#'    P(X <= t0) = pexp(t0, rate = lambda)
#'    This gives the probability that the event occurs by time t0.
#'
#' 3. The probability of the event occurring between the start and end times is given by:
#'    P(t0 < X <= t) = P(X <= t) - P(X <= t0)
#'    This represents the difference between the two cumulative probabilities.
#'
#' @param lambda The rate parameter of the exponential distribution (numeric).
#' @param time_start The start time (numeric).
#' @param time_end The end time (numeric).
#'
#' @return The probability of the event occurring between `time_start` and `time_end` as a numeric value.
#'
#' @examples
#' calculate_exponential_probability(0.1, 2, 5)
#'
#' @export
calculate_exponential_probability <- function(lambda, time_start, time_end) {
  # Calculate the probability of the event occurring between the two times
  pexp(time_end, rate = lambda) - pexp(time_start, rate = lambda)
}

#' Calculate Probability Between Two Bounds
#'
#' This function calculates the probability of a normally distributed variable falling
#' between two specified bounds. It uses the properties of the normal distribution to compute
#' this probability.
#'
#' The calculation involves the following steps:
#'
#' 1. Compute the cumulative distribution function (CDF) at the upper bound:
#'    P(X <= b) = pnorm(b, mean, sd)
#'    This gives the probability that the variable is less than or equal to the upper bound b.
#'
#' 2. Compute the CDF at the lower bound:
#'    P(X <= a) = pnorm(a, mean, sd)
#'    This gives the probability that the variable is less than or equal to the lower bound a.
#'
#' 3. The probability of the variable being between the lower and upper bounds is given by:
#'    P(a < X <= b) = P(X <= b) - P(X <= a)
#'    This represents the difference between the two cumulative probabilities.
#'
#' @param lower_bound The lower bound (numeric).
#' @param upper_bound The upper bound (numeric).
#' @param mean The mean of the distribution (numeric).
#' @param sd The standard deviation of the distribution (numeric).
#'
#' @return The probability of the variable being between `lower_bound` and `upper_bound` as a numeric value.
#'
#' @examples
#' prob_between(10, 20, 15, 5)
#'
#' @export
prob_between <- function(lower_bound, upper_bound, mean, sd) {
  # Calculate the probability of the variable being between the two bounds
  pnorm(upper_bound, mean = mean, sd = sd) - pnorm(lower_bound, mean = mean, sd = sd)
}

#' Calculate Confidence Interval for Sample Mean
#'
#' This function calculates the confidence interval for a sample mean based on the sample mean,
#' sample standard deviation, sample size, and desired confidence level. The confidence interval
#' provides a range of values within which the true population mean is expected to lie with a
#' specified level of confidence.
#'
#' The calculation involves the following steps:
#'
#' 1. Calculate the significance level (alpha) from the confidence level:
#'    alpha = 1 - confidence_level
#'
#' 2. Determine the critical z-value for the specified confidence level:
#'    z = qnorm(1 - (alpha / 2)
#'    This critical value corresponds to the desired level of confidence and accounts for the
#'    tails of the normal distribution.
#'
#' 3. Calculate the standard error (SE) of the sample mean using the formula:
#'    SE = s / sqrt(n)
#'    where s is the sample standard deviation and n is the sample size.
#'
#' 4. Calculate the margin of error using the formula:
#'    ME = z * SE
#'
#' 5. The confidence interval is then given by:
#'    CI = (sample_mean - ME, sample_mean + ME)
#'
#' @param sample_mean The mean of the sample (numeric).
#' @param sample_sd The standard deviation of the sample (numeric).
#' @param sample_size The size of the sample (integer).
#' @param confidence_level The confidence level (numeric, e.g., 0.95 for 95% confidence).
#'
#' @return A list containing the lower and upper bounds of the confidence interval:
#'   - `lower_bound`: The lower bound of the confidence interval.
#'   - `upper_bound`: The upper bound of the confidence interval.
#'
#' @examples
#' calculate_confidence_interval(50, 10, 30, 0.95)
#'
#' @export
calculate_confidence_interval <- function(sample_mean, sample_sd, sample_size, confidence_level) {
  # Calculate the significance level
  alpha <- 1 - confidence_level

  # Calculate the critical z-value for the given confidence level
  z_value <- qnorm(1 - alpha / 2)

  # Calculate the standard error of the sample mean
  standard_error <- sample_sd / sqrt(sample_size)

  # Calculate the margin of error
  margin_of_error <- z_value * standard_error

  # Return the confidence interval as a list
  list(lower_bound = sample_mean - margin_of_error, upper_bound = sample_mean + margin_of_error)
}


#' Calculate Probability of Sample Mean Greater Than or Equal to Target
#'
#' This function calculates the probability of the sample mean being greater than or equal
#' to a specified target value based on the population mean, population standard deviation,
#' and sample size. This is useful for hypothesis testing and assessing how likely it is
#' to observe a sample mean at or above a certain threshold under the null hypothesis.
#'
#' The calculation involves the following steps:
#'
#' 1. Calculate the standard error (SE) of the sample mean using the formula:
#'    SE = sigma / sqrt(n)
#'    where sigma is the population standard deviation and n is the sample size.
#'
#' 2. Compute the z-score for the target sample mean using the formula:
#'    z = (x_bar - mu) / SE
#'    where x_bar is the target sample mean and mu is the population mean.
#'
#' 3. Calculate the probability of obtaining a sample mean greater than or equal to the target
#'    value using the complementary cumulative distribution function:
#'    P(X >= x_bar) = 1 - P(Z <= z) = 1 - pnorm(z)
#'
#' @param population_mean The mean of the population (numeric).
#' @param population_sd The standard deviation of the population (numeric).
#' @param sample_size The size of the sample (integer).
#' @param sample_mean The target sample mean (numeric).
#'
#' @return The probability as a numeric value representing the likelihood of the sample mean
#'         being greater than or equal to the target.
#'
#' @examples
#' calculate_probability_mean_greater_than(100, 15, 30, 105)
#'
#' @export
calculate_probability_mean_greater_than <- function(population_mean, population_sd, sample_size, sample_mean) {
  # Calculate the standard error of the sample mean
  standard_error <- population_sd / sqrt(sample_size)

  # Calculate the z-score for the target sample mean
  z_score <- (sample_mean - population_mean) / standard_error

  # Calculate and return the probability of the sample mean being greater than or equal to the target
  1 - pnorm(z_score)
}

#' Calculate Probability of Sample Mean Less Than or Equal to Target
#'
#' This function calculates the probability of the sample mean being less than or equal
#' to a specified target value based on the population mean, population standard deviation,
#' and sample size. This is useful for hypothesis testing and assessing how likely it is
#' to observe a sample mean at or below a certain threshold under the null hypothesis.
#'
#' The calculation involves the following steps:
#'
#' 1. Calculate the standard error (SE) of the sample mean using the formula:
#'    SE = sigma / sqrt(n)
#'    where sigma is the population standard deviation and n is the sample size.
#'
#' 2. Compute the z-score for the target sample mean using the formula:
#'    z = (x_bar - mu) / SE
#'    where x_bar is the target sample mean and mu is the population mean.
#'
#' 3. Calculate the probability of obtaining a sample mean less than or equal to the target
#'    value using the cumulative distribution function:
#'    P(X <= x_bar) = P(Z <= z) = pnorm(z)
#'
#' @param population_mean The mean of the population (numeric).
#' @param population_sd The standard deviation of the population (numeric).
#' @param sample_size The size of the sample (integer).
#' @param sample_mean The target sample mean (numeric).
#'
#' @return The probability as a numeric value representing the likelihood of the sample mean
#'         being less than or equal to the target.
#'
#' @examples
#' calculate_probability_mean_less_than_or_equal(100, 15, 30, 95)
#'
#' @export
calculate_probability_mean_less_than_or_equal <- function(population_mean, population_sd, sample_size, sample_mean) {
  # Calculate the standard error of the sample mean
  standard_error <- population_sd / sqrt(sample_size)

  # Calculate the z-score for the target sample mean
  z_score <- (sample_mean - population_mean) / standard_error

  # Calculate and return the probability of the sample mean being less than or equal to the target
  pnorm(z_score)
}


#' Find X Value for a Given Probability
#'
#' This function finds the value of `x` such that the probability of a normally
#' distributed random variable being greater than `x` is equal to the specified
#' probability. This is useful for determining threshold values in statistical analyses.
#'
#' The calculation involves the following steps:
#'
#' 1. Compute the z-score corresponding to the specified probability using the quantile function
#'    for the standard normal distribution:
#'    z = qnorm(1 - probability)
#'    This gives the z-value such that the area to the right of z in the standard normal
#'    distribution equals the specified probability.
#'
#' 2. Convert the z-score to the corresponding x value using the formula:
#'    x = mu + z * sigma
#'    where mu is the mean and sigma is the standard deviation of the distribution.
#'
#' @param mean The mean of the distribution (numeric).
#' @param sd The standard deviation of the distribution (numeric).
#' @param probability The target probability (numeric, must be between 0 and 1).
#'
#' @return The `x` value corresponding to the specified probability (numeric).
#'
#' @examples
#' find_x_value(100, 15, 0.05)
#'
#' @export
find_x_value <- function(mean, sd, probability) {
  # Calculate the z-score corresponding to the target probability
  z_value <- qnorm(1 - probability)

  # Calculate and return the corresponding x value
  mean + z_value * sd
}

#' Perform a Two-Tailed Hypothesis Test for a Proportion
#'
#' This function performs a two-tailed hypothesis test for a proportion to determine
#' if the sample proportion significantly differs from the null hypothesis proportion.
#' The default null hypothesis proportion is set to 0.5.
#'
#' The calculation involves the following steps:
#'
#' Perform a Two-Tailed Hypothesis Test for a Proportion
#'
#' This function performs a two-tailed hypothesis test for a proportion to determine
#' if the sample proportion significantly differs from the null hypothesis proportion.
#' The default null hypothesis proportion is set to 0.5.
#'
#' The calculation involves the following steps:
#'
#' 1. Calculate the sample proportion (p_hat):
#'    p_hat = successes / trials
#'
#' 2. Set the null hypothesis proportion (p_0):
#'    p_0 = 0.5
#'
#' 3. Calculate the standard error (SE) for the null hypothesis proportion using:
#'    SE = sqrt((p_0 * (1 - p_0)) / trials)
#'
#' 4. Compute the z-test statistic:
#'    z = (p_hat - p_0) / SE
#'
#' 5. Calculate the p-value for the two-tailed test:
#'    p-value = 2 * (1 - P(Z <= |z|)) = 2 * (1 - pnorm(|z|))
#'
#' 6. Determine the critical z-value for the specified significance level:
#'    z_critical = qnorm(1 - (alpha / 2))
#'
#' 7. Make a decision based on the comparison of the absolute z-test statistic with the critical z-value:
#'    - If |z| > z_critical, reject the null hypothesis.
#'    - Otherwise, fail to reject the null hypothesis.
#'
#' @param successes The number of successes in the sample (integer).
#' @param trials The number of trials in the sample (sample size, integer).
#' @param significance_level The significance level for the hypothesis test (numeric, e.g., 0.05 for a 5% significance level).
#'
#' @return A list containing:
#'   - `z`: The z-test statistic (numeric).
#'   - `p_value`: The calculated p-value (numeric).
#'   - `result`: A string indicating whether to "Reject the null hypothesis" or "Fail to reject the null hypothesis."
#'
#' @examples
#' hypothesis_test_proportion(55, 100, 0.05)
#'
#' @export
hypothesis_test_proportion <- function(successes, trials, significance_level) {
  # Calculate the sample proportion
  p_hat <- successes / trials

  # Null hypothesis proportion
  p_0 <- 0.5

  # Calculate the standard error
  standard_error <- sqrt((p_0 * (1 - p_0)) / trials)

  # Calculate the z-test statistic
  z <- (p_hat - p_0) / standard_error

  # Calculate the p-value for the two-tailed test
  p_value <- 2 * (1 - pnorm(abs(z)))

  # Set significance level
  alpha <- significance_level

  # Calculate the critical z-value
  z_critical <- qnorm(1 - alpha / 2)

  # Make decision based on z-test statistic and critical z-value
  result <- if (abs(z) > z_critical) {
    "Reject the null hypothesis."
  } else {
    "Fail to reject the null hypothesis."
  }

  # Return the results as a list
  list(z = z, p_value = p_value, result = result)
}


#' Perform a One-Sample T-Test and Calculate the P-Value
#'
#' This function performs a one-sample t-test to determine if the sample mean is significantly
#' greater than a specified population mean. This test is useful for assessing whether
#' the observed sample mean provides sufficient evidence to conclude that it exceeds the
#' hypothesized population mean.
#'
#' The calculation involves the following steps:
#'
#' 1. Calculate the sample mean (x_bar) and sample standard deviation (s) from the input data.
#'
#' 2. Compute the t-statistic using the formula:
#'    t = (x_bar - mu) / (s / sqrt(n))
#'    where mu is the population mean and n is the sample size.
#'
#' 3. Use the t-distribution to calculate the p-value for the test under the alternative hypothesis that
#'    the sample mean is greater than the population mean.
#'
#' 4. Interpret the result based on the specified significance level:
#'    - If the p-value is less than the significance level, conclude that the sample mean
#'      is significantly greater than the population mean.
#'    - Otherwise, conclude that it is not significantly greater.
#'
#' @param data A numeric vector of sample data (numeric).
#' @param population_mean The population mean to compare the sample mean against (numeric).
#' @param significance_level The significance level for the test (numeric, e.g., 0.05 for a 5% significance level).
#'
#' @return A list containing:
#'   - `p_value`: The calculated p-value (numeric).
#'   - `interpretation`: A string interpreting the result based on the significance level.
#'
#' @examples
#' calculate_p_value_greater_than(c(5.1, 5.5, 5.9, 6.2, 5.7), 5.0, 0.05)
#'
#' @export
calculate_p_value_greater_than <- function(data, population_mean, significance_level) {
  # Perform the one-sample t-test
  t_test_result <- t.test(data, mu = population_mean, alternative = "greater")

  # Extract the p-value from the test result
  p_value <- t_test_result$p.value

  # Interpret the result based on the significance level
  interpretation <- ifelse(
    p_value < significance_level,
    "The sample mean is significantly greater than the population mean.",
    "The sample mean is not significantly greater than the population mean."
  )

  # Return the p-value and interpretation as a list
  list(p_value = p_value, interpretation = interpretation)
}

#' Calculate Probability of Cost Being Less Than a Given Value
#'
#' This function calculates the probability that a cost is less than a specified threshold,
#' based on the mean and standard deviation of the cost distribution. This is useful for
#' understanding the likelihood of costs falling below a certain level in a normally
#' distributed context.
#'
#' The calculation involves the following steps:
#'
#' 1. Calculate the z-score corresponding to the threshold cost using the formula:
#'    z = (X - mu) / sigma
#'    where X is the threshold cost, mu is the mean cost, and sigma is the standard deviation of the cost.
#'
#' 2. Use the cumulative distribution function (CDF) of the standard normal distribution
#'    to find the probability that a cost is less than the threshold:
#'    P(X < threshold_cost) = P(Z <= z) = pnorm(z)
#'
#' @param mean_cost The mean cost of the distribution (numeric).
#' @param sd_cost The standard deviation of the cost (numeric).
#' @param threshold_cost The cost threshold to compare against (numeric).
#'
#' @return The probability of the cost being less than the threshold (numeric).
#'
#' @examples
#' calculate_probability_cost_less_than(100, 15, 90)
#'
#' @export
calculate_probability_cost_less_than <- function(mean_cost, sd_cost, threshold_cost) {
  # Calculate the z-score for the threshold cost
  z_score <- (threshold_cost - mean_cost) / sd_cost

  # Calculate and return the probability of the cost being less than the threshold
  pnorm(z_score)
}

#' Calculate Probability of Cost Being Greater Than a Given Value
#'
#' This function calculates the probability that a cost is greater than a specified threshold,
#' based on the mean and standard deviation of the cost distribution. This is useful for
#' understanding the likelihood of costs exceeding a certain level in a normally distributed
#' context.
#'
#' The calculation involves the following steps:
#'
#' 1. Calculate the z-score corresponding to the threshold cost using the formula:
#'    z = (X - mu) / sigma
#'    where X is the threshold cost, mu is the mean cost, and sigma is the standard deviation of the cost.
#'
#' 2. Use the complementary cumulative distribution function (CDF) of the standard normal distribution
#'    to find the probability that a cost is greater than the threshold:
#'    P(X > threshold_cost) = 1 - P(Z <= z) = 1 - pnorm(z)
#'
#' @param mean_cost The mean cost of the distribution (numeric).
#' @param sd_cost The standard deviation of the cost (numeric).
#' @param threshold_cost The cost threshold to compare against (numeric).
#'
#' @return The probability of the cost being greater than the threshold (numeric).
#'
#' @examples
#' calculate_probability_cost_greater_than(100, 15, 110)
#'
#' @export
calculate_probability_cost_greater_than <- function(mean_cost, sd_cost, threshold_cost) {
  # Calculate the z-score for the threshold cost
  z_score <- (threshold_cost - mean_cost) / sd_cost

  # Calculate and return the probability of the cost being greater than the threshold
  1 - pnorm(z_score)
}

#' Perform a One-Sample T-Test and Calculate the P-Value
#'
#' This function performs a one-sample t-test to determine if the sample mean is significantly
#' greater than a specified population mean. The t-test assesses whether the difference
#' between the sample mean and the population mean is statistically significant.
#'
#' The calculation involves the following steps:
#'
#' 1. Calculate the sample mean (x_bar) and sample standard deviation (s) from the input data.
#'
#' 2. Compute the t-statistic using the formula:
#'    t = (x_bar - mu) / (s / sqrt(n))
#'    where mu is the population mean and n is the sample size.
#'
#' 3. Calculate the p-value for the one-tailed t-test using the t-distribution.
#'
#' 4. Interpret the result based on the specified significance level:
#'    - If the p-value is less than the significance level, conclude that the sample mean
#'      is significantly greater than the population mean.
#'    - Otherwise, conclude that it is not significantly greater.
#'
#' @param data A numeric vector of sample data (numeric).
#' @param population_mean The population mean to compare the sample mean against (numeric).
#' @param significance_level The significance level for the test (numeric, e.g., 0.05 for a 5% significance level).
#'
#' @return A list containing:
#'   - `p_value`: The calculated p-value (numeric).
#'   - `interpretation`: A string interpreting the result based on the significance level.
#'
#' @examples
#' calculate_p_value_greater_than(c(5.1, 5.5, 5.9, 6.2, 5.7), 5.0, 0.05)
#'
#' @export
calculate_p_value_greater_than <- function(data, population_mean, significance_level) {
  # Perform the one-sample t-test
  t_test_result <- t.test(data, mu = population_mean, alternative = "greater")

  # Extract the p-value from the test result
  p_value <- t_test_result$p.value

  # Interpret the result based on the significance level
  interpretation <- ifelse(
    p_value < significance_level,
    "The sample mean is significantly greater than the population mean.",
    "The sample mean is not significantly greater than the population mean."
  )

  # Return the p-value and interpretation as a list
  list(p_value = p_value, interpretation = interpretation)
}

#' cor3var: Calculate the correlation coefficient from standard deviations and slope
#'
#' This function calculates the correlation coefficient (r) using the standard
#' deviation of x (\( s_x \)), the standard deviation of y (\( s_y \)), and the slope of the
#' regression line. The correlation coefficient indicates the strength and direction
#' of the linear relationship between two variables.
#'
#' The calculation is based on the following formula:
#' r = (slope * s_x) / s_y
#'
#' @param s_x Standard deviation of x (numeric).
#' @param s_y Standard deviation of y (numeric).
#' @param slope The slope of the regression line (numeric).
#'
#' @return The correlation coefficient (r) (numeric).
#'
#' @examples
#' cor3var(1.5, 2.5, 0.8)
#'
#' @export
cor3var <- function(s_x, s_y, slope) {
  # Calculate the correlation coefficient (r)
  r <- (slope * s_x) / s_y
  return(r)
}


#' list_functions: List all available functions in the specified package
#'
#' This function lists all the available functions in the specified R package.
#' By default, it lists the functions in the 'allenprob' package, but you can
#' modify it to specify any package you want to inspect.
#'
#' @param package_name A string representing the name of the package to list functions from (default is "allenprob").
#' @importFrom utils lsf.str
#' @return A character vector of function names available in the specified package.
#'
#' @examples
#' list_functions()  # Lists functions in the 'allenprob' package
#' list_functions("stats")  # Lists functions in the 'stats' package
#'
#' @export
list_functions <- function(package_name = "allenprob") {
  # List all functions available in the specified package
  return(lsf.str(paste0("package:", package_name)))
}


#' r_squared: Calculate the coefficient of determination (R^2) from the correlation coefficient (r)
#'
#' This function calculates the coefficient of determination (R^2) from the given
#' correlation coefficient (r). The coefficient of determination indicates the
#' proportion of the variance in the dependent variable that is predictable
#' from the independent variable.
#'
#' The calculation is based on the formula:
#' R^2 = r^2
#'
#' @param r Correlation coefficient (numeric). Should be between -1 and 1.
#'
#' @return The coefficient of determination (R^2) (numeric). A value between 0 and 1,
#'         indicating the proportion of variance explained.
#'
#' @examples
#' r_squared(0.8)  # Returns 0.64
#'
#' @export
r_squared <- function(r) {
  # Calculate the coefficient of determination (R^2)
  r_squared <- r^2
  return(r_squared)
}


#' lm_summary: Perform linear regression and show the summary
#'
#' This function takes two data vectors x and y and fits a linear model to them.
#' The linear model is of the form:
#' y = beta_0 + beta_1 * x + epsilon
#' where beta_0 is the intercept, beta_1 is the slope, and epsilon is the error term.
#'
#' @param x A numeric vector of independent variable values (numeric).
#' @param y A numeric vector of dependent variable values (numeric).
#' @importFrom stats lm
#' @return A summary of the linear model, which includes:
#'   - Coefficients: Estimates of the model parameters (intercept and slope).
#'   - Residuals: Summary statistics of the residuals.
#'   - Statistical significance: p-values for the coefficients.
#'   - R-squared value: Proportion of variance explained by the model.
#'
#' @examples
#' lm_summary(c(1, 2, 3, 4), c(2, 4, 6, 8))  # Linear relationship should yield a slope of 2
#'
#' @export
lm_summary <- function(x, y) {
  # Fit a linear model lm(y ~ x)
  model <- lm(y ~ x)

  # Return the summary of the model
  return(summary(model))
}

#' least_squares_regression: Perform least squares regression using means, standard deviations, and correlation
#'
#' This function calculates the least squares regression line using the means of the dependent variable
#' (ybar) and the independent variable (xbar), the standard deviations of x (s_x) and y (s_y),
#' and the correlation coefficient (R) between the two variables.
#'
#' The regression line is expressed in the form:
#' y = intercept + slope * x
#'
#' The slope is calculated using the formula:
#' slope = R * (s_y / s_x)
#'
#' The intercept is calculated as:
#' intercept = ybar - slope * xbar
#'
#' @param ybar Mean of y values (numeric).
#' @param xbar Mean of x values (numeric).
#' @param s_x Standard deviation of x values (numeric).
#' @param s_y Standard deviation of y values (numeric).
#' @param R Correlation coefficient between x and y (numeric, should be between -1 and 1).
#'
#' @return A list containing:
#'   - `slope`: The slope of the regression line (numeric).
#'   - `intercept`: The intercept of the regression line (numeric).
#'
#' @examples
#' least_squares_regression(5, 3, 1.2, 2.1, 0.8)  # Example calculation
#'
#' @export
least_squares_regression <- function(ybar, xbar, s_x, s_y, R) {
  # Calculate the slope of the regression line
  slope <- R * (s_y / s_x)

  # Calculate the intercept of the regression line
  intercept <- ybar - slope * xbar

  # Return a list containing the slope and intercept
  return(list(slope = slope, intercept = intercept))
}


#' sample_mean_variance: Calculate the mean and variance of a sample distribution
#'
#' This function calculates the mean and variance of a sample distribution given two vectors:
#' a vector of values (X) and a vector of probabilities (P(X)) associated with those values.
#' The function assumes that the probabilities are normalized (i.e., they sum to 1).
#'
#' The calculations are performed as follows:
#'
#' 1. The mean is calculated using the formula:
#'    mean = sum(X * P(X))
#'
#' 2. The variance is calculated using the formula:
#'    variance = sum((X - mean)^2 * P(X))
#'
#' @param X A numeric vector of values (numeric).
#' @param P_X A numeric vector of probabilities associated with X (numeric).
#'             Length of P_X must be equal to length of X, and the probabilities should sum to 1.
#'
#' @return A list containing:
#'   - `mean`: The mean of the sample distribution (numeric).
#'   - `variance`: The variance of the sample distribution (numeric).
#'
#' @examples
#' sample_mean_variance(c(1, 2, 3), c(0.2, 0.5, 0.3))  # Example calculation
#'
#' @export
sample_mean_variance <- function(X, P_X) {
  # Calculate the mean of the distribution
  mean_value <- sum(X * P_X)

  # Calculate the variance of the distribution
  variance_value <- sum(((X - mean_value)^2) * P_X)

  # Return a list containing the mean and variance
  return(list(mean = mean_value, variance = variance_value))
}


#' expected_value: Calculate the expected value given win/lose amounts and probabilities
#'
#' This function calculates the expected value of an event based on the win amount,
#' lose amount, probability of winning, and probability of losing. The expected value
#' provides a measure of the average outcome of a probabilistic event.
#'
#' The expected value is calculated using the formula:
#' expected_value = (win_amount * prob_win) + (lose_amount * prob_loss)
#'
#' @param win_amount The amount won if the event is successful (numeric).
#' @param lose_amount The amount lost if the event is unsuccessful (numeric).
#' @param prob_win The probability of winning the event (numeric, should be between 0 and 1).
#' @param prob_loss The probability of losing the event (numeric, should be between 0 and 1).
#'
#' @return The expected value of the event (numeric).
#'
#' @examples
#' expected_value(100, -50, 0.6, 0.4)  # Expected value calculation
#'
#' @export
expected_value <- function(win_amount, lose_amount, prob_win, prob_loss) {
  # Calculate the expected value
  expected_value <- (win_amount * prob_win) + (lose_amount * prob_loss)

  return(expected_value)
}


#' Calculate Mean and Standard Deviation of the Sampling Distribution
#'
#' This function calculates the mean and standard deviation (standard error) of the
#' sampling distribution for a given population proportion and sample size.
#' The mean of the sampling distribution for proportions is equal to the population
#' proportion, and the standard deviation (standard error) is derived from the formula
#' for the standard deviation of a proportion.
#'
#' The calculations are performed as follows:
#'
#' 1. The mean of the sampling distribution is equal to the population proportion:
#'    mean = p
#'
#' 2. The standard deviation (standard error) of the sampling distribution is calculated using the formula:
#'    standard deviation = sqrt(p * (1 - p) / n)
#'
#' @param p The population proportion (numeric, should be between 0 and 1).
#' @param n The sample size (integer, must be greater than 0).
#'
#' @return A list containing:
#'   - `mean`: The mean of the sampling distribution (numeric).
#'   - `standard_deviation`: The standard deviation (standard error) of the sampling distribution (numeric).
#'
#' @examples
#' calculate_sampling_distribution(0.67, 64)  # Example calculation
#'
#' @export
calculate_sampling_distribution <- function(p, n) {
  # Calculate the mean of the sampling distribution
  mean_sampling_distribution <- p

  # Calculate the standard deviation (standard error) of the sampling distribution
  standard_deviation_sampling_distribution <- sqrt((p * (1 - p)) / n)

  # Return a list containing the mean and standard deviation
  list(mean = mean_sampling_distribution, standard_deviation = standard_deviation_sampling_distribution)
}


#' Calculate Required Sample Size for Confidence Interval of a Proportion
#'
#' This function calculates the minimum sample size needed to estimate the proportion
#' of a certain outcome in a population with a specified confidence level and margin of error.
#' This is essential for designing surveys and experiments to ensure that the estimates
#' fall within the desired accuracy.
#'
#' The sample size is calculated using the following formula:
#' n = (z * sqrt(p * (1 - p)) / E)^2
#'
#' Where:
#' - z is the z-score corresponding to the desired confidence level,
#' - p is the estimated proportion of the outcome,
#' - E is the acceptable margin of error.
#'
#' @param proportion The proportion of the outcome in the population (numeric, should be between 0 and 1).
#' @param confidence_level The desired confidence level (numeric, e.g., 0.99 for 99% confidence).
#' @param margin_of_error The acceptable margin of error for the confidence interval (numeric, should be greater than 0).
#'
#' @return The calculated sample size as an integer (numeric).
#'
#' @examples
#' calculate_sample_size_for_proportion(0.10, 0.99, 0.13)  # Example calculation
#'
#' @export
calculate_sample_size_for_proportion <- function(proportion, confidence_level, margin_of_error) {
  # Calculate the z-score for the specified confidence level
  z_score <- qnorm(1 - (1 - confidence_level) / 2)

  # Calculate the sample size using the formula
  sample_size <- (z_score * sqrt(proportion * (1 - proportion)) / margin_of_error) ^ 2

  # Round up to the next whole number since sample size must be an integer
  ceiling(sample_size)
}

#' Perform a Hypothesis Test for a Fair Coin Flip
#'
#' This function tests the claim that a coin is fair (p = 0.5) against the alternative
#' that it is not fair (p ≠ 0.5) using a two-tailed hypothesis test. The function
#' calculates the z statistic for the observed number of heads in a series of coin flips,
#' and compares it against the critical value derived from the specified significance level.
#'
#' The null hypothesis (H0) states that the coin is fair, while the alternative hypothesis (H1) states that
#' the coin is not fair.
#'
#' @param heads An integer representing the number of heads observed in the coin flips (integer).
#' @param total_flips An integer representing the total number of coin flips (integer).
#' @param alpha A numeric value representing the significance level (default is 0.01, should be between 0 and 1).
#'
#' @return A list containing:
#'   - `z_value`: The calculated z statistic (numeric).
#'   - `critical_value`: The critical z value for the two-tailed test (numeric).
#'   - `decision`: A string indicating whether to "Reject H0" or "Fail to Reject H0".
#'
#' @examples
#' coin_test_result <- coin_hypothesis_test(heads = 44, total_flips = 100)
#' print(coin_test_result)
#'
#' @export
coin_hypothesis_test <- function(heads, total_flips, alpha = 0.01) {
  # Hypothesized proportion for a fair coin
  p_0 <- 0.5

  # Calculate the sample proportion
  p_hat <- heads / total_flips

  # Calculate the standard error of the sample proportion
  standard_error <- sqrt((p_0 * (1 - p_0)) / total_flips)

  # Calculate the z statistic
  z_value <- (p_hat - p_0) / standard_error

  # Calculate the critical value for a two-tailed test
  critical_value <- qnorm(1 - alpha / 2)

  # Make a decision based on the z value and critical value
  if (abs(z_value) > critical_value) {
    decision <- "Reject H0"
  } else {
    decision <- "Fail to Reject H0"
  }

  # Return the results as a list
  return(list(z_value = z_value,
              critical_value = critical_value,
              decision = decision))
}

#' Determine the Conclusion of a Hypothesis Test Based on the p-Value
#'
#' This function evaluates the evidence against the null hypothesis based on the
#' computed p-value from a hypothesis test. It provides a conclusion about whether
#' to reject or fail to reject the null hypothesis, with a threshold for very strong
#' evidence set at a p-value of 0.001.
#'
#' @param p_value A numeric value representing the computed p-value from the hypothesis test.
#'                It should be a number between 0 and 1.
#'
#' @return A string indicating the conclusion:
#'   - "Reject the null hypothesis" if the p-value < 0.001.
#'   - "Fail to reject the null hypothesis" if the p-value >= 0.001.
#'
#' @examples
#' conclusion <- hypothesis_conclusion(p_value = 0.0005)
#' print(conclusion)  # Expected output: "Reject the null hypothesis"
#'
#' @export
hypothesis_conclusion <- function(p_value) {
  # Check the validity of the p-value
  if (!is.numeric(p_value) || p_value < 0 || p_value > 1) {
    stop("p_value must be a numeric value between 0 and 1.")
  }

  # Determine the conclusion based on the p-value
  if (p_value < 0.001) {
    conclusion <- "Reject the null hypothesis"
  } else {
    conclusion <- "Fail to reject the null hypothesis"
  }

  return(conclusion)
}
