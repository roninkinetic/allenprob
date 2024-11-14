utils::globalVariables(c("x"))
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
#' that it is not fair (p <> 0.5) using a two-tailed hypothesis test. The function
#' calculates the z statistic for the observed number of heads in a series of coin flips,
#' and compares it against the critical value derived from the specified significance level.
#'
#' The null hypothesis (H0) states that the coin is fair, while the alternative hypothesis (H1) states that
#' the coin is not fair.
#'
#' @param heads An integer representing the number of heads observed in the coin flips (integer).
#' @param total_flips An integer representing the total number of coin flips (integer).
#' @param alpha A numeric value representing the significance level (default is 0.01).
#'
#' @return A list containing:
#'   - `z_value`: The calculated z statistic (numeric).
#'   - `critical_value`: The critical z value for the two-tailed test (numeric).
#'   - `decision`: A string indicating whether to "Reject H0" or "Fail to Reject H0".
#'
#' @examples
#' coin_test_result <- coin_hypothesis_test(heads = 56, total_flips = 100)
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


#' Calculate Probability of Costs Between Two Values
#'
#' This function calculates the probability that a cost, following a normal distribution, falls between two specified values.
#'
#' The costs are assumed to follow a normal distribution with a specified mean and standard deviation.
#'
#' @param lower_bound The lower bound of the cost range (numeric).
#' @param upper_bound The upper bound of the cost range (numeric).
#' @param mean The mean of the normal distribution (numeric).
#' @param sd The standard deviation of the normal distribution (numeric).
#' @return The probability that a cost falls between the lower and upper bounds (numeric).
#' @examples
#' probability_between_costs(41, 107, 74, 22)
#' @export
probability_between_costs <- function(lower_bound, upper_bound, mean, sd) {
  # Calculate the cumulative probabilities at the lower and upper bounds
  prob_lower <- pnorm(lower_bound, mean = mean, sd = sd)
  prob_upper <- pnorm(upper_bound, mean = mean, sd = sd)

  # Calculate the probability that the cost is between the two bounds
  probability <- prob_upper - prob_lower

  return(probability)
}

#' Calculate Mean and Standard Deviation of the Sampling Distribution of Proportions
#'
#' This function calculates the mean and standard deviation of the sampling distribution
#' for a given population proportion and sample size.
#'
#' The mean of the sampling distribution is equal to the population proportion, and the
#' standard deviation (standard error) is calculated based on the population proportion
#' and the sample size.
#'
#' @param population_proportion The proportion of the population with the characteristic of interest (numeric).
#' @param sample_size The size of the sample (integer).
#' @return A list containing the mean and standard deviation of the sampling distribution.
#' @examples
#' sampling_distribution_parameters(0.74, 100)
#' @export
sampling_distribution_parameters <- function(population_proportion, sample_size) {
  # Calculate the mean of the sampling distribution
  mean_sampling_distribution <- population_proportion

  # Calculate the standard deviation (standard error) of the sampling distribution
  standard_deviation_sampling_distribution <- sqrt((population_proportion * (1 - population_proportion)) / sample_size)

  # Return the results as a list
  return(list(mean = mean_sampling_distribution, standard_deviation = standard_deviation_sampling_distribution))
}

#' Solve Uniform Distribution Problems
#'
#' This function solves problems related to a uniform distribution, where the density curve is defined
#' between a lower bound (x = 1) and an upper bound (x = 8). It calculates the height of the uniform distribution curve,
#' the density function value at a given point, the percentage of observations in a specified range,
#' and the probability of specific values (including below, above, and equal).
#'
#' @param x A numeric value or a vector of numeric values specifying the bounds for the calculations.
#' @param lower_bound The lower bound of the uniform distribution (default is 1).
#' @param upper_bound The upper bound of the uniform distribution (default is 8).
#'
#' @return A list containing:
#'   - `height`: The height of the uniform distribution curve.
#'   - `density_function`: The value of the uniform density function at `x`.
#'   - `percent_in_range`: The percent of observations falling within the specified range.
#'   - `percent_below`: The percent of observations below the specified value.
#'   - `percent_above`: The percent of observations above the specified value.
#'   - `percent_equal`: The percent of observations equal to a specified value (always 0 for continuous distributions).
#'
#' @examples
#' uniform_distribution(2, 1, 8)
#' uniform_distribution(c(2, 5), 1, 8)  # Percent between 2 and 5
#' uniform_distribution(4, 1, 8)  # Percent below 4
#' uniform_distribution(6, 1, 8)  # Percent above 6
#' uniform_distribution(7, 1, 8)  # Percent equal to 7 (always 0)
#'
#' @export
uniform_distribution <- function(x, lower_bound, upper_bound) {
  # Calculate the height of the uniform distribution curve
  height <- 1 / (upper_bound - lower_bound)

  # Calculate density function value at x
  density_function <- ifelse(x >= lower_bound & x <= upper_bound, height, 0)

  # Calculate percent between two points if two values are provided
  if (length(x) == 2) {
    percent_in_range <- (x[2] - x[1]) * height * 100
  } else {
    percent_in_range <- NA
  }

  # Calculate percent below x
  if (!is.null(x)) {
    percent_below <- (x - lower_bound) * height * 100
  }

  # Calculate percent above x
  percent_above <- (upper_bound - x) * height * 100

  # Percent equal to x (always 0 for continuous distributions)
  percent_equal <- 0

  return(list(height = height,
              density_function = density_function,
              percent_in_range = percent_in_range,
              percent_below = percent_below,
              percent_above = percent_above,
              percent_equal = percent_equal))
}

#' Calculate Probability Between Two Values for a Uniform Distribution
#'
#' This function calculates the probability (percent) that an observation falls between two given values
#' within the range of a uniform distribution defined between a lower and upper bound.
#'
#' @param lower_bound The lower bound of the uniform distribution (default is 1).
#' @param upper_bound The upper bound of the uniform distribution (default is 8).
#' @param x1 The first value in the range.
#' @param x2 The second value in the range.
#'
#' @return The percentage of observations falling between `x1` and `x2`.
#'
#' @examples
#' uniform_probability_between(1, 8, 2, 5)  # Probability between 2 and 5
#' @export
uniform_probability_between <- function(lower_bound, upper_bound, x1, x2) {
  # Calculate the height of the uniform distribution curve
  height <- 1 / (upper_bound - lower_bound)

  # Calculate the probability (percent) between x1 and x2
  percent_in_range <- (x2 - x1) * height * 100

  return(percent_in_range)
}

#' Compute CDF, Probabilities, and Expectations for the Breakage Point Distribution
#'
#' This function calculates the cumulative distribution function (CDF), probabilities,
#' expectations, and variances for the distribution of the breakage point in a 12-inch bar.
#' The breakage point, denoted as Y, has a probability density function (PDF) defined by the user.
#'
#' The function can compute the following:
#' - The CDF of Y.
#' - Specific probabilities such as P(Y <equalsto x1), P(Y > x2), and P(x1 <equalsto Y <equalsto x2).
#' - The expected value E(Y), the second moment E(Y^2), and the variance Var(Y).
#' - The probability that the break occurs more than a specified distance from the expected break point.
#' - The expected length of the shorter segment when the break occurs.
#'
#' @param pdf_function A function that calculates the PDF of Y.
#' @param a Lower bound for the probability calculations (default is 0).
#' @param b Upper bound for the probability calculations (default is 12, for the length of the bar).
#' @param expected_break The expected break point (default is the mean of the distribution).
#' @param x1 First threshold for probability calculations (default is 4).
#' @param x2 Second threshold for probability calculations (default is 6).
#' @param distance_from_expected The distance from the expected break point for probability calculation (default is 2).
#' @return A list containing:
#'   - `cdf` The cumulative distribution function.
#'   - `probabilities` A list of probabilities for specific ranges.
#'   - `expectations` A list containing the expected value, second moment, and variance.
#'   - `shorter_segment_expectation` The expected length of the shorter segment when the break occurs.
#' @importFrom stats integrate
#'
#' @examples
#' # Define a simple PDF function for a uniform distribution over the range [0, 12]
#' pdf_function <- function(y) { ifelse(y >= 0 & y <= 12, 1 / 12, 0) }
#' breakage_analysis(pdf_function, x1 = 4, x2 = 6, distance_from_expected = 2)
#'
#' @export
breakage_analysis <- function(pdf_function, a = 0, b = 12, expected_break = 6, x1 = 4, x2 = 6, distance_from_expected = 2) {
  # Function to calculate the CDF from the given PDF
  cdf_function <- function(y) {
    integrate(pdf_function, a, y)$value
  }

  # Compute the CDF of Y
  cdf_values <- sapply(seq(a, b, length.out = 100), cdf_function)

  # Calculate probabilities
  P_Y_less_than_x1 <- cdf_function(x1)
  P_Y_greater_than_x2 <- 1 - cdf_function(x2)
  P_Y_between_x1_and_x2 <- cdf_function(x2) - cdf_function(x1)

  probabilities <- list(P_Y_less_than_x1 = P_Y_less_than_x1,
                        P_Y_greater_than_x2 = P_Y_greater_than_x2,
                        P_Y_between_x1_and_x2 = P_Y_between_x1_and_x2)

  # Compute Expectations and Variance
  E_Y <- integrate(function(y) y * pdf_function(y), a, b)$value
  E_Y2 <- integrate(function(y) y^2 * pdf_function(y), a, b)$value
  variance_Y <- E_Y2 - E_Y^2

  expectations <- list(E_Y = E_Y,
                       E_Y2 = E_Y2,
                       Var_Y = variance_Y)

  # Probability that the break occurs more than `distance_from_expected` from expected break point
  P_more_than_distance_from_expected <- 1 - cdf_function(expected_break + distance_from_expected) + cdf_function(expected_break - distance_from_expected)

  # Expected length of the shorter segment
  shorter_segment_expectation <- expected_break / 2  # Assuming uniform break distribution

  return(list(cdf = cdf_values,
              probabilities = probabilities,
              expectations = expectations,
              shorter_segment_expectation = shorter_segment_expectation))
}


#' Compute Parameters and Probability for Gamma Distribution with Mean and Variance
#'
#' This function calculates the parameters (alpha and beta) of the gamma distribution
#' based on the provided mean and variance. It also computes the probability that the
#' daily power consumption exceeds a specified threshold (12 million kilowatt-hours).
#'
#' The gamma distribution has the following characteristics:
#' - The mean is equal to the product of alpha and beta.
#' - The variance is equal to alpha times the square of beta.
#'
#' The function solves for the values of alpha (shape parameter) and beta (scale parameter),
#' and then calculates the probability that the power consumption exceeds the threshold using the cumulative
#' distribution function (CDF) of the gamma distribution.
#'
#' @param mu The mean of the distribution (in millions of kilowatt-hours).
#' @param sigma2 The variance of the distribution (in millions of kilowatt-hours squared).
#' @param threshold The threshold value for which the probability of exceeding is calculated (default is 12).
#' @return A list containing:
#'   - `alpha`: The shape parameter of the gamma distribution.
#'   - `beta`: The scale parameter of the gamma distribution.
#'   - `prob_exceed_threshold`: The probability that the power consumption exceeds the threshold.
#'
#' @importFrom stats pgamma
#' @examples
#' result <- gamma_with_mean_and_var(mu = 6, sigma2 = 12)
#' print(result)
#'
#' @export
gamma_with_mean_and_var <- function(mu, sigma2, threshold = 12) {
  # Calculate alpha and beta from the given mean (mu) and variance (sigma2)
  beta <- sigma2 / mu
  alpha <- mu^2 / sigma2

  # Calculate the probability that the power consumption exceeds the threshold
  prob_exceed_threshold <- 1 - pgamma(threshold, shape = alpha, scale = beta)

  # Return the results as a list
  return(list(alpha = alpha, beta = beta, prob_exceed_threshold = prob_exceed_threshold))
}


#' Compute CDF, Probabilities, and Expectations for the Breakage Point Distribution
#'
#' This function calculates the cumulative distribution function (CDF), probabilities,
#' expectations, and variances for the distribution of the breakage point in a 12-inch bar.
#' The breakage point, denoted as Y, has a probability density function (PDF) defined by the user.
#'
#' The function can compute the following:
#' - The CDF of Y.
#' - Specific probabilities such as P(Y <equals x1), P(Y > x2), and P(x1 <equals Y <equals x2).
#' - The expected value E(Y), the second moment E(Y^2), and the variance Var(Y).
#' - The probability that the break occurs more than a specified distance from the expected break point.
#' - The expected length of the shorter segment when the break occurs.
#'
#' @param pdf_function A function that calculates the PDF of Y.
#' @param a Lower bound for the probability calculations (default is 0).
#' @param b Upper bound for the probability calculations (default is 12, for the length of the bar).
#' @param expected_break The expected break point (default is the mean of the distribution).
#' @param x1 First threshold for probability calculations (default is 4).
#' @param x2 Second threshold for probability calculations (default is 6).
#' @param distance_from_expected The distance from the expected break point for probability calculation (default is 2).
#' @return A list containing:
#'   - `cdf` The cumulative distribution function.
#'   - `probabilities` A list of probabilities for specific ranges.
#'   - `expectations` A list containing the expected value, second moment, and variance.
#'   - `shorter_segment_expectation` The expected length of the shorter segment when the break occurs.
#'
#' @examples
#' # Define a simple PDF function for a uniform distribution over the range [0, 12]
#' pdf_function <- function(y) { ifelse(y >= 0 & y <= 12, 1 / 12, 0) }
#' breakage_analysis(pdf_function, x1 = 4, x2 = 6, distance_from_expected = 2)
#'
#' @export
y_cdf_function <- function(pdf_function, a = 0, b = 12, expected_break = 6, x1 = 4, x2 = 6, distance_from_expected = 2) {
  # Function to calculate the CDF from the given PDF
  cdf_function <- function(y) {
    integrate(pdf_function, a, y)$value
  }

  # Compute the CDF of Y
  cdf_values <- sapply(seq(a, b, length.out = 100), cdf_function)

  # Calculate probabilities
  P_Y_less_than_x1 <- cdf_function(x1)
  P_Y_greater_than_x2 <- 1 - cdf_function(x2)
  P_Y_between_x1_and_x2 <- cdf_function(x2) - cdf_function(x1)

  probabilities <- list(P_Y_less_than_x1 = P_Y_less_than_x1,
                        P_Y_greater_than_x2 = P_Y_greater_than_x2,
                        P_Y_between_x1_and_x2 = P_Y_between_x1_and_x2)

  # Compute Expectations and Variance
  E_Y <- integrate(function(y) y * pdf_function(y), a, b)$value
  E_Y2 <- integrate(function(y) y^2 * pdf_function(y), a, b)$value
  variance_Y <- E_Y2 - E_Y^2

  expectations <- list(E_Y = E_Y,
                       E_Y2 = E_Y2,
                       Var_Y = variance_Y)

  # Probability that the break occurs more than `distance_from_expected` from expected break point
  P_more_than_distance_from_expected <- 1 - cdf_function(expected_break + distance_from_expected) + cdf_function(expected_break - distance_from_expected)

  # Expected length of the shorter segment
  shorter_segment_expectation <- expected_break / 2  # Assuming uniform break distribution

  return(list(cdf = cdf_values,
              probabilities = probabilities,
              expectations = expectations,
              shorter_segment_expectation = shorter_segment_expectation))
}
#' Compute CDF, Probabilities, and Expectations for the Breakage Point Distribution
#'
#' This function calculates the cumulative distribution function (CDF), probabilities,
#' expectations, and variances for the distribution of the breakage point in a 12-inch bar.
#' The breakage point, denoted as Y, has a probability density function (PDF) defined by the user.
#'
#' The function can compute the following:
#' - The CDF of Y.
#' - Specific probabilities such as P(Y <equals x1), P(Y > x2), and P(x1 <equals Y <equals x2).
#' - The expected value E(Y), the second moment E(Y^2), and the variance Var(Y).
#' - The probability that the break occurs more than a specified distance from the expected break point.
#' - The expected length of the shorter segment when the break occurs.
#'
#' @param pdf_function A function that calculates the PDF of Y.
#' @param a Lower bound for the probability calculations (default is 0).
#' @param b Upper bound for the probability calculations (default is 12, for the length of the bar).
#' @param expected_break The expected break point (default is the mean of the distribution).
#' @param x1 First threshold for probability calculations (default is 4).
#' @param x2 Second threshold for probability calculations (default is 6).
#' @param distance_from_expected The distance from the expected break point for probability calculation (default is 2).
#' @return A list containing:
#'   - `cdf` The cumulative distribution function.
#'   - `probabilities` A list of probabilities for specific ranges.
#'   - `expectations` A list containing the expected value, second moment, and variance.
#'   - `shorter_segment_expectation` The expected length of the shorter segment when the break occurs.
#'
#' @examples
#' # Define a simple PDF function for a uniform distribution over the range [0, 12]
#' pdf_function <- function(y) { ifelse(y >= 0 & y <= 12, 1 / 12, 0) }
#' breakage_analysis(pdf_function, x1 = 4, x2 = 6, distance_from_expected = 2)
#'
#' @export
y_pdf_cdf <- function(pdf_function, a = 0, b = 12, expected_break = 6, x1 = 4, x2 = 6, distance_from_expected = 2) {
  # Function to calculate the CDF from the given PDF
  cdf_function <- function(y) {
    integrate(pdf_function, a, y)$value
  }

  # Compute the CDF of Y
  cdf_values <- sapply(seq(a, b, length.out = 100), cdf_function)

  # Calculate probabilities
  P_Y_less_than_x1 <- cdf_function(x1)
  P_Y_greater_than_x2 <- 1 - cdf_function(x2)
  P_Y_between_x1_and_x2 <- cdf_function(x2) - cdf_function(x1)

  probabilities <- list(P_Y_less_than_x1 = P_Y_less_than_x1,
                        P_Y_greater_than_x2 = P_Y_greater_than_x2,
                        P_Y_between_x1_and_x2 = P_Y_between_x1_and_x2)

  # Compute Expectations and Variance
  E_Y <- integrate(function(y) y * pdf_function(y), a, b)$value
  E_Y2 <- integrate(function(y) y^2 * pdf_function(y), a, b)$value
  variance_Y <- E_Y2 - E_Y^2

  expectations <- list(E_Y = E_Y,
                       E_Y2 = E_Y2,
                       Var_Y = variance_Y)

  # Probability that the break occurs more than `distance_from_expected` from expected break point
  P_more_than_distance_from_expected <- 1 - cdf_function(expected_break + distance_from_expected) + cdf_function(expected_break - distance_from_expected)

  # Expected length of the shorter segment
  shorter_segment_expectation <- expected_break / 2  # Assuming uniform break distribution

  return(list(cdf = cdf_values,
              probabilities = probabilities,
              expectations = expectations,
              shorter_segment_expectation = shorter_segment_expectation))
}

#' This function plots the normal distribution for a given mean and standard deviation.
#' It also answers the question: According to the Empirical Rule, what are the bounds for
#' the middle 68% of the data?
#'
#' @param mu The mean of the distribution.
#' @param sigma The standard deviation of the distribution.
#' @param x_range The range of x-axis values for plotting.
#'
#' @return A plot of the normal distribution and the middle 68% range.
#'
#' @importFrom graphics curve abline
#' @importFrom stats dnorm
#' @examples
#' plot_normal_and_empirical(mu = 82, sigma = 4)
#'
#' @export
plot_normal_and_empirical <- function(mu, sigma, x_range = c(mu - 4*sigma, mu + 4*sigma)) {
  # Plotting the normal distribution
  curve(dnorm(x, mean = mu, sd = sigma), from = x_range[1], to = x_range[2],
        col = "blue", lwd = 2, ylab = "Density", xlab = "X", main = "Normal Distribution")

  # Marking the middle 68% according to the Empirical Rule
  x68_low <- mu - sigma
  x68_high <- mu + sigma
  abline(v = c(x68_low, x68_high), col = "red", lwd = 2, lty = 2)

  # Displaying bounds for middle 68% in the console
  cat("The middle 68% of the data falls between:", x68_low, "and", x68_high, "\n")
}

#' Find Probability P(X < 83)
#'
#' This function calculates the probability that a normal random variable X is less than a given value.
#'
#' @param x The value to calculate the probability for.
#' @param mu The mean of the distribution.
#' @param sigma The standard deviation of the distribution.
#'
#' @return The probability P(X < x).
#'
#' @examples
#' prob_less_than_83(mu = 82, sigma = 4, x = 83)
#'
#' @export
prob_less_than_83 <- function(mu, sigma, x) {
  pnorm(x, mean = mu, sd = sigma)
}


#' Find Probability P(X > 79)
#'
#' This function calculates the probability that a normal random variable X is greater than a given value.
#'
#' @param x The value to calculate the probability for.
#' @param mu The mean of the distribution.
#' @param sigma The standard deviation of the distribution.
#'
#' @return The probability P(X > x).
#'
#' @examples
#' prob_greater_than_79(mu = 82, sigma = 4, x = 79)
#'
#' @export
prob_greater_than_79 <- function(mu, sigma, x) {
  1 - pnorm(x, mean = mu, sd = sigma)
}


#' Find Probability P(73 < X < 84)
#'
#' This function calculates the probability that a normal random variable X lies between two values.
#'
#' @param lower The lower bound of the interval.
#' @param upper The upper bound of the interval.
#' @param mu The mean of the distribution.
#' @param sigma The standard deviation of the distribution.
#'
#' @return The probability P(lower < X < upper).
#'
#' @examples
#' prob_between_73_and_84(mu = 82, sigma = 4, lower = 73, upper = 84)
#'
#' @export
prob_between_73_and_84 <- function(mu, sigma, lower, upper) {
  pnorm(upper, mean = mu, sd = sigma) - pnorm(lower, mean = mu, sd = sigma)
}


#' Find Probability P(X <equals x)
#'
#' This function calculates the probability that a normal random variable X is less than or equal to a given value.
#'
#' @param x The value to calculate the probability for.
#' @param mu The mean of the distribution.
#' @param sigma The standard deviation of the distribution.
#'
#' @return The probability P(X <equals x).
#'
#' @examples
#' prob_less_than_or_equal_x(mu = 82, sigma = 4, x = 83)
#'
#' @export
prob_less_than_or_equal_x <- function(mu, sigma, x) {
  pnorm(x, mean = mu, sd = sigma)
}


#' Find Probability P(X >equals x)
#'
#' This function calculates the probability that a normal random variable X is greater than or equal to a given value.
#'
#' @param x The value to calculate the probability for.
#' @param mu The mean of the distribution.
#' @param sigma The standard deviation of the distribution.
#'
#' @return The probability P(X >equals x).
#'
#' @examples
#' prob_greater_than_or_equal_x(mu = 82, sigma = 4, x = 79)
#'
#' @export
prob_greater_than_or_equal_x <- function(mu, sigma, x) {
  1 - pnorm(x, mean = mu, sd = sigma)
}


#' Find Probability P(X = upper and lower bounds)
#'
#' This function calculates the probability that a normal random variable X is exactly at the lower and upper bounds.
#'
#' @param lower The lower bound.
#' @param upper The upper bound.
#' @param mu The mean of the distribution.
#' @param sigma The standard deviation of the distribution.
#'
#' @return The probability P(X = lower and upper).
#'
#' @examples
#' prob_exact_bounds(mu = 82, sigma = 4, lower = 73, upper = 84)
#'
#' @export
prob_exact_bounds <- function(mu, sigma, lower, upper) {
  pnorm(upper, mean = mu, sd = sigma) - pnorm(lower, mean = mu, sd = sigma)
}


#' Find x such that P(X < x) = 0.97725
#'
#' This function calculates the value of x such that the cumulative probability P(X < x) equals a given probability.
#'
#' @param probability The cumulative probability for which to find x (e.g., 0.97725).
#' @param mu The mean of the distribution.
#' @param sigma The standard deviation of the distribution.
#'
#' @return The value of x corresponding to the given cumulative probability.
#'
#' @examples
#' find_x_given_probability(mu = 82, sigma = 4, probability = 0.97725)
#'
#' @export
find_x_given_probability <- function(mu, sigma, probability) {
  qnorm(probability, mean = mu, sd = sigma)
}

#' Mean and Standard Deviation for Z
#'
#' This function returns the mean and standard deviation for the standard normal random variable Z,
#' which is always 0 for the mean and 1 for the standard deviation.
#'
#' @return A list containing the mean and standard deviation of Z.
#'
#' @examples
#' mean_and_sd_Z()
#'
#' @export
mean_and_sd_Z <- function() {
  list(mean = 0, sd = 1)
}

#' Plot the Distribution of Z
#'
#' This function plots the standard normal distribution with a mean of 0 and a standard deviation of 1.
#'
#' @param x_range The range of x values to plot the distribution.
#'
#' @return A plot of the standard normal distribution.
#'
#' @examples
#' plot_standard_normal_distribution()
#'
#' @export
plot_standard_normal_distribution <- function(x_range = c(-4, 4)) {
  curve(dnorm(x), from = x_range[1], to = x_range[2], col = "blue",
        main = "Standard Normal Distribution", xlab = "Z", ylab = "Density")
}

#' Find Probability P(Z < target)
#'
#' This function calculates the probability that the standard normal variable Z is less than a given target.
#'
#' @param target The value of Z for which to calculate the probability.
#'
#' @return The probability P(Z < target).
#'
#' @examples
#' prob_Z_less_than_target(1.2)
#'
#' @export
prob_Z_less_than_target <- function(target) {
  pnorm(target)
}

#' Find Probability P(Z > target)
#'
#' This function calculates the probability that the standard normal variable Z is greater than a given target.
#'
#' @param target The value of Z for which to calculate the probability.
#'
#' @return The probability P(Z > target).
#'
#' @examples
#' prob_Z_greater_than_target(1.2)
#'
#' @export
prob_Z_greater_than_target <- function(target) {
  1 - pnorm(target)
}

#' Find Probability P(lower < Z < upper)
#'
#' This function calculates the probability that the standard normal variable Z lies between two values, lower and upper.
#'
#' @param lower The lower bound for Z.
#' @param upper The upper bound for Z.
#'
#' @return The probability P(lower < Z < upper).
#'
#' @examples
#' prob_Z_between_target_values(-0.45, 1.96)
#'
#' @export
prob_Z_between_target_values <- function(lower, upper) {
  pnorm(upper) - pnorm(lower)
}

#' Find c such that P(Z < c) = target_probability
#'
#' This function finds the value of c such that P(Z < c) equals a given probability.
#'
#' @param target_probability The target probability for which to find the value of c.
#'
#' @return The value of c such that P(Z < c) = target_probability.
#'
#' @examples
#' find_c_for_probability_target(0.845)
#'
#' @export
find_c_for_probability_target <- function(target_probability) {
  qnorm(target_probability)
}

#' Find c such that P(Z > c) = target_probability
#'
#' This function finds the value of c such that P(Z > c) equals a given probability.
#'
#' @param target_probability The target probability for the upper tail.
#'
#' @return The value of c such that P(Z > c) = target_probability.
#'
#' @examples
#' find_c_for_probability_greater_than_target(0.845)
#'
#' @export
find_c_for_probability_greater_than_target <- function(target_probability) {
  qnorm(1 - target_probability)
}

#' Find c such that P(-c < Z < c) = target_probability
#'
#' This function finds the value of c such that P(-c < Z < c) equals a given probability.
#'
#' @param target_probability The desired cumulative probability for the two-tailed probability.
#'
#' @return The value of c such that P(-c < Z < c) = target_probability.
#'
#' @examples
#' find_c_for_two_tailed_probability_target(0.845)
#'
#' @export
find_c_for_two_tailed_probability_target <- function(target_probability) {
  qnorm((1 + target_probability) / 2)
}

#' Calculate Mean and Standard Error of the Sampling Distribution
#'
#' This function calculates the mean and standard error of the sampling distribution of the sample mean.
#' Given the population mean and population variance, it computes the mean and standard error for a sample of size n.
#'
#' @param population_mean The population mean (numeric).
#' @param population_variance The population variance (numeric).
#' @param sample_size The size of the sample (numeric).
#'
#' @return A list containing:
#'   - `sampling_mean`: The mean of the sampling distribution.
#'   - `standard_error`: The standard error of the sampling distribution.
#'
#' @examples
#' sampling_distribution_mean_and_se(67, 36, 100)
#'
#' @export
sampling_distribution_mean_and_se <- function(population_mean, population_variance, sample_size) {
  sampling_mean <- population_mean
  standard_error <- sqrt(population_variance / sample_size)

  return(list(sampling_mean = sampling_mean, standard_error = standard_error))
}

#' Find P(Xbar < bound)
#'
#' This function calculates the probability that the sample mean Xbar is less than a specified bound,
#' assuming the sample mean follows a normal distribution.
#'
#' @param bound The value of the bound for Xbar (numeric).
#' @param population_mean The population mean (numeric).
#' @param population_variance The population variance (numeric).
#' @param sample_size The sample size (numeric).
#'
#' @return The probability P(Xbar) < bound).
#' @examples
#' prob_Xbar_less_than_bound(70, 67, 36, 100)
#'
#' @export
prob_Xbar_less_than_bound <- function(bound, population_mean, population_variance, sample_size) {
  standard_error <- sqrt(population_variance / sample_size)
  z_score <- (bound - population_mean) / standard_error
  pnorm(z_score)
}

#' Find P(Xbar > bound)
#'
#' This function calculates the probability that the sample mean Xbar is greater than a specified bound,
#' assuming the sample mean follows a normal distribution.
#'
#' @param bound The value of the bound for Xbar (numeric).
#' @param population_mean The population mean (numeric).
#' @param population_variance The population variance (numeric).
#' @param sample_size The sample size (numeric).
#'
#' @return The probability P(Xbar) > bound).
#'
#' @examples
#' prob_Xbar_greater_than_bound(70, 67, 36, 100)
#'
#' @export
prob_Xbar_greater_than_bound <- function(bound, population_mean, population_variance, sample_size) {
  standard_error <- sqrt(population_variance / sample_size)
  z_score <- (bound - population_mean) / standard_error
  1 - pnorm(z_score)
}

#' Find P(Xbar <equals bound)
#'
#' This function calculates the probability that the sample mean Xbar is less than or equal to a specified bound,
#' assuming the sample mean follows a normal distribution.
#'
#' @param bound The value of the bound for Xbar (numeric).
#' @param population_mean The population mean (numeric).
#' @param population_variance The population variance (numeric).
#' @param sample_size The sample size (numeric).
#'
#' @return The probability P(Xbar <equals bound).
#'
#' @examples
#' prob_Xbar_less_than_or_equals_bound(70, 67, 36, 100)
#'
#' @export
prob_Xbar_less_than_or_equals_bound <- function(bound, population_mean, population_variance, sample_size) {
  standard_error <- sqrt(population_variance / sample_size)
  z_score <- (bound - population_mean) / standard_error
  pnorm(z_score)
}

#' Find P(Xbar >equals bound)
#'
#' This function calculates the probability that the sample mean Xbar is greater than or equal to a specified bound,
#' assuming the sample mean follows a normal distribution.
#'
#' @param bound The value of the bound for Xbar (numeric).
#' @param population_mean The population mean (numeric).
#' @param population_variance The population variance (numeric).
#' @param sample_size The sample size (numeric).
#'
#' @return The probability P(Xbar >equals bound).
#'
#' @examples
#' prob_Xbar_greater_than_or_equals_bound(70, 67, 36, 100)
#'
#' @export
prob_Xbar_greater_than_or_equals_bound <- function(bound, population_mean, population_variance, sample_size) {
  standard_error <- sqrt(population_variance / sample_size)
  z_score <- (bound - population_mean) / standard_error
  1 - pnorm(z_score)
}

#' Find P(lower <equals Xbar <equals upper)
#'
#' This function calculates the probability that the sample mean (Xbar) lies between two bounds, lower and upper,
#' assuming the sample mean follows a normal distribution.
#'
#' @param lower The lower bound for (Xbar) (numeric).
#' @param upper The upper bound for (Xbar) (numeric).
#' @param population_mean The population mean (numeric).
#' @param population_variance The population variance (numeric).
#' @param sample_size The sample size (numeric).
#'
#' @return The probability P(lower <equals (Xbar) <equals upper).
#'
#' @examples
#' prob_Xbar_between_bounds(45, 74, 67, 36, 100)
#'
#' @export
prob_Xbar_between_bounds <- function(lower, upper, population_mean, population_variance, sample_size) {
  standard_error <- sqrt(population_variance / sample_size)
  z_lower <- (lower - population_mean) / standard_error
  z_upper <- (upper - population_mean) / standard_error
  pnorm(z_upper) - pnorm(z_lower)
}

#' Compare the Distributions for X and Xbar
#'
#' This function compares the distributions for X (a random variable) and Xbar (the sample mean)
#' by calculating the mean and standard deviation for both. It also returns the difference in the standard deviations.
#' Additionally, the function returns an explanation of the key differences between the two distributions.
#'
#' @param population_mean The population mean (numeric).
#' @param population_sd The population standard deviation (numeric).
#' @param sample_size The sample size (numeric).
#'
#' @return A list containing:
#'   - `mean_X`: The mean of the distribution for X.
#'   - `sd_X`: The standard deviation of the distribution for X.
#'   - `mean_Xbar`: The mean of the distribution for Xbar.
#'   - `sd_Xbar`: The standard deviation of the distribution for Xbar.
#'   - `difference_in_sd`: The difference in the standard deviations between X and Xbar.
#'   - `explanation`: A string explanation of the key differences between the two distributions.
#'
#' @examples
#' compare_distributions(67, 6, 30)
#'
#' @export
compare_distributions <- function(population_mean, population_sd, sample_size) {
  # Mean of X is the population mean
  mean_X <- population_mean

  # Standard deviation of X is the population standard deviation
  sd_X <- population_sd

  # Mean of Xbar is also the population mean
  mean_Xbar <- population_mean

  # Standard deviation of Xbar (standard error) is the population standard deviation divided by the square root of the sample size
  sd_Xbar <- population_sd / sqrt(sample_size)

  # Calculate the difference in standard deviations
  difference_in_sd <- sd_X - sd_Xbar

  # Explanation of key differences
  explanation <- paste(
    "Key differences between the distribution for X and the distribution for Xbar:\n",
    "1. The mean of X is the population mean, and the mean of Xbar is also the population mean.\n",
    "2. The standard deviation of X is the population standard deviation.\n",
    "3. The standard deviation of Xbar (also called the standard error) is smaller than the standard deviation of X.\n",
    "4. The standard deviation of Xbar is equal to the population standard deviation divided by the square root of the sample size.\n",
    "5. As the sample size increases, the standard deviation of Xbar decreases, making the sampling distribution narrower.\n"
  )

  # Return the results as a list
  return(list(mean_X = mean_X, sd_X = sd_X, mean_Xbar = mean_Xbar, sd_Xbar = sd_Xbar,
              difference_in_sd = difference_in_sd, explanation = explanation))
}


#' Calculate the Confidence Interval for a New Confidence Level with Detailed Information
#'
#' This function calculates the confidence interval for a new confidence level based on a given original
#' confidence interval, sample size, and critical t-values for both the original and new confidence levels.
#' It uses the t-distribution to calculate the new margin of error and confidence interval.
#'
#' @param lower_original The lower bound of the original confidence interval.
#' @param upper_original The upper bound of the original confidence interval.
#' @param n The sample size.
#' @param confidence_original The original confidence level (e.g., 0.90 for 90%).
#' @param confidence_new The new confidence level (e.g., 0.99 for 99%).
#' @return A list containing the following:
#'   - `new_confidence_interval`: The new confidence interval as a vector.
#'   - `sample_mean`: The sample mean calculated from the original confidence interval.
#'   - `margin_error_original`: The margin of error for the original confidence interval.
#'   - `t_original`: The t-value for the original confidence level.
#'   - `t_new`: The t-value for the new confidence level.
#'   - `sample_standard_deviation`: The sample standard deviation.
#'   - `margin_error_new`: The margin of error for the new confidence interval.
#'
#' @examples
#' calculate_new_ci(15.34, 36.66, 5, 0.90, 0.99)
#' @export
calculate_new_ci <- function(lower_original, upper_original, n, confidence_original, confidence_new) {

  # Calculate the sample mean from the original confidence interval
  sample_mean <- (lower_original + upper_original) / 2

  # Calculate the margin of error from the original confidence interval
  margin_error_original <- (upper_original - lower_original) / 2

  # Calculate sample standard deviation (s) based on the margin of error and sample size
  t_original <- qt(1 - (1 - confidence_original) / 2, df = n - 1)  # t value for the original confidence level
  s <- margin_error_original * sqrt(n) / t_original

  # Calculate the margin of error for the new confidence interval
  t_new <- qt(1 - (1 - confidence_new) / 2, df = n - 1)  # t value for the new confidence level
  margin_error_new <- t_new * (s / sqrt(n))

  # Calculate the new confidence interval
  lower_new <- sample_mean - margin_error_new
  upper_new <- sample_mean + margin_error_new

  # Return the detailed results as a list
  return(list(
    new_confidence_interval = c(lower_new, upper_new),
    sample_mean = sample_mean,
    margin_error_original = margin_error_original,
    t_original = t_original,
    t_new = t_new,
    sample_standard_deviation = s,
    margin_error_new = margin_error_new
  ))
}

#' Explanation of Type I and Type II Errors
#'
#' This function explains the concept of Type I and Type II errors in the context of hypothesis testing.
#' It provides an interpretation based on a given scenario, where the null hypothesis and alternative hypothesis are defined.
#'
#' @param null_hypothesis The null hypothesis.
#' @param alternative_hypothesis The alternative hypothesis.
#' @param error_type The type of error to explain: "Type I" or "Type II".
#' @return A string explaining the selected type of error.
#'
#' @examples
#' error_explanation(
#'   "The team will not get the first down",
#'   "The team will get the first down",
#'   "Type I"
#' )
#' @export
error_explanation <- function(null_hypothesis, alternative_hypothesis, error_type) {
  if (error_type == "Type I") {
    return(paste("Type I Error (False Positive): You reject the null hypothesis ('", null_hypothesis,
                 "') when in fact it is true. In this case, you would decide to go for the first down, but the team will not succeed.", sep = ""))
  } else if (error_type == "Type II") {
    return(paste("Type II Error (False Negative): You fail to reject the null hypothesis ('", null_hypothesis,
                 "') when in fact it is false. In this case, you would decide to punt, but the team would have succeeded in getting the first down.", sep = ""))
  } else {
    return("Invalid error type. Please choose either 'Type I' or 'Type II'.")
  }
}

#' Perform a Two-Tailed t-Test
#'
#' This function performs a two-tailed hypothesis test to determine if the mean gas mileage of a sample is
#' different from the published mean.
#'
#' @param sample_mean The sample mean (e.g., 31.6 mpg).
#' @param population_mean The population mean (e.g., 33.5 mpg).
#' @param sample_sd The sample standard deviation (e.g., 3.4 mpg).
#' @param sample_size The sample size (e.g., 12).
#' @param alpha The significance level (default is 0.05).
#'
#' @return A list containing the t statistic, p-value, and decision ("Reject H0" or "Fail to reject H0").
#'
#' @importFrom stats pt
#' @examples
#' two_tailed__test(31.6, 33.5, 3.4, 12, 0.05)
#' @export
two_tailed__test <- function(sample_mean, population_mean, sample_sd, sample_size, alpha = 0.05) {
  # Calculate the standard error
  standard_error <- sample_sd / sqrt(sample_size)

  # Calculate the t statistic
  t_statistic <- (sample_mean - population_mean) / standard_error

  # Calculate the degrees of freedom
  df <- sample_size - 1

  # Two-tailed p-value calculation
  p_value <- 2 * pt(abs(t_statistic), df, lower.tail = FALSE)  # Multiply by 2 for two-tailed test

  # Make the decision
  if (p_value < alpha) {
    decision <- "Reject H0"
  } else {
    decision <- "Fail to reject H0"
  }

  # Return the result
  return(list(t_statistic = t_statistic, p_value = p_value, decision = decision))
}

#' Calculate the required sample size for a 90% confidence interval with a specified margin of error
#'
#' This function calculates the required sample size for a confidence interval of a population proportion
#' with a given margin of error and confidence level.
#'
#' @param margin_of_error The desired margin of error for the confidence interval (numeric).
#' @param confidence_level The desired confidence level (numeric, between 0 and 1).
#' @param p_estimate The estimated proportion (default is 0.5 for a fair coin).
#'
#' @return The required sample size (numeric).
#'
#' @examples
#' required_sample_size(0.1, 0.9)
#' @export
required_sample_size <- function(margin_of_error, confidence_level, p_estimate = 0.5) {
  z_score <- qnorm(1 - (1 - confidence_level) / 2)
  n <- (z_score^2 * p_estimate * (1 - p_estimate)) / (margin_of_error^2)
  return(ceiling(n))
}

#' Perform a one-tailed hypothesis test for a population proportion
#'
#' This function tests whether the proportion of people who pray for world peace is less than
#' 80% based on a sample proportion and a specified significance level.
#'
#' @param successes The number of people in the sample who pray for world peace.
#' @param sample_size The total number of people in the sample.
#' @param hypothesized_proportion The hypothesized population proportion (default is 0.80).
#' @param alpha The significance level (default is 0.10).
#'
#' @return A list containing:
#'   - `z_value`: The calculated z statistic (numeric).
#'   - `critical_value`: The critical z value for the one-tailed test (numeric).
#'   - `decision`: A string indicating whether to "Reject H0" or "Fail to Reject H0".
#'
#' @examples
#' hypothesis_test_proportion_2(77, 110, 0.80, 0.10)
#' @export
hypothesis_test_proportion_2 <- function(successes, sample_size, hypothesized_proportion = 0.80, alpha = 0.10) {
  # Calculate the sample proportion
  p_hat <- successes / sample_size

  # Calculate the standard error
  standard_error <- sqrt((hypothesized_proportion * (1 - hypothesized_proportion)) / sample_size)

  # Calculate the z statistic
  z_value <- (p_hat - hypothesized_proportion) / standard_error

  # Critical z value for a one-tailed test at alpha significance level
  critical_value <- qnorm(1 - alpha)

  # Decision based on the z value
  if (z_value < critical_value) {
    decision <- "Reject H0"
  } else {
    decision <- "Fail to Reject H0"
  }

  # Return results as a list
  return(list(z_value = z_value, critical_value = critical_value, decision = decision))
}

#' Construct a Confidence Interval for the Variance of a Normally Distributed Population
#'
#' This function constructs a 95% confidence interval for the variance of a normally distributed population based on a sample of data.
#' It uses the chi-squared distribution to calculate the confidence interval for the population variance (sigma^2) and compares the
#' sample variance to the claimed population variance.
#'
#' The formula for the confidence interval is based on the chi-squared distribution:
#'
#' The confidence interval for the variance is given by:
#'
#' ( (n-1) * s^2 ) / chi_squared_critical_upper <= sigma^2 <= ( (n-1) * s^2 ) / chi_squared_critical_lower
#'
#' Where:
#' - s^2 is the sample variance
#' - n is the sample size
#' - chi_squared_critical_upper is the upper chi-squared critical value for (1 - alpha/2) with n-1 degrees of freedom
#' - chi_squared_critical_lower is the lower chi-squared critical value for alpha/2 with n-1 degrees of freedom
#'
#' @param data A numeric vector containing the sample data.
#' @param alpha The significance level (default is 0.05 for a 95% confidence interval).
#' @return A list containing:
#'   - lower_bound: The lower bound of the confidence interval for the population variance.
#'   - upper_bound: The upper bound of the confidence interval for the population variance.
#'   - sample_variance: The sample variance.
#'   - sample_size: The sample size.
#'   - degrees_of_freedom: The degrees of freedom used in the chi-squared distribution.
#'
#' @importFrom stats qt qchisq var
#' @examples
#' battery_lifetimes <- c(1.9, 2.4, 3.0, 3.5, 4.2)
#' result <- confidence_interval_variance(battery_lifetimes)
#' print(result)
#'
#' @export
confidence_interval_variance <- function(data, alpha = 0.05) {
  # Sample size and variance
  n <- length(data)
  sample_variance <- var(data)

  # Degrees of freedom
  degrees_of_freedom <- n - 1

  # Chi-squared critical values
  chi_squared_critical_lower <- qchisq(alpha / 2, degrees_of_freedom)
  chi_squared_critical_upper <- qchisq(1 - alpha / 2, degrees_of_freedom)

  # Confidence interval for the population variance
  lower_bound <- (degrees_of_freedom * sample_variance) / chi_squared_critical_upper
  upper_bound <- (degrees_of_freedom * sample_variance) / chi_squared_critical_lower

  # Return the results as a list
  return(list(
    lower_bound = lower_bound,
    upper_bound = upper_bound,
    sample_variance = sample_variance,
    sample_size = n,
    degrees_of_freedom = degrees_of_freedom
  ))
}

#' Perform a Paired T-Test for MPG Differences
#'
#' This function conducts a paired t-test to test whether there is a significant difference
#' between the professor's calculated miles per gallon (MPG) and the car's computer estimate.
#' It performs a two-tailed hypothesis test for the mean difference.
#'
#' The hypotheses are:
#' - H0: The mean difference in MPG = 0 (no difference)
#' - Ha: The mean difference in MPG <> 0 (there is a difference)
#'
#' @param computer A numeric vector representing the MPG values from the car's computer.
#' @param driver A numeric vector representing the MPG values calculated by the professor.
#' @return A list containing:
#'   - `t_statistic`: The calculated t statistic.
#'   - `p_value`: The p-value for the test.
#'   - `decision`: A string indicating whether to "Reject H0" or "Fail to reject H0".
#'
#' @examples
#' computer <- c(41.5, 45, 43.2, 43.2, 48.4, 46.8, 39.2, 43.5, 44.3, 43.3)
#' driver <- c(36.5, 40.5, 41, 38.8, 45.4, 45.7, 34.2, 39.8, 44.9, 47.5)
#' result <- paired_t_test_mpg(computer, driver)
#' print(result)
#'
#' @export
paired_t_test_mpg <- function(computer, driver) {
  # Check if the input vectors have the same length
  if (length(computer) != length(driver)) {
    stop("The computer and driver data must have the same length.")
  }

  # Calculate the differences
  differences <- computer - driver

  # Perform the t-test
  t_test_result <- t.test(differences)

  # Extract the t-statistic and p-value
  t_statistic <- t_test_result$statistic
  p_value <- t_test_result$p.value

  # Make a decision based on the p-value
  decision <- ifelse(p_value < 0.05, "Reject H0", "Fail to reject H0")

  # Return the results
  return(list(
    t_statistic = t_statistic,
    p_value = p_value,
    decision = decision
  ))
}

#' Compare Two Proportions
#'
#' This function performs a hypothesis test and calculates a confidence interval for the difference between
#' two population proportions (e.g., proportion of freshmen and seniors opposing a plan). It also determines
#' whether to reject or fail to reject the null hypothesis and provides a simple explanation of the result.
#'
#' @param p1 The proportion of freshmen who oppose the plan.
#' @param n1 The sample size of freshmen.
#' @param p2 The proportion of seniors who oppose the plan.
#' @param n2 The sample size of seniors.
#' @param confidence_level The confidence level for the confidence interval (e.g., 0.95 for 95% confidence).
#' @param alpha The significance level (e.g., 0.05 for a 5% significance level).
#'
#' @return A list containing:
#'   - `test_statistic`: The z-test statistic for comparing the two proportions.
#'   - `p_value`: The p-value of the hypothesis test.
#'   - `reject_H0`: A message indicating whether to reject or fail to reject the null hypothesis.
#'   - `confidence_interval`: The confidence interval for the difference in proportions.
#'   - `explanation`: A simple explanation of the result.
#'
#' @examples
#' result <- compare_two_proportions(160/200, 200, 20/100, 100, 0.95, 0.05)
#' print(result)
#'
#' @export
compare_two_proportions <- function(p1, n1, p2, n2, confidence_level = 0.95, alpha = 0.05) {

  # Step 1: Hypothesis Test
  # Null hypothesis: p1 = p2
  # Alternative hypothesis: p1 != p2

  # Calculate the pooled proportion
  pooled_proportion <- (p1 * n1 + p2 * n2) / (n1 + n2)

  # Standard error of the difference in proportions
  standard_error <- sqrt(pooled_proportion * (1 - pooled_proportion) * (1/n1 + 1/n2))

  # Test statistic (z-value)
  test_statistic <- (p1 - p2) / standard_error

  # Step 2: Calculate p-value
  p_value <- 2 * (1 - pnorm(abs(test_statistic)))  # two-tailed test

  # Step 3: Determine whether to reject the null hypothesis
  reject_H0 <- ifelse(p_value < alpha, "Reject H0", "Fail to reject H0")

  # Step 4: Confidence Interval for the difference in proportions
  margin_of_error <- qnorm(1 - (1 - confidence_level) / 2) * standard_error
  lower_bound <- (p1 - p2) - margin_of_error
  upper_bound <- (p1 - p2) + margin_of_error
  confidence_interval <- c(lower_bound, upper_bound)

  # Step 5: Explanation for a non-statistical audience
  explanation <- ifelse(reject_H0 == "Reject H0",
                        "There is a statistically significant difference in the proportions of freshmen and seniors who oppose the plan.",
                        "There is no statistically significant difference in the proportions of freshmen and seniors who oppose the plan.")

  # Return the results as a list
  return(list(
    test_statistic = test_statistic,
    p_value = p_value,
    reject_H0 = reject_H0,
    confidence_interval = confidence_interval,
    explanation = explanation
  ))
}

#' Perform a Two-Sample T-Test for the Difference in Means
#'
#' This function performs a two-sample t-test to determine if the mean of one group is greater than the mean of another group.
#' The function calculates the test statistic, the p-value, and provides a conclusion based on the significance level.
#'
#' @param mean_1 The mean of the first sample (numeric).
#' @param sd_1 The standard deviation of the first sample (numeric).
#' @param n_1 The sample size of the first sample (integer).
#' @param mean_2 The mean of the second sample (numeric).
#' @param sd_2 The standard deviation of the second sample (numeric).
#' @param n_2 The sample size of the second sample (integer).
#' @param alpha The significance level for the test (default is 0.05).
#'
#' @return A list containing:
#'   - `t_statistic`: The t-test statistic (numeric).
#'   - `p_value`: The p-value for the hypothesis test.
#'   - `reject_H0`: A string indicating whether to reject or fail to reject the null hypothesis.
#'
#' @importFrom stats pt
#' @examples
#' result <- two_sample_t_test(10.40, 4.83, 97, 9.26, 4.68, 148, 0.05)
#' print(result)
#'
#' @export
two_sample_t_test <- function(mean_1, sd_1, n_1,
                              mean_2, sd_2, n_2,
                              alpha = 0.05) {

  # Calculate the standard error of the difference in means
  se_diff <- sqrt((sd_1^2 / n_1) + (sd_2^2 / n_2))

  # Calculate the t-test statistic for a one-tailed test (first group > second group)
  t_statistic <- (mean_1 - mean_2) / se_diff

  # Calculate the degrees of freedom using the formula for unequal variances
  df <- ((sd_1^2 / n_1 + sd_2^2 / n_2)^2) /
    (( (sd_1^2 / n_1)^2 / (n_1 - 1)) +
       ((sd_2^2 / n_2)^2 / (n_2 - 1)))

  # Calculate the p-value for the one-tailed test
  p_value <- 1 - pt(t_statistic, df)  # For a one-tailed test

  # Make a decision based on the p-value
  reject_H0 <- ifelse(p_value < alpha, "Reject H0: The mean of group 1 is greater than group 2",
                      "Fail to reject H0: No significant difference between the means")

  # Return the results as a list
  return(list(t_statistic = t_statistic,
              p_value = p_value,
              reject_H0 = reject_H0))
}
#' @title Hypothesis Test for Candidate Approval Rating
#' @description Performs a hypothesis test to determine if Candidate A's approval rating is greater than 50%.
#' @param n Sample size, the number of voters surveyed.
#' @param x The number of voters who said they might vote for Candidate A.
#' @param p0 The hypothesized population proportion (0.50 in this case).
#' @param alpha The significance level for the test (default is 0.05).
#' @return A list containing the detailed steps, including the assumptions, critical z-value, test statistic,
#' p-value, decision, and conclusion, with all work shown.
#' @details
#' This function calculates the test statistic and p-value for a one-sample proportion hypothesis test.
#' It checks the assumptions, sets up the null and alternative hypotheses, determines the rejection region,
#' and interprets the results based on the significance level.
#' @examples
#' hypothesis_test_approval(n = 100, x = 54, p0 = 0.50, alpha = 0.05)
#' @export
hypothesis_test_approval <- function(n, x, p0 = 0.50, alpha = 0.05) {
  # Step (a): Check the assumption
  np0 <- n * p0
  n1_minus_p0 <- n * (1 - p0)
  assumption_check <- np0 >= 10 && n1_minus_p0 >= 10
  assumption_details <- paste(
    "Check np0 >= 10 and n(1 - p0) >= 10:",
    "np0 =", np0,
    ", n(1 - p0) =", n1_minus_p0,
    "=> Normal approximation is valid:", assumption_check
  )

  # Step (b): State the hypotheses
  hypotheses <- list(
    H0 = "p = 0.50 (Candidate A's approval rating is 50%)",
    Ha = "p > 0.50 (Candidate A's approval rating is greater than 50%)"
  )

  # Step (c): Determine the rejection region
  z_critical <- qnorm(1 - alpha)
  rejection_region <- paste("Reject H0 if z > z_critical, where z_critical =", round(z_critical, 4))

  # Step (d): Calculate the test statistic
  p_hat <- x / n
  standard_error <- sqrt((p0 * (1 - p0)) / n)
  z <- (p_hat - p0) / standard_error
  test_statistic_details <- paste(
    "Calculate z = (p_hat - p0) / SE:",
    "p_hat =", p_hat,
    ", p0 =", p0,
    ", SE =", round(standard_error, 4),
    ", z =", round(z, 4)
  )

  # Step (e): Obtain the p-value
  p_value <- 1 - pnorm(z)
  p_value_details <- paste("One-tailed p-value for z =", round(z, 4), "=> p-value =", round(p_value, 4))

  # Step (f): Decision and Conclusion
  decision <- if (z > z_critical) "Reject H0" else "Do not reject H0"
  conclusion <- if (z > z_critical) {
    "There is significant evidence that Candidate A's approval rating is greater than 50%."
  } else {
    "There is not enough evidence to conclude that Candidate A's approval rating is greater than 50%."
  }

  # Return all details in the output
  return(list(
    assumption_details = assumption_details,
    hypotheses = hypotheses,
    rejection_region = rejection_region,
    test_statistic_details = test_statistic_details,
    p_value_details = p_value_details,
    decision = decision,
    conclusion = conclusion
  ))
}

#' @title Hypothesis Test with Rejection Region for Average Radiation Level
#' @description Performs a hypothesis test to determine if the average radiation level in a laboratory differs from a specified mean (e.g., 420).
#' @param n Sample size, the number of radiation level measurements.
#' @param sample_mean The sample mean of radiation levels.
#' @param mu The hypothesized population mean (default is 420).
#' @param sigma The population standard deviation of the radiation levels (default is 10).
#' @param alpha The significance level for the test (default is 0.01).
#' @return A list containing the test type, rejection region, test statistic, p-value, decision, and conclusion, with all work shown.
#' @details
#' This function calculates the test statistic and p-value for a two-tailed hypothesis test on the average radiation level.
#' It provides step-by-step output including the rejection region, test statistic, and decision criteria.
#' @examples
#' hypothesis_test_with_rejection_region(n = 49, sample_mean = 415.7, mu = 420, sigma = 10, alpha = 0.01)
#' @export
hypothesis_test_with_rejection_region <- function(n, sample_mean, mu = 420, sigma = 10, alpha = 0.01) {
  # Step (a): Determine the type of test
  test_type <- "Two-tailed test, since Ha: μ ≠ 420"

  # Step (b): Determine the rejection region
  z_critical <- qnorm(1 - alpha / 2)  # two-tailed test, so we use alpha/2
  rejection_region <- paste("Reject H0 if z < -", round(z_critical, 4), " or z > ", round(z_critical, 4))

  # Step (c): Calculate the test statistic
  standard_error <- sigma / sqrt(n)
  z <- (sample_mean - mu) / standard_error
  test_statistic_details <- paste(
    "Calculate z = (sample_mean - mu) / SE:",
    "sample_mean =", sample_mean,
    ", mu =", mu,
    ", SE =", round(standard_error, 4),
    ", z =", round(z, 4)
  )

  # Step (d): Obtain the p-value
  p_value <- 2 * (1 - pnorm(abs(z)))  # two-tailed p-value
  p_value_details <- paste("Two-tailed p-value for z =", round(z, 4), "=> p-value =", round(p_value, 4))

  # Step (c): Decision based on rejection region
  decision <- if (abs(z) > z_critical) "Reject H0" else "Do not reject H0"
  conclusion <- if (abs(z) > z_critical) {
    "There is significant evidence that the average radiation level differs from 420."
  } else {
    "There is not enough evidence to conclude that the average radiation level differs from 420."
  }

  # Return all details in the output
  return(list(
    test_type = test_type,
    rejection_region = rejection_region,
    test_statistic_details = test_statistic_details,
    p_value_details = p_value_details,
    decision = decision,
    conclusion = conclusion
  ))
}

#' @title Generic Hypothesis Test for Population Mean
#' @description Performs a hypothesis test to determine if the population mean differs from a specified value.
#' @param data A numeric vector containing the sample data.
#' @param mu The hypothesized population mean to test against.
#' @param alpha The significance level for the test (default is 0.05).
#' @return A list containing the test type, test statistic, p-value, confidence interval, decision, and conclusion, with all work shown.
#' @details
#' This function calculates the test statistic, p-value, and confidence interval for a two-tailed hypothesis test on the population mean.
#' It assumes the data are normally distributed and uses a t-test when the population standard deviation is unknown.
#' @examples
#' data <- c(34,54,73,38,89,52,75,33,50,39,42,42,40,66,72,85,28,71,52,47,41,36,33,38,49,51,55,63,72,78)
#' hypothesis_test_mean(data = data, mu = 50, alpha = 0.05)
#' @export
hypothesis_test_mean <- function(data, mu, alpha = 0.05) {
  # Step (a): Determine the type of test
  test_type <- "Two-tailed test, since Ha: μ ≠ hypothesized mean"

  # Step (b): Calculate sample statistics
  n <- length(data)
  sample_mean <- mean(data)
  sample_sd <- sd(data)
  standard_error <- sample_sd / sqrt(n)

  # Step (c): Calculate the test statistic
  t_stat <- (sample_mean - mu) / standard_error
  test_statistic_details <- paste(
    "Calculate t = (sample_mean - mu) / SE:",
    "sample_mean =", round(sample_mean, 4),
    ", mu =", mu,
    ", SE =", round(standard_error, 4),
    ", t =", round(t_stat, 4)
  )

  # Step (d): Determine the rejection region
  t_critical <- qt(1 - alpha / 2, df = n - 1)
  rejection_region <- paste("Reject H0 if t < -", round(t_critical, 4), " or t > ", round(t_critical, 4))

  # Step (e): Obtain the p-value
  p_value <- 2 * (1 - pt(abs(t_stat), df = n - 1))  # two-tailed p-value
  p_value_details <- paste("Two-tailed p-value for t =", round(t_stat, 4), "=> p-value =", round(p_value, 4))

  # Step (f): Decision and Conclusion
  decision <- if (abs(t_stat) > t_critical) "Reject H0" else "Do not reject H0"
  conclusion <- if (abs(t_stat) > t_critical) {
    "There is significant evidence that the population mean differs from the hypothesized mean."
  } else {
    "There is not enough evidence to conclude that the population mean differs from the hypothesized mean."
  }

  # Step (g): Calculate the confidence interval
  margin_of_error <- t_critical * standard_error
  ci_lower <- sample_mean - margin_of_error
  ci_upper <- sample_mean + margin_of_error
  confidence_interval <- paste("Confidence interval: (", round(ci_lower, 4), ", ", round(ci_upper, 4), ")")

  # Return all details in the output
  return(list(
    test_type = test_type,
    test_statistic_details = test_statistic_details,
    rejection_region = rejection_region,
    p_value_details = p_value_details,
    decision = decision,
    conclusion = conclusion,
    confidence_interval = confidence_interval
  ))
}

#' @title Hypothesis Test for Population Proportion
#' @description Performs a hypothesis test to determine if the population proportion differs from a specified value.
#' @param n Sample size, the number of patients or cases.
#' @param x The number of patients or cases with the outcome of interest (e.g., adverse symptoms).
#' @param p0 The hypothesized population proportion (e.g., 0.10 for 10%).
#' @param alpha The significance level for the test (default is 0.05).
#' @return A list containing the sampling distribution, test statistic, p-value, rejection region, decision, and conclusion, with all work shown.
#' @details
#' This function checks assumptions, calculates the test statistic, p-value, and provides a one-tailed or two-tailed hypothesis test on the population proportion.
#' @examples
#' hypothesis_test_proportion(n = 440, x = 23, p0 = 0.10, alpha = 0.05)
#' @export
hypothesis_test_proportion <- function(n, x, p0, alpha = 0.05) {
  # Step (a): Sampling distribution of the sample proportion
  sample_proportion <- x / n
  standard_error <- sqrt((p0 * (1 - p0)) / n)
  sampling_distribution <- paste("Sampling distribution: N(mean =", round(p0, 4),
                                 ", standard deviation =", round(standard_error, 4), ")")

  # Step (b1): Assumption check
  np0 <- n * p0
  n1_minus_p0 <- n * (1 - p0)
  assumption_check <- np0 >= 10 && n1_minus_p0 >= 10
  assumption_details <- paste(
    "Check np0 >= 10 and n(1 - p0) >= 10:",
    "np0 =", np0,
    ", n(1 - p0) =", n1_minus_p0,
    "=> Normal approximation is valid:", assumption_check
  )

  # Step (b2): State the hypotheses
  hypotheses <- list(
    H0 = paste("p =", p0, "(The proportion of patients with adverse symptoms is", p0, ")"),
    Ha = paste("p <", p0, "(The proportion of patients with adverse symptoms is less than", p0, ")")
  )

  # Step (c): Determine the rejection region for a left-tailed test
  z_critical <- qnorm(alpha)
  rejection_region <- paste("Reject H0 if z <", round(z_critical, 4))

  # Step (d): Calculate the test statistic
  z <- (sample_proportion - p0) / standard_error
  test_statistic_details <- paste(
    "Calculate z = (sample_proportion - p0) / SE:",
    "sample_proportion =", round(sample_proportion, 4),
    ", p0 =", p0,
    ", SE =", round(standard_error, 4),
    ", z =", round(z, 4)
  )

  # Step (e): Obtain the p-value for a one-tailed test
  p_value <- pnorm(z)
  p_value_details <- paste("One-tailed p-value for z =", round(z, 4), "=> p-value =", round(p_value, 4))

  # Step (f): Decision and Conclusion
  decision <- if (z < z_critical) "Reject H0" else "Do not reject H0"
  conclusion <- if (z < z_critical) {
    "There is strong evidence that fewer than 10% of patients taking this medication have adverse symptoms."
  } else {
    "There is not enough evidence to conclude that fewer than 10% of patients taking this medication have adverse symptoms."
  }

  # Return all details in the output
  return(list(
    sampling_distribution = sampling_distribution,
    assumption_details = assumption_details,
    hypotheses = hypotheses,
    rejection_region = rejection_region,
    test_statistic_details = test_statistic_details,
    p_value_details = p_value_details,
    decision = decision,
    conclusion = conclusion
  ))
}

#' @title Hypothesis Test for Population Proportion (Right-Tailed)
#' @description Performs a hypothesis test to determine if the population proportion is greater than a specified value.
#' @param n Sample size, the number of cases (e.g., voters).
#' @param x The number of cases with the outcome of interest (e.g., voters supporting a candidate).
#' @param p0 The hypothesized population proportion (e.g., 0.50 for 50%).
#' @param alpha The significance level for the test (default is 0.05).
#' @return A list containing the sampling distribution, assumptions check, hypotheses, test statistic, p-value, rejection region, decision, and conclusion.
#' @details
#' This function checks assumptions, calculates the test statistic, and provides a one-tailed hypothesis test on the population proportion.
#' @examples
#' hypothesis_test_proportion_right_tailed(n = 100, x = 54, p0 = 0.50, alpha = 0.05)
#' @export
hypothesis_test_proportion_right_tailed <- function(n, x, p0, alpha = 0.05) {
  # Step (a): Check the assumptions
  sample_proportion <- x / n
  np0 <- n * p0
  n1_minus_p0 <- n * (1 - p0)
  assumption_check <- np0 >= 10 && n1_minus_p0 >= 10
  assumption_details <- paste(
    "Check np0 >= 10 and n(1 - p0) >= 10:",
    "np0 =", np0,
    ", n(1 - p0) =", n1_minus_p0,
    "=> Normal approximation is valid:", assumption_check
  )

  # Step (b): State the hypotheses
  hypotheses <- list(
    H0 = paste("p =", p0, "(The proportion is", p0, ")"),
    Ha = paste("p >", p0, "(The proportion is greater than", p0, ")")
  )

  # Step (c): Determine the rejection region for a right-tailed test
  z_critical <- qnorm(1 - alpha)
  rejection_region <- paste("Reject H0 if z >", round(z_critical, 4))

  # Step (d): Calculate the test statistic
  standard_error <- sqrt((p0 * (1 - p0)) / n)
  z <- (sample_proportion - p0) / standard_error
  test_statistic_details <- paste(
    "Calculate z = (sample_proportion - p0) / SE:",
    "sample_proportion =", round(sample_proportion, 4),
    ", p0 =", p0,
    ", SE =", round(standard_error, 4),
    ", z =", round(z, 4)
  )

  # Step (e): Obtain the p-value for a right-tailed test
  p_value <- 1 - pnorm(z)
  p_value_details <- paste("One-tailed p-value for z =", round(z, 4), "=> p-value =", round(p_value, 4))

  # Step (f): Decision and Conclusion
  decision <- if (z > z_critical) "Reject H0" else "Do not reject H0"
  conclusion <- if (z > z_critical) {
    "There is significant evidence that the proportion is greater than the hypothesized value."
  } else {
    "There is not enough evidence to conclude that the proportion is greater than the hypothesized value."
  }

  # Return all details in the output
  return(list(
    sampling_distribution = paste("Sampling distribution: N(mean =", p0,
                                  ", standard deviation =", round(standard_error, 4), ")"),
    assumption_details = assumption_details,
    hypotheses = hypotheses,
    rejection_region = rejection_region,
    test_statistic_details = test_statistic_details,
    p_value_details = p_value_details,
    decision = decision,
    conclusion = conclusion
  ))
}

#' @title Hypothesis Test for Population Proportion (Left-Tailed)
#' @description Performs a hypothesis test to determine if the population proportion is less than a specified value.
#' @param n Sample size, the number of cases (e.g., patients).
#' @param x The number of cases with the outcome of interest (e.g., patients with adverse symptoms).
#' @param p0 The hypothesized population proportion (e.g., 0.10 for 10%).
#' @param alpha The significance level for the test (default is 0.05).
#' @return A list containing the sampling distribution, assumptions check, hypotheses, test statistic, p-value, rejection region, decision, and conclusion.
#' @details
#' This function checks assumptions, calculates the test statistic, and performs a one-tailed (left-tailed) hypothesis test on the population proportion.
#' @examples
#' hypothesis_test_proportion_left_tailed(n = 100, x = 54, p0 = 0.50, alpha = 0.05)
#' @export
hypothesis_test_proportion_left_tailed <- function(n, x, p0, alpha = 0.05) {
  # Step (a): Check the assumptions
  sample_proportion <- x / n
  np0 <- n * p0
  n1_minus_p0 <- n * (1 - p0)
  assumption_check <- np0 >= 10 && n1_minus_p0 >= 10
  assumption_details <- paste(
    "Check np0 >= 10 and n(1 - p0) >= 10:",
    "np0 =", np0,
    ", n(1 - p0) =", n1_minus_p0,
    "=> Normal approximation is valid:", assumption_check
  )

  # Step (b): State the hypotheses
  hypotheses <- list(
    H0 = paste("p =", p0, "(The proportion is", p0, ")"),
    Ha = paste("p <", p0, "(The proportion is less than", p0, ")")
  )

  # Step (c): Determine the rejection region for a left-tailed test
  z_critical <- qnorm(alpha)
  rejection_region <- paste("Reject H0 if z <", round(z_critical, 4))

  # Step (d): Calculate the test statistic
  standard_error <- sqrt((p0 * (1 - p0)) / n)
  z <- (sample_proportion - p0) / standard_error
  test_statistic_details <- paste(
    "Calculate z = (sample_proportion - p0) / SE:",
    "sample_proportion =", round(sample_proportion, 4),
    ", p0 =", p0,
    ", SE =", round(standard_error, 4),
    ", z =", round(z, 4)
  )

  # Step (e): Obtain the p-value for a left-tailed test
  p_value <- pnorm(z)
  p_value_details <- paste("One-tailed p-value for z =", round(z, 4), "=> p-value =", round(p_value, 4))

  # Step (f): Decision and Conclusion
  decision <- if (z < z_critical) "Reject H0" else "Do not reject H0"
  conclusion <- if (z < z_critical) {
    "There is significant evidence that the proportion is less than the hypothesized value."
  } else {
    "There is not enough evidence to conclude that the proportion is less than the hypothesized value."
  }

  # Return all details in the output
  return(list(
    sampling_distribution = paste("Sampling distribution: N(mean =", p0,
                                  ", standard deviation =", round(standard_error, 4), ")"),
    assumption_details = assumption_details,
    hypotheses = hypotheses,
    rejection_region = rejection_region,
    test_statistic_details = test_statistic_details,
    p_value_details = p_value_details,
    decision = decision,
    conclusion = conclusion
  ))
}


#' @title Paired t-Test for Mean Difference
#' @description Performs a paired t-test to determine if there is a significant difference between two paired samples.
#' @param sample1 A numeric vector representing the first set of paired observations (e.g., campus bookstore prices).
#' @param sample2 A numeric vector representing the second set of paired observations (e.g., Amazon prices).
#' @param alpha The significance level for the test (default is 0.05).
#' @return A list containing the test type, test statistic, p-value, confidence interval, decision, and conclusion, with all work shown.
#' @details
#' This function calculates the test statistic, p-value, and confidence interval for a paired t-test.
#' It assumes that the data are approximately normally distributed.
#' @examples
#' campus_prices <- c(99.34, 51.53, 20.45, 97.22, 61.89, 58.17, 61.63, 44.63, 96.69, 48.88)
#' amazon_prices <- c(113.94, 61.44, 31.59, 108.29, 78.44, 65.74, 63.49, 40.39, 117.99, 58.94)
#' paired_t_test(sample1 = campus_prices, sample2 = amazon_prices, alpha = 0.05)
#' @export
paired_t_test <- function(sample1, sample2, alpha = 0.05) {
  # Check that the two samples are paired and have the same length
  if (length(sample1) != length(sample2)) {
    stop("The two samples must be of the same length for a paired t-test.")
  }

  # Calculate the differences
  differences <- sample1 - sample2

  # Step (a): State the test type
  test_type <- "Paired t-test to check if the mean difference is significantly different from zero"

  # Step (b): Calculate the sample statistics
  n <- length(differences)
  mean_diff <- mean(differences)
  sd_diff <- sd(differences)
  standard_error <- sd_diff / sqrt(n)

  # Step (c): Calculate the test statistic
  t_stat <- mean_diff / standard_error
  test_statistic_details <- paste(
    "Calculate t = mean_diff / SE:",
    "mean_diff =", round(mean_diff, 4),
    ", SE =", round(standard_error, 4),
    ", t =", round(t_stat, 4)
  )

  # Step (d): Obtain the p-value for a two-tailed test
  p_value <- 2 * (1 - pt(abs(t_stat), df = n - 1))
  p_value_details <- paste("Two-tailed p-value for t =", round(t_stat, 4), "=> p-value =", round(p_value, 4))

  # Step (e): Determine the confidence interval for the mean difference
  t_critical <- qt(1 - alpha / 2, df = n - 1)
  margin_of_error <- t_critical * standard_error
  ci_lower <- mean_diff - margin_of_error
  ci_upper <- mean_diff + margin_of_error
  confidence_interval <- paste("Confidence interval for the mean difference: (",
                               round(ci_lower, 4), ", ", round(ci_upper, 4), ")")

  # Step (f): Decision and Conclusion
  decision <- if (p_value < alpha) "Reject H0" else "Do not reject H0"
  conclusion <- if (p_value < alpha) {
    "There is significant evidence of a difference in the mean of paired samples."
  } else {
    "There is not enough evidence to conclude a difference in the mean of paired samples."
  }

  # Return all details in the output
  return(list(
    test_type = test_type,
    test_statistic_details = test_statistic_details,
    p_value_details = p_value_details,
    confidence_interval = confidence_interval,
    decision = decision,
    conclusion = conclusion
  ))
}


#' @title One-Tailed Two-Sample t-Test for Difference in Means
#' @description Performs a one-tailed t-test to determine if there is a significant difference in the means of two independent samples.
#' @param n1 Sample size for the first sample.
#' @param mean1 Sample mean for the first sample.
#' @param sd1 Sample standard deviation for the first sample.
#' @param n2 Sample size for the second sample.
#' @param mean2 Sample mean for the second sample.
#' @param sd2 Sample standard deviation for the second sample.
#' @param alpha The significance level for the test (default is 0.05).
#' @param direction A string indicating the direction of the test: "greater" if testing if sample 1 is greater than sample 2, or "less" if testing if sample 1 is less than sample 2.
#' @return A list containing the test type, test statistic, p-value, confidence interval, decision, and conclusion, with all work shown.
#' @details
#' This function calculates the test statistic, p-value, and confidence interval for a one-tailed two-sample t-test. It assumes that the samples are independent, normally distributed, and have equal variances.
#' @examples
#' one_tailed_t_test(n1 = 15, mean1 = 540, sd1 = 21, n2 = 15, mean2 = 554, sd2 = 15, alpha = 0.05, direction = "less")
#' @export
one_tailed_t_test <- function(n1, mean1, sd1, n2, mean2, sd2, alpha = 0.05, direction = "greater") {
  # Step (a): State the test type
  test_type <- paste("One-tailed two-sample t-test for difference in means, testing if",
                     ifelse(direction == "greater", "mean1 > mean2", "mean1 < mean2"))

  # Step (b): Calculate pooled standard deviation
  sp <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))

  # Step (c): Calculate the test statistic
  standard_error <- sp * sqrt(1 / n1 + 1 / n2)
  t_stat <- (mean1 - mean2) / standard_error
  test_statistic_details <- paste(
    "Calculate t = (mean1 - mean2) / SE:",
    "mean1 =", mean1,
    ", mean2 =", mean2,
    ", SE =", round(standard_error, 4),
    ", t =", round(t_stat, 4)
  )

  # Step (d): Determine the critical t-value and rejection region for a one-tailed test
  df <- n1 + n2 - 2
  t_critical <- qt(1 - alpha, df)
  rejection_region <- if (direction == "greater") {
    paste("Reject H0 if t >", round(t_critical, 4))
  } else {
    paste("Reject H0 if t <", round(-t_critical, 4))
  }

  # Step (e): Obtain the p-value for the test
  if (direction == "greater") {
    p_value <- 1 - pt(t_stat, df)
  } else {
    p_value <- pt(t_stat, df)
  }
  p_value_details <- paste("One-tailed p-value for t =", round(t_stat, 4), "=> p-value =", round(p_value, 4))

  # Step (f): Decision and Conclusion
  decision <- if ((direction == "greater" && t_stat > t_critical) || (direction == "less" && t_stat < -t_critical)) "Reject H0" else "Do not reject H0"
  conclusion <- if (decision == "Reject H0") {
    "There is significant evidence supporting the hypothesis about the difference in means as specified by the direction of the test."
  } else {
    "There is not enough evidence to support the hypothesis about the difference in means as specified by the direction of the test."
  }

  # Step (g): Calculate the confidence interval for the difference in means
  t_critical_ci <- qt(1 - alpha / 2, df)
  margin_of_error <- t_critical_ci * standard_error
  ci_lower <- (mean1 - mean2) - margin_of_error
  ci_upper <- (mean1 - mean2) + margin_of_error
  confidence_interval <- paste("90% Confidence interval for the difference in means: (",
                               round(ci_lower, 4), ", ", round(ci_upper, 4), ")")

  # Return all details in the output
  return(list(
    test_type = test_type,
    test_statistic_details = test_statistic_details,
    p_value_details = p_value_details,
    rejection_region = rejection_region,
    confidence_interval = confidence_interval,
    decision = decision,
    conclusion = conclusion
  ))
}

#' @title Expected Value of a Transformed Random Variable with Detailed Explanation
#' @description Computes the expected value E[Y] where Y is a linear transformation of X, with specific conditions on the range of X. Returns detailed step-by-step explanation suitable for use in a free response answer.
#' @param a Coefficient of X in the linear transformation Y = a * X + b.
#' @param b Constant term in the linear transformation Y = a * X + b.
#' @param cdf Function defining the cumulative distribution function of X.
#' @param lower_bound Lower bound of X for the transformation (e.g., 0 for X >= 0).
#' @param upper_bound Upper bound of X for the transformation (e.g., 1 for X <= 1).
#' @param cutoff_point Point where Y = 0 if X < cutoff_point (e.g., -1 for X < 0).
#' @return A character vector containing a detailed, step-by-step explanation of the calculations and results for E[Y].
#' @details
#' This function integrates over the specified range of X to compute E[Y] for a linear transformation of X,
#' conditioned on X being within a specified range.
#' @examples
#' # Define the CDF of X: F(x) = (1/2)x + (1/2), for -1 <= x <= 1
#' cdf_x <- function(x) {
#'   result <- numeric(length(x))
#'   result[x < -1] <- 0
#'   result[x > 1] <- 1
#'   result[x >= -1 & x <= 1] <- 0.5 * x[x >= -1 & x <= 1] + 0.5
#'   return(result)
#' }
#' expected_value_y_detailed(a = 8, b = 6, cdf = cdf_x, lower_bound = 0, upper_bound = 1, cutoff_point = 0)
#' @export
expected_value_y_detailed <- function(a, b, cdf, lower_bound, upper_bound, cutoff_point) {

  # Step 1: Calculate the probability P(X >= cutoff_point)
  p_x_positive <- 1 - cdf(cutoff_point)
  step1 <- paste("Step 1: Calculate P(X >= cutoff_point).",
                 "This is the probability that X is greater than or equal to the cutoff point,",
                 "which can be calculated as P(X >= cutoff_point) = 1 - F(cutoff_point).",
                 "Given F(cutoff_point) =", round(cdf(cutoff_point), 4),
                 ", we find that P(X >= cutoff_point) =", round(p_x_positive, 4))

  # Step 2: Calculate the probability P(X < cutoff_point)
  p_x_negative <- cdf(cutoff_point)
  step2 <- paste("Step 2: Calculate P(X < cutoff_point).",
                 "This is the probability that X is less than the cutoff point,",
                 "calculated as P(X < cutoff_point) = F(cutoff_point).",
                 "We find that P(X < cutoff_point) =", round(p_x_negative, 4))

  # Step 3: Define the integral to compute E[Y | X >= cutoff_point]
  step3 <- paste("Step 3: Define the integral to compute E[Y | X >= cutoff_point].",
                 "This integral computes the expected value of Y given that X is in the range [",
                 lower_bound, ", ", upper_bound, "].",
                 "Since Y = a * X + b when X >= cutoff_point, we integrate (a * X + b) * f(X) over this range.")

  # Function to calculate E[Y | X >= cutoff_point]
  integrate_y_given_x_positive <- function() {
    integrate(function(x) (a * x + b) * (cdf(x + 0.0001) - cdf(x - 0.0001)) / 0.0002,
              lower = lower_bound, upper = upper_bound)$value
  }

  # Step 4: Compute E[Y | X >= cutoff_point]
  e_y_positive <- integrate_y_given_x_positive()
  step4 <- paste("Step 4: Calculate E[Y | X >= cutoff_point].",
                 "Using the integral defined in Step 3, we find that E[Y | X >= cutoff_point] =",
                 round(e_y_positive, 4))

  # Step 5: Calculate E[Y] using the law of total expectation
  e_y <- e_y_positive * p_x_positive + 0 * p_x_negative
  step5 <- paste("Step 5: Calculate E[Y] using the law of total expectation.",
                 "Since Y = 0 when X < cutoff_point, we have E[Y] = E[Y | X >= cutoff_point] * P(X >= cutoff_point) + E[Y | X < cutoff_point] * P(X < cutoff_point).",
                 "Substituting values from Steps 1, 2, and 4, we find E[Y] =", round(e_y, 4))

  # Return the detailed explanation as a single character vector
  return(c(step1, step2, step3, step4, step5, paste("Final Result: E[Y] =", round(e_y, 4))))
}


#' @title Paired t-Test for Pre-Test and Post-Test Scores with Detailed Calculations
#' @description This function performs a paired t-test to evaluate whether there is a significant improvement in scores from pre-test to post-test.
#' It also provides a confidence interval for the mean difference and detailed steps for each calculation.
#' @param pre_scores A numeric vector representing the pre-test scores.
#' @param post_scores A numeric vector representing the post-test scores.
#' @param alpha The significance level for the test (default is 0.05).
#' @param conf_level The confidence level for the confidence interval (default is 0.90).
#' @return A list containing the hypothesis statements, test statistic, p-value, confidence interval, detailed calculations, and interpretations.
#' @details
#' This function calculates the mean difference between post-test and pre-test scores and uses a paired t-test to test for a significant improvement.
#' It provides interpretations of the hypothesis test and the confidence interval for the mean difference, with step-by-step calculations.
#' @examples
#' pre_scores <- c(30, 30, 29, 23, 29, 30, 24, 32, 28, 16)
#' post_scores <- c(28, 30, 31, 18, 31, 31, 32, 33, 28, 20)
#' paired_t_test_pre_post(pre_scores, post_scores, alpha = 0.05, conf_level = 0.90)
#' @export
paired_t_test_pre_post <- function(pre_scores, post_scores, alpha = 0.05, conf_level = 0.90) {
  # Step (a): State the hypotheses
  hypotheses <- list(
    H0 = "The mean difference in scores (post-test - pre-test) is 0.",
    Ha = "The mean difference in scores (post-test - pre-test) is greater than 0, indicating an improvement."
  )

  # Check that the lengths of the two samples are the same
  if (length(pre_scores) != length(post_scores)) {
    stop("The length of pre-test scores must equal the length of post-test scores.")
  }

  # Calculate the differences
  differences <- post_scores - pre_scores

  # Step-by-step calculations
  n <- length(differences)
  mean_diff <- mean(differences)
  sd_diff <- sd(differences)
  standard_error <- sd_diff / sqrt(n)
  t_stat <- mean_diff / standard_error
  p_value <- 1 - pt(t_stat, df = n - 1)  # one-tailed test

  # Step-by-step explanation of calculations
  calculations <- list(
    "Calculate differences (post-test - pre-test) for each individual" = differences,
    paste("Step 1: Calculate mean of differences =", round(mean_diff, 4)),
    paste("Step 2: Calculate standard deviation of differences =", round(sd_diff, 4)),
    paste("Step 3: Calculate standard error of the mean difference (SE = sd / sqrt(n)) =", round(standard_error, 4)),
    paste("Step 4: Calculate t-statistic (t = mean_diff / SE) =", round(t_stat, 4)),
    paste("Step 5: Calculate p-value for one-tailed test (P(T > t) =", round(p_value, 4))
  )

  # Confidence Interval Calculations
  t_critical <- qt(1 - (1 - conf_level) / 2, df = n - 1)
  margin_of_error <- t_critical * standard_error
  ci_lower <- mean_diff - margin_of_error
  ci_upper <- mean_diff + margin_of_error
  confidence_interval <- paste("(", round(ci_lower, 4), ", ", round(ci_upper, 4), ")", sep = "")

  # Interpretations
  hypothesis_interpretation <- if (p_value < alpha) {
    "Since the p-value is less than the significance level, we reject the null hypothesis. There is significant evidence that attending the summer institute improves listening skills."
  } else {
    "Since the p-value is greater than the significance level, we do not reject the null hypothesis. There is not enough evidence to conclude that attending the summer institute improves listening skills."
  }

  confidence_interval_interpretation <- paste(
    "With", conf_level * 100, "% confidence, the mean increase in listening scores due to attending the summer institute is between",
    round(ci_lower, 4), "and", round(ci_upper, 4)
  )

  # Return all details in the output
  return(list(
    hypotheses = hypotheses,
    calculations = calculations,
    test_statistic = paste("t =", round(t_stat, 4)),
    p_value = paste("p-value =", round(p_value, 4)),
    conclusion = hypothesis_interpretation,
    confidence_interval = confidence_interval,
    confidence_interval_interpretation = confidence_interval_interpretation
  ))
}
