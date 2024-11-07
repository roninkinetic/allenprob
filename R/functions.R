#' Calculate Probability for Sample Proportions
#'
#' Calculates the probability for sample proportions based on a specified population proportion, sample size, and target sample proportion.
#'
#' @param population_proportion The proportion in the population (e.g., 0.5).
#' @param sample_size The size of the sample (an integer).
#' @param sample_proportion_target The target sample proportion.
#' @return The probability as a numeric value.
#' @examples
#' sample_proportion_prob(0.5, 100, 0.55)
#' @export
sample_proportion_prob <- function(population_proportion, sample_size, sample_proportion_target) {
  std_error <- sqrt(population_proportion * (1 - population_proportion) / sample_size)
  z_score <- (sample_proportion_target - population_proportion) / std_error
  pnorm(z_score)
}

#' Calculate Confidence Interval for a Proportion
#'
#' Calculates the confidence interval for a proportion based on the number of successes, sample size, and confidence level.
#'
#' @param successes The number of successes in the sample.
#' @param sample_size The size of the sample (an integer).
#' @param confidence_level The confidence level (e.g., 0.95 for 95% confidence).
#' @return A list with the lower and upper bounds of the confidence interval.
#' @examples
#' calculate_proportion_confidence_interval(50, 100, 0.95)
#' @export
calculate_proportion_confidence_interval <- function(successes, sample_size, confidence_level) {
  p_hat <- successes / sample_size
  standard_error <- sqrt((p_hat * (1 - p_hat)) / sample_size)
  alpha <- 1 - confidence_level
  z_value <- qnorm(1 - alpha / 2)
  margin_of_error <- z_value * standard_error
  list(lower_bound = p_hat - margin_of_error, upper_bound = p_hat + margin_of_error)
}

#' Calculate Sample Mean Probability
#'
#' Calculates the probability of obtaining a sample mean at least as extreme as a target value given a population mean, standard deviation, and sample size.
#'
#' @param population_mean The mean of the population.
#' @param population_std_dev The standard deviation of the population.
#' @param sample_size The size of the sample.
#' @param sample_mean_target The target sample mean.
#' @return The probability as a numeric value.
#' @examples
#' calculate_sample_mean_probability(100, 15, 30, 105)
#' @export
calculate_sample_mean_probability <- function(population_mean, population_std_dev, sample_size, sample_mean_target) {
  std_error <- population_std_dev / sqrt(sample_size)
  z_score <- (sample_mean_target - population_mean) / std_error
  1 - pnorm(z_score)
}

#' Calculate Probability Between Two Times for an Exponential Distribution
#'
#' Calculates the probability of an event occurring between two times for an exponential distribution with a specified rate.
#'
#' @param lambda The rate parameter of the exponential distribution.
#' @param time_start The start time.
#' @param time_end The end time.
#' @return The probability of the event occurring between time_start and time_end.
#' @examples
#' calculate_exponential_probability(0.1, 2, 5)
#' @export
calculate_exponential_probability <- function(lambda, time_start, time_end) {
  pexp(time_end, rate = lambda) - pexp(time_start, rate = lambda)
}

#' Calculate Probability Between Two Bounds
#'
#' Calculates the probability of a normally distributed variable being between two specified bounds.
#'
#' @param lower_bound The lower bound.
#' @param upper_bound The upper bound.
#' @param mean The mean of the distribution.
#' @param sd The standard deviation of the distribution.
#' @return The probability of the variable being between lower_bound and upper_bound.
#' @examples
#' prob_between(10, 20, 15, 5)
#' @export
prob_between <- function(lower_bound, upper_bound, mean, sd) {
  pnorm(upper_bound, mean = mean, sd = sd) - pnorm(lower_bound, mean = mean, sd = sd)
}

#' Calculate Confidence Interval for Sample Mean
#'
#' Calculates the confidence interval for a sample mean based on sample mean, sample standard deviation, sample size, and confidence level.
#'
#' @param sample_mean The mean of the sample.
#' @param sample_sd The standard deviation of the sample.
#' @param sample_size The size of the sample.
#' @param confidence_level The confidence level (e.g., 0.95 for 95% confidence).
#' @return A list with the lower and upper bounds of the confidence interval.
#' @examples
#' calculate_confidence_interval(50, 10, 30, 0.95)
#' @export
calculate_confidence_interval <- function(sample_mean, sample_sd, sample_size, confidence_level) {
  alpha <- 1 - confidence_level
  z_value <- qnorm(1 - alpha / 2)
  standard_error <- sample_sd / sqrt(sample_size)
  margin_of_error <- z_value * standard_error
  list(lower_bound = sample_mean - margin_of_error, upper_bound = sample_mean + margin_of_error)
}

#' Calculate Probability of Sample Mean Greater Than or Equal to Target
#'
#' Calculates the probability of the sample mean being greater than or equal to a target value.
#'
#' @param population_mean The mean of the population.
#' @param population_sd The standard deviation of the population.
#' @param sample_size The size of the sample.
#' @param sample_mean The target sample mean.
#' @return The probability as a numeric value.
#' @examples
#' calculate_probability_mean_greater_than(100, 15, 30, 105)
#' @export
calculate_probability_mean_greater_than <- function(population_mean, population_sd, sample_size, sample_mean) {
  standard_error <- population_sd / sqrt(sample_size)
  z_score <- (sample_mean - population_mean) / standard_error
  1 - pnorm(z_score)
}

#' Calculate Probability of Sample Mean Less Than or Equal to Target
#'
#' Calculates the probability of the sample mean being less than or equal to a target value.
#'
#' @param population_mean The mean of the population.
#' @param population_sd The standard deviation of the population.
#' @param sample_size The size of the sample.
#' @param sample_mean The target sample mean.
#' @return The probability as a numeric value.
#' @examples
#' calculate_probability_mean_less_than_or_equal(100, 15, 30, 95)
#' @export
calculate_probability_mean_less_than_or_equal <- function(population_mean, population_sd, sample_size, sample_mean) {
  standard_error <- population_sd / sqrt(sample_size)
  z_score <- (sample_mean - population_mean) / standard_error
  pnorm(z_score)
}

#' Find X Value for a Given Probability
#'
#' Finds the value of `x` such that the probability of a normally distributed random variable being greater than `x` is equal to the specified probability.
#'
#' @param mean The mean of the distribution.
#' @param sd The standard deviation of the distribution.
#' @param probability The target probability.
#' @return The `x` value corresponding to the specified probability.
#' @examples
#' find_x_value(100, 15, 0.05)
#' @export
find_x_value <- function(mean, sd, probability) {
  z_value <- qnorm(1 - probability)
  mean + z_value * sd
}

#' Perform a Two-Tailed Hypothesis Test for a Proportion
#'
#' Performs a two-tailed hypothesis test for a proportion to determine if the sample proportion significantly differs from the null hypothesis proportion (default is 0.5).
#'
#' @param successes The number of successes in the sample.
#' @param trials The number of trials in the sample (sample size).
#' @param significance_level The significance level for the hypothesis test (e.g., 0.05 for a 5% significance level).
#' @return A list containing the z-test statistic, p-value, and the decision result ("Reject the null hypothesis" or "Fail to reject the null hypothesis").
#' @examples
#' hypothesis_test_proportion(55, 100, 0.05)
#' @export
hypothesis_test_proportion <- function(successes, trials, significance_level) {
  p_hat <- successes / trials
  p_0 <- 0.5
  standard_error <- sqrt((p_0 * (1 - p_0)) / trials)
  z <- (p_hat - p_0) / standard_error
  p_value <- 2 * (1 - pnorm(abs(z)))
  alpha <- significance_level
  z_critical <- qnorm(1 - alpha / 2)
  result <- if (abs(z) > z_critical) {
    "Reject the null hypothesis."
  } else {
    "Fail to reject the null hypothesis."
  }
  list(z = z, p_value = p_value, result = result)
}

#' Perform a One-Sample T-Test and Calculate the P-Value
#'
#' Performs a one-sample t-test to determine if the sample mean is significantly greater than a specified population mean.
#'
#' @param data A numeric vector of sample data.
#' @param population_mean The population mean to compare the sample mean against.
#' @param significance_level The significance level for the test (e.g., 0.05 for a 5% significance level).
#' @return A list containing the p-value and an interpretation of the result based on the significance level.
#' @examples
#' calculate_p_value_greater_than(c(5.1, 5.5, 5.9, 6.2, 5.7), 5.0, 0.05)
#' @export
calculate_p_value_greater_than <- function(data, population_mean, significance_level) {
  t_test_result <- t.test(data, mu = population_mean, alternative = "greater")
  p_value <- t_test_result$p.value
  interpretation <- ifelse(
    p_value < significance_level,
    "The sample mean is significantly greater than the population mean.",
    "The sample mean is not significantly greater than the population mean."
  )
  list(p_value = p_value, interpretation = interpretation)
}

#' Calculate Probability of Cost Being Less Than a Given Value
#'
#' Calculates the probability that a cost is less than a specified threshold, based on the mean and standard deviation of the cost.
#'
#' @param mean_cost The mean cost of the distribution.
#' @param sd_cost The standard deviation of the cost.
#' @param threshold_cost The cost threshold to compare against.
#' @return The probability of the cost being less than the threshold.
#' @examples
#' calculate_probability_cost_less_than(100, 15, 90)
#' @export
calculate_probability_cost_less_than <- function(mean_cost, sd_cost, threshold_cost) {
  z_score <- (threshold_cost - mean_cost) / sd_cost
  pnorm(z_score)
}

#' Calculate Probability of Cost Being Greater Than a Given Value
#'
#' Calculates the probability that a cost is greater than a specified threshold, based on the mean and standard deviation of the cost.
#'
#' @param mean_cost The mean cost of the distribution.
#' @param sd_cost The standard deviation of the cost.
#' @param threshold_cost The cost threshold to compare against.
#' @return The probability of the cost being greater than the threshold.
#' @examples
#' calculate_probability_cost_greater_than(100, 15, 110)
#' @export
calculate_probability_cost_greater_than <- function(mean_cost, sd_cost, threshold_cost) {
  z_score <- (threshold_cost - mean_cost) / sd_cost
  1 - pnorm(z_score)
}

#' Perform a One-Sample T-Test and Calculate the P-Value
#'
#' Performs a one-sample t-test to determine if the sample mean is significantly greater than a specified population mean.
#'
#' @param data A numeric vector of sample data.
#' @param population_mean The population mean to compare the sample mean against.
#' @param significance_level The significance level for the test (e.g., 0.05 for a 5% significance level).
#' @return A list containing the p-value and an interpretation of the result based on the significance level.
#' @examples
#' calculate_p_value_greater_than(c(5.1, 5.5, 5.9, 6.2, 5.7), 5.0, 0.05)
#' @export
calculate_p_value_greater_than <- function(data, population_mean, significance_level) {
  t_test_result <- t.test(data, mu = population_mean, alternative = "greater")
  p_value <- t_test_result$p.value
  interpretation <- ifelse(
    p_value < significance_level,
    "The sample mean is significantly greater than the population mean.",
    "The sample mean is not significantly greater than the population mean."
  )
  list(p_value = p_value, interpretation = interpretation)
}

#' cor3var: Calculate the correlation coefficient from standard deviations and slope
#'
#' This function calculates the correlation coefficient (r) using the standard
#' deviation of x (s_x), the standard deviation of y (s_y), and the slope of the
#' regression line.
#'
#' @param s_x Standard deviation of x.
#' @param s_y Standard deviation of y.
#' @param slope The slope of the regression line.
#'
#' @return The correlation coefficient (r).
#' @examples
#' cor3var(1.5, 2.5, 0.8)
#'
#' @export
cor3var <- function(s_x, s_y, slope) {
  # Calculate the correlation coefficient (r)
  r <- (slope * s_x) / s_y
  return(r)
}

#' list_functions: List all available functions in the package
#'
#' This function lists all the available functions in the 'statHelper' package.
#'
#' @return A character vector of function names.
#' @examples
#' list_functions()
#'
#' @export
list_functions <- function() {
  # List all functions available in the statHelper package
  return(lsf.str("package:statHelper"))
}

#' r_squared: Calculate the coefficient of determination (R^2) given the correlation coefficient (r)
#'
#' This function calculates the coefficient of determination (R^2) from the correlation coefficient (r).
#'
#' @param r Correlation coefficient.
#'
#' @return The coefficient of determination (R^2).
#' @examples
#' r_squared(0.8)
#'
#' @export
r_squared <- function(r) {
  # Calculate the coefficient of determination (R^2)
  r_squared <- r^2
  return(r_squared)
}

#' lm_summary: Perform linear regression and show the summary
#'
#' This function takes two data vectors x and y and fits a linear model (lm) to them. It then returns the summary of the model.
#'
#' @param x A numeric vector of independent variable values.
#' @param y A numeric vector of dependent variable values.
#'
#' @return The summary of the linear model.
#' @examples
#' lm_summary(c(1, 2, 3, 4), c(2, 4, 6, 8))
#'
#' @export
lm_summary <- function(x, y) {
  # Fit a linear model lm(y ~ x)
  model <- lm(y ~ x)

  # Return the summary of the model
  return(summary(model))
}

#' least_squares_regression: Perform least squares regression given ybar, xbar, s_x, s_y, and R
#'
#' This function calculates the least squares regression line using the mean of y (ybar),
#' the mean of x (xbar), the standard deviations of x (s_x), and y (s_y), and the correlation coefficient (R).
#'
#' The regression line is of the form: y = intercept + slope * x
#'
#' @param ybar Mean of y values.
#' @param xbar Mean of x values.
#' @param s_x Standard deviation of x values.
#' @param s_y Standard deviation of y values.
#' @param R Correlation coefficient between x and y.
#'
#' @return A list containing the slope and intercept of the regression line.
#' @examples
#' least_squares_regression(5, 3, 1.2, 2.1, 0.8)
#'
#' @export
least_squares_regression <- function(ybar, xbar, s_x, s_y, R) {
  # Calculate the slope
  slope <- R * (s_y / s_x)

  # Calculate the intercept
  intercept <- ybar - slope * xbar

  # Return a list with the slope and intercept
  return(list(slope = slope, intercept = intercept))
}

#' sample_mean_variance: Calculate the mean and variance of a sample distribution
#'
#' This function takes two vectors, X (values) and P(X) (probabilities), and calculates
#' the mean and variance of the sample distribution.
#'
#' @param X A numeric vector of values.
#' @param P_X A numeric vector of probabilities associated with X.
#'
#' @return A list containing the mean and variance of the sample distribution.
#' @examples
#' sample_mean_variance(c(1, 2, 3), c(0.2, 0.5, 0.3))
#'
#' @export
sample_mean_variance <- function(X, P_X) {
  # Calculate the mean
  mean_value <- sum(X * P_X)

  # Calculate the variance
  variance_value <- sum(((X - mean_value)^2) * P_X)

  return(list(mean = mean_value, variance = variance_value))
}

#' expected_value: Calculate the expected value given win/lose amounts and probabilities
#'
#' This function calculates the expected value of an event given the win amount, lose amount,
#' probability of win, and probability of loss.
#'
#' @param win_amount The amount won if the event is successful.
#' @param lose_amount The amount lost if the event is unsuccessful.
#' @param prob_win Probability of winning the event.
#' @param prob_loss Probability of losing the event.
#'
#' @return The expected value of the event.
#' @examples
#' expected_value(100, -50, 0.6, 0.4)
#'
#' @export
expected_value <- function(win_amount, lose_amount, prob_win, prob_loss) {
  # Calculate the expected value
  expected_value <- (win_amount * prob_win) + (lose_amount * prob_loss)

  return(expected_value)
}
