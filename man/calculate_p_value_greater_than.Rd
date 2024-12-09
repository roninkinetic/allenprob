% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{calculate_p_value_greater_than}
\alias{calculate_p_value_greater_than}
\title{Perform a One-Sample T-Test and Calculate the P-Value}
\usage{
calculate_p_value_greater_than(data, population_mean, significance_level)

calculate_p_value_greater_than(data, population_mean, significance_level)
}
\arguments{
\item{data}{A numeric vector of sample data (numeric).}

\item{population_mean}{The population mean to compare the sample mean against (numeric).}

\item{significance_level}{The significance level for the test (numeric, e.g., 0.05 for a 5\% significance level).}
}
\value{
A list containing:
\itemize{
\item \code{p_value}: The calculated p-value (numeric).
\item \code{interpretation}: A string interpreting the result based on the significance level.
}

A list containing:
\itemize{
\item \code{p_value}: The calculated p-value (numeric).
\item \code{interpretation}: A string interpreting the result based on the significance level.
}
}
\description{
This function performs a one-sample t-test to determine if the sample mean is significantly
greater than a specified population mean. This test is useful for assessing whether
the observed sample mean provides sufficient evidence to conclude that it exceeds the
hypothesized population mean.

This function performs a one-sample t-test to determine if the sample mean is significantly
greater than a specified population mean. The t-test assesses whether the difference
between the sample mean and the population mean is statistically significant.
}
\details{
The calculation involves the following steps:
\enumerate{
\item Calculate the sample mean (x_bar) and sample standard deviation (s) from the input data.
\item Compute the t-statistic using the formula:
t = (x_bar - mu) / (s / sqrt(n))
where mu is the population mean and n is the sample size.
\item Use the t-distribution to calculate the p-value for the test under the alternative hypothesis that
the sample mean is greater than the population mean.
\item Interpret the result based on the specified significance level:
\itemize{
\item If the p-value is less than the significance level, conclude that the sample mean
is significantly greater than the population mean.
\item Otherwise, conclude that it is not significantly greater.
}
}

The calculation involves the following steps:
\enumerate{
\item Calculate the sample mean (x_bar) and sample standard deviation (s) from the input data.
\item Compute the t-statistic using the formula:
t = (x_bar - mu) / (s / sqrt(n))
where mu is the population mean and n is the sample size.
\item Calculate the p-value for the one-tailed t-test using the t-distribution.
\item Interpret the result based on the specified significance level:
\itemize{
\item If the p-value is less than the significance level, conclude that the sample mean
is significantly greater than the population mean.
\item Otherwise, conclude that it is not significantly greater.
}
}
}
\examples{
calculate_p_value_greater_than(c(5.1, 5.5, 5.9, 6.2, 5.7), 5.0, 0.05)

calculate_p_value_greater_than(c(5.1, 5.5, 5.9, 6.2, 5.7), 5.0, 0.05)

}
