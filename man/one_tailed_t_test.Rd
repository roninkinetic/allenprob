% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{one_tailed_t_test}
\alias{one_tailed_t_test}
\title{One-Tailed Two-Sample t-Test for Difference in Means}
\usage{
one_tailed_t_test(
  n1,
  mean1,
  sd1,
  n2,
  mean2,
  sd2,
  alpha = 0.05,
  direction = "greater"
)
}
\arguments{
\item{n1}{Sample size for the first sample.}

\item{mean1}{Sample mean for the first sample.}

\item{sd1}{Sample standard deviation for the first sample.}

\item{n2}{Sample size for the second sample.}

\item{mean2}{Sample mean for the second sample.}

\item{sd2}{Sample standard deviation for the second sample.}

\item{alpha}{The significance level for the test (default is 0.05).}

\item{direction}{A string indicating the direction of the test: "greater" if testing if sample 1 is greater than sample 2, or "less" if testing if sample 1 is less than sample 2.}
}
\value{
A list containing the test type, test statistic, p-value, confidence interval, decision, and conclusion, with all work shown.
}
\description{
Performs a one-tailed t-test to determine if there is a significant difference in the means of two independent samples.
}
\details{
This function calculates the test statistic, p-value, and confidence interval for a one-tailed two-sample t-test. It assumes that the samples are independent, normally distributed, and have equal variances.
}
\examples{
one_tailed_t_test(n1 = 15, mean1 = 540, sd1 = 21, n2 = 15, mean2 = 554, sd2 = 15, alpha = 0.05, direction = "less")
}
