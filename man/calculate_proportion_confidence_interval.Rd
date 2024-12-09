% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{calculate_proportion_confidence_interval}
\alias{calculate_proportion_confidence_interval}
\title{Calculate Confidence Interval for a Proportion}
\usage{
calculate_proportion_confidence_interval(
  successes,
  sample_size,
  confidence_level
)
}
\arguments{
\item{successes}{The number of successes in the sample (an integer).}

\item{sample_size}{The size of the sample (an integer).}

\item{confidence_level}{The confidence level (e.g., 0.95 for 95\% confidence).}
}
\value{
A list containing the lower and upper bounds of the confidence interval:
\itemize{
\item \code{lower_bound}: The lower bound of the confidence interval.
\item \code{upper_bound}: The upper bound of the confidence interval.
}
}
\description{
This function calculates the confidence interval for a population proportion
based on the number of successes observed in a sample, the sample size,
and the desired confidence level.
}
\details{
The calculation involves the following steps:
\enumerate{
\item Compute the sample proportion (p_hat) using the formula:
p_hat = successes / sample_size
\item Calculate the standard error (SE) of the sample proportion using the formula:
SE = sqrt((p_hat * (1 - p_hat)) / n)
where n is the sample size.
\item Determine the critical z-value for the specified confidence level:
z = qnorm(1 - (alpha / 2))
where alpha is the significance level (e.g., alpha = 1 - confidence_level).
\item Calculate the margin of error using the formula:
ME = z * SE
\item The confidence interval is then given by:
CI = (p_hat - ME, p_hat + ME)
}
}
\examples{
calculate_proportion_confidence_interval(50, 100, 0.95)

}
