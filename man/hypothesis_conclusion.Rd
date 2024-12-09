% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{hypothesis_conclusion}
\alias{hypothesis_conclusion}
\title{Determine the Conclusion of a Hypothesis Test Based on the p-Value}
\usage{
hypothesis_conclusion(p_value)
}
\arguments{
\item{p_value}{A numeric value representing the computed p-value from the hypothesis test.
It should be a number between 0 and 1.}
}
\value{
A string indicating the conclusion:
\itemize{
\item "Reject the null hypothesis" if the p-value < 0.001.
\item "Fail to reject the null hypothesis" if the p-value >= 0.001.
}
}
\description{
This function evaluates the evidence against the null hypothesis based on the
computed p-value from a hypothesis test. It provides a conclusion about whether
to reject or fail to reject the null hypothesis, with a threshold for very strong
evidence set at a p-value of 0.001.
}
\examples{
conclusion <- hypothesis_conclusion(p_value = 0.0005)
print(conclusion)  # Expected output: "Reject the null hypothesis"

}
