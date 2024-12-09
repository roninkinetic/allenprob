% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{paired_t_test_pre_post}
\alias{paired_t_test_pre_post}
\title{Paired t-Test for Pre-Test and Post-Test Scores with Detailed Calculations}
\usage{
paired_t_test_pre_post(pre_scores, post_scores, alpha = 0.05, conf_level = 0.9)
}
\arguments{
\item{pre_scores}{A numeric vector representing the pre-test scores.}

\item{post_scores}{A numeric vector representing the post-test scores.}

\item{alpha}{The significance level for the test (default is 0.05).}

\item{conf_level}{The confidence level for the confidence interval (default is 0.90).}
}
\value{
A list containing the hypothesis statements, test statistic, p-value, confidence interval, detailed calculations, and interpretations.
}
\description{
This function performs a paired t-test to evaluate whether there is a significant improvement in scores from pre-test to post-test.
It also provides a confidence interval for the mean difference and detailed steps for each calculation.
}
\details{
This function calculates the mean difference between post-test and pre-test scores and uses a paired t-test to test for a significant improvement.
It provides interpretations of the hypothesis test and the confidence interval for the mean difference, with step-by-step calculations.
}
\examples{
pre_scores <- c(30, 30, 29, 23, 29, 30, 24, 32, 28, 16)
post_scores <- c(28, 30, 31, 18, 31, 31, 32, 33, 28, 20)
paired_t_test_pre_post(pre_scores, post_scores, alpha = 0.05, conf_level = 0.90)
}
