% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{sampling_distribution_parameters}
\alias{sampling_distribution_parameters}
\title{Calculate Mean and Standard Deviation of the Sampling Distribution of Proportions}
\usage{
sampling_distribution_parameters(population_proportion, sample_size)
}
\arguments{
\item{population_proportion}{The proportion of the population with the characteristic of interest (numeric).}

\item{sample_size}{The size of the sample (integer).}
}
\value{
A list containing the mean and standard deviation of the sampling distribution.
}
\description{
This function calculates the mean and standard deviation of the sampling distribution
for a given population proportion and sample size.
}
\details{
The mean of the sampling distribution is equal to the population proportion, and the
standard deviation (standard error) is calculated based on the population proportion
and the sample size.
}
\examples{
sampling_distribution_parameters(0.74, 100)
}
