% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{compare_distributions}
\alias{compare_distributions}
\title{Compare the Distributions for X and Xbar}
\usage{
compare_distributions(population_mean, population_sd, sample_size)
}
\arguments{
\item{population_mean}{The population mean (numeric).}

\item{population_sd}{The population standard deviation (numeric).}

\item{sample_size}{The sample size (numeric).}
}
\value{
A list containing:
\itemize{
\item \code{mean_X}: The mean of the distribution for X.
\item \code{sd_X}: The standard deviation of the distribution for X.
\item \code{mean_Xbar}: The mean of the distribution for Xbar.
\item \code{sd_Xbar}: The standard deviation of the distribution for Xbar.
\item \code{difference_in_sd}: The difference in the standard deviations between X and Xbar.
\item \code{explanation}: A string explanation of the key differences between the two distributions.
}
}
\description{
This function compares the distributions for X (a random variable) and Xbar (the sample mean)
by calculating the mean and standard deviation for both. It also returns the difference in the standard deviations.
Additionally, the function returns an explanation of the key differences between the two distributions.
}
\examples{
compare_distributions(67, 6, 30)

}
