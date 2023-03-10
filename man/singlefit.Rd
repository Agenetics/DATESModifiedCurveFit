% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runfits.R
\name{singlefit}
\alias{singlefit}
\title{Run single exponential decay fit}
\usage{
singlefit(prefix, aa = 0.03, tt = 10, distance = distance)
}
\arguments{
\item{prefix}{Location of *out file}
}
\value{
A dataframe with 1 Admixture time estimates, standard error, Shapiro-Wilk p-value of fit
}
\description{
Estimate admixture timing for 1 admixtures
}
\examples{
out1 <- singlefit("/home/user/Downloads/DATESexample/example/50_40_0.1_0.1/Simulation.out:");
}
