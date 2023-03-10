% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runfits.R
\name{try100twofit}
\alias{try100twofit}
\title{Wrapper around trytwofit function}
\usage{
try100twofit(folder, i = 0, jj = 50, distance = distance)
}
\arguments{
\item{folder}{Location of *out file}
}
\value{
A dataframe with 2 Admixture time estimates, standard errors, Shapiro-Wilk p-value of fit; NULL output if no convergent solution
}
\description{
Run the trytwofit function 5 times and choose the outpur with the highest Z scores
}
\examples{
out1 <- try100twofit("/home/user/Downloads/DATESexample/example/50_40_0.1_0.1/Simulation.out:");
}
