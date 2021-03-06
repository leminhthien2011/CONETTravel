% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/populationdynamic_func.R
\name{populationdynamicfunc}
\alias{populationdynamicfunc}
\title{This function gives the estimation for the dynamic population.}
\usage{
populationdynamicfunc(country)
}
\arguments{
\item{country}{list include total travel in: total_in, total travel out each day,
total_out, durationtravel, population P}
}
\value{
the dynamic population of the coutry with travel take in to account
}
\description{
This function gives the estimation for the dynamic population.
}
\examples{
\dontrun{
library(CONETTravel)
P= 10^6
durationtravel = 84
total_in = round(runif(84,300,400),digits=0)
total_out = round(runif(84,300,400),digits=0)
country = list(P =P, total_in = total_in, total_out= total_out,durationtravel=durationtravel )
populationdynamicfunc(country)
}
}
