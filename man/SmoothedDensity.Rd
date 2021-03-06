% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AuxiliaryFunctions.R
\name{SmoothedDensity}
\alias{SmoothedDensity}
\title{Derivatves of the smoothed log-concave density}
\usage{
SmoothedDensity(z, x, phi, gam)
}
\arguments{
\item{z}{Where the function is evaluated.}

\item{x}{Vector [m] of knots.}

\item{phi}{Vector [m] of values of the log-concave density corresponding to
the knots \code{x}.}

\item{gam}{Bandwidth paramter of the Gaussian smoothing.}
}
\value{
A list with the value, 1st  and second order derivative of the smoothed
  log-concave density at \code{z}.
}
\description{
Derivatves of the smoothed log-concave density
}
