% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AuxiliaryFunctions.R
\name{GetDerivativeSRLC}
\alias{GetDerivativeSRLC}
\title{Derivatives of smooth (re-weighted) density}
\usage{
GetDerivativeSRLC(y, logcon, h)
}
\arguments{
\item{y}{Where the function should be evaluated.}

\item{logcon}{Fitted log-concave density of class \code{dlc} from package
\code{logcondens}.}

\item{h}{Bandwith used for localization.}
}
\value{
Vector of length3 containing \eqn{\hat{f}^*_h(z), \hat{f}^*_h(z)',
  \hat{f}^*(z)''}.
}
\description{
Calculation of the smooth (re-weighted) density and the first and second order
derivative. To be precise, the density is estimated without normalization and
the derivatives thereof ar calculated. However, calculation of the modes are
not effected.
}
