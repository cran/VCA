\name{getV}
\alias{getV}
\title{Determine V-Matrix for a 'VCA' Object.}
\usage{
getV(obj)
}
\arguments{
  \item{obj}{(VCA) object}
}
\value{
(VCA) object with additional elements in the 'Matrices'
element, including matrix \eqn{V}.
}
\description{
Determine the estimated variance-covariance matrix of
observations \eqn{y}.
}
\details{
A linear mixed model can be written as \eqn{y = Xb + Zg +
e}, where \eqn{y} is the column vector of observations,
\eqn{X} and \eqn{Z} are design matrices assigning fixed
(\eqn{b}), respectively, random (\eqn{g}) effects to
observations, and \eqn{e} is the column vector of residual
errors. The variance-covariance matrix of \eqn{y} is equal
to \eqn{Var(y) = ZGZ^{-T} + R}{Var(y) = ZGZ' + R}, where
\eqn{R} is the variance-covariance matrix of \eqn{e}.
}
\author{
Andre Schuetzenmeister
\email{andre.schuetzenmeister@roche.com}
}

