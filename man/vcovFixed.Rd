\name{vcovFixed}
\alias{vcovFixed}
\title{Calculate Variance-Covariance Matrix and Standard Errors of Fixed Effects for an 'VCA' Object.}
\usage{
vcovFixed(obj)
}
\arguments{
  \item{obj}{... (VCA) object for which the
  variance-covariance matrix of fixed effects shall be
  calculated}
}
\value{
(matrix) corresponding to the variance-covariance matrix of
fixed effects
}
\description{
The variance-covariance matrix of fixed effects for the
linear mixed model in 'obj' is calculated.
}
\details{
The variance-covariance matrix of fixed effects for a
linear mixed model corresponds to matrix
\eqn{(X^{T}V^{-1}X)^{-}}{(X'V"X)`}, where >\eqn{^{T}}{'}<
denotes the transpose operator, >\eqn{^{-1}}{"}< the
regular matrix inverse, and >\eqn{^{-}}{`}< the generalized
(Moore-Penrose) inverse of a matrix.
}
\examples{
data(dataEP05A2_1)
fit <- anovaVCA(y~day/(run), dataEP05A2_1)
vcov(fit)
}

