\name{MPinv}
\alias{MPinv}
\title{Moore-Penrose Inverse of a Matrix.}
\usage{
MPinv(X, tol = sqrt(.Machine$double.eps))
}
\arguments{
  \item{X}{(object) two-dimensional, for which a
  Moore-Penrose inverse has to be computed}

  \item{tol}{(numeric) tolerance value to be used in
  comparisons}
}
\value{
(object) A Moore-Penrose inverse of X.
}
\description{
This function is originally defined in package 'MASS'. It
was adapted to be able to deal with matrices from the
'Matrix' package, e.g. sparse matrices.
}
\author{
Authors of the 'MASS' package.
}

