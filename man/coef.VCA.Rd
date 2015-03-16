\name{coef.VCA}
\alias{coef.VCA}
\title{Extract Fixed Effects from 'VCA' Object.}
\usage{
\method{coef}{VCA}(object, ...)
}
\arguments{
  \item{object}{(VCA) object where fixed effects shall be
  extracted}

  \item{...}{additional parameters}
}
\description{
For coneniently using objects of class 'VCA' with other
packages expecting this function, e.g. the 'multcomp'
package for general linear hypotheses for parametric
models.
}
\examples{
data(dataEP05A2_1)
fit1 <- anovaMM(y~day/(run), dataEP05A2_1)
coef(fit1)
fit2 <- anovaVCA(y~day/run, dataEP05A2_1)
coef(fit2)
}

