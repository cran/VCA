% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{ranef.VCA}
\alias{ranef.VCA}
\title{Extract Random Effects from 'VCA' Object}
\usage{
\method{ranef}{VCA}(
  object,
  term = NULL,
  mode = c("raw", "student", "standard"),
  quiet = FALSE,
  ...
)
}
\arguments{
\item{object}{(VCA) object from which random effects shall be extracted}

\item{term}{(character) string specifying a term (factor) for which random effects 
should be extracted, one can also specify an integer which is interpreted
as i-th element of 'obj$res.assign$terms'}

\item{mode}{(character) string or abbreviation specifying whether "raw" residuals
should be returned or a transformed version c("student" or "standard")}

\item{quiet}{(logical) TRUE = will suppress any warning, which will be issued otherwise}

\item{...}{additional parameters}
}
\description{
Extract random effects and possibly apply a transformation to them (standardization,
studentization).
}
\details{
Extracting the 'RandomEffects' element of an 'VCA' object if this exists and applying
standardization (mean 0, sd 1) or studentization. For studentized random effects 
the i-th random effects is divided by the i-th main diagonal element of matrix \eqn{O = GZ^{T}QZG}{O = GZ'QZG},
where \eqn{G} is the covariance-matrix of random effects, \eqn{Z} is a design matrix assigning 
random effects to observations and matrix \eqn{Q = V^{-1}(I - H)}{Q = V"(I - H)} (see \code{\link{residuals.VCA}} for further details).
}
\examples{
\dontrun{
data(dataEP05A2_1)
fit <- anovaVCA(y~day/run, dataEP05A2_1)
ranef(fit)

# get variable-specific random effects (REs)
# both extract the same REs
ranef(fit, "day")
ranef(fit, 1)

# get standardized REs
ranef(fit, "day:run", "standard")

# or studentized REs
ranef(fit, 2, "stu")
}
}
\references{
Searle, S.R, Casella, G., McCulloch, C.E. (1992), Variance Components, Wiley New York	

Laird, N.M., Ware, J.H., 1982. Random effects models for longitudinal data. Biometrics 38, 963-974.

Schuetzenmeister, A. and Piepho, H.P. (2012). Residual analysis of linear mixed models using a simulation approach.
Computational Statistics and Data Analysis, 56, 1405-1416
}
