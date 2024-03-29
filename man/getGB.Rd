% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{getGB}
\alias{getGB}
\title{Giesbrecht & Burns Approximation of the Variance-Covariance Matrix of Variance Components}
\usage{
getGB(obj, tol = 1e-12)
}
\arguments{
\item{obj}{(object) with list-type structure, e.g. \code{VCA} object fitted by ANOVA
or a premature \code{VCA} object fitted by REML}

\item{tol}{(numeric) values < 'tol' will be considered being equal to zero}
}
\value{
(matrix) corresponding to the Giesbrecht & Burns approximation
of the variance-covariance matrix of variance components
}
\description{
Compute variance covariance matrix of variance components of a linear mixed model
via the method stated in Giesbrecht and Burns (1985).
}
\details{
This function is not intended to be called by users and therefore not exported.
}
\examples{
\dontrun{
data(dataEP05A2_3)
fit <- anovaVCA(y~day/run, dataEP05A2_3)
fit <- solveMME(fit)		# some additional matrices required
getGB(fit)
}
}
\references{
Searle, S.R, Casella, G., McCulloch, C.E. (1992), Variance Components, Wiley New York

Giesbrecht, F.G. and Burns, J.C. (1985), Two-Stage Analysis Based on a Mixed Model: Large-Sample
Asymptotic Theory and Small-Sample Simulation Results, Biometrics 41, p. 477-486
}
\seealso{
\code{\link{vcovVC}}, \code{\link{remlVCA}}, \code{\link{remlMM}}
}
\author{
Andre Schuetzenmeister \email{andre.schuetzenmeister@roche.com},
Florian Dufey \email{florian.dufey@roche.com}
}
