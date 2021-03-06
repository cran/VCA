% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{as.matrix.VCA}
\alias{as.matrix.VCA}
\title{Standard 'as.matrix' Method for 'VCA' S3-Objects}
\usage{
\method{as.matrix}{VCA}(x, ...)
}
\arguments{
\item{x}{(VCA) object}

\item{...}{additional arguments to be passed to or from methods.}
}
\value{
(matrix) equal to x$aov.tab with additional attributes "Mean" and "Nobs"
}
\description{
Standard 'as.matrix' Method for 'VCA' S3-Objects
}
\examples{
\dontrun{
data(dataEP05A2_1)
fit <- anovaVCA(y~day/run, dataEP05A2_1)
as.matrix(fit)
}
}
\seealso{
\code{\link{as.matrix.VCAinference}}
}
\author{
Andre Schuetzenmeister \email{andre.schuetzenmeister@roche.com}
}
