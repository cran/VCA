\name{solveMM}
\alias{solveMM}
\title{Solve Mixed Model Equations.}
\usage{
solveMM(obj)
}
\arguments{
  \item{obj}{... (VCA) object}
}
\value{
(VCA) object, which has additional attributes 'ranef'
corresponding to the column vector of estimated random
effects, "fixef" being the column vector of estimated fixed
effects.  Attribute "Matrices" has additional elements
referring to the elements of the MMEs.
}
\description{
Function solves the Mixed Model Equations (MMEs) to
estimate fixed and random effects. It is for internal use
only, thus, not exported.
}
\examples{
\dontrun{
data(dataEP05A2_1)
fit <- anovaVCA(y~day/run, dataEP05A2_1, NegVC=TRUE)
fit <- VCA:::solveMM(fit)
ranef(fit)
}
}
\author{
Andre Schuetzenmeister
\email{andre.schuetzenmeister@roche.com}
}

