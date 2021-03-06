% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{getMat}
\alias{getMat}
\title{Extract a Specific Matrix from a 'VCA' Object}
\usage{
getMat(obj, mat)
}
\arguments{
\item{obj}{... (VCA) object}

\item{mat}{... (character) string specifying the matrix to be extracted}
}
\value{
(matrix) as requested by the user
}
\description{
For convenience only, extracting a specific matrix from the 
"Matrices" element of a 'VCA' object if this matrix exists.
}
\details{
When 'mat="Z"' the design matrix of random effects will be returned.
If one is interested in the design matrix of random effects for a specific
variance component use a name like "Z" + NAME, where NAME has to be equal to
the name of the VC in the 'VCA' object (see examples). The same applies to 
the A-matrices in the quadratic forms, use "A" + NAME for extracting a specific 
A-matrix.
}
\examples{
\dontrun{
data(dataEP05A2_1)
fit <- anovaVCA(y~day/run, dataEP05A2_1)
getMat(fit, "Z")
getMat(fit, "Zday")
getMat(fit, "Zday:run")
getMat(fit, "Zerror")
fit2 <- anovaMM(y~day/(run), dataEP05A2_1)
getMat(fit2, "V")			 	# Var(y)
getMat(fit2, "G")				# Var(re)
}
}
