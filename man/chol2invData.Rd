\name{chol2invData}
\docType{data}
\alias{chol2invData}
\title{Dataset Generating Error in Function 'chol2inv'}
\description{
	This dataset was added because it generates an error in function 'chol2inv' when
	trying to invert the variance-covariance matrix 'V' of the mixed model 'value~ID+(Site)'.
	This dataset and the associated mixed model are part of the unit-test collection of the
	package.
}
\usage{chol2invData}
\format{A data frame with 158 observations on the following 3 variables.
\itemize{
\item value\cr\cr
	The response variable.\cr\cr
\item ID\cr\cr
	Variable with 6 levels corresponding to samples.\cr\cr
\item Site\cr\cr
	Variable with 3 levels corresponding to sites/devices.\cr
}}
\keyword{datasets}