% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vca.R
\name{print.VCAinference}
\alias{print.VCAinference}
\title{Standard Print Method for Objects of Class 'VCAinference'}
\usage{
\method{print}{VCAinference}(x, digits = 4L, what = c("all", "VC", "SD", "CV", "VCA"), ...)
}
\arguments{
\item{x}{(VCAinference) object}

\item{digits}{(integer) number of decimal digits.}

\item{what}{(character) one of "all", "VC", "SD", "CV", "VCA" specifying which part of the 'VCA'-object is to be printed.}

\item{...}{additional arguments to be passed to or from methods.}
}
\value{
invisibly returns sub-elements of 'x' specified via 'what'
}
\description{
Prints the list-type 'VCAinference'-object as tabulated output.
}
\details{
Formats the list-type objects of class 'VCAinference' for a more comprehensive
presentation of results, which are easier to grasp. The default is to show the complete
object (VCA ANOVA-table, VC-, SD-, and CV-CIs). Using parameter 'what' allows to
restrict the printed output to certain parts. Print-function invisibly returns a matrix
or a list of matrices, depending on the values of 'what', i.e. it can be used as for
packing the inference-information in one or multiple matrix-objects and extracting it/them.
}
\examples{
\dontrun{
# load data (CLSI EP05-A2 Within-Lab Precision Experiment) 
data(dataEP05A2_1)

# perform ANOVA-estimation of variance components for a nested design
res <- anovaVCA(y~day/run, Data=dataEP05A2_1)
res
inf <- VCAinference(res)
inf

# show certain parts and extract them invisibly
CVmat <- print(inf, what="CV")
CVmat

# show numerical values with more digits
print(inf, digit=12)
}
}
\seealso{
\code{\link{VCAinference}}, \code{\link{anovaVCA}}
}
\author{
Andre Schuetzenmeister \email{andre.schuetzenmeister@roche.com}
}
