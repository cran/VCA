\name{getL}
\alias{getL}
\title{Construct Linear Contrast Matrix for Hypothesis Tests.}
\usage{
getL(obj, s, what = c("fixef", "lsmeans"))
}
\arguments{
  \item{obj}{(VCA) object}

  \item{s}{(character) string or vector of strings,
  denoting one or multiple linear contrasts}

  \item{what}{(character) string specifying whether to
  construct contrast matrices of fixed effects ("fixed") or
  LS Means ("lsmeans"), abbreviations are allowed.}
}
\value{
(matrix) representing one linear hypothesis of fixed
effects or LS Means per row
}
\description{
Function constructs coefficient/contrast matrices from a
string-representation of linear hypotheses.
}
\details{
Function constructs matrices expressing custom linear
hypotheses of fixed effects or LS Means. The user has to
specify a string denoting this contrast which is then
transformed into a coefficient/contrast matrix. This string
may contain names of fixed effects belonging to same same
fixed term, numeric coefficients and mathematical operators
"+" and "-" (see examples).
}
\examples{
data(dataEP05A2_2)
fit <- anovaMM(y~day/(run), dataEP05A2_2)
L <- getL(fit, c("day1-day2", "day5-day10"), what="fixef")
L
test.fixef(fit, L=L)

# another custom hypothesis
L2 <- getL(fit, "0.25*day1+0.25*day2+0.5*day3-0.5*day4-0.5*day5")

# more complex model
data(VCAdata1)
dataS2 <- VCAdata1[VCAdata1$sample==2,]
fit.S2 <- anovaMM(y~(lot+device)/day/(run), dataS2)
L3 <- getL(fit.S2, c("lot1-lot2", "lot1:device3:day19-lot1:device3:day20",
						"lot1:device1:day1-lot1:device1:day2"))
test.fixef(fit.S2, L3)
}
\author{
Andre Schuetzenmeister
\email{andre.schuetzenmeister@roche.com}
}

