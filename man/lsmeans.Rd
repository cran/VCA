\name{lsmeans}
\alias{lsmeans}
\title{Least Squares Means of Fixed Effects.}
\usage{
lsmeans(obj, var = NULL, type = c("simple", "complex"),
  ddfm = c("contain", "residual", "satterthwaite"), quiet = FALSE)
}
\arguments{
  \item{obj}{(VCA) object having at least one fixed effect}

  \item{var}{(character) string specifying a fixed effects
  variable for which LS Means should be computed, defaults
  to all fixed effects, i.e. for each level of a fixed
  effects variable ls means will be computed}

  \item{type}{(character) "simple" = fast version of
  computing LS means}

  \item{ddfm}{(character) string specifying the method used
  for computing the degrees of freedom of the t-statistic.
  Only used when type="complex".  Available methods are
  "contain", "residual", and "satterthwaite".}

  \item{quiet}{(logical) TRUE = suppress warning messages,
  e.g. for non-estimable contrasts}
}
\value{
(matrix) with LS Means of fixed effects and respective
standard errors, in case of 'type="complex"'
}
\description{
Computes Least Squares Means (LS Means) of fixed effects
for fitted mixed models of class 'VCA'.
}
\details{
Function computes LS Means of fixed effects and their
corresponding standard errors. In case of setting argument
'type' equal to "complex" (or any abbreviation) a
\eqn{t}-test is performed on each LS Mean, returning
degrees of freedom, t-statistic and corresponding p-values.
One can choose from one of three denominator degrees of
freedom ('ddfm')-methods.

Actually, function \code{\link{test.fixef}} is called with
the "no intercept" version of the fitted model. The
"complex" option is significantly slower for unbalanced
designs (see \code{\link{test.fixef}} for details). In case
that the 'VarCov' element of the 'VCA' object already
exists (calling \code{\link{vcovVC}}), which is the most
time consuming part, results can be obtained in less amount
of time.

Standard Errors of LS Means are computed as \eqn{TPT^{T}}{T
* P * T'}, where \eqn{T} is the LS Means generating
contrast matrix and \eqn{P} is the variance-covariance
matrix of fixed effects.
}
\examples{
\dontrun{
data(dataEP05A2_2)
fit1 <- anovaMM(y~day/(run), dataEP05A2_2)
lsmeans(fit1)
lsmeans(fit1,, "complex")

# a more complex model
data(VCAdata1)
fit2 <- anovaMM(y~(lot+device)/(day)/(run), VCAdata1[VCAdata1$sample==2,])
lsmeans(fit2, "lot")
lsmeans(fit2, "device", "complex")

# pre-computed 'VarCov' element saves time
system.time(lsm1 <- lsmeans(fit2, "device", "complex"))
fit2$VarCov <- vcovVC(fit2)
system.time(lsm2 <- lsmeans(fit2, "device", "complex"))
lsm1
lsm2
}}
\author{
Andre Schuetzenmeister
\email{andre.schuetzenmeister@roche.com}
}

