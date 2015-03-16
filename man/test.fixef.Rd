\name{test.fixef}
\alias{test.fixef}
\title{Perform t-Tests for Linear Contrasts on Fixed Effects.}
\usage{
test.fixef(obj, L, ddfm = c("contain", "residual", "satterthwaite"),
  method.grad = "simple", tol = 1e-12, quiet = FALSE, opt = TRUE,
  onlyDF = FALSE, ...)
}
\arguments{
  \item{obj}{(VCA) object}

  \item{L}{(numeric) vector or matrix, specifying linear
  combinations of the fixed effects, in the latter case,
  each line represents a disctinct linear contrast}

  \item{ddfm}{(character) string specifying the method used
  for computing the denominator degrees of freedom for
  tests of fixed effects or LS Means. Available methods are
  "contain", "residual", and "satterthwaite".}

  \item{method.grad}{(character) string specifying the
  method to be used for approximating the gradient of the
  variance-covariance matrix of fixed effects at the
  estimated covariance parameter estimates (see function
  'grad' (numDeriv) for details)}

  \item{tol}{(numeric) value specifying the numeric
  tolerance for testing equality to zero}

  \item{quiet}{(logical) TRUE = suppress warning messages,
  e.g. for non-estimable contrasts}

  \item{opt}{(logical) TRUE = tries to optimize computation
  time by avoiding unnecessary computations for balanced
  datasets (see details).}

  \item{onlyDF}{(logical) TRUE = only the specified type of
  degrees of freedom are determined without carrying out
  the actual hypothesis test(s)}

  \item{...}{further parameters (for internal use
  actually)}
}
\value{
(numeric) vector or matrix with 4 elements/columns
corresponding to "Estimate", "t Value", "DF", and "Pr >
|t|".
}
\description{
This function performs t-Tests for one or multiple linear
combinations (contrasts) of estimated fixed effects.
}
\details{
Here, the same procedure as in SAS PROC MIXED
ddfm=satterthwaite (sat) is implemented. This
implementation was inspired by the code of function
'calcSatterth' of R-package 'lmerTest'. Thanks to the
authors for this nice implementation. \cr Note, that
approximated Satterthwaite degrees of freedom might differ
from 'lmerTest' and SAS PROC MIXED. Both use the inverse
Fisher-information matrix as approximation of the
variance-covariance matrix of variance components
(covariance parameters). Here, either the exact algorithm
for ANOVA-estimators of variance components, described in
Searle et. al (1992) p. 176, or the approximation presented
in Giesbrecht and Burns (19985) are used. For balanced
designs their will be no differences, usually. In case of
balanced designs, the Satterthwaite approximation is equal
to the degrees of freedom of the highest order random term
in the model (see examples).
}
\examples{
\dontrun{
data(dataEP05A2_2)
ub.dat <- dataEP05A2_2[-c(11,12,23,32,40,41,42),]
fit <- anovaMM(y~day/(run), ub.dat)
fe <- fixef(fit)
fe
lc.mat <- getL( fit, c("day1-day2", "day3-day6"))
lc.mat
test.fixef(fit, lc.mat)

# some inferential statistics about fixed effects estimates
L <- diag(nrow(fe))
rownames(L) <- colnames(L) <- rownames(fe)
test.fixef(fit, L)

# using different "residual" method determining DFs
test.fixef(fit, L, ddfm="res")

# having 'opt=TRUE' is a good idea to save time
# (in case of balanced designs)
data(VCAdata1)
datS3 <- VCAdata1[VCAdata1$sample==3,]
fit <- anovaMM(y~(lot+device)/(day)/run, datS3)
fit$VarCov <- vcovVC(fit)
fe <- fixef(fit)
L <- diag(nrow(fe))
rownames(L) <- colnames(L) <- rownames(fe)
system.time(tst1 <- test.fixef(fit, L))
system.time(tst2 <- test.fixef(fit, L, opt=FALSE))
tst1
tst2
}}
\author{
Andre Schuetzenmeister
\email{andre.schuetzenmeister@roche.com} inspired by
authors of R-package 'lmerTest'
}
\references{
Searle, S.R, Casella, G., McCulloch, C.E. (1992), Variance
Components, Wiley New York

Giesbrecht, F.G. and Burns, J.C. (1985), Two-Stage Analysis
Based on a Mixed Model: Large-Sample Asymptotic Theory and
Small-Sample Simulation Results, Biometrics 41, p. 477-486

SAS Help and Documentation PROC MIXED (MODEL-statement,
Option 'ddfm'), SAS Institute Inc., Cary, NC, USA
}
\seealso{
\code{\link{test.lsmeans}}, \code{\link{getL}}
}

