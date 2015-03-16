\name{fixef.VCA}
\alias{fixef.VCA}
\title{Extract Fixed Effects from 'VCA' Object.}
\usage{
\method{fixef}{VCA}(object, type = c("simple", "complex"),
  ddfm = c("contain", "residual", "satterthwaite"), tol = 1e-12,
  quiet = FALSE, ...)
}
\arguments{
  \item{object}{(VCA) object where fixed effects shall be
  extracted}

  \item{type}{(character) string or partial string,
  specifying whether to return "simple" (reduced) or a
  rather "complex" (more detailed) information about fixed
  effects}

  \item{ddfm}{(character) string specifying the method used
  for computing the degrees of freedom of the t-statistic.
  Only used when type="complex".  Available methods are
  "contain", "residual", and "satterthwaite".}

  \item{tol}{(numeric) value representing the numeric
  tolerance use in comparisons, values smaller than 'tol'
  will be considered equal to 0}

  \item{quiet}{(logical) TRUE = suppress warning messages,
  e.g. for non-estimable contrasts}

  \item{...}{additional parameters}
}
\description{
Conviniently extracting the 'FixedEffects' element of an
'VCA' object.
}
\details{
The default is to return the fixed effects estimates
together with their standard errors. If setting
'type="complex"' or to an abbreviation (e.g. "c")
additional inferential statistics on these estimates will
be returned, i.e. "t Value", "DF" and respective p-value
"Pr > |t|". One can choose one of three denominator degrees
of freedom ('ddfm')-methods.
}
\examples{
data(dataEP05A2_1)
fit <- anovaVCA(y~day/(run), dataEP05A2_1)
fixef(fit)

\dontrun{
# for complex models it might take some time computing complex output
data(VCAdata1)
fit <- anovaMM(y~(lot+device)/(day)/(run), VCAdata1[VCAdata1$sample==2,])
fixef(fit, "c")
}}

