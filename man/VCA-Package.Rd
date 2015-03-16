\docType{package}
\name{VCA-package}
\alias{VCA}
\alias{VCA-package}
\title{(V)ariance (C)omponent (A)nalysis.}
\description{This package implements ANOVA-type estimation of variance components for linear mixed models.
	     It uses a method of moments approach to estimate (predict) variance components (VC) implemented in functions
	     \code{\link{anovaMM}} and \code{\link{anovaVCA}}. The former represents a general, unrestricted approach to fitting
	     linear mixed models, whereas the latter is tailored for performing a Variance Component Analysis (VCA) on random models.
	     Experiments of this type frequently occur in performance evaluation analyses of diagnostic tests or analyzers (devices)
	     quantifying various types of precision (see e.g. guideline EP05-A2/A3 of the Clinical and Laboratory Standards Institute - CLSI). 
	     
	     The general Satterthwaite approximation of denominator degrees of freedom for tests of fixed effects is implemented as used in
	     SAS PROC MIXED. Results differ for unbalanced designs because of the different approaches to estimating the covariance matrix of
	     variance components. Here, two algorithms are implemented, \eqn{1^{st}}{1st} the "exact" method described in Searle et. al (1992), 
	     \eqn{2^{nd}}{2nd} an approximation described in Giesbrecht and Burns (1985). See \code{\link{test.fixef}} for details on this topic.
	     
         Furthermore, the Satterthwaite approximation of degrees of freedom for individual VCs and total variance is implemented. These are employed
         in Chi-Squared tests of estimated variances against a claimed value (total, error), as well as in Chi-Squared based confidence intervals (CI) (see \code{\link{VCAinference}}). 
	     Whenever ANOVA-type estimated VCs become negative, the default is to set them to 0. ANOVA mean squares used within the Satterthwaite approximation
	     will be adapted to this situation by re-computing ANOVA mean squares (\eqn{s_{MS}}{MS}) as \eqn{s_{MS} = C * \sigma^{2 }}{MS = C * VC},
	     where \eqn{C} is a coefficient matrix and a function of the design matrix and \eqn{\sigma^{2}}{VC} is the column-vector of adapted variance components.
	     Additionally to fitting linear mixed models and performing VCA-analyses, various plotting methods are implemented, e.g. a variability chart 
	     visualizing the variability in sub-classes emerging from an experimental design (\code{\link{varPlot}}).\cr
	     
	     ANOVA sums of squares are constructed employing quadratic forms in \eqn{y}, the column vector of observations. For sufficiently complex
	     designs, i.e. large design matrices, unbalanced data, ..., this approach is significantly faster then using \code{anova.lm} instead. These
	     matrices are also needed in computing the covariance matrix of VCs following the "scm" method (see Searle et al. 1992).
	     These quadratic forms are such that \eqn{s^{SS}_{i} = y^{T}A_{i}y}{SS_i = y' A_i y}, where \eqn{s^{SS}_{i}}{SS_i} is the sum of squares of the i-th factor, \eqn{A_i} is
	     the \eqn{(N \times N)}{(N x N)} matrix generating the quadratic form corresponding to the i-th factor, and \eqn{y} is the column-vector of observations. For these quadratic forms
	     holds \eqn{X^{T}A_{i}X = 0}{X' A_i X = 0} for mixed models, and \eqn{1^{T}A_{i}1 = 0}{1' A_i 1 = 0} for random models (\eqn{X} is the design matrix of fixed effects,
	     reducing to a column-vector of 1s in the random model case). All this can be read in Searle et. al (1992), chapter 5. 
	     
	     Note: The 'UnitTests' directory within the package-directory contains a pre-defined test-suite which can be run by sourcing
             'RunAllTests.R' for user side testing (installation verification). It requires the 'RUnit' package and checks the numerical
             equivalence to reference results (SAS PROC MIXED method=type1, SAS PROC VARCOMP) for balanced and unbalanced data and different experimental designs.}
\details{
  \tabular{ll}{
    Package: \tab VCA\cr
    Type: \tab Package \cr
    Version: \tab 1.1 \cr
    Date: \tab 2015-03-16\cr
    License:  \tab GPL (>=3) \cr
    LazyLoad: \tab yes \cr
  }
}
\author{
	Andre Schuetzenmeister 
	\email{andre.schuetzenmeister@roche.com}
}
\keyword{package}
