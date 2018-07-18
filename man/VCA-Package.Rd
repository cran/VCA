\docType{package}
\name{VCA-package}
\alias{VCA}
\alias{VCA-package}
\title{(V)ariance (C)omponent (A)nalysis.}
\description{
 This package implements ANOVA-type estimation of variance components (VC) for linear mixed models (LMM), and provides 
 Restricted Maximum Likelihood (REML) estimation incorporating functionality of the \code{lme4} package. For models fitted by REML
 the typical VCA-table is derived, also containing the variances of VC, which are approximated by the method outlined in Giesbrecht & Burns (1985). 
 REML-estimation is available via functions \code{\link{remlVCA}} for variance component analysis (VCA) and \code{\link{remlMM}} for fitting
 general LMM.
 ANOVA-methodology is a special method of moments approach for estimating (predicting) variance components implemented in functions
 \code{\link{anovaMM}} and \code{\link{anovaVCA}}. The former represents a general, unrestricted approach to fitting
 linear mixed models, whereas the latter is tailored for performing a VCA on random models.
 Experiments of this type frequently occur in performance evaluation analyses of diagnostic tests or analyzers (devices)
 quantifying various types of precision (see e.g. guideline EP05-A2/A3 of the Clinical and Laboratory Standards Institute - CLSI). 
 
 The general Satterthwaite approximation of denominator degrees of freedom for tests of fixed effects (\code{\link{test.fixef}}) and
 LS Means (\code{\link{test.lsmeans}}) is implemented as used in SAS PROC MIXED.	Results differ for unbalanced designs because of the
 different approaches to estimating the covariance matrix of variance components. Here, two algorithms are implemented for models fitted
 via ANOVA, \eqn{1^{st}}{1st} the "exact" method described in Searle et. al (1992), \eqn{2^{nd}}{2nd} an approximation described in Giesbrecht & Burns (1985).
 The latter is also used for models fitted by REML. See \code{\link{test.fixef}} and \code{\link{getGB}} for details on this topic.
 
 Furthermore, the Satterthwaite approximation of degrees of freedom for individual VCs and total variance is implemented. These are employed
 in Chi-Squared tests of estimated variances against a claimed value (total, error), as well as in Chi-Squared based confidence intervals (CI) (see \code{\link{VCAinference}}). 
 Whenever ANOVA-type estimated VCs become negative, the default is to set them equal to 0. ANOVA mean squares used within the Satterthwaite approximation
 will be adapted to this situation by re-computing ANOVA mean squares (\eqn{s_{MS}}{MS}) as \eqn{s_{MS} = C * \sigma^{2 }}{MS = C * VC},
 where \eqn{C} is a coefficient matrix and a function of the design matrix and \eqn{\sigma^{2}}{VC} is the column-vector of adapted variance components.
 Total variance corresponds to a conservative estimate of the total variability in these cases, i.e. it will be larger than e.g. the total variance of the
 same model fitted by REML, because the negative VC will not contribute to total variance. See the documentation \code{\link{anovaVCA}} and \code{\link{anovaMM}}
 for details, specifically argument \code{NegVC}.
 
 Additionally to fitting linear mixed models and performing VCA-analyses, various plotting methods are implemented, e.g. a variability chart 
 visualizing the variability in sub-classes emerging from an experimental design (\code{\link{varPlot}}). Random effects and residuals can be transformed 
 and plotted using function \code{\link{plotRandVar}}. Standardization and studentization are generally available, Pearson-type transformation is only
 available for residuals. Plotting (studentized) random variates of a LMM should always be done to reveal potential problems of the fitted model, e.g. 
 violation of model assumptions and/or whether there are outlying observations.\cr
 
 There are two approaches to estimating ANOVA sums (SSQ) of squares. Method 1) constructs SSQ employing quadratic forms in \eqn{y}, 
 the column vector of observations. For sufficiently complex designs, i.e. large design matrices, unbalanced data, etc., this approach is significantly
 faster than using \code{anova.lm} instead. These matrices are also required in computing the covariance matrix of VCs following the "scm" method (see Searle et al. 1992).
 These quadratic forms are such that \eqn{s^{SS}_{i} = y^{T}A_{i}y}{SS_i = y' A_i y}, where \eqn{s^{SS}_{i}}{SS_i} is the sum of squares of the i-th factor, \eqn{A_i} is
 the \eqn{(N \times N)}{(N x N)} matrix generating the quadratic form corresponding to the i-th factor, and \eqn{y} is the column-vector of observations. For these quadratic forms
 holds \eqn{X^{T}A_{i}X = 0}{X' A_i X = 0} for mixed models, and \eqn{1^{T}A_{i}1 = 0}{1' A_i 1 = 0} for random models (\eqn{X} is the design matrix of fixed effects,
 reducing to a column-vector of 1s in the random model case). All this can be read in Searle et. al (1992), chapter 5. 
 Method 2) for constructing ANOVA-SSQ makes use of the SWEEP-operator (Goodnight 1979) reducing the computational time dramatically compared to method 1). It
 is the default method for ANOVA-type estimation. See the examples in \code{\link{anovaVCA}} for a comparison of methods. 
 
 Further reduction of the computation time can be achieved using Intel's \code{Math Kernel Library (MKL)}. When the package is loaded it will be automatically
 checked whether this is the case or not.
 
 In LS Means computation of fitted LMM it is possible to compute LS Means using specific values of covariables, which is equivalent to using option 'AT' in 
 the 'lsmeans'-statement of SAS PROC MIXED. It is also possible to apply other than the default weighting scheme for (fixed) factor-variables. See the details section
 in \code{\link{lsmeans}} and the description of argument \code{at}.
 
 Note: The 'UnitTests' directory within the package-directory contains a pre-defined test-suite which can be run by sourcing
       'RunAllTests.R' for user side testing (installation verification). It requires the 'RUnit' package and checks the numerical
        equivalence to reference results (SAS PROC MIXED method=type1/reml, SAS PROC VARCOMP) for balanced and unbalanced data and different experimental designs.}
\details{
  \tabular{ll}{
    Package: \tab VCA\cr
    Type: \tab Package \cr
    Version: \tab 1.3.4 \cr
    Date: \tab 2018-07-18\cr
    License:  \tab GPL (>=3) \cr
    LazyLoad: \tab yes \cr
  }
}
\author{
	Andre Schuetzenmeister \email{andre.schuetzenmeister@roche.com},
	Florian Dufey \email{florian.dufey@contractors.roche.com}
}
\references{
	Searle, S.R, Casella, G., McCulloch, C.E. (1992), Variance Components, Wiley New York
	
	Goodnight, J.H. (1979), A Tutorial on the SWEEP Operator, The American Statistician, 33:3, 149-158
	
	Giesbrecht, F.G. and Burns, J.C. (1985), Two-Stage Analysis Based on a Mixed Model: Large-Sample Asymptotic Theory and Small-Sample Simulation Results, Biometrics 41, p. 477-486
		
	Satterthwaite, F.E. (1946),  An Approximate Distribution of Estimates of Variance Components., Biometrics Bulletin 2, 110-114
		
	Gaylor,D.W., Lucas,H.L., Anderson,R.L. (1970), Calculation of Expected Mean Squares by the Abbreviated Doolittle and Square Root Methods., Biometrics 26 (4): 641-655
	
	SAS Help and Documentation PROC MIXED, SAS Institute Inc., Cary, NC, USA
		
}
\keyword{package}
