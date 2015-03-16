\name{anovaVCA}
\alias{anovaVCA}
\title{ANOVA-Type Estimation of Variance Components for Random Models.}
\usage{
anovaVCA(form, Data, NegVC = FALSE, VarVC.method = c("scm", "gb"))
}
\arguments{
  \item{form}{(formula) specifying the model to be fit, a
  response variable left of the '~' is mandatory}

  \item{Data}{(data.frame) storing all variables referenced
  in 'form'}

  \item{NegVC}{(logical) FALSE = negative variance
  component estimates (VC) will be set to 0 and they will
  not contribute to the total variance (as done in SAS PROC
  NESTED, conservative estimate of total variance). The
  original ANOVA estimates can be found in attribute
  'VCoriginal'.  The degrees of freedom of the total
  variance are based on adapted mean squares (MS), i.e.
  adapted MS are computed as \eqn{D * VC}, where VC is the
  column vector with negative VCs set to 0. \cr TRUE =
  negative variance component estimates will not be set to
  0 and they will contribute to the total variance
  (original definition of the total variance).}

  \item{VarVC.method}{(character) string specifying whether
  to use the algorithm given in the \eqn{1^{st}}{1st}
  reference ("scm") or in the \eqn{2^{nd}}{2nd} refernce
  ("gb"). Method "scm" (Searle, Casella, McCulloch) is the
  exact algorithm but slower, "gb" (Giesbrecht, Burns) is
  termed "rough approximation" by the authors, but
  sufficiently exact compared to e.g. SAS PROC MIXED
  (method=type1) which uses the inverse of the
  Fisher-Information matrix as approximation. For balanced
  designs all methods give identical results, only in
  unbalanced designs differences occur.}
}
\value{
(object) of class 'VCA'
}
\description{
This function equates observed ANOVA sums of squares
(\eqn{SS}) to their expected values and solves the
resulting system of linear equations for variance
components. ANOVA SS are computed via matrices
\eqn{A_i}{A_i} expressing them as quadratic forms in
\eqn{y} as \eqn{ss_i = y^{T}A_{i}y}{ss_i = y' * A_i * y}.
Matrices \eqn{A_i} are also used to compute the
variance-coveriance matrix of variance components (VC).
}
\details{
Function \code{anovaVCA} is tailored for performing
Variance Component Analyses (VCA) for random models,
assuming all VCs as factor variables, i.e. their levels
correspond to distinct columns in the design matrix (dummy
variables). Any predictor variables are automatically
converted to factor variables, since continuous variables
may not be used on the right side of the formula 'form'.

For diagnostics, a key parameter is "precision", i.e. the
accuracy of a quantification method influenced by varying
sources of random error. This type of experiments is
requested by regulatory authorities to proof the quality of
diagnostic test, e.g. quantifying intermediate precision
according to CLSI guideline EP5-A2/A3. No, fixed effects
are allowed besides the intercept. Whenever fixed effects
are part of the model to be analyzed, use function
\code{\link{anovaMM}} instead.

Function \code{anovaVCA} represents a special form of the
"method of moments" approach applicable to arbitrary random
models either balanced or unbalanced. The system of linear
equations, which is built from the ANOVA-Type 1 sums of
squares, is closely related to the method used by SAS PROC
VARCOMP, where ANOVA mean squares (\eqn{MS}) are used (see
\code{\link{getCmatrix}}). The former can be written as
\eqn{ss = C * s} and the latter as \eqn{ms = D * s}, where
\eqn{C} and \eqn{D} denote the respective coefficient
matrices, \eqn{s} the column-vector of variance components
(VC) to be estimated/predicted, and \eqn{ss} and \eqn{ms}
the column vector of ANOVA sum of squares, respectively,
mean squares. Mutliplying element \eqn{d_ij} of matrix
\eqn{D} by element \eqn{c_in} of matrix \eqn{C} (\eqn{i,j =
1,...,n}), results in matrix \eqn{C}. Thus, \eqn{C} can
easily be converted to \eqn{D} by the inverse operation.
Matrix \eqn{D} is used to estimate total degrees of freedom
(DF) according to Satterthwaite (1946).

If computing matrices \eqn{A_i} generating ANOVA \eqn{SS}
as quadratic forms in \eqn{y}, i.e. \eqn{ss_i = y' * A_i *
y} has two benefits. It allows obtaing ANOVA \eqn{SS} as
well as computing the variance-covariance of VCs.
Furthermore, it is much faster than fitting the linear
model via \code{\link{lm}} and calling function
\code{\link{anova}} on the 'lm' object for complex models,
where complex refers to the number of columns of the design
matrix and the degree of unbalancedness. Degrees of freedom
for the \eqn{i}-th term are obtained by function
\code{\link{anovaDF}}. Besides the algorithm for obtaining
the variance-covariance matrix of VCs described in Searle
et. al (1992), which is termed "exact", a \eqn{2^{nd}}{2nd}
option for computing this matrix is implemented, which is
described in Giesbrecht and Burns (1985). It can be chosen
by setting 'VarVC.method="gb"'.
}
\examples{
# load data (CLSI EP05-A2 Within-Lab Precision Experiment)
data(dataEP05A2_2)

# perform ANOVA-estimation of variance components
res <- anovaVCA(y~day/run, dataEP05A2_2)
res

# desing with two main effects (ignoring the hierarchical structure of the design)
anovaVCA(y~day+run, dataEP05A2_2)

# compute confidence intervals, perform F- and Chi-Squared tests
INF <- VCAinference(res, total.claim=3.5, error.claim=2)
INF

### load data from package
data(VCAdata1)

data_sample1 <- VCAdata1[VCAdata1$sample==1,]

### plot data for visual inspection (there is no variance between runs on a day)
varPlot(y~lot/day/run, data_sample1)

### estimate VCs for 4-level hierarchical design (error counted) for sample_1 data
anovaVCA(y~lot/day/run, data_sample1)

### using different model (ignoring the hierarchical structure of the design)
anovaVCA(y~lot+day+lot:day:run, data_sample1)

### same model with unbalanced data
anovaVCA(y~lot+day+lot:day:run, data_sample1[-c(1,11,15),])

### use the numerical example from the CLSI EP05-A2 guideline (p.25)
data(dataEP05A2_example)
res.ex <- anovaVCA(result~day/run, dataEP05A2_example)

### also perform Chi-Squared tests
### Note: in guideline claimed SD-values are used, here, claimed variances are used
VCAinference(res.ex, total.claim=3.4^2, error.claim=2.5^2)

#' # load another example dataset and extract the "sample==1" subset
data(VCAdata1)
sample1 <- VCAdata1[which(VCAdata1$sample==1),]

# generate an additional factor variable and random errors according to its levels
sample1$device <- gl(3,28,252)
set.seed(505)
sample1$y <- sample1$y + rep(rep(rnorm(3,,.25), c(28,28,28)),3)

# fit a crossed-nested design with main factors 'lot' and 'device'
# and nested factors 'day' and 'run' nested below
res1 <- anovaVCA(y~(lot+device)/day/run, sample1)
res1
}
\author{
Andre Schuetzenmeister
\email{andre.schuetzenmeister@roche.com}
}
\references{
Searle, S.R, Casella, G., McCulloch, C.E. (1992), Variance
Components, Wiley New York

Giesbrecht, F.G. and Burns, J.C. (1985), Two-Stage Analysis
Based on a Mixed Model: Large-Sample Asymptotic Theory and
Small-Sample Simulation Results, Biometrics 41, p. 477-486

Satterthwaite, F.E. (1946), An Approximate Distribution of
Estimates of Variance Components., Biometrics Bulletin 2,
110-114

Gaylor,D.W., Lucas,H.L., Anderson,R.L. (1970), Calculation
of Expected Mean Squares by the Abbreviated Doolittle and
Square Root Methods., Biometrics 26 (4): 641-655

SAS Help and Documentation PROC MIXED, SAS Institute Inc.,
Cary, NC, USA

SAS Help and Documentation PROC VARCOMP, SAS Institute
Inc., Cary, NC, USA
}
\seealso{
\code{\link{anovaMM}}, \code{\link{stepwiseVCA}},
\code{\link{print.VCA}}, \code{\link{VCAinference}},
\code{\link{getCmatrix}}, \code{\link{ranef}},
\code{\link{plotRandVar}}
}

