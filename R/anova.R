# TODO: Add comment
# 
# Author: schueta6
###############################################################################


#'ANOVA-Type Estimation of Variance Components for Random Models
#'
#'This function equates observed ANOVA Type-I sums of squares (\eqn{SS}) to their expected values and solves the resulting system of linear equations
#'for variance components. 
#'
#'For diagnostics, a key parameter is "precision", i.e. the accuracy of a quantification method influenced by varying sources of random error. 
#'This type of experiments is requested by regulatory authorities to proof the quality of diagnostic tests, e.g. quantifying intermediate
#'precision according to CLSI guideline EP5-A2/A3. No, fixed effects are allowed besides the intercept. 
#'Whenever fixed effects are part of the model to be analyzed, use function \code{\link{anovaMM}} instead.
#'
#'Function \code{anovaVCA} is tailored for performing Variance Component Analyses (VCA) for random models, assuming all VCs as factor variables, i.e. their levels
#'correspond to distinct columns in the design matrix (dummy variables). Any predictor variables are automatically converted to factor variables, since continuous
#'variables may not be used on the right side of the formula 'form'.
#'
#'ANOVA \eqn{SS} are computed employing the SWEEP-operator (Goodnight 1979, default).
#'according to Searle et al. (1992) which corresponds to \code{VarVC.method="scm"}.
#'
#'Function \code{anovaVCA} represents a special form of the "method of moments" approach applicable to arbitrary random models either balanced or unbalanced.
#'The system of linear equations, which is built from the ANOVA Type-I sums of squares, is closely related to the method used 
#'by SAS PROC VARCOMP, where ANOVA mean squares (\eqn{MS}) are used. The former can be written as \eqn{ss = C * s}
#'and the latter as \eqn{ms = D * s}, where \eqn{C} and \eqn{D} denote the respective coefficient matrices, \eqn{s} the column-vector
#'of variance components (VC) to be estimated/predicted, and \eqn{ss} and \eqn{ms} the column vector of ANOVA sum of squares, respectively, mean squares. 
#'Mutliplying element \eqn{d_ij} of matrix \eqn{D} by element \eqn{c_in} of matrix \eqn{C} (\eqn{i,j = 1,...,n}), results in 
#'matrix \eqn{C}. Thus, \eqn{C} can easily be converted to \eqn{D} by the inverse operation. Matrix \eqn{D} is used to estimate
#'total degrees of freedom (DF) according to Satterthwaite (1946).
#'
#'The method for computing ANOVA Type-I \eqn{SS} is much faster than fitting the linear model via \code{\link{lm}} and calling function \code{\link{anova}} on the 'lm' object
#'for complex models, where complex refers to the number of columns of the design matrix and the degree of unbalancedness. \eqn{DF} are directly derived from the SWEEP-operator as the number of linearly independent
#'columns of the partial design matrix corresponding to a specific \eqn{VC}.
#'
#'@param form			(formula) specifying the model to be fit, a response variable left of the '~' is mandatory
#'@param Data			(data.frame) containing all variables referenced in 'form'
#'@param by				(factor, character) variable specifying groups for which the analysis should be performed individually,
#'						i.e. by-processing
#'@param NegVC			(logical) FALSE = negative variance component estimates (VC) will be set to 0 and they will not contribute to the total variance 
#'						(as done in SAS PROC NESTED, conservative estimate of total variance). The original ANOVA estimates can be found in element 'VCoriginal'. 
#'						The degrees of freedom of the total variance are based on adapted mean squares (MS), i.e. adapted MS are computed as \eqn{D * VC}, where VC is 
#'						the column vector with negative VCs set to 0. \cr
#'						TRUE = negative variance component estimates will not be set to 0 and they will contribute to the total variance (original definition of the total variance).
#'@param VarVC.method	(character) string specifying whether to use the algorithm given in Searle et al. (1992) which corresponds to \code{VarVC.method="scm"} or in
#'						Giesbrecht and Burns (1985) which can be specified via "gb". Method "scm" (Searle, Casella, McCulloch)
#'						is the exact algorithm, "gb" (Giesbrecht, Burns) is termed "rough approximation"
#'						by the authors, but sufficiently exact compared to e.g. SAS PROC MIXED (method=type1) which
#'						uses the inverse of the Fisher-Information matrix as approximation. For balanced designs all
#'						methods give identical results, only in unbalanced designs differences occur. 
#'@param MME			(logical) TRUE = (M)ixed (M)odel (E)quations will be solved, i.e. 'VCA' object will have additional elements
#'						"RandomEffects", "FixedEffects", "VarFixed" (variance-covariance matrix of fixed effects) and the "Matrices"
#'						element has addional elements corresponding to intermediate results of solving MMEs.
#'						FALSE = do not solve MMEs, which reduces the computation time for very complex models significantly.
#'@param quiet			(logical) TRUE = will suppress any warning, which will be issued otherwise 
#'@param order.data		(logical) TRUE = class-variables will be ordered increasingly, FALSE = ordering of class-variables
#'						will remain as is
#'
#'@return (object) of class 'VCA'
#'
#'@aliases anovaVCA
#'
#'@seealso \code{\link{anovaMM}}, \code{\link{remlVCA}}, \code{\link{remlMM}}, \code{\link{print.VCA}}, \code{\link{VCAinference}}, 
#'\code{\link{ranef}}, \code{\link{plotRandVar}}, \code{\link{stepwiseVCA}}
#'
#'@references 
#'
#'Searle, S.R, Casella, G., McCulloch, C.E. (1992), Variance Components, Wiley New York
#'
#'Goodnight, J.H. (1979), A Tutorial on the SWEEP Operator, The American Statistician, 33:3, 149-158
#'
#'Giesbrecht, F.G. and Burns, J.C. (1985), Two-Stage Analysis Based on a Mixed Model: Large-Sample
#'Asymptotic Theory and Small-Sample Simulation Results, Biometrics 41, p. 477-486
#'
#'Satterthwaite, F.E. (1946),  An Approximate Distribution of Estimates of Variance Components., 
#'Biometrics Bulletin 2, 110-114
#'
#'Gaylor,D.W., Lucas,H.L., Anderson,R.L. (1970), Calculation of Expected Mean Squares by the Abbreviated Doolittle and Square Root Methods., 
#'Biometrics 26 (4): 641-655
#'
#'SAS Help and Documentation PROC MIXED, SAS Institute Inc., Cary, NC, USA
#'
#'SAS Help and Documentation PROC VARCOMP, SAS Institute Inc., Cary, NC, USA
#'
#'@author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#'
#'@examples 
#'\dontrun{
#'
#'# load data (CLSI EP05-A2 Within-Lab Precision Experiment) 
#'data(dataEP05A2_2)
#'
#'# perform ANOVA-estimation of variance components
#'res <- anovaVCA(y~day/run, dataEP05A2_2)
#'res
#'
#'# design with two main effects (ignoring the hierarchical structure of the design)
#'anovaVCA(y~day+run, dataEP05A2_2)
#'
#'# compute confidence intervals, perform F- and Chi-Squared tests
#'INF <- VCAinference(res, total.claim=3.5, error.claim=2)
#'INF
#'
#'### load data from package
#'data(VCAdata1)
#'
#'data_sample1 <- VCAdata1[VCAdata1$sample==1,]
#'
#'### plot data for visual inspection
#'varPlot(y~lot/day/run, data_sample1)
#'
#'### estimate VCs for 4-level hierarchical design (error counted) for sample_1 data
#'anovaVCA(y~lot/day/run, data_sample1)
#'
#'### using different model (ignoring the hierarchical structure of the design)
#'anovaVCA(y~lot+day+lot:day:run, data_sample1)
#'
#'### same model with unbalanced data
#'anovaVCA(y~lot+day+lot:day:run, data_sample1[-c(1,11,15),])
#'
#'### use the numerical example from the CLSI EP05-A2 guideline (p.25)
#'data(Glucose,package="VCA")
#'res.ex <- anovaVCA(result~day/run, Glucose)
#'
#'### also perform Chi-Squared tests
#'### Note: in guideline claimed SD-values are used, here, claimed variances are used
#'VCAinference(res.ex, total.claim=3.4^2, error.claim=2.5^2)
#'
#'### now use the six sample reproducibility data from CLSI EP5-A3
#'### and fit per sample reproducibility model
#'data(CA19_9)
#'fit.all <- anovaVCA(result~site/day, CA19_9, by="sample")
#'
#'reproMat <- data.frame(
#'Sample=c("P1", "P2", "Q3", "Q4", "P5", "Q6"),
#'Mean= c(fit.all[[1]]$Mean, fit.all[[2]]$Mean, fit.all[[3]]$Mean, 
#'fit.all[[4]]$Mean, fit.all[[5]]$Mean, fit.all[[6]]$Mean),
#'Rep_SD=c(fit.all[[1]]$aov.tab["error","SD"], fit.all[[2]]$aov.tab["error","SD"],
#'fit.all[[3]]$aov.tab["error","SD"], fit.all[[4]]$aov.tab["error","SD"],
#'fit.all[[5]]$aov.tab["error","SD"], fit.all[[6]]$aov.tab["error","SD"]),
#'Rep_CV=c(fit.all[[1]]$aov.tab["error","CV[%]"],fit.all[[2]]$aov.tab["error","CV[%]"],
#'fit.all[[3]]$aov.tab["error","CV[%]"],fit.all[[4]]$aov.tab["error","CV[%]"],
#'fit.all[[5]]$aov.tab["error","CV[%]"],fit.all[[6]]$aov.tab["error","CV[%]"]),
#'WLP_SD=c(sqrt(sum(fit.all[[1]]$aov.tab[3:4,"VC"])),sqrt(sum(fit.all[[2]]$aov.tab[3:4, "VC"])),
#'sqrt(sum(fit.all[[3]]$aov.tab[3:4,"VC"])),sqrt(sum(fit.all[[4]]$aov.tab[3:4, "VC"])),
#'sqrt(sum(fit.all[[5]]$aov.tab[3:4,"VC"])),sqrt(sum(fit.all[[6]]$aov.tab[3:4, "VC"]))),
#'WLP_CV=c(sqrt(sum(fit.all[[1]]$aov.tab[3:4,"VC"]))/fit.all[[1]]$Mean*100,
#'sqrt(sum(fit.all[[2]]$aov.tab[3:4,"VC"]))/fit.all[[2]]$Mean*100,
#'sqrt(sum(fit.all[[3]]$aov.tab[3:4,"VC"]))/fit.all[[3]]$Mean*100,
#'sqrt(sum(fit.all[[4]]$aov.tab[3:4,"VC"]))/fit.all[[4]]$Mean*100,
#'sqrt(sum(fit.all[[5]]$aov.tab[3:4,"VC"]))/fit.all[[5]]$Mean*100,
#'sqrt(sum(fit.all[[6]]$aov.tab[3:4,"VC"]))/fit.all[[6]]$Mean*100),
#'Repro_SD=c(fit.all[[1]]$aov.tab["total","SD"],fit.all[[2]]$aov.tab["total","SD"],
#'fit.all[[3]]$aov.tab["total","SD"],fit.all[[4]]$aov.tab["total","SD"],
#'fit.all[[5]]$aov.tab["total","SD"],fit.all[[6]]$aov.tab["total","SD"]),
#'Repro_CV=c(fit.all[[1]]$aov.tab["total","CV[%]"],fit.all[[2]]$aov.tab["total","CV[%]"],
#'fit.all[[3]]$aov.tab["total","CV[%]"],fit.all[[4]]$aov.tab["total","CV[%]"],
#'fit.all[[5]]$aov.tab["total","CV[%]"],fit.all[[6]]$aov.tab["total","CV[%]"]))
#'
#'for(i in 3:8) reproMat[,i] <- round(reproMat[,i],digits=ifelse(i%%2==0,1,3))
#'reproMat
#'
#'# now plot the precision profile over all samples
#'plot(reproMat[,"Mean"], reproMat[,"Rep_CV"], type="l", main="Precision Profile CA19-9",
#'xlab="Mean CA19-9 Value", ylab="CV[%]")
#'grid()
#'points(reproMat[,"Mean"], reproMat[,"Rep_CV"], pch=16)
#'
#'# load another example dataset and extract the "sample==1" subset
#'data(VCAdata1)
#'sample1 <- VCAdata1[which(VCAdata1$sample==1),]
#'
#'# generate an additional factor variable and random errors according to its levels
#'sample1$device <- gl(3,28,252)                                      
#'set.seed(505)
#'sample1$y <- sample1$y + rep(rep(rnorm(3,,.25), c(28,28,28)),3)     
#'
#'# fit a crossed-nested design with main factors 'lot' and 'device' 
#'# and nested factors 'day' and 'run' nested below 
#'res1 <- anovaVCA(y~(lot+device)/day/run, sample1) 
#'res1
#'
#'# fit same model for each sample using by-processing
#'lst <- anovaVCA(y~(lot+device)/day/run, VCAdata1, by="sample")
#'lst
#'
#'# now fitting a nonsense model on the complete dataset "VCAdata1" 
#'# where the SWEEP-operator is the new default since package version 1.2
#'# takes ~5s
#'system.time(res.sw <- anovaVCA(y~(sample+lot+device)/day/run, VCAdata1))
#'# applying functions 'anova' and 'lm' in the same manner takes ~ 265s
#'system.time(res.lm <- anova(lm(y~(sample+lot+device)/day/run, VCAdata1)))
#'res.sw
#'res.lm
#'}


anovaVCA <- function(	form, Data, by=NULL, NegVC=FALSE,
		VarVC.method=c("scm","gb"), MME=FALSE, quiet=FALSE, order.data=TRUE)
{
	if(!is.null(by))
	{
		stopifnot(is.character(by))
		stopifnot(by %in% colnames(Data))
		stopifnot(is.factor(by) || is.character(by))
		
		levels  <- unique(Data[,by])
		res <- lapply(levels, function(x) {
					tmp.res <- try(anovaVCA(form=form, Data[Data[,by] == x,], NegVC=NegVC, VarVC.method=VarVC.method, MME=MME, quiet=quiet), silent=TRUE)
					if(is(tmp.res, "try-error") && !quiet)
						warning(paste0("Error for '", by, ".", x, "':\n\t", attr(tmp.res, "condition")$message))
					tmp.res
		})
		names(res) <- paste(by, levels, sep=".")
		return(res)
	}
	
	stopifnot(class(form) == "formula")
	stopifnot(identical(class(Data),"data.frame"))
	stopifnot(nrow(Data) > 2)                                               # at least 2 observations for estimating a variance
	stopifnot(is.logical(NegVC))
	
	if(is.null(.GlobalEnv$msgEnv))											# may be removed after loading the package
		msgEnv <<- new.env(parent=emptyenv())
	
	VarVC.method <- match.arg(VarVC.method)
	tobj <- terms(form, simplify=TRUE, keep.order=TRUE)						# expand nested factors if necessary, retain ordering of terms in the formula
	if(!attr(tobj, "response"))
		stop("You need to include a response variable in the formula!")
	form <- formula(tobj)
	Data <- orderData(Data, tobj, quiet=quiet,
			exclude.numeric=FALSE, order.data=order.data)					# convert all variables to factors
	resp <- as.character(tobj[[2]])											
	res       	<- list()
	res$call  	<- match.call()
	res$Type  	<- "Random Model"
	if(is.null(rownames(Data)))												# assign rownames if missing to record results of re-ordering taking place before the analysis
		rownames(Data) <- 1:nrow(Data)
	res$data  	<- Data
	res$terms 	<- tobj
	res$VCnames <- c(attr(tobj, "term.labels"), "error")
	res$terms.classes <- sapply(Data[,rownames(attr(tobj, "factors"))[-1]], class)
	int <- res$intercept <- attr(tobj, "intercept") == 1	
	
	res$response <- resp
	
	stopifnot(resp %in% colnames(Data))
	stopifnot(is.numeric(Data[,resp]))
	
	Ndata <- nrow(Data)
	
	rmInd <- integer()
	
	resp.NA <- is.na(Data[,resp])
	
	if(any(resp.NA))
	{    
		rmInd <- c(rmInd, which(resp.NA))
		if(!quiet)
			message("There are ", length(which(resp.NA))," missing values for the response variable (obs: ", paste(which(resp.NA), collapse=", "), ")!")
	}    
	
	fac  <- attr(tobj, "term.labels")
	vars <- rownames(attr(tobj, "factors"))[-1]                             # remove response
	Nvc  <- length(fac) + 1 
	
	if(!is.null(vars))
	{
		for(i in vars)                                                          # convert all nested factors as factor-objects
		{
			if(length(unique(Data[,i])) == Ndata)
				stop("Variable '", i,"': number of levels of each grouping factor must be < number of observations!")
			
			if( any(is.na(Data[,i])))
			{
				NAind <- which(is.na(Data[,i]))
				rmInd <- c(rmInd, NAind)
				if(!quiet)
					message("Variable '", i,"' has ",length(NAind)," missing values (obs: ", paste(NAind, collapse=", "), ")!" )
			}
			if(length(levels(Data[,i])) >= 0.75*nrow(Data) && !quiet)
				message("Variable >>> ", i," <<< has at least 0.75 * nrow(Data) levels!")
		}
		rmInd <- unique(rmInd)
	}	
	
	if(length(rmInd) > 0)
		Data <- Data[-rmInd,]
	
	attr(res$data, "analysis.order") <- rownames(Data)						# actually record the ordering of the data 
	
	Mean <- mean(Data[,resp], na.rm=TRUE)                                   # mean computed after removing incomplete observations
	vcol <- rownames(attr(tobj, "factors"))
	if(is.null(vcol))
		vcol <- resp
	Data <- na.omit(Data[,vcol, drop=F])
	Nobs <- N <- nrow(Data)
	
	allObsEqual <- FALSE
	if(all(Data[,resp] == Data[1,resp]))										# no variance detectable?
	{
		allObsEqual <- TRUE
		Data.org <- Data
		Data[,resp] <- Data[,resp] + rnorm(nrow(Data))
		
		if(!quiet)
			warning("All values of response variable ", paste0("'", resp, "'"), " are equal!")
	}
	
	tmp.res <- getSSQsweep(Data, tobj)									# determine ANOVA Type-1 sum of squares using sweeping

	gc(verbose=FALSE)
	
	Lmat    <- tmp.res$Lmat
	CVC     <- tmp.res$CVC
	aov.tab <- tmp.res$aov.tab												# basic ANOVA-table
	rownames(aov.tab) <- res$VCnames
	DF <- aov.tab[,"DF"]
	SS <- aov.tab[,"SS"]
	VC <- aov.tab[,"VC"]
	C <-CVC$C
	Ci<-CVC$Ci
	if(allObsEqual)
	{
		aov.tab[,c("SS", "MS","VC","SD")] <- 0
		Data <- Data.org
		SS[1:length(SS)] <- 0
	}
	if(allObsEqual)															# prevent singlular C-matrix in case of all obs equal situations
	{
		ind <- which(upper.tri(C, TRUE))
		ind <- ind[which(C[ind] %in% c(0, NaN))]
		if(length(ind) > 0)
			C[ind] <- 1
	}
	
	C2  <- apply(C, 2, function(x) x/DF)                                    # coefficient matrix for mean squares (MS)
	Ci2 <- solve(C2)
#	colnames(aov.tab) <- c("DF", "SS", "MS", "VC")
	
	IndNegVC <- 1:nrow(aov.tab)                                             # should negative VC-estimates contribute to total variance 
	NegVCmsg <- ""                                                          # message text in case that VC is negative
	VCoriginal <- aov.tab[,"VC"]                                            # keep original ANOVA-estimates of VCs
	if(!NegVC)
	{
		IndNegVC <- which(aov.tab[,"VC"] < 0)  
		if(length(IndNegVC) > 0)                                            # there are negative VC
		{
			aov.tab[IndNegVC, "VC"] <- 0                                    # set negative VC to 0
			NegVCmsg <- "* VC set to 0"
		}
	}    
	totVC   <- sum(aov.tab$VC)
	aov.tab <- rbind(total=c(NA, NA, NA, totVC, sqrt(totVC)), aov.tab)    
	
	aov.tab["total", "DF"] <- SattDF(c(C2 %*% aov.tab[-1,"VC"]), Ci=Ci2, DF=DF)   	# will automatically adapt ANOVA-MS if any VCs were set to 0 
	
#	suppressWarnings(aov.tab <- cbind(aov.tab, SD=sqrt(aov.tab[,"VC"])))    		# warnings suppressed because sqrt of negative numbers doese not exists
	aov.tab <- cbind(aov.tab, "CV[%]"=aov.tab[,"SD"]*100/Mean)
	aov.tab <- cbind(aov.tab, "%Total"=aov.tab[,"VC"]*100/totVC)
	aov.tab <- aov.tab[,c("DF", "SS", "MS", "VC", "%Total", "SD", "CV[%]")]    
	aov.tab <- as.matrix(aov.tab)
	aov.tab <- apply(aov.tab, 1:2, function(x) ifelse(is.nan(x), NA, x))
	
	res$aov.tab <- aov.tab
	res$rmInd   <- rmInd
	
	res$Nvc			 <- Nvc
	res$VarVC.method <- VarVC.method
	res$Mean         <- Mean
	res$Nobs         <- Nobs
	res$EstMethod    <- "ANOVA"
	res$VCoriginal   <- VCoriginal
	res$NegVCmsg     <- NegVCmsg
	res$NegVC        <- NegVC
	if(Nobs != Ndata)
		res$Nrm <- Ndata - Nobs                                	# save number of observations that were removed due to missing data
	
	Lmat$C.SS	<- C
	Lmat$C.MS	<- C2
	Lmat$Ci.SS  <- Ci
	Lmat$Ci.MS  <- Ci2
	Lmat$VCvar  <- CVC$VCvar
	Lmat$y      <- matrix(Data[, resp], ncol=1)								# column vector of observations
	Lmat$X      <- matrix(ifelse(int, 1, 0), ncol=1, nrow=nrow(Data))		# design matrix of fixed effects
	colnames(Lmat$X) <- "int"
	if(int)
	{
		res$fe.assign <- 0
		attr(res$fe.assign, "terms") <- "int"
	}
	Lmat$rf.ind 	<- 1:(nrow(aov.tab)-1)
	Lmat$VCall  	<- aov.tab[-1, "VC"]									# results based on scaled values
	res$Matrices 	<- Lmat                                   	# needed for variance-covariance matrix of VCs
	
	res$formula  <- form
	res$balanced <- ifelse(isBalanced(form, Data), "balanced","unbalanced" )
	class(res) <- "VCA"    
	
	if(MME)
		res <- solveMME(res)
	
	if(allObsEqual)
		res$aov.tab[,"%Total"] <- 0	
	
	gc(verbose=FALSE)										# trigger garbage collection
	return(res)
}





#'ANOVA-Type Estimation of Mixed Models
#'
#'Estimate/Predict random effects employing ANOVA-type estimation and obtain generalized least squares estimates
#'of fixed effects for any linear mixed model including random models and linear models.
#'
#'A Linear Mixed Model, noted in standard matrix notation, can be written as {y = Xb + Zg + e}, where
#'\eqn{y} is the column vector of observations, \eqn{X} and \eqn{Z}{Z} are design matrices assigning fixed (\eqn{b}),
#'respectively, random (\eqn{g}) effects to observations, and \eqn{e} is the column vector of residual errors.
#'Whenever there is an intercept in the model, i.e. the substring "-1" is not part of the model formula, the same
#'restriction as in SAS PROC MIXED is introduced setting the last fixed effect equal to zero. Note, that the results
#'of an linear contrasts are not affected by using an intercept or not, except that constrained fixed effects cannot
#'be part of such contrasts (one could use the intercept estimated instead).
#'
#'Here, no further restrictions on the type of model are made. One can fit mixed models as well as random models, which
#'constitute a sub-set of mixed models (intercept being the only fixed effect). Variables must be either of type "numeric"
#'or "factor". "character" variables are automatically converted to factors and the response variable has to be numeric, of course.
#'In case that 'class(Data[,i])' is neither one of these three options, an error is issued.
#'Even simple linear models can be fitted, i.e. models without a random part (without \eqn{Zg}{Zg}) besides the
#'residual errors. In this case, an Analysis of Variance (ANOVA) table is computed in the same way as done by function 'anova.lm'.
#'
#'One drawback of using ANOVA-type estimation of random effects is, that random effects are independent, i.e they have
#'zero covariance by definition \eqn{cov(g_i,g_j) = 0}. Another one is that estimated variance components may become negative
#'under certain conditions. The latter situation is addressed by setting negative variance estimates equal to zero and adapting
#'ANOVA mean squares (MS) as \eqn{MS = C * VC}, where \eqn{C} is a coefficient matrix and a function of the design matrix \eqn{[X Z]}
#'and \eqn{VC} is the column-vector of adapted variance components. The Satterthwaite approximation of total degrees of freedom 
#'(DF for total variance) will use adapted \eqn{MS}-values.
#'
#'Note, that setting negative VCs equal to zero results in a conservative estimate of the total variance, i.e. it will be larger than
#'the estimate including the negative VC(s). Use parameter 'NegVC=TRUE' to explicitly allow negative variance estimates. 
#'
#'For further details on ANOVA Type-I estimation methods see \code{\link{anovaVCA}}.
#'
#'@param form			(formula) specifying the linear mixed model (fixed and random part of the model),
#'						all random terms need to be enclosed by round brackets. Any variable not being bracketed
#'						will be considered as fixed. Interaction terms containing at least one random factor
#'						will automatically be random (Piepho et al. 2003). All terms appearing in the model 
#'						(fixed or random) need to be compliant with the regular expression "^[^[\\.]]?[[:alnum:]_\\.]*$",
#'						i.e. they may not start with a dot and may then only consist of alpha-numeric characters, 
#'						dot and underscore. Otherwise, an error will be issued.
#'@param Data			(data.frame) containing all variables referenced in 'form', note that variables can only be
#'						of type "numeric", "factor" or "character". The latter will be automatically converted to "factor".
#'@param by				(factor, character) variable specifying groups for which the analysis should be performed individually,
#'						i.e. by-processing
#'@param VarVC.method	(character) string specifying whether to use the algorithm given in Searle et al. (1992) which corresponds to \code{VarVC.method="scm"} or in
#'						Giesbrecht and Burns (1985) which can be specified via "gb". Method "scm" (Searle, Casella, McCulloch)
#'						is the exact algorithm, "gb" (Giesbrecht, Burns) is termed "rough approximation"
#'						by the authors, but sufficiently exact compared to e.g. SAS PROC MIXED (method=type1) which
#'						uses the inverse of the Fisher-Information matrix as approximation. For balanced designs all
#'						methods give identical results, only in unbalanced designs differences occur. 
#'@param NegVC			(logical) FALSE = negative variance component estimates (VC) will be set to 0 and they will not 
#'						contribute to the total variance (as done e.g. in SAS PROC NESTED, conservative estimate of total variance). 
#'						The original ANOVA estimates can be found in element 'VCoriginal'. 
#'						The degrees of freedom of the total variance are based on adapted mean squares (MS) (see details).
#'						TRUE = negative variance component estimates will not be set to 0 and they will contribute to the total 
#'						variance (original definition of the total variance).
#'@param quiet			(logical) TRUE = will suppress any warning, which will be issued otherwise 
#'@param order.data		(logical) TRUE = class-variables will be ordered increasingly, FALSE = ordering of class-variables
#'						will remain as is
#'@return (VCA) object
#'
#'@author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#'
#'@references	
#'
#'Searle, S.R, Casella, G., McCulloch, C.E. (1992), Variance Components, Wiley New York	
#'
#'Goodnight, J.H. (1979), A Tutorial on the SWEEP Operator, The American Statistician, 33:3, 149-158
#'
#'Giesbrecht, F.G. and Burns, J.C. (1985), Two-Stage Analysis Based on a Mixed Model: Large-Sample
#'Asymptotic Theory and Small-Sample Simulation Results, Biometrics 41, p. 477-486
#'
#'H.P.Piepho, A.Buechse and K.Emrich (2003), A Hitchhiker's Guide to Mixed Models for Randomized Experiments,
#'J.Agronomy & Crop Science 189, p. 310-322
#'
#'Gaylor,D.W., Lucas,H.L., Anderson,R.L. (1970), Calculation of Expected Mean Squares by the Abbreviated Doolittle and Square Root Methods., 
#'Biometrics 26 (4): 641-655 
#'
#'SAS Help and Documentation PROC MIXED, SAS Institute Inc., Cary, NC, USA
#'
#'@seealso \code{\link{anovaVCA}}
#'
#'@examples 
#'\dontrun{
#'
#'data(dataEP05A2_2)
#'
#'# assuming 'day' as fixed, 'run' as random
#'anovaMM(y~day/(run), dataEP05A2_2)
#'
#'# assuming both as random leads to same results as
#'# calling anovaVCA
#'anovaMM(y~(day)/(run), dataEP05A2_2)
#'anovaVCA(y~day/run, dataEP05A2_2)
#'
#'# use different approaches to estimating the covariance of 
#'# variance components (covariance parameters)
#'dat.ub <- dataEP05A2_2[-c(11,12,23,32,40,41,42),]			# get unbalanced data
#'m1.ub <- anovaMM(y~day/(run), dat.ub, VarVC.method="scm")
#'m2.ub <- anovaMM(y~day/(run), dat.ub, VarVC.method="gb")		
#'V1.ub <- round(vcovVC(m1.ub), 12)
#'V2.ub <- round(vcovVC(m2.ub), 12)
#'all(V1.ub == V2.ub)
#'
#'# fit a larger random model
#'data(VCAdata1)
#'fitMM1 <- anovaMM(y~((lot)+(device))/(day)/(run), VCAdata1[VCAdata1$sample==1,])
#'fitMM1
#'# now use function tailored for random models
#'fitRM1 <- anovaVCA(y~(lot+device)/day/run, VCAdata1[VCAdata1$sample==1,])
#'fitRM1
#'
#'# there are only 3 lots, take 'lot' as fixed 
#'fitMM2 <- anovaMM(y~(lot+(device))/(day)/(run), VCAdata1[VCAdata1$sample==2,])
#'
#'# the following model definition is equivalent to the one above,
#'# since a single random term in an interaction makes the interaction
#'# random (see the 3rd reference for details on this topic)
#'fitMM3 <- anovaMM(y~(lot+(device))/day/run, VCAdata1[VCAdata1$sample==2,])
#'
#'# fit same model for each sample using by-processing
#'lst <- anovaMM(y~(lot+(device))/day/run, VCAdata1, by="sample")
#'lst
#'
#'# fit mixed model originally from 'nlme' package
#'
#'library(nlme)
#'data(Orthodont)
#'fit.lme <- lme(distance~Sex*I(age-11), random=~I(age-11)|Subject, Orthodont) 
#'
#'# re-organize data for using 'anovaMM'
#'Ortho <- Orthodont
#'Ortho$age2 <- Ortho$age - 11
#'Ortho$Subject <- factor(as.character(Ortho$Subject))
#'fit.anovaMM1 <- anovaMM(distance~Sex*age2+(Subject)*age2, Ortho)
#'
#'# use simplified formula avoiding unnecessary terms
#'fit.anovaMM2 <- anovaMM(distance~Sex+Sex:age2+(Subject)+(Subject):age2, Ortho)
#'
#'# and exclude intercept
#'fit.anovaMM3 <- anovaMM(distance~Sex+Sex:age2+(Subject)+(Subject):age2-1, Ortho)
#'
#'# compare results
#'fit.lme
#'fit.anovaMM1
#'fit.anovaMM2
#'fit.anovaMM3
#'
#'# are there a sex-specific differences?
#'cmat <- getL(fit.anovaMM3, c("SexMale-SexFemale", "SexMale:age2-SexFemale:age2")) 
#'cmat
#'
#'test.fixef(fit.anovaMM3, L=cmat)
#'
#'# former versions of the package used R-function 'lm' and 'anova',
#'# which is significantly slower for sufficiently large/complex models
#'data(realData)
#'datP1 <- realData[realData$PID==1,]
#'system.time(anova.lm.Tab <- anova(lm(y~lot/calibration/day/run, datP1)))
#'# Using the sweeping approach for estimating ANOVA Type-1 sums of squares
#'# this is now the default setting. 
#'system.time(anovaMM.Tab1  <- anovaMM(y~lot/calibration/day/run, datP1))
#'
#'# compare results, note that the latter corresponds to a linear model,
#'# i.e. without random effects. Various matrices have already been computed,
#'# e.g. "R", "V" (which are identical in this case).
#'anova.lm.Tab
#'anovaMM.Tab1
#'}
#'
#'@seealso \code{\link{anovaVCA}}, \code{\link{VCAinference}}, \code{\link{remlVCA}}, \code{\link{remlMM}}
#'\code{\link{ranef}}, \code{\link{fixef}}, \code{\link{vcov}}, \code{\link{vcovVC}}, 
#'\code{\link{test.fixef}}, \code{\link{test.lsmeans}}, \code{\link{plotRandVar}}

anovaMM <- function(form, Data, by=NULL, VarVC.method=c( "scm","gb"),
					NegVC=FALSE, quiet=FALSE, order.data=TRUE)
{
	if(!is.null(by))
	{
		stopifnot(is.character(by))
		stopifnot(by %in% colnames(Data))
		stopifnot(is.factor(by) || is.character(by))
		
		levels  <- unique(Data[,by])
		res <- lapply(levels, function(x) {
					tmp.res <- try(anovaMM(form=form, Data[Data[,by] == x,], NegVC=NegVC, VarVC.method=VarVC.method, quiet=quiet), silent=TRUE)
					if(is(tmp.res, "try-error") && !quiet)
						warning(paste0("Error for '", by, ".", x, "':\n\t", attr(tmp.res, "condition")$message))
					tmp.res
		})
		names(res) <- paste(by, levels, sep=".")
		return(res)
	}
	
	stopifnot(class(form) == "formula")
	stopifnot(identical(class(Data),"data.frame"))
	stopifnot(nrow(Data) > 2)                                               	# at least 2 observations for estimating a variance
	
	if(is.null(.GlobalEnv$msgEnv))												# may removed after loading the package
		msgEnv <<- new.env(parent=emptyenv())
	
	VarVC.method <- match.arg(VarVC.method)
	
	res <- list()
	res$call <- match.call()
	if(is.null(rownames(Data)))													# assign rownames if missing to record results of re-ordering taking place before the analysis
		rownames(Data) <- 1:nrow(Data)
	tobj <- terms(form, simplify=TRUE, keep.order=TRUE)                    		# expand nested factors if necessary retain order of terms in the formula
	Data <- orderData(Data, tobj, quiet=quiet, order.data=order.data)
	res$data <- Data															# as provided 
	org.form <- form
	
	if(length(attr(tobj, "term.labels")) == 0)									# handle pure-error models with 'anovaVCA'
		return(anovaVCA(form, Data))
	
	int   <- res$intercept <- attr(tobj, "intercept") == 1						# has intercept
	form  <- formula(tobj)
	res$terms <- tobj
	
	if(any(!grepl("^[^[\\.]]?[[:alnum:]_\\.]*$", rownames(attr(tobj, "factors")))))
		stop("There are terms in the model formula where regular expression '^[^[\\.]]?[[:alnum:]_\\.]*$' does not fit!")
	
	if(!attr(tobj, "response"))
		stop("You need to include a response variable in the fixed effects formula!")
	resp <- as.character(form)[2]
	res$response <- resp
	
	stopifnot(resp %in% colnames(Data))
	stopifnot(is.numeric(Data[,resp]))
	
	rf <- gregexpr("\\([[:alnum:]_\\.]*\\)", as.character(org.form)[3])			# check for random effects
	if(rf[[1]][1] != -1)														# identify random variables
	{
		len <- attr(rf[[1]], "match.length")
		pos <- rf[[1]]
		tmp <- NULL
		for(i in 1:length(len))
		{
			tmp <- c(tmp, substr(as.character(org.form)[3], pos[i]+1, pos[i]+len[i]-2))				# remember random factors, exclude brackets
		}
		rf <- tmp
	}
	
	Ndata <- nrow(Data)	
	rmInd <- integer()	
	resp.NA <- is.na(Data[,resp])
	
	if(any(resp.NA))
	{    
		rmInd <- c(rmInd, which(resp.NA))
		if(!quiet)
			message("There are ", length(which(resp.NA))," missing values for the response variable (obs: ", paste(which(resp.NA), collapse=", "), ")!")
		res$resp.NA <- rmInd
	}    
	
	fac  	<- attr(tobj, "term.labels")
	if(length(fac) > 1)
	{
		rf.ind  <- which(apply(sapply(rf, function(x) regexpr(x, fac)), 1, function(x) any(x>0)))		
	} else {
		if(length(rf) > 0)
		{
			if(rf == fac)
				rf.ind <- 1
			else
				rf.ind <- numeric(0)
		}
	}
	vars    <- rownames(attr(tobj, "factors"))[-1]	                        # remove response
	Nvc     <- length(fac) + 1 
	if(length(rf.ind) > 0)													# at least one random term in 'form'
	{
		res$random <- fac[rf.ind]
		res$fixed  <- fac[-rf.ind]
	}
	else
	{
		res$random <- character(0)											# only fixed effects
		res$fixed  <- fac
	}
	res$VCnamesall <-c(fac,"error")
	res$VCnames <- c(res$random, "error")
	res$Nvc  	<- length(res$VCnames)											# error is one additional VC
	res$Type 	<- if(length(res$fixed) == 0)
				"Random Model"
			else
				"Mixed Model"	
	
	for(i in rev(vars))														# check Data for consistency
	{
		if(length(unique(Data[,i])) == Ndata && !is.numeric(Data[,i]))
			stop("Variable '", i,"': number of levels of each grouping factor must be < number of observations!")
		
		if( any(is.na(Data[,i])))
		{
			NAind <- which(is.na(Data[,i]))
			rmInd <- c(rmInd, NAind)
			if(!quiet)
				message("Variable '", i,"' has ",length(NAind)," missing values (obs: ", paste(NAind, collapse=", "), ")!" )
		}
	}
	attr(res$data, "analysis.order") <- rownames(Data)						# actually record the ordering of the data 
	
	rmInd <- unique(rmInd)
	
	if(length(rmInd) > 0)
		Data <- Data[-rmInd,]												
	
	Data <- na.omit(Data[,rownames(attr(tobj, "factors"))])					# get rid of incomplete observations
	Mean <- mean(Data[,resp], na.rm=TRUE)                                   # mean computed after removing incomplete observations
	Nobs <- N <- nrow(Data)
	y    <- matrix(Data[,resp], ncol=1)										# vector of observations
	
	allObsEqual <- FALSE
	if(all(Data[,resp] == Data[1,resp]))									# no variance detectable?
	{
		allObsEqual <- TRUE
		Data.org <- Data
		Data[,resp] <- Data[,resp] + rnorm(nrow(Data))
		
		if(!quiet)
			warning("All values of response variable ", paste0("'", resp, "'"), " are equal!")
	}
	
	tmp.res <- getSSQsweep(Data, tobj, res$random)
	
	
	Lmat    <- tmp.res$Lmat
	CVC     <- tmp.res$CVC
	aov.tab <- tmp.res$aov.tab												# basic ANOVA-table
	rownames(aov.tab) <- res$VCnamesall
	DF <- aov.tab[,"DF"]
	SS <- aov.tab[,"SS"]
	VC <- aov.tab[,"VC"]
	C <-CVC$C
	Ci<-CVC$Ci
#	colnames(aov.tab) <- c("DF", "SS", "MS", "VC")
	
	if(allObsEqual)
	{
		aov.tab[,c("SS", "MS")] <- 0
		Data <- Data.org
		SS[1:length(SS)] <- 0
	}
	rownames(aov.tab) <- c(attr(tobj, "term.labels"), "error")
	
	res$Mean <- Mean
	res$formula <- org.form
	res$Nobs <- Nobs
	res$aov.org <- aov.tab
	
	rf.ind  <- c(rf.ind, nrow(aov.tab))
	
	
	# at this point Zre comprises fixed and random effects
	C2  <- apply(C, 2, function(x) x/DF)                                    # coefficient matrix for mean squares (MS)
	Ci2 <- solve(C2)
	VCorg <- VC
	VC  <- VCorg <- as.matrix( Ci %*% SS)                                   # solve for VCs (p.173)
	VCnam <- rownames(aov.tab)
	VCnam[length(VCnam)] <- "error"
	rownames(VC) <- rownames(VCorg) <- VCnam
	colnames(VC) <- colnames(VCorg) <- "Estimate"
	VCvar <- CVC$VCvar
#    VCvar <- t(C)%*%VCvar%*%C
#    VCvar <- Ci[rf.ind,]%*%VCvar%*%t(Ci[rf.ind,])
	VCvar <- VCvar[rf.ind,rf.ind]
	
	VC  <- VC[rf.ind] 														# only use VC for random terms, inter-fixed effects variation not in scope
	aov.tab <- aov.tab[rf.ind,,drop=F]
	
	rownames(aov.tab)[nrow(aov.tab)] <- "error"
	aov.tab$VC <- VC
	if (nrow(aov.tab)>1) rownames(VCvar)  <- colnames(VCvar) <- rownames(aov.tab)
	
	res$NegVCmsg <- ""
	res$VCoriginal <- aov.tab[, "VC"]
	
	if(!NegVC)
	{
		IndNegVC <- which(aov.tab[,"VC"] < 0)  
		if(length(IndNegVC) > 0)                                            # there are negative VC
		{
			aov.tab[IndNegVC, "VC"] <- 0                                    # set negative VC to 0
			res$NegVCmsg <- "* VC set to 0"
		}
	} 
	
	totVC   <- sum(aov.tab$VC)
	aov.tab <- rbind(total=c(NA, NA, NA, totVC, sqrt(totVC)), aov.tab) 	
	aov.tab["total", "DF"] <- SattDF(c(C2[rf.ind, rf.ind] %*% aov.tab[-1, "VC"]), 	# will automatically adapt ANOVA-MS if any VCs were set to 0 
			Ci=Ci2[rf.ind, rf.ind, drop=F], DF=DF[rf.ind])  
	
#	suppressWarnings(aov.tab <- cbind(aov.tab, SD=sqrt(aov.tab[,"VC"])))    		# warnings suppressed because sqrt of negative numbers doese not exists
	aov.tab <- cbind(aov.tab, "CV[%]"=aov.tab[,"SD"]*100/Mean)
	aov.tab <- cbind(aov.tab, "%Total"=aov.tab[,"VC"]*100/totVC)
	aov.tab <- aov.tab[,c("DF", "SS", "MS", "VC", "%Total", "SD", "CV[%]")]    
	aov.tab <- as.matrix(aov.tab)
	aov.tab <- apply(aov.tab, 1:2, function(x) ifelse(is.nan(x), NA, x))
	
	res$EstMethod <- "ANOVA"
	
	if(Nobs != Ndata)
		res$Nrm <- Ndata - Nobs                            				# save number of observations that were removed due to missing data
	
	Lmat$y <- y															# column vector of observations
	if(int)
	{
		Lmat$X    <- matrix(1, ncol=1, nrow=nrow(Data))					# design matrix of fixed effects: include intercept --> needs a restriction
		colnames(Lmat$X) <- "int"
		fe.assign <- 0													# '0' indicates intercept
	}
	else
		fe.assign <- NULL
	
	fe <- fac[-rf.ind]
	
	INT <- ifelse(int, "1", "-1")
	if(length(fe) > 0)
		INT <- paste(INT, "+", sep="")
	
	fixed.form <- paste(resp, "~", INT, paste(fe, collapse="+"), sep="")
	fixed.form <- as.formula(fixed.form)
	
	old.opt <- options(contrasts=c("contr.SAS", "contr.poly"))
	
	suppressWarnings({
				X1 	<- model.matrix(fixed.form, data=Data)						# with SAS contrasts but without the columns where restrictions apply
				X10	<- apply(X1, 2, function(x) all(x==0))
				X2 <- model.matrix(	fixed.form, data=Data,						# full model matrix with re-setted contrasts
						contrasts.arg=lapply(Data[,sapply(Data, is.factor), drop=FALSE],
								contrasts, contrasts=FALSE))
				X20	<- apply(X2, 2, function(x) all(x==0))
				
				X2.asgn 	<- attr(X2, "assign")									# keep info		
				
				if(any(X10))
					X1 <- X1[,-which(X10),drop=FALSE]
				
				if(any(X20))
				{
					X2 <- X2[,-which(X20),drop=FALSE]
					X2.asgn 	<- X2.asgn[-which(X20)]
				}				
			})
	
	options(old.opt)													# reset contrasts option
	
	fixed.terms <- terms(fixed.form)
	fe.terms 	<- attr(fixed.terms, "term.labels")
	
	X2[,!colnames(X2) %in% colnames(X1)] <- 0							# transfer restriction from X1 to X2
	if(int)
	{
		colnames(X2)[1] <- "int"
		fe.terms <- c("int", fe.terms)
	}
	
	fe.assign <- X2.asgn
	attr(fe.assign, "terms") <- fe.terms
	
	Lmat$X <- X2
	
	res$fe.assign 	<- fe.assign										# mapping columns of X to fixed terms in the model formula
	res$fixed.terms <- fixed.terms
	
	Lmat$rf.ind <- rf.ind												# indices of random effects
	Lmat$VCall  <- VCorg												# VC-estimates as if all factors were random
	Lmat$C.SS   <- C
	Lmat$C.MS   <- C2 
	Lmat$Ci.SS  <- Ci
	Lmat$Ci.MS  <- Ci2
	Lmat$VCvar   <- VCvar
	
	res$VarVC.method <- VarVC.method
	res$aov.tab  <- aov.tab
	res$Matrices <- Lmat
	res$balanced <- if(isBalanced(form, Data)) 
				"balanced"  
			else 
				"unbalanced"
	
	class(res)     <- "VCA"
	
	if(length(Lmat$rf.ind) == 1)										# contains only the residual error
	{
		res$Type <- "Linear Model"
		tmp      <- res$aov.org
		tmp		 <- tmp[,-which(colnames(tmp) %in% c("VC", "SD"))]
		tmp 	 <- cbind(tmp, "F value" = tmp[,"MS"]/tmp[nrow(tmp),"MS"])
		tmp		 <- cbind(tmp, "Pr(>F)" = pf(tmp[,"F value"]/tmp[nrow(tmp), "F value"], tmp[, "DF"], tmp[nrow(tmp), "DF"], lower.tail=FALSE))
		tmp[nrow(tmp), c("F value", "Pr(>F)")] <- NA	
		res$aov.org <- tmp
	}
	res <- solveMME(res)
	
	if(allObsEqual)
		res$aov.tab[,"%Total"] <- 0	
	
	gc(verbose=FALSE)													# trigger garbage collection
	return(res)
}




#'Bottom-Up Step-Wise VCA-Analysis of the Complete Dataset
#'
#'Function performs step-wise VCA-analysis on a fitted VCA-object by leaving out N-1 to 0
#'top-level variance components (VC).
#'
#'This function uses the complete data to quantify sub-sets of variance components.
#'In each step the current total variance is estimated by subtracting the sum of all left-out VCs
#'from the total variance of the initial VCA object. Doing this guarantees that the contribution to the total 
#'variance which is due to left-out VCs is accounted for, i.e. it is estimated but not included/reported.
#'The degrees of freedom (DFs) of the emerging total variances of sub-sets are determined using the Satterthwaite
#'approximation. This is achieved by extracting the corresponding sub-matrix from the coefficient matrix \eqn{C} of
#'the 'VCA' object, the sub-vector of ANOVA mean squares, and the sub-vector of degrees of freedom and calling
#'function \code{\link{SattDF}} method="total".
#'
#'This step-wise procedure starts one-level above error (repeatability) and ends at the level of the upper-most VC.
#'It can only be used on models fitted by ANOVA Type-1, i.e. by function \code{\link{anovaVCA}}.
#'
#'@param obj			(VCA) object representing the complete analysis
#'@param VarVC.method	(character) string specifying the algorithm to be used for estimating variance-covariance matrix
#'						of VCs (see \code{\link{anovaMM}} for details).
#'
#'@return (list) of (simplified) 'VCA' objects representing analysis-result of sub-models
#'
#'@author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#'
#'@examples 
#'\dontrun{
#'data(VCAdata1)
#'datS7L1 <- VCAdata1[VCAdata1$sample == 7 & VCAdata1$lot == 1, ]
#'fit0 <- anovaVCA(y~device/day/run, datS7L1, MME=TRUE)
#'
#'# complete VCA-analysis result
#'fit0
#'
#'# perform step-wise (bottom-up) VCA-analyses
#'sw.res <- stepwiseVCA(fit0)
#'sw.res
#'
#'# get CIs on intermediate precision 
#'VCAinference(sw.res[["device:day"]])
#'}

stepwiseVCA <- function(obj, VarVC.method=c("scm", "gb"))
{
	stopifnot(obj$Type == "Random Model")					# only random models correspond to true VCA-analyses
	stopifnot(obj$EstMethod == "ANOVA")
	VarVC.method <- match.arg(VarVC.method)
	tab <- obj$aov.tab
	Nstp  <- nrow(tab)-2
	res <- vector("list", length=Nstp)
	names(res) <- rev(rownames(tab)[2:(Nstp+1)])
	Ci <- getMat(obj, "Ci.MS")
	N  <- nrow(tab)
	VarCov <- vcovVC(obj, method=VarVC.method)
	
	for(i in 1:Nstp)
	{
		C <- Ci[(N-i-1):(N-1), (N-i-1):(N-1)]
		M <- tab[(N-i):(N), "MS"]
		D <- tab[(N-i):(N), "DF"]
		
		DF <- SattDF(M, C, D)
		mat <- tab[(N-i):N, ]	
		mat <- rbind(total=c(DF=DF, SS=NA, MS=NA, VC=sum(mat[,"VC"]), "%Total"=100, SD=NA, "CV[%]"=NA), mat)
		mat[, "%Total"] <- 100 * mat[,"VC"] / mat["total", "VC"]
		mat["total", "SD"] <- sqrt(mat["total", "VC"])
		mat["total", "CV[%]"] <- 100 * mat["total", "SD"] / obj$Mean
		ele <- list(aov.tab=mat)
		ele$NegVCmsg <- obj$NegVCmsg
		class(ele) <- "VCA"
		ele$Type <- "Random Model"
		ele$EstMethod <- obj$EstMethod 
		ele$Mean <- obj$Mean
		ele$Nobs <- obj$Nobs
		ele$balanced <- obj$balanced
		ele$Matrices <- list(Ci.MS=C, rf.ind=1:(nrow(mat)-1))

		ele$VarCov <- VarCov[(N-i-1):(N-1), (N-i-1):(N-1)]
		ele$VarVC.method <- VarVC.method
		
		res[[i]] <- ele
	}
	return(res)
}


#'Calling F90-implementation of the SWEEP-Operator
#'
#'Function calls a fast Fortran90-implementation of the SWEEP operator using the
#'transpose of the original augmented matrix \eqn{X'X} (see \code{\link{getSSQsweep}}).
#'In the sweeping step, also the C matrix, needed to obtain the variance estimates from 
#'the sum of squares and the Covariance matrix of the estimates are calculated. 
#'
#'This is an utility-function not intended to be called directly.
#'
#'@param M			(matrix) matrix, representing the augmented matrix \eqn{X'X}
#'@param asgn		(integer) vector, identifying columns in \eqn{M} corresponding to variables, 
#'					respectively, to their coefficients
#'@param thresh		(numeric) value used to check whether the influence of the a coefficient
#'					to reducing the error sum of squares is small enough to conclude that the
#'					corresponding column in \eqn{X'X} is a linear combination of preceding 
#'					columns
#'@param tol		(numeric) value used to check numerical equivalence to zero
#'@param Ncpu		(integer) number of cores to be used for parallel processing
#'(not yet used)
#'
#'@author Florian Dufey \email{florian.dufey@@roche.com}
#'
#'@return (list) with eight elements:\cr
#'\itemize{
#'\item{SSQ}{(numeric) vector of ANOVA sum of squares}
#'\item{LC}{(integer) vector indicating linear dependence of each column}
#'\item{DF}{(integer) degrees of freedom}
#'\item{C}{(double precision) Matrix relating the sums of squares to the variances}
#'\item{Ci}{(double precision) inverse of matrix relating the sums of squares to the variances}
#'\item{VC}{(double precision)  variance}
#'\item{SD}{(double precision) standard deviations}
#'\item{Var}{(double precision) covariance matrix of the estimated variances}
#'}
#'
#'@references 
#'Goodnight, J.H. (1979), A Tutorial on the SWEEP Operator, The American Statistician, 33:3, 149-158
#'

Fsweep <- function(M, asgn,  thresh=1e-10, tol=1e-10, Ncpu=1)
{
	#dyn.load('vca.dll')
	nr <- nrow(M)
	nSSQ <- length(unique(asgn))
	nSSQ2 <- nSSQ*nSSQ
	stopifnot(nr == ncol(M))
	LC <- ZeroK <- rep(0, length(asgn))
	
	VC <- rep(0.0,nSSQ)
	SD <- rep(0.0,nSSQ)
	C <- matrix(rep(0,nSSQ2),ncol=nSSQ,nrow=nSSQ)
	Ci <- matrix(rep(0,nSSQ2),ncol=nSSQ,nrow=nSSQ)
	Var <- matrix(rep(0,nSSQ2),ncol=nSSQ,nrow=nSSQ)
	ind <- is.nan(M) | is.na(M) | M == Inf | M == -Inf;
	if(any(ind)) M[which(ind)] <- 0
	info <- 0
	
	
	swept <- .Fortran("Gsweep", M=as.double(M), NumK=as.integer(length(asgn)), k=as.integer(asgn), thresh=as.double(thresh), 
			nr=as.integer(nr), LC=as.integer(LC), 
			tol=as.double(tol),nSSQ=as.integer(nSSQ), SSQ=as.double(rep(0.0, nSSQ)),DF=as.integer(rep(0,nSSQ)),VC=as.double(VC), SD=as.double(SD),C=as.double(C),Ci=as.double(Ci),Var=as.double(Var),info=as.integer(info),PACKAGE="VCA")

	info <- swept$info
	
	if(info!=0)
		stop("Inverting the C Matrix failed! Double check the specified model! Are there replicates?")
		
	res <- list(SSQ=swept$SSQ,sweptLC = swept$LC, DF = swept$DF,
				LC=tapply(swept$LC, asgn, function(x) length(which(x==1))),
				VC=swept$VC, SD=swept$SD,C=matrix(swept$C,ncol=nSSQ,nrow=nSSQ),
				Ci=matrix(swept$Ci,ncol=nSSQ,nrow=nSSQ),
				Var=matrix(swept$Var,ncol=nSSQ,nrow=nSSQ))
	
	#dyn.unload('vca.dll')
	return(res)
}

#'ANOVA Sum of Squares via Sweeping
#'
#'Compute ANOVA Type-1 sum of squares for linear models.
#'
#'This function performs estimation of ANOVA Type-1 sum of squares
#'using the SWEEP-operator (see reference), operating on the augmented
#'matrix \eqn{X'X}, where \eqn{X} represents the design matrix not differentiating
#'between fixed and random factors. See the numerical example in \code{\link{Fsweep}}
#'exemplifying the type of augmentation of \eqn{X'X} on which sweeping is carried out.
#'
#'This is an utility function not intended to be called directly.
#'For each term in the formula the design-matrix \eqn{Z} is constructed.
#'Matrix \eqn{X} corresponds to binding all these \eqn{Z}-matrices together column-wise.
#'
#'Degrees of freedom for each term are determined by subtracting the number of
#'linearly dependent columns from the total number of column in X asigned to a
#'specific term.
#'
#'@param Data			(data.frame) with the data
#'@param tobj			(terms) object derived from original formula object
#'@param random		(character) vector, optionally containing information about each
#'model term, whether it is random or fixed (only used in mixed models)
#'
#'@return (list) representing the  with variables:\cr
#'\itemize{
#'\item{aov.tab}{basic ANOVA-table with degrees of freedom (DF), SS and MS}
#'\item{Lmat}{(list) with components 'Z' and 'A'}
#'}
#'
#'@author 	Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com},
#' 			Florian Dufey \email{florian.dufey@@roche.com}
#'
#'@seealso \code{\link{Fsweep}}
#'
#'@references 
#'
#'Goodnight, J.H., (1979), A Tutorial on the SWEEP Operator, The American Statistician, 33:3, p.149-158
#'
#'@examples 
#'\dontrun{
#'data(dataEP05A2_1)
#'res <- VCA:::getSSQsweep(dataEP05A2_1, terms(y~day/run))
#'str(res)
#'}

getSSQsweep <- function(Data, tobj, random=NULL)
{	
	form     <- formula(tobj)
	resp     <- as.character(form)[2]
	fac      <- attr(tobj, "term.labels")
	int      <- attr(tobj, "intercept") == 1 
	
	N        <- nrow(Data)															
	SS       <- numeric()
	DF       <- numeric()
	Lmat     <- list()                                                      	# compute A-matrices, used in constructing covariance-matrix of VCs
	Lmat$Z   <- list()
	Lmat$Zre <- Matrix(nrow=N, ncol=0)
	y        <- Matrix(Data[,resp], ncol=1)
	
	ord <- attr(tobj, "order")
	vars <- rownames(attr(tobj, "factors"))
	fac <- c(fac, "error")
	DF  <- rep(0, length(fac))
	done <- rep(FALSE, length(fac))					# keep track which DFs have been computed
	names(DF) <- fac
	vars <- vars[-1]
	Nvc  <- length(fac) 							# error not included
	N <- nrow(Data)
	asgn <- NULL
	Lmat <- list()
	Lmat$Z <- list()
	CVC <- list()
	
	if(int)
	{
		Lmat$Zre <- Matrix(1, nrow=N, ncol=1)
		colnames(Lmat$Zre) <- "int"
		asgn     <- 0
	}else{
		Lmat$Zre <- Matrix(0, nrow=N, ncol=1)
		colnames(Lmat$Zre) <- "dummy"
		asgn     <- -1
	}
	
	NumK <- c(rep(0, Nvc-1), N)
	
	for(i in 1:Nvc)                 				# construct design-matrices Z for each term, and matrix X (here Zre)                       
	{
		if(i < Nvc)
		{
			tmpMM <- model.matrix(as.formula(paste(resp, "~", fac[i], "-1", sep="")), Data)
			Lmat$Z[[i]] <- Matrix(tmpMM)
			all0 <- apply(Lmat$Z[[i]], 2, function(x) all(x==0))
			if(any(all0))
			{
				ZeroCol <- which(all0)
				Lmat$Z[[i]] <- Lmat$Z[[i]][,-ZeroCol]
			}	
			if(!is.null(random))
				attr(Lmat$Z[[i]], "type") <- ifelse(fac[i] %in% random, "random", "fixed")
			
			attr(Lmat$Z[[i]], "term") <- fac[i]
			
			NumK[i] <- ncol(Lmat$Z[[i]])
			Lmat$Zre    <- Matrix(cbind(as.matrix(Lmat$Zre), as.matrix(Lmat$Z[[i]])))
			asgn <- c(asgn, rep(i, ncol(Lmat$Z[[i]])))
		}
		else
		{
			Lmat$Z[[Nvc]] <- Diagonal(N)							# include error design matrix
			attr(Lmat$Z[[Nvc]], "type") <- "random"
			attr(Lmat$Z[[Nvc]], "term") <- "error"
		}
	}
	
	attr(Lmat$Zre, "assign") <- asgn
	X <- Lmat$Zre
	Xt <- t(X)
	y  <- Matrix(Data[,resp], ncol=1)
	yt <- t(y)
	
	M <- rbind(	cbind(as.matrix(Xt%*%X), as.matrix(Xt%*%y)), 
			cbind(as.matrix(yt%*%X), as.matrix(yt%*%y)))	
	
	uind <- unique(asgn)							# all factors
	SS <- LC <- NULL
	nr <- nrow(M)
	if (!int) M[1,1]=N
	V 		<- rep(1, ncol(M))
	swept 	<- Fsweep(M, asgn=asgn)					# sweep matrix M
	LC    	<- swept$LC
	SS	  	<- swept$SSQ
	SSQ   	<- swept$SSQ
	DF    	<- swept$DF
	C     	<- swept$C
	Ci    	<- swept$Ci
	VC    	<- swept$VC
	SD    	<- swept$SD
	VCvar 	<- swept$Var
	CVC$C 	<- C
	CVC$Ci 	<- Ci
	CVC$VCvar <- VCvar
	aov.tab <- data.frame(DF=DF, SS=SSQ, MS=SSQ/DF, VC=VC, SD=SD)	# basic ANOVA-table
	
	return(list(aov.tab=aov.tab, CVC=CVC,Lmat=Lmat))
}

