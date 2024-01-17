
#'Perform (V)ariance (C)omponent (A)nalysis via REML-Estimation
#'
#'Function performs a Variance Component Analysis (VCA) using Restricted Maximum Likelihood (REML)
#'to fit the random model, i.e. a linear mixed model (LMM) where the intercept is the only fixed effect.
#'
#'Here, a variance component model is fitted by REML using the \code{\link{lmer}} function of the
#'\code{lme4}-package. For all models the Giesbrechnt & Burns (1985) approximation of the variance-covariance
#'matrix of variance components (VC) is applied. A Satterthwaite approximation of the degrees of freedom
#'for all VC and total variance is based on this approximated matrix using \eqn{df=2Z^2}, where
#'\eqn{Z} is the Wald statistic \eqn{Z=\sigma^2/se(\sigma^2)}, and \eqn{\sigma^2} is here used for an
#'estimated variance. The variance of total variability, i.e. the sum of all VC is computed via summing
#'up all elements of the variance-covariance matrix of the VC.
#'Note, that for large datasets approximating the variance-covariance matrix of VC is computationally expensive
#'and may take very long. There is no Fisher-information matrix available for 'merMod' objects, which can
#'serve as approximation. To avoid this time-consuming step, use argument 'VarVC=FALSE' but remember,
#'that no confidence intervals for any VC will be available. If you use Microsoft's R Open, formerly known
#'as Revolution-R, which comes with Intel's Math Kernel Library (MKL), this will be automatically detected
#'and an environment-optimized version will be used, reducing the computational time very much (see examples).
#'
#'@param form          (formula) specifying the model to be fit, a response variable left of the '~' is mandatory
#'@param Data          (data.frame) containing all variables referenced in 'form'
#'@param by			(factor, character) variable specifying groups for which the analysis should be performed individually,
#'i.e. by-processing
#'@param VarVC			(logical) TRUE = the variance-covariance matrix of variance components will be approximated using 
#'the method found in Giesbrecht & Burns (1985), which also serves as basis for applying a Satterthwaite
#'approximation of the degrees of freedom for each variance component, FALSE = leaves out this step, 
#'no confidence intervals for VC will be available
#'@param quiet			(logical) TRUE = will suppress any messages or warnings, which will be issued otherwise
#'@param order.data	(logical) TRUE = class-variables will be ordered increasingly, FALSE = ordering of class-variables
#'will remain as is 
#'
#'@seealso \code{\link{remlMM}}, \code{\link{VCAinference}}, \code{\link{ranef.VCA}}, \code{\link{residuals.VCA}},
#'\code{\link{anovaVCA}}, \code{\link{anovaMM}}, \code{\link{plotRandVar}}, \code{\link{lmer}}
#'
#'@author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#'
#'@examples
#'\dontrun{
#'
#'# a VCA standard example
#'data(dataEP05A2_3)
#'
#'# fit it by ANOVA first, then by REML
#'fit0 <- anovaVCA(y~day/run, dataEP05A2_3) 
#'fit1 <- remlVCA(y~day/run, dataEP05A2_3)
#'fit0
#'fit1
#'
#'# make example unbalanced
#'set.seed(107)
#'dat.ub <- dataEP05A2_3[-sample(1:80, 7),]
#'fit0ub <- anovaVCA(y~day/run, dat.ub) 
#'fit1ub <- remlVCA(y~day/run, dat.ub) 
#'
#'# not that ANOVA- and REML-results now differ
#'fit0ub
#'fit1ub
#'
#'### Use the six sample reproducibility data from CLSI EP5-A3
#'### and fit per sample reproducibility model
#'data(CA19_9)
#'fit.all <- remlVCA(result~site/day, CA19_9, by="sample")
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
#'# REML-estimation not yes optimzed to the same degree as
#'# ANOVA-estimation. Note, that no variance-covariance matrix
#'# for the REML-fit is computed (VarVC=FALSE)!
#'# Note: A correct analysis would be done per-sample, this is just
#'#       for illustration.
#'data(VCAdata1)
#'# with complete sweeping implemented as FORTRAN-routine fit 
#'system.time(fit0 <- anovaVCA(y~sample+(device+lot)/day/run, VCAdata1))
#'system.time(fit1 <- remlVCA(y~sample+(device+lot)/day/run, VCAdata1, VarVC=FALSE))
#'
#'# The previous example will also be interesting for environments using MKL.
#'# Run it once in a GNU-R environment and once in a MKL-environment
#'# and compare computational time of both. Note, that 'VarVC' is now set to TRUE
#'# and variable "sample" is put into the brackets increasing the number of random
#'# effects by factor 10. On my Intel Xeon E5-2687W 3.1 GHz workstation it takes
#'# ~ 400s with GNU-R and ~25s with MKL support (MRO) both run under Windows.
#'system.time(fit2 <- remlVCA(y~(sample+device+lot)/day/run, VCAdata1, VarVC=TRUE))
#'
#'# using the SWEEP-Operator is even faster 
#'system.time(fit3 <- anovaVCA(y~(sample+device+lot)/day/run, VCAdata1))
#'fit2
#'fit3
#'}

remlVCA <- function(form, Data, by=NULL, VarVC=TRUE, quiet=FALSE, order.data=TRUE)
{
	Call <- match.call()
	
	
	
	if(!is.null(by))
	{
		stopifnot(is.character(by))
		stopifnot(by %in% colnames(Data))
		stopifnot(is.factor(by) || is.character(by))
		
		levels  <- unique(Data[,by])
		res <- lapply(levels, function(x) {
					tmp.res <- try(remlVCA(form=form, Data[Data[,by] == x,], VarVC=VarVC, quiet=quiet), silent=TRUE)
					if(is(tmp.res, "try-error") && !quiet)
						warning(paste0("Error for '", by, ".", x, "':\n\t", attr(tmp.res, "condition")$message))
					tmp.res
				})
		names(res) <- paste(by, levels, sep=".")
		return(res)
	}
	
	stopifnot(class(form) == "formula")
	stopifnot(is.logical(VarVC))
	stopifnot(is.logical(quiet))
	stopifnot(identical(class(Data),"data.frame"))
	stopifnot(nrow(Data) > 2)                                               # at least 2 observations for estimating a variance
	
	if(is.null(.GlobalEnv$msgEnv))											# may removed after loading the package
		msgEnv <<- new.env(parent=emptyenv())
	
	org.form <- form
	trms <- terms(form)														# convert VCA-formula to valid lmer-formula
	stopifnot(attr(trms, "response") == 1)
	lab  <- attr(trms, "term.labels")
	Data <- orderData(Data, trms, quiet=quiet,
			exclude.numeric=FALSE, order.data=order.data)		# convert all variables to factors
	
	allObsEqual <- FALSE
	
	if(length(lab) == 0)
	{
		allObsEqual <- TRUE
		if(!quiet)
			warning("No random effects specified! Call function 'anovaVCA' instead!")
		return(anovaVCA(form, Data))
	}
	
	lab  <- paste("(1|", lab, ")", sep="")
	resp <- rownames(attr(trms, "factors"))[1]
	vars <- rownames(attr(trms, "factors"))[-1]
	
	allObsEqual <- FALSE
	
	if(all(Data[,resp] == Data[1,resp]))										# no variance detectable?
	{
		allObsEqual <- TRUE
		
		if(!quiet)
			warning("All values of response variable ", paste0("'", resp, "'"), " are equal!")
	}
	
	form <- as.formula(paste(resp, "~", paste(lab, collapse="+"), sep=""))
	
	stopifnot(resp %in% colnames(Data))
	stopifnot(is.numeric(Data[,resp]))
	
	rmInd 	<- integer()	
	resp.NA <- is.na(Data[,resp])
	
	if(any(resp.NA))
	{    
		rmInd <- c(rmInd, which(resp.NA))
		if(!quiet)
			message("There are ", length(which(resp.NA))," missing values for the response variable (obs: ", paste(which(resp.NA), collapse=", "), ")!")
	}    
	
	if(!is.null(vars))
	{
		for(i in vars)                                                  # convert all nested factors as factor-objects
		{
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
	
	vcol <- rownames(attr(trms, "factors"))
	if(is.null(vcol))
		vcol <- resp
	Ndata <- nrow(Data)
	Data  <- na.omit(Data[,vcol, drop=F])
	Nobs <- nrow(Data)
	
	if(quiet || allObsEqual)
		suppressWarnings(fit <- lmer(form, Data, control=lmerControl(optimizer="bobyqa")))						# fit via 'lmer'
	else
		fit <- lmer(form, Data, control=lmerControl(optimizer="bobyqa"))
	
	res <- list(call=Call, Type="Random Model", EstMethod="REML", data=Data, terms=trms,
			intercept=as.logical(attr(trms, "intercept")), response=resp)
	
	res$Mean 	 	 <- mean(Data[,resp], na.rm=TRUE)
	res$formula		 <- org.form												# as specified by the user
	res$form 	 	 <- form
	res$Nvc  	 	 <- length(lab) + 1											# error is one additional VC
	res$VCnames  	 <- c(attr(trms, "term.labels", "error"))
	res$NegVCmsg 	 <- ""														# there will never be anything to report
	res$VarVC.method <- "gb"
	res$balanced <- if(isBalanced(as.formula(trms), Data)) 
				"balanced"  
			else 
				"unbalanced" 
	
	res$terms.classes <- sapply(Data[,rownames(attr(trms, "factors"))[-1]], class)
	
	if(Nobs != Ndata)
		res$Nrm <- Ndata - Nobs                         # save number of observations that were removed due to missing data
	
	res$Nobs <- Nobs
	
	tmp <- lmerSummary(	obj=fit, VarVC=VarVC, 			# construct table similar to aov-table and approximate vcovVC
			terms=attr(trms, "term.labels"),
			Mean=res$Mean)	
	res <- c(res, tmp)
	class(res) <- "VCA"
	
	if(allObsEqual)
	{
		res$aov.tab[,"DF"] <- NA
		res$aov.tab[,"%Total"] <- 0	
	}
	res
}



#'Fit Linear Mixed Models via REML
#'
#'Function fits Linear Mixed Models (LMM) using Restricted Maximum Likelihood (REML).
#'
#'The model is formulated exactly as in function \code{\link{anovaMM}}, i.e. random terms need be enclosed by round brackets.
#'All terms appearing in the model (fixed or random) need to be compliant with the regular expression "^[^[\\.]]?[[:alnum:]_\\.]*$",
#'i.e. they may not start with a dot and may then only consist of alpha-numeric characters, 
#'dot and underscore. Otherwise, an error will be issued.
#'
#'Here, a LMM is fitted by REML using the \code{\link{lmer}} function of the \code{lme4}-package. 
#'For all models the Giesbrechnt & Burns (1985) approximation of the variance-covariance
#'matrix of variance components (VC) can be applied ('VarVC=TRUE'). A Satterthwaite approximation of the degrees of freedom
#'for all VC and total variance is based on this approximated matrix using \eqn{df=2Z^2}, where
#'\eqn{Z} is the Wald statistic \eqn{Z=\sigma^2/se(\sigma^2)}, and \eqn{\sigma^2} is here used for an
#'estimated variance. The variance of total variability, i.e. the sum of all VC is computed via summing
#'up all elements of the variance-covariance matrix of the VC.
#'One can constrain the variance-covariance matrix of random effects \eqn{G} to be either diagonal ('cov=FALSE'), i.e.
#'all random effects are indpendent of each other (covariance is 0). If 'cov=TRUE' (the default) matrix \eqn{G} will be
#'constructed as implied by the model returned by function \code{\link{lmer}}. 
#'
#'As for objects returned by function \code{\link{anovaMM}} linear hypotheses of fixed effects or LS Means can be
#'tested with functions \code{\link{test.fixef}} and \code{\link{test.lsmeans}}. Note, that option "contain" does
#'not work for LMM fitted via REML.
#'
#'Note, that for large datasets approximating the variance-covariance matrix of VC is computationally expensive
#'and may take very long. There is no Fisher-information matrix available for 'merMod' objects, which can
#'serve as approximation. To avoid this time-consuming step, use argument 'VarVC=FALSE' but remember,
#'that no confidence intervals for any VC will be available. If you use Microsoft's R Open, formerly known
#'as Revolution-R, which comes with Intel's Math Kernel Library (MKL), this will be automatically detected
#'and an environment-optimized version will be used, reducing the computational time considerably (see examples).
#'
#'@param form          (formula) specifying the model to be fit, a response variable left of the '~' is mandatory, random terms
#'have to be enclosed in brackets (see details for definition of valid model terms)
#'@param Data          (data.frame) containing all variables referenced in 'form'
#'@param by			(factor, character) variable specifying groups for which the analysis should be performed individually,
#'i.e. by-processing
#'@param VarVC			(logical) TRUE = the variance-covariance matrix of variance components will be approximated using 
#'the method found in Giesbrecht & Burns (1985), which also serves as basis for applying a Satterthwaite
#'approximation of the degrees of freedom for each variance component, FALSE = leaves out this step, 
#'no confidence intervals for VC will be available
#'@param cov			(logical) TRUE = in case of non-zero covariances a block diagonal matrix will be constructed,
#'FALSE = a diagonal matrix with all off-diagonal element being equal to zero will be contructed
#'@param quiet			(logical) TRUE = will suppress any messages or warning, which will be issued otherwise 
#'@param order.data	(logical) TRUE = class-variables will be ordered increasingly, FALSE = ordering of class-variables
#'will remain as is
#'
#'@seealso \code{\link{remlVCA}}, \code{\link{VCAinference}}, \code{\link{ranef.VCA}}, \code{\link{residuals.VCA}},
#'\code{\link{anovaVCA}}, \code{\link{anovaMM}}, \code{\link{plotRandVar}},  \code{\link{test.fixef}},  
#'\code{\link{test.lsmeans}}, \code{\link{lmer}}
#'
#'@author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#'
#'@examples
#'\dontrun{
#'data(dataEP05A2_2)
#'
#'# assuming 'day' as fixed, 'run' as random
#'remlMM(y~day/(run), dataEP05A2_2)
#'
#'# assuming both as random leads to same results as
#'# calling anovaVCA
#'remlMM(y~(day)/(run), dataEP05A2_2)
#'anovaVCA(y~day/run, dataEP05A2_2)
#'remlVCA(y~day/run, dataEP05A2_2)
#'
#'# fit a larger random model
#'data(VCAdata1)
#'fitMM1 <- remlMM(y~((lot)+(device))/(day)/(run), VCAdata1[VCAdata1$sample==1,])
#'fitMM1
#'# now use function tailored for random models
#'fitRM1 <- anovaVCA(y~(lot+device)/day/run, VCAdata1[VCAdata1$sample==1,])
#'fitRM1
#'
#'# there are only 3 lots, take 'lot' as fixed 
#'fitMM2 <- remlMM(y~(lot+(device))/(day)/(run), VCAdata1[VCAdata1$sample==2,])
#'
#'# the following model definition is equivalent to the one above,
#'# since a single random term in an interaction makes the interaction
#'# random (see the 3rd reference for details on this topic)
#'fitMM3 <- remlMM(y~(lot+(device))/day/run, VCAdata1[VCAdata1$sample==2,])
#'
#'# fit same model for each sample using by-processing
#'lst <- remlMM(y~(lot+(device))/day/run, VCAdata1, by="sample")
#'lst
#'
#'# fit mixed model originally from 'nlme' package
#'
#'library(nlme)
#'data(Orthodont)
#'fit.lme <- lme(distance~Sex*I(age-11), random=~I(age-11)|Subject, Orthodont) 
#'
#'# re-organize data for using 'remlMM'
#'Ortho <- Orthodont
#'Ortho$age2 <- Ortho$age - 11
#'Ortho$Subject <- factor(as.character(Ortho$Subject))
#'fit.remlMM1 <- remlMM(distance~Sex*age2+(Subject)*age2, Ortho)
#'
#'# use simplified formula avoiding unnecessary terms
#'fit.remlMM2 <- remlMM(distance~Sex+age2+Sex:age2+(Subject)+age2:(Subject), Ortho)
#'
#'# and exclude intercept
#'fit.remlMM3 <- remlMM(distance~Sex+Sex:age2+(Subject)+(Subject):age2-1, Ortho)
#'
#'# now use exclude covariance of per-subject intercept and slope
#'# as for models fitted by function 'anovaMM'
#'fit.remlMM4 <- remlMM(distance~Sex+Sex:age2+(Subject)+(Subject):age2-1, Ortho, cov=FALSE)
#'
#'# compare results
#'fit.lme
#'fit.remlMM1
#'fit.remlMM2
#'fit.remlMM3
#'fit.remlMM4
#'
#'# are there a sex-specific differences?
#'cmat <- getL(fit.remlMM3, c("SexMale-SexFemale", "SexMale:age2-SexFemale:age2")) 
#'cmat
#'
#'test.fixef(fit.remlMM3, L=cmat)
#'}

remlMM <- function(form, Data, by=NULL, VarVC=TRUE, cov=TRUE, quiet=FALSE, order.data=TRUE)
{
	Call <- match.call()
	
	if(!is.null(by))
	{
		stopifnot(is.character(by))
		stopifnot(by %in% colnames(Data))
		stopifnot(is.factor(by) || is.character(by))
		
		levels  <- unique(Data[,by])
		res <- lapply(levels, function(x) {
					tmp.res <- try(remlMM(form=form, Data[Data[,by] == x,], VarVC=VarVC, cov=cov, quiet=quiet), silent=TRUE)
					if(is(tmp.res, "try-error") && !quiet)
						warning(paste0("Error for '", by, ".", x, "':\n\t", attr(tmp.res, "condition")$message))
					tmp.res
				})
		names(res) <- paste(by, levels, sep=".")
		return(res)
	}
	
	stopifnot(class(form) == "formula")
	stopifnot(is.logical(VarVC))
	stopifnot(is.logical(quiet))
	stopifnot(identical(class(Data),"data.frame"))
	stopifnot(nrow(Data) > 2)                                               # at least 2 observations for estimating a variance
	
	if(is.null(.GlobalEnv$msgEnv))												# may removed after loading the package
		msgEnv <<- new.env(parent=emptyenv())
	
	org.form <- form
	trms <- terms(form, simplify=TRUE, keep.order=TRUE)						# convert VCA-formula to valid lmer-formula
	stopifnot(attr(trms, "response") == 1)
	resp <- rownames(attr(trms, "factors"))[1]
	org.form <- form
	Data <- orderData(Data, trms, quiet=quiet,
			order.data=order.data)
	
	if(any(!grepl("^[^[\\.]]?[[:alnum:]_\\.]*$", rownames(attr(trms, "factors")))))
		stop("There are terms in the model formula where regular expression '^[^[\\.]]?[[:alnum:]_\\.]*$' does not fit!")
	
	stopifnot(resp %in% colnames(Data))
	stopifnot(is.numeric(Data[,resp]))
	
	res <- list(call=Call,  EstMethod="REML", data=Data, terms=trms,
			response=resp)
	
	allObsEqual <- FALSE
	
	if(all(Data[,resp] == Data[1,resp]))										# no variance detectable?
	{
		allObsEqual <- TRUE
		if(!quiet)
			warning("All values of response variable ", paste0("'", resp, "'"), " are equal!")
	}
	
	int <- res$intercept <- attr(trms, "intercept") == 1						# has intercept	
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
	
	if(any(resp.NA))								# remove missing data and warn
	{    
		rmInd <- c(rmInd, which(resp.NA))
		if(!quiet)
			message("There are ", length(which(resp.NA))," missing values for the response variable (obs: ", paste(which(resp.NA), collapse=", "), ")!")
		res$resp.NA <- rmInd
	}    
	
	fac  	<- attr(trms, "term.labels")
	if(length(fac) > 1)
	{
		rf.ind  <- which(apply(sapply(rf, function(x) regexpr(x, fac)), 1, function(x) any(x>0)))		
	}
	else
	{
		if(length(rf) > 0)
		{
			if(rf == fac)
				rf.ind <- 1
			else
				rf.ind <- numeric(0)
		}
	}
	vars    <- rownames(attr(trms, "factors"))[-1]	                        # remove response
	Nvc     <- length(fac) + 1 
	Order   <- NULL
	
	for(i in rev(vars))														# check Data for consistency
	{
		if( any(is.na(Data[,i])))
		{
			NAind <- which(is.na(Data[,i]))
			rmInd <- c(rmInd, NAind)
			if(!quiet)
				message("Variable '", i,"' has ",length(NAind)," missing values (obs: ", paste(NAind, collapse=", "), ")!" )
		}
	}
	
	if(length(rf.ind) > 0)													# at least one random term in 'form'
	{
		res$random <- fac[rf.ind]
		res$fixed  <- fac[-rf.ind]
	}
	else
	{
		if(!quiet)
			warning("No random terms in the model! Call 'anovaMM' insead!")
		
		return(anovaMM(form, Data))
	}
	
	res$terms.classes <- sapply(Data[,rownames(attr(trms, "factors"))[-1]], class)
	
	res$Type <- if(length(res$fixed) == 0)
				"Random Model"
			else
				"Mixed Model"			
	
	lab       <- attr(trms, "term.labels")	
	fixed     <- paste(res$fixed, collapse="+")	
	var.num   <- names(res$terms.classes[which(res$terms.classes == "numeric")])
	fac 	  <- attr(trms, "factors")[var.num,res$fixed,drop=FALSE]		# restrict to columns of fixed effects variables and rows with numeric variables
	fe.num 	  <- character()
	
	fac 	  <- fac[which(apply(fac, 1, any)),,drop=FALSE]
	fac 	  <- fac[,which(apply(fac, 2, any)),drop=FALSE]
	fe.num 	  <- colnames(fac)
	iacs      <- which(grepl(":", res$random))										# interaction terms in the formula
	num.rand  <- sapply(var.num, function(x) grepl(x, res$random)) 
	num.rand  <- as.matrix(num.rand)
	num.fixed <- sapply(var.num, function(x) grepl(x, res$fixed))
	random 	  <- ""
	trmMat    <- matrix(nrow=length(res$random), ncol=2, dimnames=list(NULL, c("form", "subject")))
	
	for(i in 1:length(res$random))			# process each random term
	{
		if(nchar(random) == 0)
			sep <- ""
		else
			sep <- "+"
		
		if(length(num.rand) > 0 && any(num.rand[i,]))				# random term with numeric variable found
		{
			if(i %in% iacs)
			{
				splt 	<- unlist(strsplit(res$random[i], ":"))
				tmp.num <- which(splt %in% var.num)
				numVar 	<- splt[ tmp.num]
				others  <- paste(splt[-tmp.num], collapse=":")				
				ind 	<- which(trmMat[,2] == others)				# subject terms may only occur once
				
				if(length(ind) > 0)
					trmMat[ind,2] <- NA
				
				trmMat[i,] <- c(numVar, others)	
			}
			else
				stop("Numeric variables may only occur in random interaction terms!")
		}
		else
			trmMat[i,] <- c("1", res$random[i])
	}
	trmMat 	<- na.omit(trmMat)
	random 	<- paste( apply(trmMat, 1, function(x) paste("(", x[1], "|", x[2], ")", sep="")), collapse="+")	
	form 	<- paste(  resp, "~", fixed, "+", random, sep="")
	form 	<- as.formula(form)
	vcol 	<- rownames(attr(trms, "factors"))
	
	if(is.null(vcol))
		vcol <- resp
	
	Ndata <- nrow(Data)
	Data  <- na.omit(Data[,vcol, drop=F])
	Nobs  <- nrow(Data)
	
	if(quiet || allObsEqual)
		suppressWarnings(fit <- lmer(form, Data, control=lmerControl(optimizer="bobyqa")))		# fit via 'lmer'
	else
		fit <- lmer(form, Data, control=lmerControl(optimizer="bobyqa"))
	
	
	res$Mean 		 <- mean(Data[,resp], na.rm=TRUE)
	res$formula		 <- org.form						# as.specified by the user
	res$form 	 	 <- form							# as used in the lmer-call
	res$NegVCmsg 	 <- ""
	res$VarVC.method <- "gb"
	
	if(int)
	{
		X <- matrix(1, ncol=1, nrow=nrow(Data))			# design matrix of fixed effects: include intercept --> needs a restriction
		colnames(X) <- "int"
		fe.assign <- 0									# '0' indicates intercept
	}
	else
	{
		fe.assign <- NULL
		X <- matrix(1, ncol=0, nrow=nrow(Data))	
	}
	
	fixed.form <- as.formula(paste(	resp, "~", 
					paste(	ifelse(int, "1", "-1"), 
							fixed, sep=ifelse(fixed=="", "", "+")),
					sep=""))
	
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
	
	res$fe.assign 	<- fe.assign										# mapping columns of X to fixed terms in the model formula
	res$fixed.terms <- fixed.terms
	res$balanced <- if(isBalanced(as.formula(trms), Data)) 
				"balanced"  
			else 
				"unbalanced" 
	
	if(Nobs != Ndata)
		res$Nrm <- Ndata - Nobs                         # save number of observations that were removed due to missing data
	
	res$Nobs <- Nobs
	
	tmp <- lmerSummary(	obj=fit, VarVC=VarVC, 			# construct table similar to aov-table and approximate vcovVC
			terms=res$random,
			Mean=res$Mean,
			cov=cov, X=X2)				
	
	tmp$Matrices$y <- Matrix(Data[,resp], ncol=1)
	res <- c(res, tmp)
	class(res) <- "VCA"
	
	if(allObsEqual)
	{
		res$aov.tab[,"DF"] <- NA
		res$aov.tab[,"%Total"] <- 0	
	}
	res
}



#'Construct Variance-Covariance Matrix of Random Effects for Models Fitted by Function 'lmer'
#'
#'This function restricts the variance-covariance matrix of random effects \eqn{G} to be either
#'diagonal ('cov=FALSE') or to take any non-zero covariances into account (default, 'cov=TRUE').
#'
#'This function is not intended to be called directly by users and therefore not exported!
#'
#'@param obj		(object) inheriting from class 'lmerMod'
#'@param cov		(logical) TRUE = in case of non-zero covariances a block diagonal matrix will be constructed,
#'FALSE = a diagonal matrix with all off-diagonal element being equal to zero will be contructed
#'
#'@return (Matrix) representing the variance-covariance structure of random effects \eqn{G}
#'
#'@author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#'
#'@examples 
#'\dontrun{
#'library(lme4)
#'data(Orthodont)
#'Ortho <- Orthodont
#'Ortho$age2 <- Ortho$age - 11
#'Ortho$Subject <- factor(as.character(Ortho$Subject))
#'fit <-lmer(distance~Sex+Sex:age2+(age2|Subject), Ortho) 
#'G1 <- VCA:::lmerG(fit, cov=FALSE)
#'G2 <- VCA:::lmerG(fit, cov=TRUE)
#'G1[1:10,1:10]
#'G2[1:10,1:10]
#'}

lmerG <- function(obj, cov=FALSE)
{
	stopifnot(inherits(obj, "lmerMod"))
	
	vc  <- VarCorr(obj)
	Nre <- unlist(lapply(lme4::ranef(obj), nrow))					# number of random effects per variance component 
	
	lst <- list()
	for(i in 1:length(Nre))
	{
		tmp.vc <- vc[[i]]
		
		if(!cov && nrow(tmp.vc) > 1)
		{
			tmp.vc <- diag(diag(tmp.vc))
		}
		
		lst <- c(lst, replicate(Nre[i], tmp.vc, simplify=FALSE))	# remove covariances from off-diagonal
		
	}
	G <- bdiag(lst)
	G
}


#'Derive and Compute Matrices for Objects Fitted by Function 'lmer'
#'
#'Function derives and computes all matrices required for down-stream
#'analyses of VCA-objects fitted with REML via function \code{\link{lmer}}.
#'
#'Mixed Model Equations (MME) are solved for fixed and random effects applying the same
#'constraints as in \code{\link{anovaMM}}. 
#'The most elaborate and therefore time consuming part is to prepare all matrices required for 
#'approximating the variance-covariance matrix of variance components (see \code{\link{getGB}}).
#'To reduce the computational time, this function tries to optimize object-classes depending
#'on whether Intel's (M)ath (K)ernel (L)ibrary could be loaded or not. MKL appears to be more
#'performant with ordinary matrix-objects, whereas all other computations are perfomred using
#'matrix-representations of the \code{Matrix}-package.
#'
#'This function is not intended to be called directly by users and therefore not exported.
#'
#'@param obj		(object) inheriting from 'lmerMod'
#'@param tab		(data.frame) representing the basic VCA-table
#'@param terms		(character) vector used for ordering variance components
#'@param cov		(logical) take non-zero covariances among random effects into account (TRUE) or
#'not (FALSE), the latter is the default in this package and also implemented in
#'\code{\link{remlVCA}}, \code{\link{anovaVCA}}, and \code{\link{anovaMM}}.
#'@param X			(matrix) design matrix of fixed effects as constructed to meet VCA-package requirements
#'
#'@return (list), a premature 'VCA' object
#'
#'@seealso \code{\link{remlVCA}}, \code{\link{remlMM}}
#'
#'@author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}

lmerMatrices <- function(obj, tab=NULL, terms=NULL, cov=FALSE, X=NULL)
{
	stopifnot(inherits(obj, "lmerMod"))
	
	re.org <- re <- lme4::ranef(obj)	# use lme4's ranef S3-method
	
	VCnam  <- NULL
	reInd  <- list()
	last   <- 0
	count  <- 1
	REnam  <- names(re)
	REnamZ <- NULL
	
	for(i in 1:length(re))				# transform random effects into VCA-compatible structure and derive					
	{									# column-index in random effects design matrix Z for all random effects
		if(ncol(re[[i]]) > 1)			# regression model, i.e. random effects in multi-column matrix
		{
			REi 	<- re[[i]]
			NRi 	<- nrow(REi)
			NCi 	<- ncol(REi)
			trm 	<- names(re)[i]
			tmp.re 	<- ind <- NULL			
			
			for(j in 1:NCi)
			{
				reInd[[count]] <- seq(last+1, last+NRi*NCi, by=NCi)		# columns in Z corresponding to a specific random effect
				
				if(j == 1)
				{
					tmp.nam <- paste(trm, rownames(re[[i]]), sep="")
					names(reInd)[count] <- REnam[count]
				}
				else
				{
					tmp.nam <- c(tmp.nam, paste(colnames(REi)[j], tmp.nam, sep=":"))
					names(reInd)[count] <- paste(colnames(REi)[j], REnam[i], sep=":")
				}
				last  <- last + 1
				count <- count + 1
				
				if(j == 1)
					ind <- eval(parse(text=paste("((",j,"-1)*",NRi,"+1):(", j,"*", NRi, ")", sep="")))
				else
					ind <- paste(ind, eval(parse(text=paste("((",j,"-1)*",NRi,"+1):(", j,"*", NRi, ")", sep=""))), sep=", ")
			}
			ind 	<- paste("c(", paste(ind, collapse=","), ")", sep="")
			ind 	<- eval(parse(text=ind))
			REnamZ 	<- c(REnamZ, tmp.nam[ind])							# column names in Z-matrix			
			cn 		<- colnames(REi)
			cn 		<- paste(names(re)[i], cn, sep=":")
			cn 		<- gsub(":\\(Intercept\\)", "", cn)					# get rif of ()
			VCnam 	<- c(VCnam, cn)										# names of Variance components
		}
		else		# single column random effects matrix
		{
			reInd[[count]] <- (last+1):(last+nrow(re[[i]]))
			names(reInd)[count] <- REnam[count]
			last   <- last + nrow(re[[i]])
			count  <- count + 1
			nami   <- unlist(strsplit(REnam[i], ":"))
			rni    <- rownames(re[[i]])
			VCnam  <- c(VCnam, REnam[i])			
			REnamZ <- c(REnamZ, sapply(rni, function(x) paste(paste(nami, unlist(strsplit(x, ":")), sep=""), collapse=":")))
		}
	}
	REnamZ 		<- as.character(REnamZ)						# names for the random effects and the columns in Z
	reInd  		<- reInd[terms]								# order according to order of variables in user-specified formula
	reInd  		<- reInd[which(!is.na(names(reInd)))]
	Nre	   		<- as.numeric(unlist(lapply(re.org, nrow)))
	sNre   		<- sum(Nre)
	re.assign 	<- list(ind=integer(sNre), terms=names(reInd))
	VC	  	  	<- tab[-c(1, nrow(tab)), "VC"]
	
	for(i in 1:length(reInd))
		re.assign$ind[reInd[[i]]] <- i
	
	Nvc	  <- nrow(tab)-1							# minus total and error								
	Zt 	  <- obj@pp$Zt
	
	if(check4MKL())
		Zt <- as.matrix(Zt)
	
	if(is.null(X) || !is(X, "matrix")) 
		X <- model.matrix(obj, type="fixed")
	
	Xt <- t(X)
	Z  <- t(Zt)
	colnames(Z) <- REnamZ
	G  <- lmerG(obj, cov=cov)						# construct G-matrix
	
	if(check4MKL())
		R	<- diag(nrow(obj@frame))* tab[nrow(tab),"VC"]
	else
		R	<- Diagonal(nrow(obj@frame))* tab[nrow(tab),"VC"]
	
	V	 <- Z %*% G %*% Zt + R						# variance-covariance matrix of observations
	
	Vi	 <- solve(V)
	ViX  <- Vi %*% X
	XtVi <- Xt %*% Vi
	P 	 <- MPinv(Xt %*% ViX)   					# re-compute since 'vcov(obj)' differs from e.g. SAS PROC MIXED
	Q	 <- Vi - ViX %*% P %*% XtVi
	T    <- P %*% XtVi
	
	res 		 <- list()
	res$Nvc 	 <- Nvc
	res$VCnames  <- c(terms, "error") 
	y  <- obj@resp$y
	fe <- P %*% XtVi %*% y						# solve mixed model equations for fixed effects
	
	colnames(fe) <- "Estimate"
	rownames(fe) <- colnames(X)
	iind <- which(names(fe) == "(Intercept)")
	
	if(length(iind) > 0)
		names(fe)[iind] <- "int"
	
	res$Matrices 		<- list(Zre=Z, G=G, R=R, V=V, Vi=Vi, Q=Q, X=X, T=T, y=y)
	res$FixedEffects  	<- fe	
	re  				<- G %*% Zt %*% Vi %*% (y - X %*% fe) 		# estimate random effects solving Mixed Model Equations
	re 					<- Matrix(re)								# re-order random effects accordingly
	rownames(re) 		<- REnamZ
	res$RandomEffects 	<- re
	res$re.assign 		<- re.assign
	res$VarFixed  		<- P
	res$ColOrderZ 		<- reInd			# keep information about original col-indexing
	res
}

#'Derive VCA-Summary Table from an Object Fitted via Function \code{\link{lmer}}
#'
#'This function builds a variance components analysis (VCA) table
#'from an object representing a model fitted by function \code{\link{lmer}}
#'of the \code{lme4} R-package. 
#'
#'It applies the approximation of the variance-covariance
#'matrix of variance components according to Giesbrecht & Burns (1985) and uses this
#'information to approximate the degrees of freedom according to Satterthwaite
#'(see SAS PROC MIXED documentation option 'CL').
#'
#'This function can be used to create a VCA-results table from almost any fitted 'lmerMod'-object, i.e. one can
#'apply it to a model fitted via function \code{\link{lmer}} of the \code{lme4}-package. The only 
#'additional argument that needs to be used is 'tab.only' (see examples).
#'
#'@param obj		(lmerMod) object as returned by function \code{\link{lmer}}
#'@param VarVC		(logical) TRUE = the variance-covariance matrix of variance components will be approximated
#'					following the Giesbrecht & Burns approach, FALSE = it will not be approximated	
#'@param terms		(character) vector, optionally defining the order of variance terms to be used
#'@param Mean		(numeric) mean value used for CV-calculation
#'@param cov		(logical) TRUE = in case of non-zero covariances a block diagonal matrix will be constructed,
#'					FALSE = a diagonal matrix with all off-diagonal elements being equal to zero will be contructed
#'@param X			(matrix) design matrix of fixed effects as constructed to meet VCA-package requirements
#'@param tab.only	(logical) TRUE = will return only the VCA-results table as 'data.frame', argument 'VarVC' will 
#'					be automatically set to 'FALSE' (see details)
#'
#'@return (list) still a premature 'VCA'-object but close to a "complete" 'VCA'-object
#'
#'@seealso \code{\link{remlVCA}}, \code{\link{remlMM}}
#'
#'@author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#'
#'@references
#'Searle, S.R, Casella, G., McCulloch, C.E. (1992), Variance Components, Wiley New York
#'
#'Giesbrecht, F.G. and Burns, J.C. (1985), Two-Stage Analysis Based on a Mixed Model: Large-Sample
#'Asymptotic Theory and Small-Sample Simulation Results, Biometrics 41, p. 477-486 
#'
#'@examples 
#'\dontrun{
#'# fit a model with a VCA-function first
#'data(VCAdata1)
#'fit0 <- remlVCA(y~(device+lot)/day/run, subset(VCAdata1, sample==5))
#'fit0
#'
#'# fit the same model with function 'lmer' of the 'lme4'-package
#'library(lme4)
#'fit1 <- lmer(y~(1|device)+(1|lot)+(1|device:lot:day)+(1|device:lot:day:run),
#'subset(VCAdata1, sample==5))
#'lmerSummary(fit1, tab.only=TRUE)
#'}

lmerSummary <- function(obj, VarVC=TRUE, terms=NULL, Mean=NULL, cov=FALSE, X=NULL, tab.only=FALSE)
{
	stopifnot(inherits(obj, "lmerMod"))	
	if(tab.only)
	{
		VarVC <- FALSE
		
		if(is.null(Mean))
			Mean <- mean(obj@resp$y, na.rm=TRUE)
	}
	
	Sum  <- as.data.frame(summary(obj)$varcor)
	
	Sum[which(Sum[,"var1"] %in% c(NA, "(Intercept)")), "var1"] <- ""
	
	Sum[,"grp"] <- apply(Sum[,c("grp", "var1")], 1, 
			function(x)
			{
				if(x[2] == "")
					return(x[1])
				else
					return(paste(rev(x), collapse=":"))
			}
	)
	re.cor <- NULL
	
	if(any(!is.na(Sum[,"var2"])))
	{
		ind.cor <- which(!is.na(Sum[,"var2"]))
		re.cor  <- Sum[ind.cor,,drop=FALSE]
		Sum <- Sum[-ind.cor,]
	}
	
	if(tab.only)
	{
		if(is.null(terms))
		{
			terms <- Sum[,"grp"]
			terms <- terms[-length(terms)]
		}
	}
	
	Sum  <- Sum[,-c(2,3)]
	rownames(Sum) <- Sum[,"grp"]
	Sum <- Sum[c(terms, "Residual"), ]
	colnames(Sum) <- c("Name", "VC", "SD")
	Sum <- rbind(c(Name="total", VC=sum(Sum[,"VC"]), SD=sqrt(sum(Sum[,"VC"]))), Sum)
	rownames(Sum) <- Sum[,"Name"]
	Sum$VC <- as.numeric(Sum$VC)
	Sum$SD <- as.numeric(Sum$SD)
	Sum$Perc <- 100*Sum$VC/Sum$VC[1]
	Sum$CV	 <- 100*Sum$SD/Mean
	
	if(!tab.only)										# complete VCA-object shall be created
	{
		obj <- lmerMatrices(obj, tab=Sum, terms=terms,	# compute some required matrices
				cov=cov, X=X)							
		obj$aov.tab <- Sum								# required for Giesbrecht & Burns approximation
	}
	
	if(VarVC)
	{
		varVC <- obj$VarCov <- getGB(obj)				# apply Giesbrecht & Burns approximation optimized for MKL
		varVC <- diag(varVC)
		varVC <- c(sum(obj$VarCov), varVC)				# variance of total is sum of all elements of the variance-covariance matrix
		seVC  <- sqrt(varVC)
		Sum$varVC <- varVC
		Sum$Wald <- Sum$VC/seVC							# Wald-statistic
		Sum$DF	 <- 2*Sum$Wald^2						# see SAS-Doc of PROC MIXED
		Sum 	 <- Sum[,c(8,2,4,3,5,6)]
		colnames(Sum)[c(3,5,6)] <- c("%Total", "CV[%]", "Var(VC)")
	}
	else
	{	
		Sum <- Sum[,c(2,4,3,5)]
		colnames(Sum)[c(2,4)] <- c("%Total", "CV[%]")
	}
	
	if(check4MKL())			# remaining part of the package computes with Matrix-package
	{
		obj$Matrices$Zre 	<- Matrix(obj$Matrices$Zre)
		obj$Matrices$G   	<- Matrix(obj$Matrices$G)
		obj$Matrices$R 		<- Matrix(obj$Matrices$R)
		obj$Matrices$V 		<- Matrix(obj$Matrices$V)
		obj$Matrices$Vi 	<- Matrix(obj$Matrices$Vi)
		obj$Matrices$Q 		<- Matrix(obj$Matrices$Q)
		obj$Matrices$X	    <- Matrix(obj$Matrices$X)
		obj$Matrices$T	    <- Matrix(obj$Matrices$T)
	}	
	
	rownames(Sum)[nrow(Sum)] <- "error"
	if(tab.only)							# return just the VCA-table
	{
		Sum <- as.matrix(Sum)
		attr(Sum, "Mean") <- Mean
		return(Sum)
	}
	obj$aov.tab <- Sum
	obj$re.cor <- re.cor					# save correlation among random terms
	obj
}



#' Intermediate Precision for remlVCA-fitted objects of class 'VCA'
#' 
#' Intermediate precision in this context here means any sum of variances
#' below the full model originally fitted. A typical use case could be 
#' reproducibility-experiments with a single lot or multiple lots, where
#' a pooled version of within-lab precision shall be determined.
#' 
#' @param obj		(object) of class 'VCA' fitted by 'remlVCA'
#' @param vc		(character) string specifying the variance component
#' 					up to which an intermediate precision shall be derived
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @examples 
#' data(dataEP05A2_3)
#' res <- remlVCA(y~day/run, dataEP05A2_3)
#' IPday <- getIP.remlVCA(res, "day:run")
#' VCAinference(IPday)

getIP.remlVCA <- function(obj, vc){
	stopifnot(class(obj) == "VCA")
	stopifnot(obj$EstMethod == "REML")
	stopifnot(vc %in% rownames(obj$aov.tab)[-1])
	
	vVC <- vcovVC(obj)
	idx <- which(rownames(vVC) == vc)
	tab <- obj$aov.tab[(idx+1):nrow(obj$aov.tab),]
	vVC <- vVC[idx:nrow(vVC), idx:nrow(vVC)]			# restricted variance-covariance matrix of VC
	tVC <- sum(tab[,"VC"])								# total variability
	tSE <- sqrt(sum(vVC))								# standard error of total VC
	wld <- tVC/tSE 
	DF	<- 2*wld^2
	tab[,"%Total"] <- 100 * tab[,"VC"]/tVC
	tab <- rbind(total=c(DF, tVC, 100, sqrt(tVC), 100*sqrt(tVC)/obj$Mean, tSE^2), tab)
	res <- list()
	res$aov.tab <- tab
	res$Mean <- obj$Mean
	res$EstMethod <- obj$EstMethod
	res$balanced  <- obj$balanced
	res$Type <- "Random Model"
	res$Nobs <- obj$Nobs
	res$NegVCmsg <- ""
	class(res) <- "VCA"

	res
}