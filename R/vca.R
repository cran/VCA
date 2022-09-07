
# Load/unload C-lib.

.onLoad <- function(lib, pkg)
{
	#library.dynam(chname="VCA", package=pkg, lib.loc=lib)
	# create VCA message environment
	msgEnv <<- new.env(parent=emptyenv())
	check4MKL()
}

.onUnload <- function(lib)
{
	#library.dynam.unload(chname="VCA", libpath=lib)
}	


# TODO: Add comment
# 
# Author: schueta6
###############################################################################


# TODO: Add comment
# 
# Author: schueta6
###############################################################################




#'Fit Linear Mixed Model by ANOVA or REML
#'
#'Function serves as interface to functions \code{\link{anovaMM}} and \code{\link{remlMM}}
#'for fitting a linear mixed model (LMM) either by ANOVA or REML. All arguments applicable
#'to either one of these functions can be specified (see \code{\link{anovaMM}} or \code{\link{remlMM}} for details).
#'
#'Besides offering a convenient interface to both functions for fitting a LMM, this function also provides all elements
#'required for standard task of fitted models, e.g. prediction, testing general linear hypotheses via R-package \code{multcomp},
#'etc. (see examples).
#'
#'@param form		(formula) specifiying the linear mixed model, random effects need to be identified by enclosing
#'them in round brackets, i.e. ~a/(b) will model factor 'a' as fixed and 'b' as random
#'@param Data		(data.frame)  containing all variables referenced in 'form', note that variables can only be of type
#'"numeric", "factor" or "character". The latter will be automatically converted to "factor"
#'@param method	(character) either "anova" to use ANOVA Type-I estimation of variance components or "reml" to use 
#'restricted maximum likelihood (REML) estimation of variance component		
#'@param scale		(logical) TRUE = scale values of the response aiming to avoid numerical problems
#'when numbers are either very small or very large, FALSE = use original scale
#'@param VarVC		(logical) TRUE = variance-covariance matrix of variance components will be computed, FALSE = it will not
#'be computed
#'@param ...		additional arguments to be passed to function \code{\link{anovaMM}} or function \code{\link{remlMM}}.
#'
#'@author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#'
#'@seealso \code{\link{fitVCA}}, \code{\link{anovaMM}}, \code{\link{remlMM}}
#'
#'@examples
#'\dontrun{
#'data(dataEP05A2_2)
#'
#'# assuming 'day' as fixed, 'run' as random
#'# Note: default method is "anova"
#'fitLMM(y~day/(run), dataEP05A2_2)
#'
#'# explicitly request "reml"
#'fitLMM(y~day/(run), dataEP05A2_2, method="reml")
#'
#'# assuming both as random leads to same results as
#'# calling anovaVCA (ANOVA is the default)
#'fitLMM(y~(day)/(run), dataEP05A2_2)
#'anovaVCA(y~day/run, dataEP05A2_2)
#'
#'# now using REML-estimation
#'fitLMM(y~(day)/(run), dataEP05A2_2, "reml")
#'remlVCA(y~day/run, dataEP05A2_2)
#'
#'# use different approaches to estimating the covariance of 
#'# variance components (covariance parameters)
#'# create unbalanced data
#'dat.ub <- dataEP05A2_2[-c(11,12,23,32,40,41,42),]
#'m1.ub <- fitLMM(y~day/(run), dat.ub, VarVC.method="scm")
#'# VarVC.method="gb" is an approximation not relying on quadratic forms
#'m2.ub <- fitLMM(y~day/(run), dat.ub, VarVC.method="gb")		
#'# REML-estimated variance components usually differ from ANOVA-estimates
#'# and so do the variance-covariance matrices
#'m3.ub <- fitLMM(y~day/(run), dat.ub, "reml", VarVC=TRUE)		
#'V1.ub <- round(vcovVC(m1.ub), 12)
#'V2.ub <- round(vcovVC(m2.ub), 12)
#'V3.ub <- round(vcovVC(m3.ub), 12)
#'
#'# fit a larger random model
#'data(VCAdata1)
#'fitMM1 <- fitLMM(y~((lot)+(device))/(day)/(run), VCAdata1[VCAdata1$sample==1,])
#'fitMM1
#'# now use function tailored for random models
#'fitRM1 <- anovaVCA(y~(lot+device)/day/run, VCAdata1[VCAdata1$sample==1,])
#'fitRM1
#'
#'# there are only 3 lots, take 'lot' as fixed 
#'fitMM2 <- fitLMM(y~(lot+(device))/(day)/(run), VCAdata1[VCAdata1$sample==2,])
#'# use REML on this (balanced) data
#'fitMM2.2 <- fitLMM(y~(lot+(device))/(day)/(run), VCAdata1[VCAdata1$sample==2,], "reml")
#'
#'# the following model definition is equivalent to the one above,
#'# since a single random term in an interaction makes the interaction
#'# random (see the 3rd reference for details on this topic)
#'fitMM3 <- fitLMM(y~(lot+(device))/day/run, VCAdata1[VCAdata1$sample==2,])
#'
#'# fit same model for each sample using by-processing
#'lst <- fitLMM(y~(lot+(device))/day/run, VCAdata1, by="sample")
#'lst
#'
#'# fit mixed model originally from 'nlme' package
#'
#'library(nlme)
#'data(Orthodont)
#'fit.lme <- lme(distance~Sex*I(age-11), random=~I(age-11)|Subject, Orthodont) 
#'
#'# re-organize data
#'Ortho <- Orthodont
#'Ortho$age2 <- Ortho$age - 11
#'Ortho$Subject <- factor(as.character(Ortho$Subject))
#'fit.anovaMM1 <- fitLMM(distance~Sex*age2+(Subject)*age2, Ortho)
#'
#'# use simplified formula avoiding unnecessary terms
#'fit.anovaMM2 <- fitLMM(distance~Sex+Sex:age2+(Subject)+(Subject):age2, Ortho)
#'
#'# and exclude intercept
#'fit.anovaMM3 <- fitLMM(distance~Sex+Sex:age2+(Subject)+(Subject):age2-1, Ortho)
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
#'# fit LMM with fixed lot and device effects and test for lot-differences
#'data(VCAdata1)
#'fitS5 <- fitLMM(y~(lot+device)/(day)/(run), subset(VCAdata1, sample==5), "reml")
#'fitS5
#'
#'# apply Tukey-HSD test to screen for lot differences
#'library(multcomp)
#'res.tuk <- glht(fitS5, linfct=mcp(lot="Tukey"))
#'summary(res.tuk)
#'
#'# compact letter display
#'res.tuk.cld <- cld(res.tuk, col=paste0("gray", c(90,60,75)))
#'plot(res.tuk.cld)
#'}

fitLMM <- function(	form, Data, method=c("anova", "reml"), scale=TRUE, VarVC=TRUE, ...)
{
	call 	<- match.call()
	method	<- match.arg(tolower(method[1]), choices=c("anova", "reml"))
	Args	<- list(...)
	
	Args$form  <- form
	Args$Data  <- Data
	
	if(scale)
	{
		if(method=="anova")
			Args$Fun  <- "anovaMM"			
		else
		{
			Args$Fun  	<- "remlMM"
			Args$VarVC 	<- VarVC
		}
		fit <- do.call("Scale", Args)
	}
	else
	{
		if(method=="anova")
			fit <- do.call("anovaMM", Args)
		else
		{
			Args$VarVC 	<- VarVC
			fit 		<- do.call("remlMM", Args)
		}
	}
	
	wasVCA <- FALSE								# record that result was VCA-object
	
	if(is(fit, "VCA"))						# if by-processing was not used 			
	{
		fit <- list(fit)
		wasVCA <- TRUE
	}
	
	for(i in 1:length(fit))								# over each list-element
	{
		if(scale)
			fit[[i]] <- reScale(fit[[i]], VarVC=VarVC)	# first re-scale
		
		fit[[i]]$call 			<- call											# prevent non-informative object-name, e.g. "form"
		fe 						<- fixef(fit[[i]])
		X  						<- getMat(fit[[i]], "X")		
		fit[[i]]$fitted.values	<- as.numeric(as.matrix(X) %*% fe[,"Estimate", drop=F])
		fixed 					<- fit[[i]]$fixed[!grepl(":", fit[[i]]$fixed)]	# no interactions
		fixed					<- sapply(fit[[i]]$data[,fixed,drop=F], class)
		fixed					<- names(fixed[which(fixed != "numeric")])
		fit[[i]]$xlevels 		<- lapply(fixed, function(x) as.character(unique(fit[[i]]$data[,x])))
		names(fit[[i]]$xlevels) <- fixed
		
		attr(fit[[i]]$Type, "fitted.by") <- "fitLMM"	# only LMM fitted by 'fitLMM' can be handled in standard ways
	}
	
	if(wasVCA)								
		fit <- fit[[1]]
	
	fit
}






#'Predictions from a Model Fitted by \code{fitLMM}
#'
#'Model returns fitted values in case \code{newdata} is NULL or evaluates
#'the fitted model employing user-specified data \code{newdata}. The default is that
#'fitted values incorporate fixed effects and random effects, leaving out the (conditional)
#'residuals only. If the interest lies in constraining predictions to the fixed effects only
#'set \code{re=NA} or incorporate just part of the random variability specifying distinct random
#'effects (see \code{re}.
#'
#'@param object			(VCA) object fitted via function \code{\link{fitLMM}}
#'@param newdata			(data.frame) with all variables required for the specified prediction,
#'i.e. the default settings require all variables of the original model, 
#'in case of \code{re=NA}, only variables corresponding to fixed effects are
#'required.
#'@param re				(character) if NULL (default) all random effects will be included, 
#'to restrict predictions to the fixed effects use \code{re=NA}, for
#'a subset of random effects included in predictions use any valid
#'random effects specification, i.e. \code{object$random}
#'@param allow.new.levels	(logical) if new levels (no part of the original fitted model) in newdata
#'are allowed. If FALSE (default), such new values in newdata will trigger
#'an error; if TRUE, then the prediction will use the unconditional (population-level)
#'values for data with previously unobserved levels (or NAs).
#'@param ...				additional arguments passdo or from other methods
#'
#'@return (numeric) vector of prediction results
#'@method predict VCA
#'@author Andre Schuetzenmeister \email{andre.schuetzeneister@@roche.com}
#'
#'@examples 
#'\dontrun{
#'# fit LMM with fixed lot and device effects and test for lot-differences
#'data(VCAdata1)
#'datS5 <- subset(VCAdata1, sample==5)
#'fitS5 <- fitLMM(y~(lot+device)/(day)/(run), datS5, "anova")
#'fitS5
#'
#'# fitted values including fixed and random effects
#'pred0 <- predict(fitS5)
#'pred0
#'# sanity check:
#'all(round(pred0 + resid(fitS5) - datS5$y, 12) == 0)
#'# restrict to fixed effects
#'predict(fitS5, re=NA)
#'# restrict to fixed effects and dayly random effects
#'# see required names
#'fitS5$random
#'predict(fitS5, re="lot:device:day")
#'
#'# check against original 'lmer'-predictions
#'# use version from VCA-package (ordinary data.frame)
#'data(Orthodont, package="VCA")
#'Ortho <- Orthodont
#'Ortho$age2 <- Ortho$age-11
#'# use exactly the same data, same ordering
#'Ortho <- orderData(Ortho, distance ~ Sex * age2 + (Subject) * age2)
#'fit.fitLMM <- fitLMM(distance ~ Sex * age2 + (Subject) * age2, Ortho, "reml")
#'library(lme4)
#'fit.lmer <- lmer(distance ~ Sex + age2 + Sex:age2 + (age2 | Subject), Ortho)
#'# check fitted value first (fixed + random effects)
#'predict(fit.lmer)
#'predict(fit.fitLMM)
#'# restrict to fixed part only
#'predict(fit.lmer, re.form=NA)  
#'predict(fit.fitLMM, re=NA)
#'# user-specified 'newdata'
#'newdata <- Ortho[45:54,]
#'newdata$age2 <- newdata$age2 + 5
#'# include fixed and random effects
#'predict(fit.lmer, newdata)
#'predict(fit.fitLMM, newdata)
#'# generate new data
#'newdata <- Ortho[45:54,]          
#'newdata$age2 <- newdata$age2 + 5
#'# predict on newdata using fixed and random effects
#'predict(fit.lmer, newdata) 
#'predict(fit.fitLMM, newdata)       
#'# restrict prediction to random Subject effects
#'predict(fit.lmer, newdata, re.form=~(1|Subject))        
#'predict(fit.fitLMM, newdata, re="Subject")
#'# restrict prediction to random per-Subject slope
#'predict(fit.lmer, newdata, re.form=~(age2-1|Subject)) 
#'predict(fit.fitLMM, newdata, re="age2:Subject")
#'}

predict.VCA <- function(object, newdata=NULL, re=NULL, allow.new.levels=FALSE, ...)
{
	if(is.null(re))								# all random effects per default
	{
		re <- object$random
		if(is.null(re))							# VCA-model
		{
			if(object$EstMethod == "REML")
				re <- object$VCnames
			else
				re <- object$VCnames[-length(object$VCnames)]
		}
	}
	else
	{
		if(!is.na(re[1]) && !all(re %in% object$random))
			stop("You can only specify random effects in 're' which are available in 'object$random'!")
	}
	
	user.nd <- FALSE							# indicate that no user-defined newdata was specified (default)
	
	if(!is.na(re[1]) && !is.null(newdata))		# all variables found in 're' in 'newdata'?
	{
		re.trms <-  unique(unlist(strsplit(object$random, ":")))
		cnd		<- colnames(newdata)
		if(any(!re.trms %in% cnd))
			stop(paste0("Variable(s): ", paste(re.trms[which(!re.trms %in% cnd)], collapse=", "), " are missing in 'newdata'!"))
	}
	
	if(is.null(attr(object$Type, "fitted.by")))
	{
		warning("Predictions can only be made on models fitted by function 'fitLMM' or 'fitVCA'!")
		return()
	}
	trms  	<- object$fixed.terms
	
	### fixed part first	
	if(is.null(newdata))
	{
		X <- getMat(object, "X")				# get design matrix of fixed effects
		newdata <- object$data	
	}
	else
	{
		user.nd <- TRUE							# indicate user-defined newdata
		Terms	<- delete.response(trms)
		mf 		<- model.frame(Terms, newdata, xlev=object$xlevels)	
		cntr    <- as.list(rep("contr.SAS", length(object$xlevels)))	# use SAS-type contrasts, i.e. last level will be constrained 
		names(cntr) <- names(object$xlevels)
		X		<- model.matrix(Terms, mf, contrasts.arg=cntr)
		asgn	<- attr(X, "assign")
		asgn0	<- object$fe.assign
		tdiff	<- table(asgn0) - table(asgn)	# missing columns in X?
		
		if(any(tdiff == 0))
		{
			X0  <- matrix(0, nrow=nrow(X), ncol=length(asgn0))		# pre-allocate matrix of right dimension
			ind <- match(asgn, asgn0) 								# add all-zero columns to match vector of fixed effects
			for(i in unique(ind))
			{
				tmp <- which(ind == i)
				if(length(tmp) > 1)
					ind[tmp] <- ind[tmp] + 0:(length(tmp)-1)
			}
			X0[,ind] <- X
			X <- X0
		}
	}
	fe 		<- fixef(object, quiet=TRUE)
	preds	<- as.numeric(as.matrix(X) %*% fe[,"Estimate", drop=F])
	
	### add random part to predictions
	if(!is.na(re[1]))
	{
		#Z 	<- getMat(object, "Z")
		ref <- ranef(object, quiet=TRUE)						# at this point all random effects would be used
		
		# generate names for random effects from 'newdata' accounting for user-specified 'newdata' 
		#fac   <- na.omit(names(object$terms.classes[re]))		# only these will get names like VarnameVarlevel 
		rtrms <- re
		#renam <- NULL		
		Z	  <- matrix(nrow=nrow(newdata), ncol=0)
		
		for(i in 1:length(rtrms))						# over all random terms
		{
			tmpZ <- getMM(as.formula(paste0("~", rtrms[i], "-1")), newdata) 	# model.matrix wo intercept
			Z  	 <- cbind(Z, tmpZ) 												# build complete Z-matrix corresponding to newdata
		}	

		ref <- ref[colnames(Z),,drop=F]						# align random effects to column-order of Z
		preds <- preds + as.numeric(Z %*% ref)		
	}	
	names(preds) <- rownames(newdata)
	return(preds)
}





#'Fit Variance Component Model by ANOVA or REML
#'
#'Function serves as interface to functions \code{\link{anovaVCA}} and \code{\link{remlVCA}}
#'for fitting a variance component models (random models) either by ANOVA or REML. All arguments applicable
#'to either one of these functions can be specified (see \code{\link{anovaVCA}} or \code{\link{remlVCA}} for details).
#'
#'@param form		(formula) specifiying the variance component model (see \code{\link{anovaVCA}} and/or \code{\link{remlVCA}})
#'@param Data		(data.frame)  containing all variables referenced in 'form'
#'@param method	(character) either "anova" to use ANOVA Type-I estimation of variance components or "reml" to use 
#'restricted maximum likelihood (REML) estimation of variance component		
#'@param scale		(logical) TRUE = scale values of the response aiming to avoid numerical problems
#'when numbers are either very small or very large, FALSE = use original scale
#'@param VarVC		(logical) TRUE = variance-covariance matrix of variance components will be computed, FALSE = it will not
#'be computed
#'@param ...		additional arguments to be passed to function \code{\link{anovaVCA}} or function \code{\link{remlVCA}}.
#'
#'@author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#'
#'@seealso \code{\link{fitLMM}}, \code{\link{anovaVCA}}, \code{\link{remlVCA}}
#'
#'@examples
#'\dontrun{
#'#load data (CLSI EP05-A2 Within-Lab Precision Experiment) 
#'data(dataEP05A2_2)
#'
#'# perform ANOVA-estimation of variance components
#'res.anova <- fitVCA(y~day/run, dataEP05A2_2, "anova")
#'# perform REML-estimation of variance components
#'res.reml <- fitVCA(y~day/run, dataEP05A2_2, "reml")
#'
#'# compare scaling vs. not scaling the response
#'fit0 <- fitVCA(y~day/run, dataEP05A2_2, "anova", scale=TRUE)
#'fit1 <- fitVCA(y~day/run, dataEP05A2_2, "anova", scale=FALSE)
#'}

fitVCA <- function(form, Data, method=c("anova", "reml"), scale=TRUE, VarVC=TRUE, ...)
{
	call <- match.call()
	method <- match.arg(tolower(method[1]), choices=c("anova", "reml"))
	
	if(scale)
	{
		if(method=="anova")
			fit <- Scale("anovaVCA", form=form, Data=Data, ...)
		else
			fit <- Scale(remlVCA, form=form, Data=Data, VarVC=VarVC, ...)
	}
	else
	{
		if(method=="anova")
			fit <- anovaVCA(form=form, Data=Data, ...)
		else
			fit <- remlVCA(form=form, Data=Data, VarVC=VarVC, ...)	
	}
	
	wasVCA <- FALSE								# record that result was VCA-object
	
	if(is(fit, "VCA"))						# if by-processing was not used 			
	{
		fit <- list(fit)
		wasVCA <- TRUE
	}
	for(i in 1:length(fit))								# over each list-element
	{
		if(scale)
			fit[[i]] <- reScale(fit[[i]], VarVC=VarVC)				# first re-scale
        
		fit[[i]]$call 			<- call								# prevent non-informative object-name, e.g. "form"
		tmp.fit 				<- fit[[i]]							# fixing an error that first occurred testing with R-4.2.1
		fe 						<- fixef(tmp.fit)#fixef(fit[[i]])
		X  						<- getMat(fit[[i]], "X")		
		fit[[i]]$fitted.values	<- as.numeric(as.matrix(X) %*% fe[,"Estimate", drop=F])
		fixed 					<- fit[[i]]$fixed[!grepl(":", fit[[i]]$fixed)]	# no interactions
		fit[[i]]$xlevels 		<- lapply(fixed, function(x) as.character(unique(fit[[i]]$data[,x])))
		names(fit[[i]]$xlevels) <- fixed
		
		attr(fit[[i]]$Type, "fitted.by") <- "fitVCA"	# only LMM fitted by 'fitLMM' can be handled in standard ways
	}
	
	if(wasVCA)								
		fit <- fit[[1]]
	
	fit
}





#'Re-Scale results of 'VCA' or 'VCAinference' 
#'
#'Function adjusts variance components (VC) and standard deviations (SD) and their respective
#'confidence intervals of 'VCAinference' objects, and the 'VCAobj' sub-element. For 'VCA' objects
#'the VC and SD values are adjusted as well as the fixed and random effects and the covariance-matrix
#'of fixed effects.
#'
#'@param obj		(object) either of class 'VCA' or 'VCAinference'
#'@param VarVC		(logical) TRUE = variance-covariance matrix of the fitted model 'obj'
#'will be computed and automatically re-scaled, FALSE = variance-covariance
#'matrix will not be computed and re-scaled. This might cause wrong results
#'in downstream analyses which require this matrix on the correct scale! Only
#'use this option if computation time really matters!
#'
#'@return (object) either of class 'VCA' or 'VCAinference', where results have been 
#'transformed back to the original scale of the response variable 
#'
#'@author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#'
#'@seealso \code{\link{Scale}}
#'
#'@examples 
#'\dontrun{
#'data(dataEP05A2_3)
#'
#'# reference values
#'fit0 <- anovaVCA(y~day/run, dataEP05A2_3, MME=TRUE)
#'inf0 <- VCAinference(fit0, VarVC=TRUE)
#'
#'fit1 <- Scale("anovaVCA", y~day/run, dataEP05A2_3, MME=TRUE)
#'inf1 <- VCAinference(fit1, VarVC=TRUE)
#'inf1 <- reScale(inf1)
#'
#'# compare to reference
#'print(inf0, what="VC")
#'print(inf1, what="VC")
#'print(inf0, what="SD")
#'print(inf1, what="SD")
#'print(inf0, what="CV")
#'print(inf1, what="CV")
#'
#'# now use REML-based estimation
#'fit0 <- remlVCA(y~day/run, dataEP05A2_3)
#'inf0 <- VCAinference(fit0)
#'
#'fit1 <- Scale("remlVCA", y~day/run, dataEP05A2_3, MME=TRUE)
#'inf1 <- VCAinference(fit1)
#'inf1 <- reScale(inf1)
#'
#'# compare to reference
#'print(inf0, what="VC")
#'print(inf1, what="VC")
#'print(inf0, what="SD")
#'print(inf1, what="SD")
#'print(inf0, what="CV")
#'print(inf1, what="CV")
#'}

reScale <- function(obj, VarVC=TRUE)
{
	if(is.list(obj) && !is(obj, "VCA") && !is(obj, "VCAinference"))
	{		
		if(!all(sapply(obj, class) %in% c("VCA", "VCAinference")))
			stop("Only lists of 'VCA' or 'VCAinference' objects are accepted!")
		
		obj.len <- length(obj)
		
		res <- lapply(obj, FUN=reScale)
		names(res) <- names(obj)
		
		if(obj.len == 1)			# mapply returns a list of length 2 in case that length(obj) was equal to 1
			res <- res[1]
		
		return(res)
	}	
	
	if(!is.null(obj$rescaled))		# re-scaling has already been done on this object
		return(obj)
	
	stopifnot(class(obj) %in% c("VCA", "VCAinference"))
	
	if(is(obj, "VCAinference"))
		VCAobj <- obj$VCAobj
	else
		VCAobj <- obj
	
	scale <- VCAobj$scale	
	
	if(VarVC || VCAobj$EstMethod == "ANOVA")						# estimate covariance-matrix of VCs before any re-scaling takes place
	{
		VCAobj <- solveMME(VCAobj)
		VCAobj$VarCov <- vcovVC(VCAobj)
	}
	
	VCAobj$aov.tab[,"VC"] 		<- VCAobj$aov.tab[,"VC"] * scale^2
	VCAobj$aov.tab[,"SD"] 		<- VCAobj$aov.tab[,"SD"] * scale
	if("Var(VC)" %in% colnames(VCAobj$aov.tab))
		VCAobj$aov.tab[,"Var(VC)"]	<- VCAobj$aov.tab[,"Var(VC)"] * scale^4
	if(VCAobj$EstMethod == "ANOVA")
	{
		VCAobj$aov.tab[,"SS"]	<- VCAobj$aov.tab[,"SS"] * scale^2
		VCAobj$aov.tab[,"MS"]	<- VCAobj$aov.tab[,"MS"] * scale^2
	}
	VCAobj$Mean 				<- VCAobj$Mean * scale
	
	VCAobj$Matrices$y <- VCAobj$Matrices$y * scale
	VCAobj$data[,VCAobj$response] <- VCAobj$data[,VCAobj$response] * scale			# re-scale response variable in the data-element

	if(!is.null(VCAobj$FixedEffects))
		VCAobj$FixedEffects 	<- try(VCAobj$FixedEffects * scale,  silent=TRUE)
	if(!is.null(VCAobj$RandomEffects))
		VCAobj$RandomEffects 	<- try(VCAobj$RandomEffects * scale, silent=TRUE)
	if(!is.null(VCAobj$VarFixed))
		VCAobj$VarFixed 		<- try(VCAobj$VarFixed * scale^2, 	 silent=TRUE)
	if(!is.null(VCAobj$VarCov))
		VCAobj$VarCov			<- try(VCAobj$VarCov * scale^4, 	 silent=TRUE)
	
	if(is(obj, "VCAinference"))
	{
		obj$VCAobj <- VCAobj
		
		# variance components
		
		obj$ConfInt$VC$OneSided$LCL <- obj$ConfInt$VC$OneSided$LCL * scale^2
		obj$ConfInt$VC$OneSided$UCL <- obj$ConfInt$VC$OneSided$UCL * scale^2
		
		obj$ConfInt$VC$TwoSided$LCL <- obj$ConfInt$VC$TwoSided$LCL * scale^2
		obj$ConfInt$VC$TwoSided$UCL <- obj$ConfInt$VC$TwoSided$UCL * scale^2
		
		# SD
		
		obj$ConfInt$SD$OneSided$LCL <- obj$ConfInt$SD$OneSided$LCL * scale
		obj$ConfInt$SD$OneSided$UCL <- obj$ConfInt$SD$OneSided$UCL * scale
		
		obj$ConfInt$SD$TwoSided$LCL <- obj$ConfInt$SD$TwoSided$LCL * scale
		obj$ConfInt$SD$TwoSided$UCL <- obj$ConfInt$SD$TwoSided$UCL * scale
	}
	else
		obj <- VCAobj
	
	obj$rescaled <- TRUE			# record that re-scaling was performed
	obj 
}





#'Solve System of Linear Equations using Inverse of Cholesky-Root
#'
#'Function solves a system of linear equations, respectively, inverts a matrix
#'by means of the inverse Cholesky-root.
#'
#'This function is intended to reduce the computational time in function
#'\code{\link{solveMME}} which computes the inverse of the square variance-
#'covariance Matrix of observations. It is considerably faster than function
#'\code{\link{solve}} (see example).
#'Whenever an error occurs, which is the case for non positive definite matrices
#''X', function \code{\link{MPinv}} is called automatically yielding a generalized
#'inverse (Moore-Penrose inverse) of 'X'.
#'
#'@param X			(matrix, Matrix) object to be inverted
#'@param quiet		(logical) TRUE = will suppress any warning, which will be issued otherwise 
#'
#'@return (matrix, Matrix) corresponding to the inverse of X
#'
#'@author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#'
#'@examples
#'\dontrun{
#'# following complex (nonsense) model takes pretty long to fit
#'system.time(res.sw <- anovaVCA(y~(sample+lot+device)/day/run, VCAdata1))
#'# solve mixed model equations (not automatically done to be more efficient)
#'system.time(res.sw <- solveMME(res.sw))
#'# extract covariance matrix of observations V
#'V1 <- getMat(res.sw, "V")
#'V2 <- as.matrix(V1)
#'system.time(V2i <- solve(V2))
#'system.time(V1i <- VCA:::Solve(V1))
#'V1i <- as.matrix(V1i)
#'dimnames(V1i) <- NULL
#'dimnames(V2i) <- NULL
#'all.equal(V1i, V2i)
#'} 

Solve <- function(X, quiet=FALSE)
{
	stopifnot(ncol(X) == nrow(X))
	
	clsMatrix <- inherits(X, "Matrix")
	if(clsMatrix)
	{
		cls <- class(X)
		X <- as.matrix(X)
	}
	Xi <- try(chol2inv(chol(X)), silent=TRUE)
	
	if(inherits(Xi, "try-error"))			# use Moore-Penrose inverse instead in case of an error
	{										# using the Cholesky-decomposition approach
		if(!quiet)
			warning("\tMatrix inversion via 'chol2inv' failed!\n\tUse generalized (Moore-Penrose) inverse (MPinv)!", sep="\n")
		Xi <- MPinv(X)
	}
	
	if(clsMatrix)
	{
		Xi <- as(Xi, "dgCMatrix")			# Solve used here for inverting V-matrix
	}
	return(Xi)
}



#'Solve Mixed Model Equations
#'
#'Function solves the Mixed Model Equations (MME) to estimate fixed and random effects.
#'
#'This function is for internal use only, thus, not exported.
#'
#'@param obj			... (VCA) object
#'
#'@return 	(VCA) object, which has additional elements "RandomEffects" corresponding to the column vector 
#'of estimated random effects, "FixedEffects" being the column vector of estimated fixed effects. 
#'Element "Matrices" has additional elements referring to the elements of the MMEs and element
#'"VarFixed" corresponds to the variance-covariance matrix of fixed effects.
#'
#'@author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#'
#'@examples
#'\dontrun{
#'data(dataEP05A2_1)
#'fit <- anovaVCA(y~day/run, dataEP05A2_1, NegVC=TRUE)
#'fit <- solveMME(fit)
#'ranef(fit)
#'}

solveMME <- function(obj)
{
	stopifnot(class(obj) == "VCA")
	mats   	<- obj$Matrices
	V      	<- mats[["V"]]
	
	if(is.null(V))						# will be the case for ANOVA-type fitted models
	{
		obj  <- getV(obj)
		mats <- obj$Matrices
		V	 <- mats[["V"]]
	}
	Z      	<- mats$Zre
	R 	   	<- mats$R
	G      	<- mats$G
	X      	<- Matrix(mats$X)
	y      	<- mats$y
	
	if(is.null(mats$Vi))
		Vi  <- Solve(V, quiet=TRUE)
	else
		Vi	<- mats$Vi
	
	K 		<- Solve(t(X) %*% Vi %*% X, quiet=TRUE)		# variance-covariance matrix of fixed effects
	T	   	<- K %*% t(X) %*% Vi
	fixed  	<- T %*% y

	mats$Vi <- Vi
	mats$T  <- T
	rownames(fixed) <- colnames(X)
	colnames(fixed) <- "Estimate"
	re		<- obj$RandomEffects
	
	if(is.null(Z))
		re <- NULL
	else
	{
		if(is.null(re))
		{
			re  <- G %*% t(Z) %*% Vi %*% (y - X %*% fixed) 
			rownames(re) <- colnames(Z)
			colnames(re) <- "Estimate"
		}
	}
	obj$RandomEffects <- re
	obj$FixedEffects  <- fixed
	obj$Matrices	  <- mats
	obj$VarFixed      <- K
	return(obj)
}





#'Extract Fixed Effects from 'VCA' Object
#'
#'For conveniently using objects of class 'VCA' with other packages expecting this
#'function, e.g. the 'multcomp' package for general linear hypotheses for parametric
#'models (currently not fully implemented).
#'
#'@param object		(VCA) object where fixed effects shall be extracted
#'@param quiet			(logical) TRUE = will suppress any warning, which will be issued otherwise 
#'@param ...			additional parameters
#'
#'@method coef VCA 
#'
#'@examples 
#'\dontrun{
#'data(dataEP05A2_1)
#'fit1 <- anovaMM(y~day/(run), dataEP05A2_1)
#'coef(fit1)
#'fit2 <- anovaVCA(y~day/run, dataEP05A2_1)
#'coef(fit2)
#'}

coef.VCA <- function(object, quiet=FALSE, ...)
{
	Call <- match.call()
	obj <- object
	stopifnot(class(obj) == "VCA")
	fe  <- fixef(obj)
	if(is.null(fe))								# solve mixed model equations first
	{
		obj  <- solveMME(obj)
		fe   <- fixef(obj)
		nam  <- as.character(as.list(Call)$object)
		if(length(nam) == 1 && nam %in% names(as.list(.GlobalEnv)))
		{
			expr <- paste(nam, "<<- obj")		# update object missing MME results
			eval(parse(text=expr))
		}
		else
		{
			if(!quiet)
				message("Some required information missing! Usually solving mixed model equations has to be done as a prerequisite!")
		}
	}
	fe  <- fe[,"Estimate", drop=F]
	nam <- rownames(fe)
	fe  <- c(fe)
	names(fe) <- nam
	return(fe)
}


#'Calculate Variance-Covariance Matrix of Fixed Effects for an 'VCA' Object
#'
#'Return the variance-covariance matrix of fixed effects for a linear mixed model
#'applicable for objects of class 'VCA'.
#'
#'Actually this function only extracts this matrix or, if not available, calls function
#'\code{\link{vcovFixed}} which performs calculations. It exists for compatibility reasons,
#'i.e. for coneniently using objects of class 'VCA' with other packages expecting this
#'function, e.g. the 'multcomp' package for general linear hypotheses for parametric
#'models.  
#'
#'@param object		 	(VCA) object for which the variance-covariance matrix of
#'fixed effects shall be calculated
#'@param quiet				(logical) TRUE = will suppress any warning, which will be issued otherwise 
#'@param ...				additional parameters
#'
#'@return (matrix) corresponding to the variance-covariance matrix of fixed effects
#'
#'@method vcov VCA
#'
#'@examples 
#'\dontrun{
#'data(dataEP05A2_1)
#'fit1 <- anovaMM(y~day/(run), dataEP05A2_1)
#'vcov(fit1)
#'
#'fit2 <- anovaVCA(y~day/run, dataEP05A2_1)
#'vcov(fit2)
#'}

vcov.VCA <- function(object, quiet=FALSE, ...)
{
	obj <- object
	stopifnot(class(obj) == "VCA")
	if(!obj$intercept && length(obj$fixed) == 0)
	{
		if(!quiet)
			warning("There is no variance-convariance matrix of fixed effects for this object!")
		return(NA)
	}
	if(is.null(obj$VarFixed))
		return(vcovFixed(obj, quiet=quiet))
	else
		return(obj$VarFixed)
}


#'Degrees of Freedom for Testing Linear Contrasts of Fixed Effects and Least Square Means
#'
#'There are three methods implemented, which are all available in SAS PROC MIXED, namely 
#'"contain", "residual", and "satterthwaite" approximations. See the documentation of SAS
#'PROC MIXED for details on this topic.
#'
#'The implementation of the Satterthwaite approximation was inspired by the code of function 
#''calcSatterth' of R-package 'lmerTest'.
#'
#'@param obj			(VCA) object
#'@param L				(numeric) vector specifying the linear combination of the fixed effect or
#'						LS Means
#'@param ddfm			(character) string specifying the method used for computing the denominator
#'						degrees of freedom for tests of fixed effects or LS Means. Available methods are
#'						"contain", "residual", and "satterthwaite".
#'@param tol			(numeric) value specifying the numeric tolerance for testing equality to zero
#'@param method.grad	(character) string specifying the method to be used for approximating the gradient
#'						of the variance-covariance matrix of fixed effects at the estimated covariance parameter
#'						estimates (see function 'grad' (numDeriv) for details)
#'@param opt			(logical) TRUE = tries to optimize computation time by avoiding unnecessary computations
#'						for balanced datasets (see \code{\link{test.fixef}}). 
#'@param items			(list) of pre-computed values
#'
#'@return (numeric) vector with the specified type of degrees of freedom
#'
#'@author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#'
#'@seealso \code{\link{test.fixef}}

getDDFM <- function(obj, L, ddfm=c("contain", "residual", "satterthwaite"), tol=1e-12, method.grad="simple", opt=TRUE, items=NULL)
{
	stopifnot(class(obj) == "VCA")
	stopifnot(!is.null(colnames(L)))
	ddfm <- match.arg(ddfm)
	
	if(ddfm == "residual")
	{
		return(obj$Nobs - rankMatrix(getMat(obj, "X")))			# as described in documentation of SAS PROC MIXED
	}
	else if(ddfm == "contain")
	{		
		cn <- colnames(L)										# fe or LS Means names
		cn <- cn[which(L[1,] != 0)]

		if(length(cn) == 1 && cn == "int")						# handle intercept fixed effect
		{
			return(min(obj$aov.org[obj$random, "DF"]))
		}
		fe <- obj$fixed

		tmp <- sapply(fe, function(x) gregexpr(x, cn))			# can names of fixed terms be found in column names of L?

		if(is(tmp, "list")) {								# there was only a single non-zero column in L --> list returned
			fe <- sapply(tmp, function(x) x != -1)
		} else {												# mulitple non-zero columns in L
			fe <- apply(tmp, 2, function(y) any(y==1))
		}		

		if(any(fe))
		{
			fe <- fe[which(fe)]
			fe <- names(fe)
			rn <- obj$random
			rn <- sapply(rn, function(x) all(sapply(fe, grepl, x)))

			if(any(rn))
			{
				DF <- min(obj$aov.org[names(rn), "DF"])
				return(DF)
			}
			else
			{
				tmpZ <-  getMat(obj, "Z")
				
				if(!is.null(tmpZ))
					return(obj$Nobs - rankMatrix(cbind(as.matrix(getMat(obj, "X")), as.matrix(tmpZ))))
				else
					return(obj$Nobs - rankMatrix(as.matrix(getMat(obj, "X"))))
			}
		}
		else
		{
			return(obj$Nobs - rankMatrix(cbind(as.matrix(getMat(obj, "X")), as.matrix(getMat(obj, "Z")))))	# no random term including the fixed effects term --> residual DFs
		}
	}
	else if(ddfm == "satterthwaite")
	{
		stopifnot(class(obj) == "VCA")		
		if(is.null(dim(L)))
			L <- matrix(L, nrow=1)		
		
		if(obj$balanced == "balanced" && opt && obj$EstMethod != "REML" && FALSE)
		{
			return(obj$aov.tab[2,"DF"])			
		}
		r       <- items$r
		SVD 	<- items$SVD
		nu.m 	<- NULL											# see SAS-help for SAS PROC MIXED -> model /ddfm=sat
		A <- obj$VarCov											# same names as in SAS Help of PROC MIXED, Model statement, option "ddfm=sat"
		if(is.null(A))
			A <- vcovVC(obj)
		
		VCs <- obj$aov.tab[-1, "VC"]							# all variance components but total
		
		for( m in 1:length(SVD$values) )
		{   
			g <- grad(function(x)  L %*% DfSattHelper(obj, x) %*% t(L),  VCs, method = method.grad)
			nu.m <- c(nu.m, 2 * (SVD$values[m])^2/(t(g) %*% A %*% g))
		}
		
		E  <- sum( (nu.m/(nu.m-2)) * as.numeric(nu.m>2))
		DF <- 2*E*as.numeric(E>r)/(E-r)
		
		return(DF)
	}
}






#'Calculate Variance-Covariance Matrix of Variance Components of 'VCA' objects
#'
#'This function computes the variance-covariance matrix of variance components (VC) either
#'applying the approach given in the \eqn{1^{st}}{1st} reference ('method="scm"') or using
#'the approximation given in the \eqn{2^{nd}}{2nd} reference ('method="gb"').
#'
#'
#'This function is called on a 'VCA' object, which can be the sole argument. In this case the value
#'assigned to element 'VarVC.method' of the 'VCA' object will be used.
#'
#'@param obj			(VCA) object
#'@param method		(character) string, optionally specifying whether to use the algorithm given in the
#'1st reference ("scm") or in the 2nd refernce ("gb"). If not not supplied, the 
#'option is used coming with the 'VCA' object.
#'@param quiet			(logical) TRUE = will suppress any warning, which will be issued otherwise 
#'
#'@return (matrix) corresponding to variance-covariance matrix of variance components
#'
#'@author 	Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com},
#'Florian Dufey \email{florian.dufey@@roche.com}
#'
#'@references
#'Searle, S.R, Casella, G., McCulloch, C.E. (1992), Variance Components, Wiley New York
#'
#'Giesbrecht, F.G. and Burns, J.C. (1985), Two-Stage Analysis Based on a Mixed Model: Large-Sample
#'Asymptotic Theory and Small-Sample Simulation Results, Biometrics 41, p. 477-486
#'
#'@examples
#'\dontrun{
#'data(realData)
#'dat1 <- realData[realData$PID==1,]
#'fit  <- anovaVCA(y~lot/calibration/day/run, dat1) 
#'vcovVC(fit)
#'vcovVC(fit, "scm")		# Searle-Casella-McCulloch method (1st reference)
#'vcovVC(fit, "gb")		# Giesbrecht and Burns method (2nd reference)
#'}

vcovVC <- function(obj, method=NULL, quiet=FALSE)
{
	Call <- match.call()
	
	stopifnot(class(obj) == "VCA")
	if(!is.null(method))
		method <- match.arg(method, c("scm", "gb"))
	else
		method <- obj$VarVC.method
	
	VCvar <- obj$VarCov
	
	if(obj$EstMethod == "REML")
	{
		if(is.null(VCvar))
		{
			if(!quiet)
				warning("When fitting a model via REML, set 'VarVC=TRUE' for computing the variance-covariance matrix of variance components!")
			return(NULL)
		}
		if(method == "scm" && !quiet)
			warning("For models fitted by REML only the Giesbrecht & Burns method is applicable!")
	}
	
	if(!is.null(VCvar))
		return(VCvar)
	
	Z  <- obj$Matrices$Z
	
	if(method == "scm")								
	{
        VCvar <-obj$Matrices$VCvar
        rownames(VCvar) <- colnames(VCvar) <- obj$VCnames
	}
	else
	{
		if(is.null(obj$VarFixed))
		{
			obj  <- solveMME(obj)
			nam  <- as.character(as.list(Call)$obj)
			if(length(nam) == 1 && nam %in% names(as.list(.GlobalEnv)))
			{
				expr <- paste(nam, "<<- obj")		# update object missing MME results
				eval(parse(text=expr))
			}
			else
			{
				if(quiet)
					message("Some required information missing! Usually solving mixed model equations has to be done as a prerequisite!")
			}
		}		
		VCvar <- getGB(obj)							# apply Giesbrecht & Burns approximation
	}
	attr(VCvar, "method") <- method					# which method was used
	
	return(VCvar)
}


#'Standard Printing Method for Objects of Class 'VCA'
#'
#'Function prints 'VCA' objects as returned e.g. by function \code{\link{anovaVCA}}.
#'
#'@param x         (VCA) object of class 'VCA' as returned by function 'anovaVCA'.
#'@param digits    (integer) number of digits numeric values are rounded to before printing.
#'@param ...       additional arguments to be passed to or from methods.
#'
#'@method print VCA

print.VCA <- function(x, digits=6L, ...)
{
	EstMeth <- x$EstMethod
	Mean <- x$Mean
	Nobs <- x$Nobs
	Nrm  <- x$Nrm
	balanced <- x$balanced
	NegVCmsg <- ifelse(is.null(x$NegVCmsg), "", x$NegVCmsg)
	
	if(x$Type == "Linear Model")
	{
		if(!"skipHeading" %in% names(list(...)))
			cat("\n\nAnalysis of Variance Table:\n---------------------------\n\n")
		tmp <- x$aov.org
		tmp$DF <- round(tmp$DF)
		printCoefmat(tmp, digits=digits, dig.tst=digits, has.Pvalue=TRUE, P.values=5, 
				na.print="", zap.ind=3, tst.ind=4, cs.ind=NULL)
		cat("\n")
		
	}
	else
	{
		if(x$EstMethod == "ANOVA")
		{
			MM <- x$Type == "Mixed Model" 
			if(!"skipHeading" %in% names(list(...)) && !MM)
				cat("\n\nResult Variance Component Analysis:\n-----------------------------------\n\n")
			if(!"skipHeading" %in% names(list(...))&& MM)
			{
				cat("\n\nANOVA-Type Estimation of Mixed Model:\n--------------------------------------\n\n")
			}
			if(MM)
			{
				fe <- fixef(x)[,"Estimate", drop=FALSE]
				nam <- rownames(fe)
				fe <- c(fe)
				names(fe) <- nam
			}
			rn <- rownames(x$aov.tab)
			cn <- colnames(x$aov.tab)
			
			#		if(!MM)
			NegVCind <- which(x$VCoriginal < 0) + 1                    	# "VCoriginal" does not contain total variance --> +1
			mat <- matrix(x$aov.tab, nrow=length(rn), dimnames=list(rn, cn))
			mat <- apply(mat, 1:2, round, digits=digits)
			mat <- cbind(Name=rn, mat)
			if("CV" %in% colnames(mat))
				colnames(mat)[which(colnames(mat) == "CV")] <- "CV[%]"
			rownames(mat) <- 1:nrow(mat)
			if(NegVCmsg != "")
			{
				mat[NegVCind, c("VC", "%Total", "SD", "CV[%]")] <- paste(mat[NegVCind, "VC"], "*", sep="")
			}
			mat <- apply(mat, 1:2, function(x) ifelse(x==NaN, NA, x))
			if(MM)
			{
				cat("\t[Fixed Effects]\n\n")
				print(round(fe, digits))
				cat("\n\n\t[Variance Components]\n\n")
			}
			print(mat, na.print="", quote=FALSE, digits=digits)
		}
		else
		{
			MM <- x$Type == "Mixed Model" 
			if(!"skipHeading" %in% names(list(...)) && !MM)
				cat("\n\nResult Variance Component Analysis:\n-----------------------------------\n\n")
			if(!"skipHeading" %in% names(list(...))&& MM)
			{
				cat("\n\nREML-Estimation of Mixed Model:\n-------------------------------\n\n")
			}
			if(MM)
			{
				fe <- fixef(x)[,"Estimate", drop=FALSE]
				nam <- rownames(fe)
				fe <- c(fe)
				names(fe) <- nam
			}
			rn <- rownames(x$aov.tab)
			cn <- colnames(x$aov.tab)
			mat <- as.matrix(x$aov.tab)
			mat <- apply(mat, 1:2, round, digits=digits)
			mat <- cbind(Name=rn, mat)
			rownames(mat) <- 1:nrow(mat)
			
			if(MM)
			{
				cat("\t[Fixed Effects]\n\n")
				print(round(fe, digits))
				cat("\n\n\t[Variance Components]\n\n")
			}
			
			print(mat, na.print="", quote=FALSE, digits=digits)
		}
	}
	
	cat("\nMean:", round(Mean, digits), 
			paste("(N = ",Nobs, ifelse(is.null(Nrm), "", paste(", ", Nrm, " observations removed due to missing data", sep="")), ")", sep=""), 
			"\n\nExperimental Design:", balanced," |  Method:", x$EstMethod)
	
	if(NegVCmsg != "" && x$Type != "Linear Model")                # there were VCs set to zero
		cat(" |", NegVCmsg, "| adapted MS used for total DF")
	cat("\n\n")
}


#'Inferential Statistics for VCA-Results
#'
#'Function \code{VCAinference} constructs one- and two-sided confidence intervals, and performs Chi-Squared tests for total
#'and error variance against claimed values for 'VCA' objects.
#'
#'This function computes confidence intervals (CI) for variance components (VC), standard deviations (SD)
#'and coefficients of variation (CV). VCs 'total' and 'error' can be tested against claimed values specifying parameters
#''total.claim' and 'error.claim'. One can also specify claim-values in terms of SD or CV (see \code{claim.type}).\cr
#'Confidence intervals for VCs are constructed either following the same rules as in SAS 9.2 PROC MIXED with option 'method=type1'
#'(ci.method="sas") or using Satterthwaite methodology throughout (ci.method="satterthwaite"). In the former approach
#'for VC total and error, which are constrained to be \eqn{>= 0}, CIs are based on the Chi-Squared distribution. Degrees of freedom
#'(DF) for total variance are approximated using the Satterthwaite approximation (which is not available in either SAS procedure).
#'For all other VCs, the CI is \eqn{[sigma^2-QNorm(alpha/2)*SE(sigma^2); sigma^2+QNorm(1-alpha/2)*SE(sigma^2)]}, where QNorm(x) indicates the x-quantile of 
#'the standard normal distribution. The second method approximates DFs for all VCs using the Satterthwaite approximation and CIs are
#'based on the corresponding Chi-Squared distribution for all VCs (see examples). 
#'Note that in the computation of the covariance-matrix of the VCs, the estimated VCs will be used. If these are requested to be set to 0 
#'(\code{NegVC=FALSE} in \code{\link{anovaVCA}}), the result might not be conformable with theory given in the first reference. 
#'The validity of this implementation was checked against SAS 9.2 PROC MIXED (method=type1), where VCs are not constrained to be >= 0. 
#'The sampling variances for VCs are obtained assuming normality throughout based on \eqn{Var(\sigma^{2} = C^{-1} * Var(m_{SS} * (C^{-1})^{T}))}{Var(sigma^2) = Ci * Var(SS) * Ci'}, 
#'where \eqn{C^{-1}}{Ci} is the inverse of the coefficient matrix equating observed Sum of Squares (SS)
#'to their expected values, and \eqn{(C^{-1})^{T}}{Ci'} indicating the transpose of \eqn{C^{-1}}{Ci} (see Searle et al. 1992, pg. 176).
#'
#'An input \code{VCA}-object can be in one of three states:\cr 
#'\cr
#'State (1) corresponds to the situation, where all VC > 0.\cr 
#'State (2) corresponds to the situation, where at least one VC < 0.\cr
#'State (3) corresponds to situations, where negative VC estimates occured but were set to 0, i.e. \code{NegVC=FALSE} - the Default.\cr
#'
#'State (2) occurs when parameter \code{NegVC} was set to TRUE in \code{\link{anovaVCA}}, state (3) represents the default-setting in 
#'function \code{\link{anovaVCA}}. If a \code{VCA}-object is in state (1), parameter \code{excludeNeg} has no effect (there are no negative VCs),
#'only parameter \code{constrainCI} is evaluated. For \code{VCA}-objects in state(2), \code{constrainCI} has no effect, because constraining
#'CIs for unconstrained VCs makes no sense. State (3) forces parameter \code{constrainCI} to be set to TRUE and one can only choose whether to
#'exclude CIs of negative VC estimates or not. Whenever VCs have to be constrained, it is straight forward to apply constraining also to any
#'CI. Note that situations outlined above only occur when parameter \code{VarVC} is set to TRUE, which causes estimation of the covariance-matrix
#'of variance components. The default is only to compute and report CIs for total and error variance, which cannot become negative.
#'
#'
#'@param obj				(object) of class 'VCA' or, alternatively, a list of 'VCA' objects, where all other argument can be 
#'							specified as vectors, where the i-th vector element applies to the i-th element of 'obj' (see examples) 
#'@param alpha				(numeric) value specifying the significance level for \eqn{100*(1-alpha)}\% confidence intervals.
#'@param total.claim		(numeric) value specifying the claim-value for the Chi-Squared test for the total variance (SD or CV, see \code{claim.type}).
#'@param error.claim		(numeric) value specifying the claim-value for the Chi-Squared test for the error variance (SD or CV, see \code{claim.type}).
#'@param claim.type			(character) one of "VC", "SD", "CV" specifying how claim-values have to be interpreted:\cr
#'							"VC" (Default) = claim-value(s) specified in terms of variance(s),\cr
#'							"SD" = claim-values specified in terms of standard deviations (SD),\cr
#'							"CV" = claim-values specified in terms of coefficient(s) of variation (CV)
#'							and are specified as percentages.\cr
#'							If set to "SD" or "CV", claim-values will be converted to variances before applying the Chi-Squared test (see examples).
#'@param VarVC				(logical) TRUE = if element "Matrices" exists (see \code{\link{anovaVCA}}), the covariance
#'							matrix of the estimated VCs will be computed (see \code{\link{vcovVC}}, which is used in CIs for 
#'							intermediate VCs if 'method.ci="sas"'. 
#'							Note, this might take very long for larger datasets, since there are many matrix operations involved. 
#'							FALSE (Default) = computing covariance matrix of VCs is omitted, as well as CIs for intermediate VCs.
#'@param excludeNeg			(logical) TRUE = confidence intervals of negative variance estimates will not be reported. \cr
#'							FALSE = confidence intervals for all VCs will be reported including those with negative VCs.\cr
#'							See the details section for a thorough explanation.
#'@param constrainCI		(logical) TRUE = CI-limits for all variance components are constrained to be >= 0.\cr
#'							FALSE = unconstrained CIs with potentially negative CI-limits will be reported.\cr
#'							which will preserve the original width of CIs.
#'							See the details section for a thorough explanation.
#'@param ci.method			(character) string or abbreviation specifying which approach to use for computing confidence intervals of variance components (VC).
#'							"sas" (default) uses Chi-Squared based CIs for total and error and normal approximation for all other VCs (Wald-limits, option "NOBOUND"
#'							in SAS PROC MIXED); "satterthwaite" will approximate DFs for each VC using the Satterthwaite approach (see \code{\link{SattDF}} for models
#'							fitted by ANOVA) and all Cis are based on the Chi-Squared distribution. This approach is conservative but avoids negative values for the lower bounds. 
#'@param quiet				(logical) TRUE = will suppress any warning, which will be issued otherwise 
#'@note						Original CIs will always be available independent of parameter-settings of \code{excludeNeg} and
#'							\code{constrainCI}. Original CIs are stored in attribute "CIoriginal" of the returned 'VCAinference'-object, e.g.
#'							'attr(obj$ConfInt$SD$OneSided, "CIoriginal")' or 'attr(obj$ConfInt$CV$TwoSided, "CIoriginal")'.
#'
#'@return  (VCAinference) object, a list with elements: \cr
#'\itemize{
#'\item{ChiSqTest}{(data.frame) with results of the Chi-Squared test}
#'\item{ConfInt}{(list) with elements \code{VC}, \code{SD}, \code{CV}, all lists themselves containing (data.frame) objects \code{OneSided} and \code{TwoSided}}
#'\item{VCAobj}{(VCA) object specified as input, if \code{VarVC=TRUE}, the 'aov.tab' element will have an extra column "Var(VC)" storing variances of VC-estimates"}
#'} 
#'
#'@author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#'
#'@references 
#'
#'Searle, S.R, Casella, G., McCulloch, C.E. (1992), Variance Components., Wiley New York
#'
#'Burdick, R., Graybill, F. (1992), Confidence Intervals on Variance Components. Marcel Dekker, Inc.
#'
#'Satterthwaite, F.E. (1946), An Approximate Distribution of Estimates of Variance Components., 
#'Biometrics Bulletin 2, 110-114
#'
#'@seealso \code{\link{print.VCAinference}}, \code{\link{anovaVCA}}
#'
#'@examples 
#'\dontrun{
#'
#'# load data (CLSI EP05-A2 Within-Lab Precision Experiment) 
#'data(dataEP05A2_1)
#'
#'# perform (V)variance (C)component (A)nalysis (also compute A-matrices)
#'res <- anovaVCA(y~day/run, dataEP05A2_1)
#'
#'# get confidence intervals for total and error (VC, SD, CV)
#'VCAinference(res)
#'
#'# additionally request CIs for all other VCs; default is to constrain 
#'# CI-limits to be >= 0
#'# first solve MME
#'res <- solveMME(res)
#'VCAinference(res, VarVC=TRUE)
#'
#'# now using Satterthwaite methodology for CIs
#'VCAinference(res, VarVC=TRUE, ci.method="satt")
#'
#'# request unconstrained CIs
#'VCAinference(res, VarVC=TRUE, constrainCI=FALSE)
#'
#'# additionally request Chi-Squared Tests of total and error, default 
#'# is that claim values are specified as variances (claim.type="VC")
#'VCAinference(res, total.claim=4.5, error.claim=3.5)
#'
#'# perform Chi-Squared Tests, where claim-values are given as SD, 
#'# compare p-values to former example
#'VCAinference(res, total.claim=sqrt(4.5), error.claim=sqrt(3.5), claim.type="SD")
#'
#'# now using Satterthwaite methodology for CIs
#'VCAinference(res, total.claim=sqrt(4.5), error.claim=sqrt(3.5), 
#'claim.type="SD", ci.method="satt")
#'
#'# now add random error to example data forcing the ANOVA-estimate of the 
#'# day-variance to be negative
#'set.seed(121)
#'tmpData <- dataEP05A2_1
#'tmpData$y <- tmpData$y + rnorm(80,,3)
#'res2 <- anovaVCA(y~day/run, tmpData)
#'
#'# call 'VCAinference' with default settings
#'VCAinference(res2)
#'
#'# extract components of the returned 'VCAinference' object
#'inf <- VCAinference(res2, total.claim=12)
#'inf$ConfInt$VC$OneSided			# one-sided CIs for variance components
#'inf$ConfInt$VC$TwoSided			# two-sided CI for variance components
#'inf$ChiSqTest
#'
#'# request CIs for all VCs, default is to exclude CIs of negative VCs (excludeNeg=TRUE) 
#'# solve MMEs first (or set MME=TRUE when calling anovaVCA)
#'res2 <- solveMME(res2)
#'VCAinference(res2, VarVC=TRUE)
#'
#'# request CIs for all VCs, including those for negative VCs, note that all CI-limits 
#'# are constrained to be >= 0
#'VCAinference(res2, VarVC=TRUE, excludeNeg=FALSE)
#'
#'# request unconstrained CIs for all VCs, including those for negative VCS
#'# one has to re-fit the model allowing the VCs to be negative
#'res3 <- anovaVCA(y~day/run, tmpData, NegVC=TRUE, MME=TRUE)
#'VCAinference(res3, VarVC=TRUE, excludeNeg=FALSE, constrainCI=FALSE)
#'
#'### use the numerical example from the CLSI EP05-A2 guideline (p.25)
#'data(Glucose,package="VCA")
#'res.ex <- anovaVCA(result~day/run, Glucose)
#'
#'### also perform Chi-Squared tests
#'### Note: in guideline claimed SD-values are used, here, claimed variances are used
#'VCAinference(res.ex, total.claim=3.4^2, error.claim=2.5^2)
#'
#'
#'# load another example dataset and extract the "sample_1" subset
#'data(VCAdata1)
#'sample1 <- VCAdata1[which(VCAdata1$sample==1),]
#'
#'# generate an additional factor variable and random errors according to its levels
#'sample1$device <- gl(3,28,252)                                      
#'set.seed(505)
#'sample1$y <- sample1$y + rep(rep(rnorm(3,,.25), c(28,28,28)),3)     
#'
#'# fit a crossed-nested design with main factors 'lot' and 'device' 
#'# and nested factors 'day' and 'run' nested below, also request A-matrices 
#'res1 <- anovaVCA(y~(lot+device)/day/run, sample1) 
#'
#'# get confidence intervals, covariance-matrix of VCs, ..., 
#'# explicitly request the covariance-matrix of variance components
#'# solve MMEs first
#'res1 <- solveMME(res1)
#'inf1 <- VCAinference(res1, VarVC=TRUE, constrainCI=FALSE)
#'inf1
#'
#'# print numerical values with more digits
#'print(inf1, digit=12)
#'
#'# print only parts of the 'VCAinference' object (see \code{\link{print.VCAinference}})
#'print(inf1, digit=12, what=c("VCA", "VC"))
#'
#'# extract complete covariance matrix of variance components 
#'# (main diagonal is part of standard output -> "Var(VC"))
#'VarCovVC <- vcovVC(inf1$VCAobj)
#'round(VarCovVC, 12)
#'
#'# use by-processing and specific argument-values for each level of the by-variable
#'data(VCAdata1)
#'fit.all <- anovaVCA(y~(device+lot)/day/run, VCAdata1, by="sample", NegVC=TRUE)
#'inf.all <- VCAinference(fit.all, total.claim=c(.1,.75,.8,1,.5,.5,2.5,20,.1,1))
#'print.VCAinference(inf.all, what="VC")
#'}

VCAinference <- function(obj, alpha=.05, total.claim=NA, error.claim=NA, claim.type="VC", 
		VarVC=FALSE, excludeNeg=TRUE, constrainCI=TRUE, ci.method="sas",
		quiet=FALSE)
{
	Call <- match.call()
	
	if(is.list(obj) && !is(obj, "VCA"))
	{
		if(!all(sapply(obj, class) == "VCA"))
			stop("Only lists of 'VCA' object are accepted!")
		
		obj.len <- length(obj)
		
		if(!"msgEnv" %in% ls(.GlobalEnv))
			msgEnv <<- new.env(parent=emptyenv())
		
		assign("VCAinference.obj.is.list", TRUE, envir=msgEnv)			# indicate that a list-type object was passed intially
		
		res <- mapply(	FUN=VCAinference, obj=obj, alpha=alpha, 
				total.claim=total.claim, error.claim=error.claim,
				claim.type=claim.type, VarVC=VarVC, excludeNeg=excludeNeg,
				constrainCI=constrainCI, ci.method=ci.method, 
				SIMPLIFY=FALSE)
		names(res) <- names(obj)
		
		if(obj.len == 1)			# mapply returns a list of length 2 in case that length(obj) was equal to 1
			res <- res[1]
		
		rm("VCAinference.obj.is.list", envir=msgEnv)
		
		return(res)
	}	
	
	stopifnot(class(obj) == "VCA")
	ci.method <- match.arg(ci.method, c("sas", "satterthwaite"))
	
	MM <- obj$Type == "Mixed Model"					# only exists for mixed models
	if(MM)
		VCs   <- obj$Matrices$VCall
	else
		VCs <- obj$aov.tab[-1,"VC"]		 			# for random models directly from ANOVA-table
	
	claim.type <- match.arg(claim.type, c("VC", "SD", "CV"))
	
	Mean <- obj$Mean
	Nvc  <- nrow(obj)
	EstMethod <- obj$EstMethod
	
	if( all(obj$aov.tab[,"VC"] > 0) )
		VCstate <- 1                               
	else if( any(obj$aov.tab[,"VC"] < 0) )
		VCstate <- 2
	else
	{
		if(obj$NegVCmsg != "")
			VCstate <- 3                            # negative VCs set to zero                               
		else                                                                
			VCstate <- 1                            # zero-VC estimated as such                                      
	}
	
	if( !is.null(obj$Matrices) && VarVC )
	{
		if(is.null(obj$VarCov))						# solve mixed model equations first
		{
			# obj  <- solveMME(obj)
			
			Lmat  <- obj$Matrices                               # different matrices needed for VCA
			
			VCvar <- vcovVC(obj, method=obj$VarVC.method) 		# get variance-covariance matrix of VCs (p.176); do not pass total VC
			
			NegVCmsg    <- obj$NegVCmsg
			VCoriginal  <- obj$VCoriginal
			Nobs        <- obj$Nobs
			Nrm         <- obj$Nrm
			balanced    <- obj$balanced
			
			obj$aov.tab <- cbind(obj$aov.tab, "Var(VC)"=c(NA, diag(VCvar)))  
			
			class(obj) <- "VCA"
			obj$NegVCmsg   <- NegVCmsg
			obj$VCoriginal <- VCoriginal
			obj$VarCov     <- VCvar                                    # store variance-covariance matrix of variance components
			obj$Mean       <- Mean
			obj$Nobs       <- Nobs
			obj$Nrm        <- Nrm
			obj$balanced   <- balanced
			
			nam0 <- deparse(Call$obj)
			nam1 <- sub("\\[.*", "", nam0)
			
			if(length(nam1) == 1 && nam1 %in% names(as.list(.GlobalEnv)))		# obj is not function call
			{
				expr <- paste(nam0, "<<- obj")						# update object missing MME results
				eval(parse(text=expr))	
			}
			else	# warning only if not called on list of VCA-objects
			{
				if( !"VCAinference.obj.is.list" %in% names(as.list(msgEnv)) && !quiet )
					message("Mixed model equations solved locally. Results could not be assigned to object!")
			}
		} else
			obj$aov.tab <- cbind(obj$aov.tab, "Var(VC)"=c(NA, diag(obj$VarCov)))  
	}
	
	if(!is.na(total.claim) && nrow(obj$aov.tab) == 1)                               # if error is the only VC no total variance exists (or is equal)
		total.claim <- NA
	
	if(nrow(obj$aov.tab) == 1)
		Nvc <- 1
	else
		Cind <- which(rownames(obj$aov.tab) %in% c("total", "error"))
	
	if(!is.na(total.claim))
	{
		if(is.numeric(total.claim))
		{
			if( (is.na(total.claim) || total.claim <= 0 ) && !quiet)
				warning("Parameter 'total.claim' is not correctly specified! Chi-Squared test for total precision is omitted!")
		}
		else
		{
			if(!quiet)
				warning("Parameter 'total.claim' is not correctly specified! Chi-Squared test for total precision is omitted!")
		}
	}
	if(!is.na(error.claim))
	{
		if(is.numeric(error.claim))
		{
			if( (is.na(error.claim) || error.claim <= 0) && !quiet)
				warning("Parameter 'error.claim' is not correctly specified! Chi-Squared test for error precision (repeatability) is omitted!")
		}
		else
		{
			if(!quiet)
				warning("Parameter 'error.claim' is not correctly specified! Chi-Squared test for error precision (repeatability) is omitted!")
		}
	}
	
	nam <- rownames(obj$aov.tab)
	Nvc <- length(nam)                                  # number of variance components including total variance
	
	# Chi-Squared tests
	
	CStest <- data.frame(Name=nam, "Claim"=rep(NA,Nvc))
	CStest$"ChiSq value" <- rep(NA, Nvc)
	CStest$"Pr (>ChiSq)" <- rep(NA, Nvc)
	
	if(!is.na(total.claim))
	{
		if(claim.type == "SD")
			total.claim <- total.claim^2                # now the Chi-Squared Test for variances can be used
		
		if(claim.type =="CV")
			total.claim <- (total.claim * Mean/100)^2   # re-caluculate variance from CV
		
		CStest$"Claim"[1] <- total.claim
		
		CStest$"ChiSq value"[1] <- obj$aov.tab[1, "VC"]*obj$aov.tab[1, "DF"]/total.claim 
		
		CStest$"Pr (>ChiSq)"[1] <- pchisq(q=CStest$"ChiSq value"[1], df=obj$aov.tab[1, "DF"], lower.tail=TRUE)
		
		if(claim.type == "SD")
			CStest$"Claim"[1] <- sqrt(total.claim)
		
		if(claim.type == "CV")
			CStest$"Claim"[1] <- sqrt(total.claim)*100/Mean
	}
	
	if(!is.na(error.claim))
	{
		if(claim.type == "SD")
			error.claim <- error.claim^2                # now the Chi-Squared Test for variances can be used
		
		if(claim.type =="CV")
			error.claim <- (error.claim * Mean/100)^2   # re-caluculate variance from CV
		
		CStest$"Claim"[Nvc]       <- error.claim
		CStest$"ChiSq value"[Nvc] <- obj$aov.tab[Nvc, "VC"]*obj$aov.tab[Nvc, "DF"]/error.claim 
		CStest$"Pr (>ChiSq)"[Nvc] <- pchisq(q=CStest$"ChiSq value"[Nvc], df=obj$aov.tab[Nvc, "DF"], lower.tail=TRUE)
		
		if(claim.type == "SD")
			CStest$"Claim"[Nvc] <- sqrt(error.claim)
		
		if(claim.type == "CV")
			CStest$"Claim"[Nvc] <- sqrt(error.claim)*100/Mean
	}
	
	rownames(CStest) <- CStest$Name
	
	# CIs on VCs, SDs and CVs
	
	indCS <- which(rownames(obj$aov.tab) %in% c("total", "error"))				# Chi-Squred dist VCs
	indN  <- which(!rownames(obj$aov.tab) %in% c("total", "error"))				# Normal dist VCs
	
	#####################
	### two-sided CIs ###
	
	CI_VC <- data.frame(Name=nam)
	CI_VC$LCL <- numeric(nrow(obj$aov.tab))
	CI_VC$UCL <- numeric(nrow(obj$aov.tab))
	
	# Chi-Squared CIs
	
	lower.qchisq <- qchisq(p=1-alpha/2, df=obj$aov.tab[indCS, "DF"])
	upper.qchisq <- qchisq(p=alpha/2, df=obj$aov.tab[indCS, "DF"])
	CI_VC$LCL[indCS] <- obj$aov.tab[indCS, "VC"]*obj$aov.tab[indCS, "DF"]/lower.qchisq
	CI_VC$UCL[indCS] <- obj$aov.tab[indCS, "VC"]*obj$aov.tab[indCS, "DF"]/upper.qchisq
	
	# CIs based on Normal-Distribution (ci.method="sas") or Chi-Squared (ci.method="satterthwaite")
	
	if(ci.method == "sas" && "Var(VC)" %in% colnames(obj$aov.tab))		# variance-covariance matrix of VCs required
	{
		CI_VC$LCL[indN] <- obj$aov.tab[indN, "VC"]+qnorm(alpha/2)*sqrt(obj$aov.tab[indN, "Var(VC)"])
		CI_VC$UCL[indN] <- obj$aov.tab[indN, "VC"]+qnorm(1-alpha/2)*sqrt(obj$aov.tab[indN, "Var(VC)"])
		DFs <- NULL
	}
	else if(ci.method == "satterthwaite")
	{
		if(obj$EstMethod == "ANOVA")
		{
			rind	<- obj$Matrices$rf.ind						# indices of random terms in the original ANOVA-table
			aov.tab <- if(obj$Type == "Random Model")			# original ANOVA table
						obj$aov.tab[-1,]					# remove row for total 
					else	
						obj$aov.org
			MS  	<- aov.tab[rind, "MS"]
			Ci  	<- getMat(obj, "Ci.MS")[rind, rind]
			DF  	<- aov.tab[rind, "DF"]
			sDF 	<- SattDF(MS, Ci, DF, type="individual") 	# satterthwaite DFs
			ind 	<- 1:(length(sDF)-1)						# without row for error 
			DFs 	<- c(obj$aov.tab[1,"DF"], sDF)				# total DF
		}
		else													# fitted by REML, only Satterthwaite DF exist
		{
			DFs <- obj$aov.tab[,"DF"]
			sDF <- DFs[-1]
			ind <- 1:(length(sDF)-1)
		}
		
		CI_VC$DF <- DFs										# add Satterthwaite DFs
		CI_VC <- CI_VC[,c(1,4,2,3)]
		lower.qchisq <- qchisq(p=1-alpha/2, df=sDF[ind])
		upper.qchisq <- qchisq(p=alpha/2, df=sDF[ind])
		CI_VC$LCL[ind+1] <- obj$aov.tab[ind+1, "VC"]*sDF[ind]/lower.qchisq		# lower limit all but total and error
		CI_VC$UCL[ind+1] <- obj$aov.tab[ind+1, "VC"]*sDF[ind]/upper.qchisq		# upper limit  - " -
	}
	else
	{
		CI_VC$LCL[indN] <- NA
		CI_VC$UCL[indN] <- NA
		DFs <- NULL
	}
	
	rownames(CI_VC) <- CI_VC$Name
	attr(CI_VC, "CIoriginal") <- CI_VC                      # keep original CI limits as attribute
	
	if(VCstate == 1)                                        # excludeNeg not evaluated, since all VCs positive
	{
		if( constrainCI && any(CI_VC$LCL < 0, na.rm=TRUE))  # any negative CI-limits?
		{
			LCLind <- which(CI_VC$LCL < 0)
			UCLind <- which(CI_VC$UCL < 0)                  # keep info about constrained LCLs
			
			attr(CI_VC, "LCLconstrained") <- LCLind
			
			CI_VC$LCL[LCLind] <- 0                          # set negative LCL to 0
			
			if(length(UCLind) > 0)                          # set negative UCL to 0 
			{
				CI_VC$UCL[UCLind] <- 0
				attr(CI_VC, "UCLconstrained") <- UCLind     # keep info about constrained UCLs
			}
		}
	}
	if(VCstate == 2)                                        # constrainCI not evaluated since negative VCs explicitly allowed
	{
		if( excludeNeg )                                    # CIs for negative VC are excluded
		{    
			NegInd <- which(obj$aov.tab[,"VC"] < 0)
			CI_VC$LCL[NegInd] <- NA                         # there has to be at least one LCL < 0
			CI_VC$UCL[NegInd] <- NA
		}
	}
	if(VCstate == 3)
	{
		if( excludeNeg )
		{
			NegInd <- which(obj$VCoriginal < 0)+1  			# which VCs were set to 0, i.e. were originally < 0
			CI_VC$LCL[NegInd] <- NA                         # there has to be at least one LCL < 0
			CI_VC$UCL[NegInd] <- NA
		}
		
		if(any(CI_VC$LCL < 0, na.rm=TRUE))                  # constrainCI implicitly set to 0, since negative VCs set to zero automatically results in constraining CIs
		{
			LCLind <- which(CI_VC$LCL < 0)
			UCLind <- which(CI_VC$UCL < 0)                  # keep info about constrained LCLs
			
			attr(CI_VC, "LCLconstrained") <- LCLind
			
			CI_VC$LCL[LCLind] <- 0                          # set negative LCL to 0
			
			if(length(UCLind) > 0)                          # set negative UCL to 0 
			{
				CI_VC$UCL[UCLind] <- 0
				attr(CI_VC, "UCLconstrained") <- UCLind     # keep info about constrained UCLs
			}
		}
	}
	
	# test whether point-estimates within CI or not
	CIwarning <- FALSE
	CIout <- NULL
	CIind <- which(obj$aov.tab[,"VC"] < CI_VC$LCL | obj$aov.tab[,"VC"] > CI_VC$UCL)
	
	if(length(CIind) > 0)
	{
		CI_VC$LCL[CIind] <- CI_VC$UCL[CIind] <- NA									# set them NA
		CIwarning <- TRUE
		CIout <- c(CIout, rownames(obj$aov.tab)[CIind])
	}
	
	
	# SD and CV computed from VC      NOTE: sqrt o.k.??????????
	
	CI_SD <- data.frame(Name=nam)
	CI_SD$DF <- DFs							# add Satterthwaite DFs if available
	Sign <- sign(CI_VC$LCL)
	CI_SD$LCL <- sqrt(abs(CI_VC$LCL))*Sign
	Sign <- sign(CI_VC$UCL)
	CI_SD$UCL <- sqrt(abs(CI_VC$UCL))*Sign
	rownames(CI_SD) <- CI_SD$Name
	
	if(any(is.na(obj$aov.tab[,"SD"])))
		CI_SD[which(is.na(obj$aov.tab[,"SD"])),2:ncol(CI_SD)] <- NA 
	
	CI_CV <- data.frame(Name=nam)
	CI_CV$DF <- DFs							# add Satterthwaite DFs if available
	CI_CV$LCL <- CI_SD$LCL * 100/Mean
	CI_CV$UCL <- CI_SD$UCL * 100/Mean
	rownames(CI_CV) <- CI_CV$Name
	
	CI_two_sided <- list(CI_VC=CI_VC, CI_SD=CI_SD, CI_CV=CI_CV)
	
	# one-sided CIs
	
	CI_VC <- data.frame(Name=nam)
	CI_VC$DF <- DFs							# add Satterthwaite DFs if available
	CI_VC$LCL <- numeric(nrow(obj$aov.tab))
	CI_VC$UCL <- numeric(nrow(obj$aov.tab))
	
	# CIs based on Chi-Squared
	
	upper.qchisq <- qchisq(p=alpha,   df=obj$aov.tab[indCS, "DF"])
	lower.qchisq <- qchisq(p=1-alpha, df=obj$aov.tab[indCS, "DF"])
	CI_VC$LCL[indCS] <- obj$aov.tab[indCS, "VC"]*obj$aov.tab[indCS, "DF"]/lower.qchisq
	CI_VC$UCL[indCS] <- obj$aov.tab[indCS, "VC"]*obj$aov.tab[indCS, "DF"]/upper.qchisq
	rownames(CI_VC) <- CI_VC$Name
	
	# CIs based on Normal Distribution or Chi-Squared 
	
	if(ci.method == "sas" && "Var(VC)" %in% colnames(obj$aov.tab))	
	{
		CI_VC$UCL[indN] <- obj$aov.tab[indN, "VC"]+qnorm(1-alpha)*sqrt(obj$aov.tab[indN, "Var(VC)"])
		CI_VC$LCL[indN] <- obj$aov.tab[indN, "VC"]+qnorm(  alpha)*sqrt(obj$aov.tab[indN, "Var(VC)"])
	}
	else if(ci.method == "satterthwaite")
	{
		lower.qchisq <- qchisq(p=1-alpha, df=sDF[ind])
		upper.qchisq <- qchisq(p=alpha, df=sDF[ind])
		CI_VC$LCL[ind+1] <- obj$aov.tab[ind+1, "VC"]*sDF[ind]/lower.qchisq
		CI_VC$UCL[ind+1] <- obj$aov.tab[ind+1, "VC"]*sDF[ind]/upper.qchisq
	}
	else
	{
		CI_VC$UCL[indN] <- NA
		CI_VC$LCL[indN] <- NA
	}
	
	attr(CI_VC, "CIoriginal") <- CI_VC                      # keep original CI limits as attribute
	
	if(VCstate == 1)                                        # excludeNeg not evaluated, since all VCs positive
	{
		if( constrainCI && any(CI_VC$LCL < 0, na.rm=TRUE))  # any negative CI-limits?
		{
			LCLind <- which(CI_VC$LCL < 0)
			UCLind <- which(CI_VC$UCL < 0)                  # keep info about constrained LCLs
			CI_VC$LCL[LCLind] <- 0 
			attr(CI_VC, "LCLconstrained") <- LCLind
			
			if(length(UCLind) > 0)                          # set negative UCL to 0 
			{
				CI_VC$UCL[UCLind] <- 0
				attr(CI_VC, "UCLconstrained") <- UCLind     # keep info about constrained UCLs
			}
		}
	}
	if(VCstate == 2)                                        # constrainCI not evaluated since negative VCs explicitly allowed
	{
		if( excludeNeg )                                    # CIs for negative VC are excluded
		{    
			NegInd <- which(obj$aov.tab[,"VC"] < 0)
			
			CI_VC$UCL[NegInd] <- NA                         # there has to be at least one UCL < 0
			CI_VC$LCL[NegInd] <- NA                         # there has to be at least one UCL < 0
		}
	}
	if(VCstate == 3)
	{
		
		if( excludeNeg )
		{
			NegInd <- which(obj$VCoriginal < 0) + 1  		# which VCs were set to 0, i.e. were originally < 0
			CI_VC$UCL[NegInd] <- NA                         # there has to be at least one UCL < 0
			CI_VC$LCL[NegInd] <- NA
		}
		
		if(any(CI_VC$LCL < 0, na.rm=TRUE))
		{
			LCLind <- which(CI_VC$LCL < 0)              # keep info about constrained LCLs
			CI_VC$LCL[LCLind] <- 0                      # set negative LCL to 0
			attr(CI_VC, "LCLconstrained") <- LCLind
		}
		
#            if(any(CI_VC$UCL < 0, na.rm=TRUE))
#            {
#                UCLind <- which(CI_VC$UCL < 0)              # keep info about constrained UCLs
#                
#                attr(CI_VC, "UCLconstrained") <- UCLind
#                
#                CI_VC$UCL[UCLind] <- 0                      # set negative UCL to 0
#            }
		
	}
	
	# test whether point-estimates within CI or not
	CIindLCL <- which(obj$aov.tab[,"VC"] < CI_VC$LCL)
	CIindUCL <- which(obj$aov.tab[,"VC"] > CI_VC$UCL)
	
	if(length(CIindLCL) > 0)
	{
		CI_VC$LCL[CIindLCL] <- NA								# set LCL to NA
		CIwarning <- TRUE
		CIout <- c(CIout, rownames(obj$aov.tab)[CIindLCL])
	}
	if(length(CIindUCL) > 0)
	{
		CI_VC$UCL[CIindUCL] <- NA										# set LCL to NA
		CIwarning <- TRUE
		CIout <- c(CIout, rownames(obj$aov.tab)[CIindUCL])
	}
	
	if(CIwarning && !quiet)
	{	
		warning("Point estimate(s) of VC(s) ", paste(paste("'", unique(CIout), "'", sep=""), sep=", ")," found outside of confidence interval, CI was set to 'NA'!")
	}
	# SD and CV computed from VC
	
	CI_SD <- data.frame(Name=nam)
	CI_SD$DF <- DFs											# add Satterthwaite DFs if available
	SignLCL  <- sign(CI_VC$LCL)
	CI_SD$LCL <- sqrt(abs(CI_VC$LCL)) * SignLCL
	SignUCL  <- sign(CI_VC$UCL)
	CI_SD$UCL <- sqrt(abs(CI_VC$UCL)) * SignUCL
	rownames(CI_SD) <- CI_SD$Name
	
	if(any(is.na(obj$aov.tab[,"SD"])))
		CI_SD[which(is.na(obj$aov.tab[,"SD"])),2:ncol(CI_SD)] <- NA 
	
	CI_CV <- data.frame(Name=nam)
	CI_CV$DF <- DFs											# add Satterthwaite DFs if available
	CI_CV$LCL <- CI_SD$LCL * 100/Mean
	CI_CV$UCL <- CI_SD$UCL * 100/Mean
	rownames(CI_CV) <- CI_CV$Name
	
	CI_one_sided <- list(CI_VC=CI_VC, CI_SD=CI_SD, CI_CV=CI_CV)
	
	CIs <- list(VC=list(OneSided=CI_one_sided$CI_VC, TwoSided=CI_two_sided$CI_VC),
			SD=list(OneSided=CI_one_sided$CI_SD, TwoSided=CI_two_sided$CI_SD),
			CV=list(OneSided=CI_one_sided$CI_CV, TwoSided=CI_two_sided$CI_CV))
	
	result <- list(ChiSqTest=CStest, ConfInt=CIs, VCAobj=obj, alpha=alpha)
	class(result) <- "VCAinference"
	attr(result, "EstMethod")   <- EstMethod
	attr(result, "excludeNeg")  <- excludeNeg
	attr(result, "constrainCI") <- constrainCI
	attr(result, "claim.type")  <- claim.type
	attr(result, "ci.method")   <- ci.method
	return(result)
}


#'Standard Print Method for Objects of Class 'VCAinference'
#'
#'Prints the list-type 'VCAinference'-object as tabulated output. 
#'
#'Formats the list-type objects of class 'VCAinference' for a more comprehensive
#'presentation of results, which are easier to grasp. The default is to show the complete
#'object (VCA ANOVA-table, VC-, SD-, and CV-CIs). Using parameter 'what' allows to
#'restrict the printed output to certain parts. Print-function invisibly returns a matrix
#'or a list of matrices, depending on the values of 'what', i.e. it can be used as for
#'packing the inference-information in one or multiple matrix-objects and extracting it/them.
#'
#'@param x         (VCAinference) object 
#'@param digits    (integer) number of decimal digits.
#'@param what      (character) one of "all", "VC", "SD", "CV", "VCA" specifying which part of the 'VCA'-object is to be printed.
#'@param ...       additional arguments to be passed to or from methods.
#'
#'@return invisibly returns sub-elements of 'x' specified via 'what'
#'
#'@author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#'
#'@method print VCAinference
#'
#'@seealso \code{\link{VCAinference}}, \code{\link{anovaVCA}}
#'
#'@examples
#'\dontrun{
#'# load data (CLSI EP05-A2 Within-Lab Precision Experiment) 
#'data(dataEP05A2_1)
#'
#'# perform ANOVA-estimation of variance components for a nested design
#'res <- anovaVCA(y~day/run, Data=dataEP05A2_1)
#'res
#'inf <- VCAinference(res)
#'inf
#'
#'# show certain parts and extract them invisibly
#'CVmat <- print(inf, what="CV")
#'CVmat
#'
#'# show numerical values with more digits
#'print(inf, digit=12)
#'}

print.VCAinference <- function(x, digits=4L, what=c("all", "VC", "SD", "CV", "VCA"), ...)
{
	if(is.list(x) && !is(x, "VCAinference"))
	{
		if(!all(sapply(x, class) == "VCAinference"))
			stop("Only lists of 'VCAinference' objects can printed!")
		
		nam <- names(x)
		lst <- list()
		
		for(i in 1:length(x))
		{
			print(nam[i])
			lst[[i]] <- print(x[[i]], digits=digits, what=what)
		}
		return()			# leave function now
	}
	
	stopifnot(class(x) == "VCAinference") 
	claim.type <- attr(x, "claim.type")
	VCAobj <- x$VCAobj
	
	ret <- list()
	
	what <- match.arg(tolower(what), tolower(c("all", "VC", "SD", "CV", "VCA")), several.ok=TRUE)
	
#	MM <- !is.null(attr(x$VCAobj, "FixedEffects"))
	
	if(VCAobj$Type == "Mixed Model")
		cat("\n\n\nInference from Mixed Model Fit\n------------------------------\n\n")
	else if(VCAobj$Type == "Linear Model")
		cat("\n\n\nInference from Linear Model Fit:\n--------------------------------\n\n")
	else
		cat("\n\n\nInference from (V)ariance (C)omponent (A)nalysis\n------------------------------------------------\n\n")
	
	if( any( c("all", "vca") %in% what) )
	{
		if(VCAobj$Type == "Linear Model")
			cat("> ANOVA Table:\n--------------\n\n")
		else
			cat("> VCA Result:\n-------------\n\n")
		print(VCAobj, digits=digits, skipHeading=TRUE)
		cat("\n")
	}    
	
	if( any( c("all", "vc") %in% what) )
	{
		cat("> VC:\n-----\n")
		VC <- as.data.frame(as.matrix(VCAobj$aov.tab)[,"VC", drop=FALSE])
		colnames(VC) <- "Estimate"
		
		if(any(!is.na(x$ChiSqTest[,"Claim"])) && claim.type == "VC")
		{
			VC$Claim <- as.numeric(rep(NA, nrow(VC)))
			VC$"ChiSq" <- as.numeric(rep(NA, nrow(VC)))
			VC$"Pr(>ChiSq)" <- as.numeric(rep(NA, nrow(VC)))
			VC[as.character(x$ChiSqTest[,"Name"]), c("Claim", "ChiSq", "Pr(>ChiSq)")] <- x$ChiSqTest[,-1]
		}
		
		VC$DF <- x$ConfInt$VC$TwoSided$DF
		
		VC$"CI LCL" <-  as.numeric(rep(NA, nrow(VC)))
		VC$"CI UCL" <-  as.numeric(rep(NA, nrow(VC)))
		
		VC[as.character(x$ConfInt$VC$TwoSided[,"Name"]), c("CI LCL", "CI UCL")] <- x$ConfInt$VC$TwoSided[,c("LCL", "UCL")]
		
		VC$"One-Sided LCL" <-  as.numeric(rep(NA, nrow(VC)))
		VC$"One-Sided UCL" <-  as.numeric(rep(NA, nrow(VC)))
		
		VC[as.character(x$ConfInt$VC$OneSided[,"Name"]), c("One-Sided LCL", "One-Sided UCL")] <-x$ConfInt$VC$OneSided[,c("LCL", "UCL")]
		
		VC <- as.matrix(VC)
		VC <- apply(VC, 1:2, round, digits)   
		
		if(!is.null(attr(x$ConfInt$VC$TwoSided, "LCLconstrained")))
			VC[attr(x$ConfInt$VC$TwoSided, "LCLconstrained"), "CI LCL"] <- paste(VC[attr(x$ConfInt$VC$TwoSided, "LCLconstrained"), "CI LCL"], "*", sep="")
		
		if(!is.null(attr(x$ConfInt$VC$TwoSided, "UCLconstrained")))
			VC[attr(x$ConfInt$VC$TwoSided, "UCLconstrained"), "CI UCL"] <- paste(VC[attr(x$ConfInt$VC$TwoSided, "UCLconstrained"), "CI UCL"], "*", sep="")
		
		if(!is.null(attr(x$ConfInt$VC$OneSided, "LCLconstrained")))
			VC[attr(x$ConfInt$VC$OneSided, "LCLconstrained"), "One-Sided LCL"] <- paste(VC[attr(x$ConfInt$VC$OneSided, "LCLconstrained"), "One-Sided LCL"], "*", sep="")
		
		if(!is.null(attr(x$ConfInt$VC$OneSided, "UCLconstrained")))
			VC[attr(x$ConfInt$VC$OneSided, "UCLconstrained"), "One-Sided UCL"] <- paste(VC[attr(x$ConfInt$VC$OneSided, "UCLconstrained"), "One-Sided UCL"], "*", sep="")
		
		print(noquote(VC), na.print="")
		ret$VC <- VC
	}
	
	if( any( c("all", "sd") %in% what) )
	{
		cat("\n> SD:\n-----\n")
		SD <- as.data.frame(as.matrix(VCAobj$aov.tab[,"SD", drop=FALSE]))
		colnames(SD) <- "Estimate"
		
		if(any(!is.na(x$ChiSqTest[,"Claim"])) && claim.type == "SD")
		{
			SD$Claim <- as.numeric(rep(NA, nrow(SD)))
			SD$"ChiSq" <- as.numeric(rep(NA, nrow(SD)))
			SD$"Pr(>ChiSq)" <- as.numeric(rep(NA, nrow(SD)))
			SD[as.character(x$ChiSqTest[,"Name"]), c("Claim", "ChiSq", "Pr(>ChiSq)")] <- x$ChiSqTest[,-1]
		}
		
		SD$DF <- x$ConfInt$SD$TwoSided$DF
		
		SD$"CI LCL" <-  as.numeric(rep(NA, nrow(SD)))
		SD$"CI UCL" <-  as.numeric(rep(NA, nrow(SD)))
		
		SD[as.character(x$ConfInt$SD$TwoSided[,"Name"]), c("CI LCL", "CI UCL")] <-x$ConfInt$SD$TwoSided[,c("LCL", "UCL")]
		
		SD$"One-Sided LCL" <-  as.numeric(rep(NA, nrow(SD)))
		SD$"One-Sided UCL" <-  as.numeric(rep(NA, nrow(SD)))
		
		SD[as.character(x$ConfInt$SD$OneSided[,"Name"]), c("One-Sided LCL", "One-Sided UCL")] <-x$ConfInt$SD$OneSided[,c("LCL", "UCL")]
		
		SD <- as.matrix(SD)
		SD <- apply(SD, 1:2, round, digits)  
		
		if(!is.null(attr(x$ConfInt$VC$TwoSided, "LCLconstrained")))                 # Note: attribute "LCLconstrained" exists only for VC since SD and CV computed from VC
			SD[attr(x$ConfInt$VC$TwoSided, "LCLconstrained"), "CI LCL"] <- paste(SD[attr(x$ConfInt$VC$TwoSided, "LCLconstrained"), "CI LCL"], "*", sep="")
		
		if(!is.null(attr(x$ConfInt$VC$TwoSided, "UCLconstrained")))                 # see above
			SD[attr(x$ConfInt$VC$TwoSided, "UCLconstrained"), "CI UCL"] <- paste(SD[attr(x$ConfInt$VC$TwoSided, "UCLconstrained"), "CI UCL"], "*", sep="")
		
		if(!is.null(attr(x$ConfInt$VC$OneSided, "LCLconstrained")))                 # see above
			SD[attr(x$ConfInt$VC$OneSided, "LCLconstrained"), "One-Sided LCL"] <- paste(SD[attr(x$ConfInt$VC$OneSided, "LCLconstrained"), "One-Sided LCL"], "*", sep="")
		
		if(!is.null(attr(x$ConfInt$VC$OneSided, "UCLconstrained")))                 # see above
			SD[attr(x$ConfInt$VC$OneSided, "UCLconstrained"), "One-Sided UCL"] <- paste(SD[attr(x$ConfInt$VC$OneSided, "UCLconstrained"), "One-Sided UCL"], "*", sep="")
		
		print(noquote(SD), na.print="")
		
		ret$SD <- SD
	}
	
	if( any( c("all", "cv") %in% what) )
	{
		cat("\n> CV[%]:\n--------\n")
		CV <- as.data.frame(as.matrix(VCAobj$aov.tab)[,"CV[%]", drop=FALSE])
		colnames(CV) <- "Estimate"
		
		if(any(!is.na(x$ChiSqTest[,"Claim"])) && claim.type == "CV")
		{
			CV$Claim <- as.numeric(rep(NA, nrow(CV)))
			CV$"ChiSq" <- as.numeric(rep(NA, nrow(CV)))
			CV$"Pr(>ChiSq)" <- as.numeric(rep(NA, nrow(CV)))
			CV[as.character(x$ChiSqTest[,"Name"]), c("Claim", "ChiSq", "Pr(>ChiSq)")] <- x$ChiSqTest[,-1]
		}
		
		CV$DF <- x$ConfInt$CV$TwoSided$DF
		
		CV$"CI LCL" <-  as.numeric(rep(NA, nrow(CV)))
		CV$"CI UCL" <-  as.numeric(rep(NA, nrow(CV)))
		
		CV[as.character(x$ConfInt$CV$TwoSided[,"Name"]), c("CI LCL", "CI UCL")] <-x$ConfInt$CV$TwoSided[,c("LCL", "UCL")]
		
		CV$"One-Sided LCL" <-  as.numeric(rep(NA, nrow(CV)))
		CV$"One-Sided UCL" <-  as.numeric(rep(NA, nrow(CV)))
		
		CV[as.character(x$ConfInt$CV$OneSided[,"Name"]), c("One-Sided LCL","One-Sided UCL")] <-x$ConfInt$CV$OneSided[,c("LCL", "UCL")]
		
		CV <- as.matrix(CV)
		CV <- apply(CV, 1:2, round, digits)   
		
		if(!is.null(attr(x$ConfInt$VC$TwoSided, "LCLconstrained")))                 # Note: attribute "LCLconstrained" exists only for VC since SD and CV computed from VC
			CV[attr(x$ConfInt$VC$TwoSided, "LCLconstrained"), "CI LCL"] <- paste(CV[attr(x$ConfInt$VC$TwoSided, "LCLconstrained"), "CI LCL"], "*", sep="")
		
		if(!is.null(attr(x$ConfInt$VC$TwoSided, "UCLconstrained")))                 # see above
			CV[attr(x$ConfInt$VC$TwoSided, "UCLconstrained"), "CI UCL"] <- paste(CV[attr(x$ConfInt$VC$TwoSided, "UCLconstrained"), "CI UCL"], "*", sep="")
		
		if(!is.null(attr(x$ConfInt$VC$OneSided, "LCLconstrained")))                 # see above
			CV[attr(x$ConfInt$VC$OneSided, "LCLconstrained"), "One-Sided LCL"] <- paste(CV[attr(x$ConfInt$VC$OneSided, "LCLconstrained"), "One-Sided LCL"], "*", sep="")
		
		if(!is.null(attr(x$ConfInt$VC$OneSided, "UCLconstrained")))                 # see above
			CV[attr(x$ConfInt$VC$OneSided, "UCLconstrained"), "One-Sided UCL"] <- paste(CV[attr(x$ConfInt$VC$OneSided, "UCLconstrained"), "One-Sided UCL"], "*", sep="")
		
		print(noquote(CV), na.print="")
		
		ret$CV <- CV
	}
	
	ci.method <- attr(x, "ci.method")
	
	if( any(c("all", "vc", "sd", "cv") %in% what) )
		cat( paste("\n\n", 100*(1-x$alpha), "% Confidence Level  ", ifelse(attr(x, "excludeNeg") && any(VCAobj$aov.tab[,"VC"] <= 0), "|  CIs for negative VCs excluded  ", ""),
						ifelse(!is.null(attr(x$ConfInt$VC$TwoSided, "LCLconstrained")), "| * CI-limits constrained to be >= 0", ""), 
						ifelse(ci.method=="sas", "\nSAS PROC MIXED method used for computing CIs", "\nSatterthwaite methodology used for computing CIs"), sep=""), "\n\n")
	
	if(length(ret) > 1)                     # extracting parts of the 'VCAinference' object
		invisible(ret)
	else if(length(ret) == 1)
		invisible(ret[[1]])
	
}


#'Extract Residuals of a 'VCA' Object
#'
#'Function extracts marginal or conditional residuals from a 'VCA' object, 
#'representing a linear mixed model.
#'
#'There are two types of residuals which can be extraced from a 'VCA' object.
#'Marginal residuals correspond to \eqn{e_m = y - \hat{y}}{e_m = y - y_hat}, where \eqn{\hat{y} = Xb}{y_hat = Xb} with \eqn{X}
#'being the design matrix of fixed effects and \eqn{b} being the column vector of fixed
#'effects parameter estimates. Conditional residuals are defined as \eqn{e_c = y - Xb - Zg},
#'where \eqn{Z} corresponds to the designs matrix of random effects \eqn{g}. 
#'Whenever 'obj' is a pure-error model, e.g. 'y~1' both options will return the same values 
#'\eqn{y - Xb} and \eqn{b} corresponds to the intercept.
#'Each type of residuals can be standardized, studentized, or transformed to pearson-type residuals. 
#'The former corresponds to a transformation of residuals to have mean 0 and variance equal to 1 (\eqn{(r - \bar{r})/\sigma_{r}}{[r - mean(r)]/sd(r)]}). 
#'Studentized residuals emerge from dividing raw residuals by the square-root of diagonal elements of the corresponding 
#'variance-covariance matrix. For conditional residuals, this is \eqn{Var(c) = P = RQR}, with \eqn{Q = V^{-1}(I - H)}{Q = V"(I - H)},
#'\eqn{H = XT} being the hat-matrix, and \eqn{T = (X^{T}V^{-1}X)^{-1}X^{T}V^{-1}}{T = (X'V"X)"X'V"}. For marginal residuals, this matrix
#'is \eqn{Var(m) = O = V - Q}. Here, >\eqn{^{T}}{'}< denotes the matrix transpose operator, 
#'and >\eqn{^{-1}}{"}< the regular matrix inverse. Pearson-type residuals are computed in the same manner as studentized, only
#'the variance-covariance matrices differ. For marginal residuals this is equal to \eqn{Var(y) = V}, for conditional residuals
#'this is \eqn{Var(c) = R} (see \code{\link{getV}} for details). 
#'
#'@param object		(VCA) object
#'@param type			(character) string specifying the type of residuals to be returned,
#'valid options are "marginal" and "conditional" or abbreviations
#'@param mode			(character) string or abbreviation specifying the specific transformation
#'applied to a certain type of residuals. There are "raw" (untransformed), 
#'"standardized", "studentized" and "pearson" (see details) residuals.
#'@param quiet			(logical) TRUE = will suppress any warning, which will be issued otherwise 
#'@param ...			additional parameters
#'
#'@author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#'
#'@method residuals VCA
#'
#'@references 
#'
#'Hilden-Minton, J. A. (1995). Multilevel diagnostics for mixed and hierarchical linear
#'models. Dissertation, University of California, Los Angeles.
#'
#'Nobre, J. S. & Singer, J. M. (2007). Residual analysis for linear mixed models. Biometrical
#'Journal, 49, 863-875.
#'
#'Schuetzenmeister, A. and Piepho, H.P. (2012). Residual analysis of linear mixed models using a simulation approach.
#'Computational Statistics and Data Analysis, 56, 1405-1416
#'
#'@examples
#'\dontrun{
#'data(VCAdata1)
#'datS1 <- VCAdata1[VCAdata1$sample==1,]
#'fit1  <- anovaVCA(y~(lot+device)/(day)/(run), datS1) 
#'
#'# default is conditional (raw) residuals
#'resid(fit1)
#'resid(fit1, "m")
#'
#'# get standardized version
#'resid(fit1, mode="stand")		# conditional residuals (default)
#'resid(fit1, "marg", "stand")		# marginal residuals
#'
#'# get studentized version, taking their 
#'# covariances into account
#'resid(fit1, mode="stud")		# conditional residuals (default)
#'resid(fit1, "marg", "stud")		# marginal residuals
#'}
#'
#'@aliases resid
#'
#'@seealso \code{\link{ranef}}, \code{\link{anovaVCA}}, \code{\link{anovaMM}}

residuals.VCA <- function(object, type=c("conditional", "marginal"), mode=c("raw", "student", "standard", "pearson"), quiet=FALSE, ...)
{		
	Call <- match.call()
	obj <- object
	
	if(is.list(obj) && !is(obj, "VCA"))
	{
		if(!all(sapply(obj, class) == "VCA"))
			stop("Only lists of 'VCA' object are accepted!")
		
		obj.len <- length(obj)
		
		if(!"msgEnv" %in% ls(.GlobalEnv))
			msgEnv <<- new.env(parent=emptyenv())
		
		assign("VCAinference.obj.is.list", TRUE, envir=msgEnv)			# indicate that a list-type object was passed intially
		
		res <- mapply(	FUN=residuals.VCA, object=obj,
				type=type[1], mode=mode[1], SIMPLIFY=FALSE)
		
		names(res) <- names(obj)
		
		if(obj.len == 1)			# mapply returns a list of length 2 in case that length(obj) was equal to 1
			res <- res[1]
		
		rm("VCAinference.obj.is.list", envir=msgEnv)
		
		return(res)
	}	
	
	stopifnot(class(obj) == "VCA")
	
	type <- match.arg(type)
	mode <- match.arg(mode)
	
	if(!is.null(obj$scale) && is.null(obj$rescaled))
		warning("The fitted model has not been re-scaled yet! Results are likely to differ from correct results!")
	
	if(is.null(obj$FixedEffects))
	{
		obj  <- solveMME(obj)
		nam0 <- deparse(Call$object)
		nam1 <- sub("\\[.*", "", nam0)
		
		if(length(nam1) == 1 && nam1 %in% names(as.list(.GlobalEnv)))
		{
			expr <- paste(nam0, "<<- obj")		# update object missing MME results
			eval(parse(text=expr))
		}
		else
		{
			if( !"VCAinference.obj.is.list" %in% names(as.list(msgEnv)) && !quiet)
				message("Some required information missing! Usually solving mixed model equations has to be done as a prerequisite!")
		}
	}
	
	Xb <- getMat(obj, "X") %*% obj$FixedEffects
	y  <- getMat(obj, "y")
	
	if(as.character(obj$terms)[3] == "1")							# pure error model
	{
		type <- "marginal"
	}
	
	if(!is.null(obj$scale) && mode != "raw")
		scale <- obj$scale
	else
		scale <- 1
	
	if(type == "marginal")											# marginal residuals
	{
		res <- y - Xb
		V   <- getMat(obj, "V")
		
		if(mode == "student")
		{			
			X  <- getMat(obj, "X")
			Vi <- getMat(obj, "Vi")
			
			#Q  <- X %*% MPinv(t(X) %*% Vi %*% X) %*% t(X)
			K   <- obj$VarFixed	/ scale^2																
			Q   <- X %*% K %*% t(X)
			res <- res/sqrt(diag(V-Q))								# apply studentization
		}
		else if(mode == "pearson")
		{
			res <- res/sqrt(diag(V))								# Pearson-type residuals
		}
	}
	else															# conditional residuals
	{
		Z   <- getMat(obj, "Z")
		if(obj$EstMethod == "REML")									# re-order Z to align with random effects
			Z <- Z[,unlist(obj$ColOrderZ)]
		
		Zg  <- Z %*% ranef(obj)										# order according to column-order of Z
		
		res <- y - Xb - Zg
		
		R  	<- getMat(obj, "R")
		
		if(mode == "student")
		{
			mats <- obj$Matrices
			
			Q	 <- mats$Q
			
			if(is.null(Q))
			{
				X  <- mats$X
				T  <- mats$T
				Vi <- mats$Vi
				mats$H <- H  <- X %*% T
				mats$Q <- Q  <- Vi %*% (diag(nrow(H))-H)
				
				nam0 <- deparse(Call$object)	
				nam1 <- sub("\\[.*", "", nam0) 
				
				if(length(nam1) == 1 && nam1 %in% names(as.list(.GlobalEnv)))	# write back to object in the calling env
				{
					obj$Matrices <- mats
					expr <- paste(nam0, "<<- obj")
					eval(parse(text=expr))
				}
			}
			P 	 <- R %*% Q %*% R
			res	 <- res / sqrt(diag(P))								# apply studentization
		}
		else if(mode == "pearson")
		{
			res <- res/sqrt(diag(R))								# Pearson-type residuals
		}
	}
	res <- c(res[,1])											
	
	if(mode=="standard")
	{
		res <- res / sd(res)										# apply standardization
	}
	names(res) <- rownames(obj$data)		
	attr(res, "type") <- type
	attr(res, "mode") <- mode
	res <- res/scale
	return(res)
}
