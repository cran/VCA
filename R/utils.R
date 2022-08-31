# TODO: Add comment
# 
# Author: schueta6
###############################################################################



#'Scale Response Variable to Ensure Robust Numerical Calculations
#'
#'Function determines scaling factor for transforming the mean of the response to a range
#'between 0.1 and 1, applies scaling of the response and binds the scaling factor to the 
#'data as attribute 'scale'.
#'
#'@param Data			(data.frame) with the data to be fitted and the response to be scaled
#'@param resp			(character) name of the (numeric) response variable
#'
#'@return (data.frame) with the response scaled according to the scaling-factor,
#'which is recorded in the attribute \code{scale} of the data set
#'
#'@author Andre Schuetzenmeister \email{andre.schuetzenmester@@roche.com}

scaleData <- function(Data=NULL, resp=NULL)
{
	Mean 	<- mean(Data[,resp], na.rm=TRUE)								# scale response variable for more robust numerical optimization
	scale 	<- 10^ceiling(log10(abs(mean(Data[,resp], na.rm=TRUE))))	
	Data[,resp] <- Data[, resp] / scale
	attr(Data, "scale") <- scale
	Data
}


#'Load 'RevoUtilsMath'-package if available
#'
#'This function is taken from the Rprofile.site file of Microsoft R Open.
#'It was added to the package namespace to avoid a NOTE during the R CMD check
#'process stating that this function is not gobally defined.
#'
#'Only change to the original version is a different bracketing scheme to match
#'the one used in the remaining source-code of the package. 
#'
#'@param package		(character) package name to load, usually this will be package
#''RevoUtilsMath' if available
#'
#'@author Authors of the Rprofile.site file in Microsoft R Open.

load_if_installed <- function(package) 
{
	if (!identical(system.file(package="RevoUtilsMath"), "")) 
	{
		do.call('library', list(package))
		return(TRUE)
	} 
	else
		return(FALSE)
}



#'Re-Order Data.Frame
#'
#'Functions attempts to standardize input data for linear mixed model analyses
#'to overcome the problem that analysis results sometimes depend on ordering of
#'the data and definition of factor-levels.
#'
#'@param Data				(data.frame) with input data intented to put into standard-order
#'@param trms				(formula, terms) object speciying a model to be fitted to \code{Data}
#'@param order.data		(logical) TRUE = variables will be increasingly ordered, FALSE = order of
#'the variables remains as is
#'@param exclude.numeric	(logical) TRUE = numeric variables will not be included in the reordering,
#'which is required whenever this variable serves as covariate in a LMM, 
#'FALSE = numeric variables will also be converted to factors, useful in 
#'VCA-analysis, where all variables are interpreted as class-variables	
#'@param quiet				(logical) TRUE = omits any (potentially) informative output regarding
#'re-ordering and type-casting of variables
#'
#'@author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#'@examples 
#'\dontrun{
#'# random ordering
#'data(dataEP05A2_1)
#'dat <- dataEP05A2_1
#'levels(dat$day) <- sample(levels(dat$day))
#'# this has direct impact e.g. on order of estimated effects
#'fit <- anovaVCA(y~day/run, dat, order.data=FALSE)
#'ranef(fit)
#'# to guarantee consistent analysis results
#'# independent of the any data orderings option
#'# 'order.data' is per default set to TRUE:
#'fit <- anovaVCA(y~day/run, dat)
#'ranef(fit)
#'# which is identical to:
#'fit2 <- anovaVCA(y~day/run, orderData(dat, y~day/run), order.data=FALSE)
#'ranef(fit2)
#'}

orderData <- function(Data, trms, order.data=TRUE, exclude.numeric=TRUE, quiet=FALSE)
{
	stopifnot(identical(class(Data),"data.frame"))
	stopifnot("formula" %in% class(trms))
	if(!is(trms, "terms"))
		trms <- terms(trms)
	vars <- rownames(attr(trms, "factors"))[-1]
	if(length(vars) == 0)
		return(Data)
	stopifnot(all(vars %in% colnames(Data)))
	ord.vars <- NULL
	
	for(i in rev(vars))
	{
		if(is.numeric(Data[,i]) && exclude.numeric)
			next
		else
		{
			if(!quiet && is.character(Data[,i]))
				message("Convert variable ", i," from \"character\" to \"factor\"!")
			ord.vars <- c(paste0("Data[,\"",i,"\"]"), ord.vars)
			cha <- as.character(Data[,i])
			num <- suppressWarnings(as.numeric(cha))						# warnings are like to occur here
			if(!any(is.na(num)) && !any(grepl("0[[:digit:]]+", cha)))		# also handle elements like "01", "001", ... those will not give NAs converting them to numeric
			{
				if(order.data)
					Data[,i] <- factor(cha, levels=as.character(sort(unique(num))))
				else
					Data[,i] <- factor(cha, levels=as.character(unique(num)))
			}
			else
			{
				if(order.data)
					Data[,i] <- factor(cha, levels=sort(unique(cha)))
				else
					Data[,i] <- factor(cha, levels=unique(cha))
			}
		}
	}
	if(order.data)
	{
		exprs <- paste0("Data <- Data[order(", paste(ord.vars, collapse=","), "),]")
		eval(parse(text=exprs))
	}
	Data
}



#'Extract the Model Frame from a 'VCA' Object
#'
#'Function returns the data-element of 'object' and
#'adds the terms-element as attribute. 
#'
#'It enables application of functions relying on the existence of 
#'this method, e.g. the functin 'glht' of the 'multcomp'
#'R-package.
#'
#'@param formula		(VCA) object
#'@param ...			additional arguments
#'@return (data.frame) with attribute 'terms'
#'@method model.frame VCA
#'@author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}

model.frame.VCA <- function(formula, ...)
{
	stopifnot(class(formula) == "VCA")
	object <- formula
	rn <- rownames(attr(object$terms, "factors"))
	MF <- object$data[,rn]
	attr(MF, "terms") <- object$terms
	MF
}



#'Model Matrix of a Fitted VCA-Object
#'
#'Function returns matrix \code{X} corresponding
#'to the design matrix of fixed effects of the fitted
#'model.
#'
#'@param object			(VCA) object
#'@param ...				further arguments
#'@method model.matrix VCA

model.matrix.VCA <- function(object, ...)
{
	getMat(object, "X")
}




#'Automatically Scale Data Calling these Functions: 'anovaVCA', 'anovaMM', 'remlVCA' or 'remlMM'
#'
#'This function scales data before fitting a linear mixed model aiming to avoid numerical problems
#'when numbers of the response variable are either very small or very large. It adds attribute "scale" 
#'to the resulting 'VCA'-object, which is used by function \code{\link{reScale}} to transform back the
#'VCA-results of a \code{VCA} or \code{VCAinference} object that was previously scaled.
#'
#'NOTE: Scaling is applied on the complete data set, without checking whether there are incomplete
#'observations or not!
#'
#'@param Fun		(expr, function, character) either a complete function call to one of "anovaVCA", "anovaMM", "remlVCA", "remlMM",
#'a character string or just the function name without quotes (see example)
#'@param form		(formula) specifying the model to fitted by 'Fun'
#'@param Data		(data.frame) with all variables specified via 'Fun'
#'@param ...		additional arguments applying to one of the four functions \code{\link{anovaVCA}},\code{\link{anovaMM}},
#'\code{\link{remlVCA}}, \code{\link{remlMM}}
#'
#'@return (object) of class 'VCA' which can be used as input for function \code{\link{VCAinference}} 
#'
#'@author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#'
#'@seealso \code{\link{reScale}}
#'
#'@examples 
#'\dontrun{
#'data(dataEP05A2_3)
#'
#'# simulate very large numbers of the response
#'dat3   <- dataEP05A2_3
#'dat3$y <- dat3$y * 1e8
#'
#'# now try to fit 21-day model to this data
#'fit <- anovaVCA(y~day/run, dat3)
#'
#'# now use 'Scale' function
#'fit1 <- Scale("anovaVCA", y~day/run, dat3)
#'fit2 <- Scale(anovaVCA, y~day/run, dat3)	# also works
#'fit3 <- Scale(anovaVCA(y~day/run, dat3)) # works as well
#'
#'# back to original scale
#'(fit1 <- reScale(fit1))
#'(fit2 <- reScale(fit2))
#'(fit3 <- reScale(fit3))
#'
#'# reference values
#'fit0 <- anovaVCA(y~day/run, dataEP05A2_3, MME=TRUE)
#'inf0 <- VCAinference(fit0, VarVC=TRUE)
#'
#'fit1 <- Scale(anovaVCA(y~day/run, dataEP05A2_3, MME=TRUE))
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
#'fit1 <- Scale("remlVCA", y~day/run, dataEP05A2_3)
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
#'
#'# scaling also works with by-processing
#'data(VCAdata1)
#'fit <- Scale(anovaVCA(y~(device+lot)/day/run, VCAdata1, by="sample"))
#'reScale(fit)
#'}

Scale <- function(Fun, form, Data, ...)
{
	call <- as.list(match.call())			# check whether a complete function call was the sole argument
	
	if(is(call$Fun, "call"))
	{
		call <- as.list(call$Fun)
		Fun  <- as.character(call[[1]])
		form <- as.formula(call[[2]])
		Data <- get(as.character(call[[3]]))
		if(length(call) > 3)
			Args <- call[4:length(call)]
		else
			Args <- NULL
	}
	else	
		Args <- list(...)	
	
	if("by" %in% names(Args))
	{
		stopifnot(is.character(Args$by))
		stopifnot(Args$by %in% colnames(Data))
		stopifnot(is.factor(Data[,Args$by]) || is.character(Data[,Args$by]))
		
		levels  	<- unique(Data[,Args$by])
		tmpArgs 	<- Args
		tmpArgs$by	<- NULL
		if(length(tmpArgs) == 0)
			tmpArgs <- NULL
		res 		<- lapply(levels, function(x) Scale(Fun=Fun, form=form, Data[Data[,Args$by] == x,], tmpArgs))
		names(res) 	<- paste0(Args$by, levels)
		return(res)
	}
	if(is.function(Fun))						# if a function was passed and not the name of it
		Fun <- deparse(substitute(Fun))
	fun  <- match.arg(Fun, choices=c("anovaVCA", "anovaMM", "remlVCA", "remlMM"))
	stopifnot(class(form) == "formula")
	stopifnot(identical(class(Data),"data.frame"))
	
	tobj <- terms(form)
	stopifnot(attr(tobj, "response")==1)
	resp		<- rownames(attr(tobj, "factors"))[1]
	
	Mean  		<- mean(Data[,resp], na.rm=TRUE)								# scale response variable for more robust numerical optimization
	scale 		<- 10^ceiling(log10(abs(mean(Data[,resp], na.rm=TRUE))))	
	Data[,resp] <- Data[,resp]/(scale)
	
	FunArgs <- list(form=form, Data=Data)
	FunArgs <- c(FunArgs, Args)
	
	obj <- do.call(Fun, FunArgs)
	
	obj$scale <- scale
	obj
}


#'Check for Availability of Intel's Math Kernel Library
#'
#'Majority of the code is borrowed from the Microsoft R Open Rprofile.site file.
#'In case MKL can be detected this information will be stored in a separate envrionment, which
#'other function know about. If so, an optimized version of function \code{\link{getGB}}
#'will be used which used ordinary matrix-objects instead of matrices defined by the
#'\code{Matrix}-package. This seems to accelerate computation time for large datasets
#'by up to factor 30.
#'
#'This function is for internal use only and therefore not exported.
#'
#'@return variable 'MKL' in envir "msgEnv" will be set to TRUE/FALSE
#'
#'@author 	Authors of the Rprofile.site file in Microsoft R Open,
#'Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}

check4MKL <- function()
{
	if("msgEnv" %in% ls(.GlobalEnv) && !is.null(msgEnv$MKL))
		return(msgEnv$MKL)
	else
	{
		msgEnv <<- new.env(parent=emptyenv())				
		
		if(Sys.info()["sysname"] == "Darwin")				# Mac OSx
		{		
			assign("MKL", TRUE, envir=msgEnv)				# set MKL to TRUE, although, it may not be installed --> no good method known to check for MKL under MacOS
		} 
		else 												# other operating systems
		{		
			MRO <- FALSE
			try(MRO <- load_if_installed("RevoUtilsMath"), silent=TRUE)			# function only exists in MRO environment
			if(!inherits(MRO, "try-error") && MRO)
				assign("MKL", TRUE, envir=msgEnv)
			else 
				assign("MKL", FALSE, envir=msgEnv)
		}
		
		return(msgEnv$MKL)
	}
}



#'Giesbrecht & Burns Approximation of the Variance-Covariance Matrix of Variance Components
#'
#'Compute variance covariance matrix of variance components of a linear mixed model
#'via the method stated in Giesbrecht and Burns (1985). 
#'
#'This function is not intended to be called by users and therefore not exported.
#'
#'@param obj		(object) with list-type structure, e.g. \code{VCA} object fitted by ANOVA
#'or a premature \code{VCA} object fitted by REML
#'@param tol		(numeric) values < 'tol' will be considered being equal to zero
#'
#'@return 	(matrix) corresponding to the Giesbrecht & Burns approximation
#'of the variance-covariance matrix of variance components
#'@author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com},
#'Florian Dufey \email{florian.dufey@@roche.com}
#'
#'@seealso \code{\link{vcovVC}}, \code{\link{remlVCA}}, \code{\link{remlMM}}
#'
#'@references
#'Searle, S.R, Casella, G., McCulloch, C.E. (1992), Variance Components, Wiley New York
#'
#'Giesbrecht, F.G. and Burns, J.C. (1985), Two-Stage Analysis Based on a Mixed Model: Large-Sample
#'Asymptotic Theory and Small-Sample Simulation Results, Biometrics 41, p. 477-486 
#'
#'@examples 
#'\dontrun{
#'data(dataEP05A2_3)
#'fit <- anovaVCA(y~day/run, dataEP05A2_3)
#'fit <- solveMME(fit)		# some additional matrices required
#'getGB(fit)
#'}

getGB <- function (obj, tol = 1e-12)
{
	Nvc <- obj$Nvc
	Z <- obj$Matrices$Zre
	Q <- obj$Matrices$Q
	if (is.null(Q)) 
	{
		X <- getMat(obj, "X")
		Vi <- getMat(obj, "Vi")
		T <- getMat(obj, "T")
		Q <- Vi - Vi %*% X %*% T
	}
	VCvar <- matrix(0, Nvc, Nvc)
	nze <- NULL
	for (i in 1:Nvc) 
	{
		VCi <- obj$VCoriginal[i]
		if (is.null(VCi))
			VCi <- obj$aov.tab[obj$re.assign$terms[i], "VC"]
		if (i < Nvc && abs(VCi) < tol) 
		{
			VCvar[i, ] <- VCvar[, i] <- 0
			next
		}
		else
			nze <- c(nze, i)
		Zi  <- Z[, which(obj$re.assign$ind == i)]
		ZiT <- t(Zi)
		
		for (j in i:Nvc) 
		{			
			Zj    <- Z[, which(obj$re.assign$ind == j)]
			ZjT   <- t(Zj)
			Zpart <- as.matrix(ZjT %*% Q %*% Zi)
			
			if(j == Nvc)
			{
				if(j == i)
					VCvar[i,j] <- sum(sum(Q*Q))
				else
				{
					ZiTQ <- as.matrix(ZiT %*% Q)
					VCvar[i,j] <- VCvar[j,i] <- sum(sum(ZiTQ*ZiTQ))		# trace A*t(B)
				}
			}
			else
				VCvar[j,i] <- VCvar[i,j] <- sum(sum(Zpart*Zpart))		# trace A*t(B)
		}
	}
	nze <- unique(nze)
	VCnames <- obj$VCnames
	if (!"error" %in% VCnames)
		VCnames <- c(VCnames, "error")
	rownames(VCvar) <- colnames(VCvar) <- VCnames
	
	VCvar[nze, nze] <- 2 * solve(VCvar[nze, nze])
	VCvar <- VCvar
	attr(VCvar, "method") <- "gb"
	VCvar
}


#'Generic Method for Extracting Random Effects from a Fitted Model
#'
#'@param object		(object)
#'@param ...			additional parameters
#'@seealso \code{\link{ranef.VCA}}

ranef <- function(object, ...)
	UseMethod("ranef")

#'Extract Random Effects from 'VCA' Object
#'
#'Extract random effects and possibly apply a transformation to them (standardization,
#'studentization).
#'
#'Extracting the 'RandomEffects' element of an 'VCA' object if this exists and applying
#'standardization (mean 0, sd 1) or studentization. For studentized random effects 
#'the i-th random effects is divided by the i-th main diagonal element of matrix \eqn{O = GZ^{T}QZG}{O = GZ'QZG},
#'where \eqn{G} is the covariance-matrix of random effects, \eqn{Z} is a design matrix assigning 
#'random effects to observations and matrix \eqn{Q = V^{-1}(I - H)}{Q = V"(I - H)} (see \code{\link{residuals.VCA}} for further details). 
#'
#'@param object		(VCA) object from which random effects shall be extracted
#'@param term			(character) string specifying a term (factor) for which random effects 
#'should be extracted, one can also specify an integer which is interpreted
#'as i-th element of 'obj$res.assign$terms'
#'@param mode			(character) string or abbreviation specifying whether "raw" residuals
#'should be returned or a transformed version c("student" or "standard")
#'@param quiet			(logical) TRUE = will suppress any warning, which will be issued otherwise 
#'@param ...			additional parameters
#'
#'@method ranef VCA
#'
#'@references 
#'
#'Searle, S.R, Casella, G., McCulloch, C.E. (1992), Variance Components, Wiley New York	
#'
#'Laird, N.M., Ware, J.H., 1982. Random effects models for longitudinal data. Biometrics 38, 963-974.
#'
#'Schuetzenmeister, A. and Piepho, H.P. (2012). Residual analysis of linear mixed models using a simulation approach.
#'Computational Statistics and Data Analysis, 56, 1405-1416
#'
#'@examples 
#'\dontrun{
#'data(dataEP05A2_1)
#'fit <- anovaVCA(y~day/run, dataEP05A2_1)
#'ranef(fit)
#'
#'# get variable-specific random effects (REs)
#'# both extract the same REs
#'ranef(fit, "day")
#'ranef(fit, 1)
#'
#'# get standardized REs
#'ranef(fit, "day:run", "standard")
#'
#'# or studentized REs
#'ranef(fit, 2, "stu")
#'}

ranef.VCA <- function(object, term=NULL, mode=c("raw", "student", "standard"), quiet=FALSE, ...)
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
		
		if(is.null(term))
		{
			res <- mapply(	FUN=ranef.VCA, object=obj,
					mode=mode[1], quiet=quiet, SIMPLIFY=FALSE)
		}
		else
		{
			res <- mapply(	FUN=ranef.VCA, object=obj, term=term,
					mode=mode[1], quiet=quiet, SIMPLIFY=FALSE)
		}
		names(res) <- names(obj)
		
		if(obj.len == 1)			# mapply returns a list of length 2 in case that length(obj) was equal to 1
			res <- res[1]
		
		rm("VCAinference.obj.is.list", envir=msgEnv)
		
		return(res)
	}	
	
	stopifnot(class(obj) == "VCA")
	mode <- match.arg(mode)
	
	if(!is.null(obj$scale) && is.null(obj$rescaled))
		warning("The fitted model has not been re-scaled yet! Results are likely to differ from correct results!")
	
	if(!is.null(obj$scale) && !mode %in% c("raw", "standard"))
		scale <- obj$scale
	else
		scale <- 1
	
	ObjNam  <- deparse(Call$object)
	ObjNam2 <- sub("\\[.*", "", ObjNam)
	
	if(is.null(obj$RandomEffects))
	{
		obj  <- solveMME(obj)
		
		if(length(ObjNam2) == 1 && ObjNam2 %in% names(as.list(.GlobalEnv)))
		{
			expr <- paste(ObjNam, "<<- obj")		# update object missing MME results
			eval(parse(text=expr))
		}
		else
		{
			if( !"VCAinference.obj.is.list" %in% names(as.list(msgEnv)) && !quiet)
				message("Mixed model equations were solved but results could not be assigned to 'VCA' object!")
		}
	}
	
	if(mode == "student" && is.null(obj$Matrices$Q))
	{
		mats <- obj$Matrices
		X  <- mats$X
		T  <- mats$T
		Vi <- mats$Vi
		mats$H <- H  <- X %*% T
		mats$Q <- Q  <- Vi %*% (diag(nrow(H))-H)
		obj$Matrices <- mats
		
		if(length(ObjNam2) == 1 && ObjNam2 %in% names(as.list(.GlobalEnv)))
		{
			expr <- paste(ObjNam, "<<- obj")		# update object missing MME results
			eval(parse(text=expr))
		}
		else
		{
			if( !"VCAinference.obj.is.list" %in% names(as.list(msgEnv)) && !quiet)
				message("Matrices 'H' and 'Q' were comuted but could not be assigned to 'VCA' object!")
		}
	}
	
	re <- obj$RandomEffects
	nam <- rownames(re)
	if(!is.null(term) && !is.character(term))
	{
		term <- obj$re.assign$terms[as.integer(term)]
	}
	
	if(mode == "standard")
	{
		tmp <- tapply(re, obj$re.assign$ind, scale)
		re  <- matrix(nrow=nrow(re))
		for(i in 1:length(obj$re.assign$terms))
			re[which(obj$re.assign$ind == i)] <- tmp[[i]]
		
		rownames(re) <- nam
	}
	else if(mode == "student")
	{
		G <- getMat(obj, "G")
		Q <- getMat(obj, "Q")
		Z <- getMat(obj, "Z")
		O <- G %*% t(Z) %*% Q %*% Z %*% G
		re <- re / sqrt(diag(O))
		re <- as.matrix(re)
		rownames(re) <- nam
	}
	
	re <- re/scale
	
	if(is.null(term) || !term %in% obj$re.assign$terms)
	{
		if(!is.null(term) && !term %in% obj$re.assign$terms && !quiet)
		{
			warning("There is no term in the random part of the formula corresponding to specified 'term'!")
		}
		
		ind <- NULL
		for(i in 1:length(obj$re.assign$terms))
			ind <- c(ind, which(obj$re.assign$ind == i))
		
		re <- re[ind,,drop=FALSE]
		
		attr(re, "mode") <- mode
		attr(re, "term") <- "all"
		
		return(re)
	}
	else
	{
		ind <- which(obj$re.assign$ind == which(obj$re.assign$terms == term)) 
		re  <- re[ind,,drop=F]
		attr(re, "mode") <- mode
		attr(re, "term") <- term
		
		return(re)
	}
}




#'Generic Method for Extracting Fixed Effects from a Fitted Model
#'
#'@param object		(object)
#'@param ...			additional parameters
#'@seealso \code{\link{fixef.VCA}}

fixef <- function(object, ...)
	UseMethod("fixef")

#'Extract Fixed Effects from 'VCA' Object
#'
#'Conveniently extracting the 'FixedEffects' element of an 'VCA' object. 
#'
#'The default is to return the fixed effects estimates together with their standard errors.
#'If setting 'type="complex"' or to an abbreviation (e.g. "c") additional inferential statistics
#'on these estimates will be returned, i.e. "t Value", "DF" and respective p-value "Pr > |t|". 
#'One can choose one of three denominator degrees of freedom ('ddfm')-methods. The implementation
#'of these methods are an attempt to align with the results of SAS PROC MIXED. See the respective
#'SAS-documentation for details.
#'
#'@param object		(VCA) object where fixed effects shall be extracted
#'@param type			(character) string or partial string, specifying whether
#'to return "simple" (reduced) or a rather "complex" (more detailed) 
#'information about fixed effects
#'@param ddfm			(character) string specifying the method used for computing the 
#'degrees of freedom of the t-statistic. Only used when type="complex".
#'Available methods are "contain", "residual", and "satterthwaite".
#'@param tol			(numeric) value representing the numeric tolerance use in comparisons, values
#'smaller than 'tol' will be considered equal to 0
#'@param quiet			(logical) TRUE = suppress warning messages, e.g. for non-estimable contrasts
#'@param ...			additional parameters
#'
#'@method fixef VCA
#'
#'@examples 
#'\dontrun{
#'data(dataEP05A2_1)
#'fit <- anovaVCA(y~day/(run), dataEP05A2_1)
#'fixef(fit)
#'
#'# for complex models it might take some time computing complex output
#'data(VCAdata1)
#'fit <- anovaMM(y~(lot+device)/(day)/(run), VCAdata1[VCAdata1$sample==2,])
#'fixef(fit, "c")
#'}

fixef.VCA <- function(object, type=c("simple", "complex"), ddfm=c("contain", "residual", "satterthwaite"), 
		tol=1e-12, quiet=FALSE, ...)
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
		
		res <- mapply(	FUN=fixef.VCA, obj=obj, type=type[1],
				ddfm=ddfm[1], tol=tol, quiet=quiet,
				SIMPLIFY=FALSE)
		
		names(res) <- names(obj)
		
		if(obj.len == 1)			# mapply returns a list of length 2 in case that length(obj) was equal to 1
			res <- res[1]
		
		rm("VCAinference.obj.is.list", envir=msgEnv)
		
		return(res)
	}	
	
	stopifnot(class(obj) == "VCA")
	type <- match.arg(type)
	if(length(ddfm) > 1 && type == "complex")
	{
		ddfm <- "satterthwaite"
		if(!quiet)
			message("Note: 'ddfm' not specified, option \"satterthwaite\" was used!")
	}
	ddfm <- match.arg(ddfm)
	
	vc <- vcov(obj, quiet=quiet)
	se <- suppressWarnings(sqrt(diag(vc)))
	fe <- obj$FixedEffects
	
	if(is.null(fe))								# solve mixed model equations first
	{
		obj  <- solveMME(obj)
		fe   <- obj$FixedEffects
		nam0 <- deparse(Call$object)
		
		nam1 <- sub("\\[.*", "", nam0)			# remove any index-operators 
		if(length(nam1) == 1 && nam1 %in% names(as.list(.GlobalEnv)))
		{
			expr <- paste(nam0, "<<- obj")		# update object missing MME results
			eval(parse(text=expr))
		}
		else			# message only if not called on list of VCA-objects
		{
			if( !"VCAinference.obj.is.list" %in% names(as.list(msgEnv)) && !quiet)
				message("Some required information missing! Usually solving mixed model equations has to be done as a prerequisite!")
		}
	}
	nam <- rownames(fe)
	
	if(type == "simple")
	{		
		fe <- matrix(fe, ncol=1)
		colnames(fe) <- "Estimate"
		fe <- cbind(fe, SE=se)		
		rownames(fe) <- nam
	}
	else
	{
		if(is.null(obj$VarCov))
			obj$VarCov <- vcovVC(obj)
		if(is.null(obj$VarFixed))
			obj$VarFixed <- vcov(obj)
		LC  <- diag(nrow(fe))	
		rownames(LC) <- rownames(fe)
		colnames(LC) <- rownames(fe)
		tst <- test.fixef(obj, LC, ddfm=ddfm, quiet=TRUE)
		tst <- cbind(tst, SE=se)
		fe  <- tst[,c(1,5,2,3,4),drop=F]
		NAs <- is.na(fe[,"Estimate"])
		if(any(NAs))
		{
			fe[which(NAs),"Estimate"] <- 0
			fe[which(NAs), 2:5] <- NA
		}
	}
	
	return(fe)
}



#'Contrast Matrix for LS Means
#'
#'Function determines appropriate contrast matrix for computing the LS Means of
#'each factor level of one or multiple fixed effects variables. 
#'
#'This functions implements the 5 rules given in the documentation of SAS PROC GLM for computing the LS Means.#'
#'The LS Means correspond to marginal means adjusted for bias introduced by unbalancedness.
#'
#'@param obj			(VCA) object
#'@param var			(character) string specifyig the fixed effects variable for which
#'the LS Means generating matrices should be computed
#'@param quiet			(logical) TRUE = will suppress any warning, which will be issued otherwise 
#'
#'@return	(matrix) where each row corresponds to a LS Means generating contrast
#'for each factor level of one or multiple fixed effects variable(s)
#'
#'@author Andre Schutzenmeister \email{andre.schuetzenmeister@@roche.com}
#'
#'@examples 
#'\dontrun{
#'data(dataEP05A2_1)
#'fit1 <- anovaMM(y~day/run, dataEP05A2_1)
#'
#'VCA:::lsmMat(fit1, "day")	# function not exported
#'VCA:::lsmMat(fit1, "run")
#'VCA:::lsmMat(fit1)			# is equal to listing all fixed terms
#'
#'# a more complex and unbalanced model
#'data(VCAdata1)
#'datS1 <- VCAdata1[VCAdata1$sample == 1, ]
#'set.seed(42)
#'datS1ub <- datS1[-sample(1:nrow(datS1))[1:25],]
#'fit2 <- anovaMM(y~(lot+device)/day/(run), datS1ub)
#'VCA:::lsmMat(fit2, c("lot", "device"))
#'}

lsmMat <- function(obj, var=NULL, quiet=FALSE)
{
	stopifnot(class(obj) == "VCA")
	if(!is.null(var))
		stopifnot(var %in% obj$fixed)
	X <- getMat(obj, "X")
	
	if(is.null(var))
		var <- obj$fixed
	terms 	  <- attr(obj$terms, "factors")
	dat.cls   <- sapply(obj$data[,rownames(terms)[-1]], class)	
	variables <- names(dat.cls)
	fe.assign <- obj$fe.assign						# assignment by integers
	fe.terms  <- attr(fe.assign, "terms")	
	
	fe.class  <- list()
	
	for(i in 1:length(fe.terms))
	{
		if(fe.assign[i] == 0)
		{
			fe.class[[i]] <- NULL
		}
		tmp <- unlist(strsplit(fe.terms[i], ":"))
		fe.class[[i]] <- dat.cls[tmp]
	}
	names(fe.class) <- fe.terms
	
	if(!all(var %in% fe.terms))
		stop("At least one element of 'var' is not a fixed term!")
	
	fe.tvec   <- fe.terms[fe.assign+ifelse(obj$intercept, 1, 0)]				# assignment by names
	
	terms.cls <- rep("factor", length(fe.terms))
	names(terms.cls) <- fe.terms
	fe.tab    <- table(fe.tvec)
	
	lsmm <- matrix(nrow=0, ncol=ncol(X))				
	cn  <- colnames(lsmm) <- colnames(X)
	rn  <- Means <- NULL 
	
	if(any(dat.cls == "numeric"))					# find non-dummy variable columns in X
	{
		num.var <- names(dat.cls[which(dat.cls == "numeric")])
		for(i in 1:length(terms.cls))
		{
			if(any(sapply(num.var, grepl, fe.terms[i])))
			{
				terms.cls[i] <- "numeric"	
			}
		}		
	}
	
	fe.cls 		<- terms.cls[fe.assign+ifelse(obj$intercept, 1, 0)]
	num.terms 	<- terms.cls == "numeric"
	
	if(any(num.terms))														# are there any numeric terms
	{
		ind <- which(num.terms)
		Means <- list()
		
		for(i in 1:length(ind))
		{
			tmp		 <- numeric()
			tmp      <- c(tmp, Nlev=as.numeric(fe.tab[fe.terms[ind[i]]]))	# number of columns in X for the current numeric term
			splt     <- unlist(strsplit(fe.terms[ind[i]], ":"))
			tmp.cls  <- dat.cls[splt]
			num.var  <- which(tmp.cls == "numeric")
			fac.var  <- which(tmp.cls == "factor")
			tmp.mean <- apply(obj$data[,names(num.var), drop=FALSE], 1, function(x) eval(parse(text=paste(x, collapse="*"))))	# multiply each value of multiple numeric variables for each row		
			tmp 	 <- c(tmp, mean=mean(tmp.mean))							# mean of current combination of numeric variables
			
			if(length(fac.var) > 0)											# at least one factor variable
			{
				for(j in 1:length(fac.var))									# determine number of levels for each factor variable
				{
					exprs <- paste("c(tmp,",names(fac.var[j]),"=length(unique(obj$data[,\"",names(fac.var[j]),"\"])))", sep="") 
					tmp   <- eval(parse(text=exprs))
				}
			}
			eval(parse(text=paste("Means[[\"", fe.terms[ind[i]], "\"]] <- tmp", sep="")))	# add to list 'Means' and use name of the numeric variable as name of the list element
		}
	}
	
	for(i in 1:length(var))									# over all fixed terms for which LS Means shall be computed
	{
		tsplt <- unlist(strsplit(var[i], ":"))				# LS Means can only be generated for pure factor variables or combinations of such, no numeric variable may interact
		
		if(any(dat.cls[tsplt] == "numeric"))
		{
			if(!quiet)
				warning("'",var[i],"' is \"numeric\", LS Means cannot be estimated,'",var[i],"' will be skipped!")
			next
		}
		lvl.ind  <- which(fe.tvec == var[i])				# indices of all columns representing levels of var[i]
		lvl.nam  <- cn[lvl.ind]								# use column names of X instead of indices
		lvl.asgn <- unique(fe.assign[lvl.ind])				# value of the fixed effects assignment indicator for the current term var[i]
		tmp.lsm  <- matrix(0, nrow=0, ncol=ncol(X))
		colnames(tmp.lsm) <- cn
		rem.ind.init  <- which(fe.tvec != var[i])			# (init)ial (rem)aining columns which need to be treated differently
		rem.nam.init  <- cn[rem.ind.init]
		rem.asgn.init <- fe.assign[rem.ind.init] 
		
		for(j in 1:length(lvl.ind))							# over levels (columns) of the current factor
		{
			rem.ind  <- rem.ind.init						# re-set for each level of the current (i-th) term
			rem.nam  <- rem.nam.init
			rem.asgn <- rem.asgn.init
			
			con.mat <- matrix(0, nrow=1, ncol=ncol(X))
			colnames(con.mat) <- cn
			splt <- unlist(strsplit(lvl.nam[j], ":"))		# get all sub-terms
			
			if(obj$intercept)
			{
				con.mat[1,"int"] <- 1
				rem.ind  <- rem.ind[ -which(cn == "int")]
				rem.nam  <- rem.nam[ -which(cn == "int")]
				rem.asgn <- rem.asgn[-which(cn == "int")]
			}
			covar <- fe.cls[rem.ind] == "numeric"			# covariates involved?
			
			if(any(covar))									# [rule 1] (SAS PROC GLM documentation contrast matrix for LS Means)
			{
				cov.ind   <- rem.ind[ which(covar)]			# column-index of current covariate in contrast matrix
				rem.ind   <- rem.ind[-which(covar)]
				rem.nam   <- rem.nam[-which(covar)]
				rem.asgn  <- rem.asgn[-which(covar)]
				cov.asgn  <- fe.assign[cov.ind]				# covariate-indices
				ucov.asgn <- unique(cov.asgn)				
				
				for(k in 1:length(ucov.asgn))				# over all terms representing covariates
				{
					tmp.term <- fe.terms[ucov.asgn[k] + ifelse(obj$intercept, 1, 0)]
					tmp.info <- Means[[tmp.term]]
					tmp.ind  <- cov.ind[which(cov.asgn == ucov.asgn[k])]
					cn.splt  <- t(sapply(cn[tmp.ind], function(x) unlist(strsplit(x, ":"))))
					
					if(nrow(cn.splt) == 1)													# atomic covariate use mean value
					{
						con.mat[,which(fe.assign == ucov.asgn[k])] <- tmp.info["mean"]
						next
					}
					
					cov.cont <- apply(cn.splt, 1, function(x) any(splt %in% x))				# any sub-terms of the current term found in the levels of the current covariate
					
					if(any(cov.cont))														# covariate is distributed according to a term that is also a sub-term  of the current factor-level 
					{
						con.mat[1,tmp.ind[which(cov.cont)]] <- tmp.info["mean"]/length(which(cov.cont))
					}
					else																	# covariate independent of the current factor-level
					{
						con.mat[1,tmp.ind] <- tmp.info["mean"]/tmp.info["Nlev"]
					}
				}
			}
			
			tmp.splt  <- unlist(strsplit(lvl.nam[j], ":"))									# [rule 2] handling all effects that are contained by the current effect
			contained <- sapply(rem.nam, function(x) 
						all(unlist(strsplit(x, ":")) %in% tmp.splt))				
			
			if(any(contained))
			{
				tmp.asgn  <- rem.asgn[which(contained)]										# this might be multiple elements and this might be different values, e.g. 
				contained <- rem.nam[which(contained)]
				rem.ind   <- rem.ind[-which(rem.asgn %in% tmp.asgn)]						# column belongig to the same term must not be hanled again --> remove them from remaining elements
				rem.nam   <- rem.nam[-which(rem.asgn %in% tmp.asgn)]
				rem.asgn  <- rem.asgn[-which(rem.asgn %in% tmp.asgn)] 
				con.mat[,contained] <- 1
			}
			
			con.mat[1,lvl.ind[j]] <- 1														# [rule 3] setting the columns corresponding to the current effect to 1, all others remain 0
			
			contain <- sapply(rem.nam, function(x) 											# [rule 4] consider effects that contain the current effect
						all(tmp.splt %in% unlist(strsplit(x, ":"))))  
			
			if(any(contain))
			{
				tmp.asgn  <- rem.asgn[which(contain)]
				utmp.asgn <- unique(tmp.asgn)
				
				for(k in 1:length(utmp.asgn))	
				{
					tmp.ind <- rem.ind[which(contain)]
					tmp.ind <- tmp.ind[which(tmp.asgn == utmp.asgn[k])]						# only those indices that belong to the k-th term (assignment value)
					con.mat[1, tmp.ind] <- 1 / length(tmp.ind)								
				}
				
				rem.ind  <- rem.ind[ -which(rem.asgn %in% utmp.asgn)]
				rem.nam  <- rem.nam[ -which(rem.asgn %in% utmp.asgn)]
				rem.asgn <- rem.asgn[-which(rem.asgn %in% utmp.asgn)]
			}
			
			if(length(rem.ind) > 0)															# [rule 5] there are still remaining, not yet treated effects
			{																				#          set these to 1/number of levels
				ufac.asgn <- unique(rem.asgn)
				
				for(k in 1:length(ufac.asgn))
				{
					con.mat[1, which(fe.assign == ufac.asgn[k])] <- 1/fe.tab[fe.terms[ufac.asgn[k] + ifelse(obj$intercept, 1, 0)] ]
				}
			}
			
			tmp.lsm <- rbind(tmp.lsm, con.mat)		
		}
		rownames(tmp.lsm) <- lvl.nam
		lsmm <- rbind(lsmm, tmp.lsm)
	}
	
	attr(lsmm, "var.class") <- dat.cls
	
	if(!is.null(Means))
		attr(lsmm, "means")	<- Means
	
	return(lsmm)
}




#'Check Whether Design Is Balanced Or Not
#'
#'Assess whether an experimental design is balanced or not.
#'
#'This function is for internal use only. Thus, it is not exported.
#'
#'The approach taken here is to check whether each cell defined by one level of a factor are all equal or
#'not. Here, data is either balanced or unbalanced, there is no concept of "planned unbalancedness" as
#'discussed e.g. in Searle et al. (1992) p.4. The expanded (simplified) formula is divided into main factors
#'and nested factors, where the latter are interaction terms. The \eqn{N}-dimensional contingency table, \eqn{N} being the
#'number of main factors, is checked for all cells containing the same number. If there are differences, the
#'dataset is classified as "unbalanced". All interaction terms are tested individually. Firstly, a single factor 
#'is generated from combining factor levels of the first \eqn{(n-1)} variables in the interaction term. The last variable
#'occuring in the interaction term is then recoded as factor-object with \eqn{M} levels. \eqn{M} is the number of factor
#'levels within each factor level defined by the first \eqn{(n-1)} variables in the interaction term. This is done to 
#'account for the independence within sub-classes emerging from the combination of the first \eqn{(n-1)} variables.
#'
#'@param form      (formula) object defining the experimental design.
#'@param Data      (data.frame) containing all variables appearing in 'form'.
#'@param na.rm     (logical) TRUE = delete rows where any is NA, FALSE = NAs are not removed, if there are NAs in the
#'response variable and all information in independent variables is available, then only the design is checked.
#'
#'
#'@return (logical) TRUE if data is balanced, FALSE if data is unbalanced (according to the definition of balance used)
#'
#'@author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#'
#'@examples 
#'
#'\dontrun{
#'data1 <- data.frame(site=gl(3,8), lot=factor(rep(c(2,3,1,2,3,1), 
#'rep(4,6))), day=rep(1:12, rep(2,12)), y=rnorm(24,25,1))
#'
#'# not all combinations of 'site' and 'lot' in 'data1'
#'
#'VCA:::isBalanced(y~site+lot+site:lot:day, data1)
#'
#'# balanced design for this model
#'
#'VCA:::isBalanced(y~lot+lot:day, data1)
#'
#'# gets unbalanced if observation is NA
#'
#'data1[1,"y"] <- NA
#'VCA:::isBalanced(y~lot+lot:day, data1)
#'VCA:::isBalanced(y~lot+lot:day, data1, FALSE)
#'}

isBalanced <- function(form, Data, na.rm=TRUE)
{
	form <- terms(form, simplify=TRUE, keep.order=TRUE)
	if(length(attr(form, "factors")) == 0)
		return(TRUE)
	
	
	
	rn   <- rownames(attr(form, "factors"))
	
	for(i in 2:length(rn))
	{
		if(!is.numeric(Data[[rn[i]]]))				# for all non-numeric variables
		{
			Data[[rn[i]]] <- factor(as.character(Data[[rn[i]]]))			# get rid of non-existing factor levels (e.g. in data sub-sets)
		}
	}
	
	if(na.rm)
	{
		stopifnot(attr(form, "response") == 1)
		resp <- rn[1]
		Data <- Data[,rn]                           # only used variable appearing in the formula
		Data <- na.omit(Data)
	}
	
	fac  <- attr(form, "term.labels")
	
	mainFac  <- character()
	balanced <- TRUE
	
	for(i in 1:length(fac))
	{
		if(!balanced)
			break
		
		tmp  <- unlist(strsplit(fac[i], ":"))		# split interaction terms into variables
		Nvar <- length(tmp) 
		
		if( Nvar == 1)                              # main factor
		{
			mainFac <- c(mainFac, tmp)       
		} 
		else                                        # nested factors ("a:b" interpreted as b nested in a, since no possible main factor "b" is taken into account)
		{
			tmp.df <- data.frame(var1=apply(Data[,tmp[-Nvar], drop=FALSE], 1, function(x){
								return(paste(paste(letters[1:length(x)], x, sep=""), collapse=""))
							}))
			var2    <- character()
			var1Lev <- unique(tmp.df$var1)
			
			for(j in 1:length(var1Lev))             # start testing each term with nested factors ("a:b")
			{
				ind  <- which(tmp.df$var1 == var1Lev[j])
				var2 <- c(var2, as.factor(as.integer(Data[ind, tmp[Nvar]]))) 
			}
			tmp.df$var2 <- var2
			tmp.tab     <- table(tmp.df)           	# generate contingency table -> ...
			
			balanced <- all(c(tmp.tab) == tmp.tab[1,1])		# ... check for equal numbers in cells of the contingency table
		}
	}
	
	if(length(mainFac) > 0 && balanced)           	# check all main factors if no evidence of unbalncedness
	{
		tmp.tab  <- table(Data[,mainFac])
		tmp.df   <- as.data.frame(tmp.tab)
		balanced <- all(tmp.df[,"Freq"] == tmp.df[1, "Freq"])
	}
	
	return(balanced)
}


#'Compute the Trace of a Matrix
#'
#'Function computes the sum of main-diagonal elements of a square matrix.
#'
#'@param x     	(matrix, Matrix) object
#'@param quiet		(logical) TRUE = will suppress any warning, which will be issued otherwise 
#'@return (numeric) value, the trace of the matrix

Trace <- function(x, quiet=FALSE)
{
	if(ncol(x) != nrow(x) && !quiet)
		warning("Matrix is not quadratic!")
	return(sum(diag(as.matrix(x))))
}



#'Standard 'as.matrix' Method for 'VCA' S3-Objects
#'
#'@param x         (VCA) object 
#'@param ...       additional arguments to be passed to or from methods.
#'
#'@return (matrix) equal to x$aov.tab with additional attributes "Mean" and "Nobs"
#'
#'@author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#'
#'@seealso \code{\link{as.matrix.VCAinference}}
#'
#'@method as.matrix VCA
#'
#'@examples 
#'\dontrun{
#'data(dataEP05A2_1)
#'fit <- anovaVCA(y~day/run, dataEP05A2_1)
#'as.matrix(fit)
#'}

as.matrix.VCA <- function(x, ...)
{
	Mean <- x$Mean
	Nobs <- x$Nobs
	mat <- as.matrix(x$aov.tab)
	attr(mat, "Mean") <- Mean
	attr(mat, "Nobs") <- Nobs
	return(mat)
}


#'Standard 'as.matrix' Method for 'VCAinference' S3-Objects
#'
#'This function makes use of the hidden feature of function \code{\link{print.VCAinference}} which invisibly returns character 
#'matrices of estimated variance components expressed as "VC" (variance component), "SD" (standard deviation) or "CV" (coefficient
#'of variation). If argument "what" is not specified, a named list will be returned with all three matrices.
#'
#'@param x         (VCAinference) object 
#'@param what		(character) one or multiple choices from "VC" (variance component), "SD" (standard deviation) or
#'"CV" (coefficient of variation)
#'@param digits	(integer) number of decimal digits to be used
#'@param ...       additional arguments to be passed to or from methods.
#'
#'@return 	(matrix) with point estimates, one- and two-sided confidence intervals and variances
#'of the estimated variance components
#'
#'@seealso \code{\link{print.VCAinference}}, \code{\link{as.matrix.VCA}}
#'
#'@author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#'
#'@method as.matrix VCAinference
#'
#'@examples 
#'\dontrun{
#'data(dataEP05A2_1)
#'fit <- anovaVCA(y~day/run, dataEP05A2_1)
#'inf <- VCAinference(fit, VarVC=TRUE)
#'as.matrix(inf, what="VC", digits=6)
#'as.matrix(inf, what="SD", digits=6)
#'as.matrix(inf, what="CV", digits=2)
#'
#'# request list of matrices
#'as.matrix(inf)
#'}

as.matrix.VCAinference <- function(x, what=c("VC", "SD", "CV"), digits=6, ...)
{
	what <- match.arg(what, choices=c("VC", "SD", "CV"), several.ok=TRUE)
	if(length(what) > 1)
	{
		obj <- x
		mat <- lapply(what, function(x) as.matrix(obj, what=x, digits=digits, ...))
		names(mat) <- what
	}
	else
	{
		out   <- capture.output(mat <- print(x, what=what, digits=digits))
		vals  <- suppressWarnings(c(mat))
		zero  <- which(vals == "0*")
		vals  <- suppressWarnings(as.numeric(vals))
		if(length(zero) > 0)
			vals[zero] <- 0
		mat  <- matrix( vals, ncol=ncol(mat), nrow=nrow(mat),
				dimnames=attr(mat, "dimnames"))
		attr(mat, "Method") <- x$VCAobj$EstMethod
		attr(mat, "conf.level") <- 1-x$alpha
		attr(mat, "unit") <- what
	}
	mat
}


#'Satterthwaite Approximation for Total Degrees of Freedom and for Single Variance Components
#'
#'This function estimates degrees of freedom of the total variance (type="total")
#'in random models or individual variance components (type="individual"). 
#'It bases on the results of the unified approach to ANOVA-type estimation 
#'of variance components as implemented in functions \code{\link{anovaVCA}} 
#'and \code{\link{anovaMM}}.
#'
#'Function is used internally, thus, it is not exported. Option 'type="total"' is used in 
#'functions \code{\link{anovaVCA}} and \code{\link{anovaMM}} for approximating total DF.
#'Option 'type="individual"' is used in function \code{\link{VCAinference}} when choosing
#''ci.method="satterthwaite"' for approximating DFs for individual variance components.
#'
#'@param MS        	(numeric) vector of sequential mean squares (ANOVA type-1).
#'@param Ci       		(matrix) where elements are numeric values representing the inverse of the coefficient
#'matrix for calculation of expected mean squares (see \code{\link{anovaVCA}}).
#'@param DF        	(numeric) vector with the degrees of freedom for each factor in a ANOVA type-1 model.
#'@param type			(character) string specifying whether "total" degrees of freedom should be approximated or those of
#'individual variance components
#'
#'@return numeric value representing the Satterthwaite DFs of the total variance.
#'
#'@author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#'
#'@examples 
#'
#'\dontrun{
#'data(dataEP05A2_2)
#'res <- anovaVCA(y~day/run, dataEP05A2_2)
#'VCA:::SattDF(res$aov.tab[-1,"MS"], getMat(res, "Ci.MS"), res$aov.tab[-1,"DF"], type="tot")
#'
#'# now approximating individual DF for variance components
#'VCA:::SattDF(res$aov.tab[-1,"MS"], getMat(res, "Ci.MS"), res$aov.tab[-1,"DF"], type="i")
#'}

SattDF <- function(MS, Ci, DF, type=c("total", "individual"))
{
	type <- match.arg(type)
	if(type == "individual")
	{
		res <- NULL
		for(i in 1:nrow(Ci))
			res <- c(res, SattDF(MS, Ci[i,,drop=FALSE], DF))
		return(res)
	}
	I <- matrix(1,nrow(Ci),1)
	t <- t(I) %*% Ci
	t <- c(t) * c(MS)
	
	r1 <- sum(t)^2 / sum(t^2/DF)        # DF total
	
	return(r1)  
}


#'Overparameterized Design Matrices
#'
#'Function \code{getMM} constructs overparameterized design matrices from a model formula and a data.frame.
#'
#'This function constructs the overparameterized design matrix for a given dataset 'Data' according to
#'the model formula 'form'. Each combination of factor-levels and or numeric variables is identified
#'and accounted for by a separate column. See examples for differences compared to function 'model.matrix' (stats).
#'This type of design matrix is used e.g. in constructing A-matrices of quadratic forms in \eqn{y} expressing
#'ANOVA sums of squares as such. This is key functionality of functions \code{\link{anovaVCA}} and \code{\link{anovaMM}}
#'used e.g. in constructing the coefficient matrix \eqn{C} whose inverse is used in solving for ANOVA Type-1 based
#'variance components.. 
#'
#'@param form			(formula) with or without response specifying the model to be fit
#'@param Data			(data.frame) with the data
#'@param keep.order		(logical) TRUE = terms in 'form' should keep their positions, otherwise
#'						main effects come first and all interactions will be put into increasing order
#'
#'@author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#'
#'@examples 
#'\dontrun{
#'# load example data (CLSI EP05-A2 Within-Lab Precision Experiment)
#'data(dataEP05A2_3)
#'tmpData <- dataEP05A2_3[1:10,] 
#'
#'# check out the differences
#'getMM(~day+day:run, tmpData)
#'model.matrix(~day+day:run, tmpData)
#'
#'# adapt factor variables in 'tmpData'
#'tmpData$day <- factor(tmpData$day)
#'
#'# check out the differences now
#'getMM(~day+day:run, tmpData)
#'model.matrix(~day+day:run, tmpData)
#'
#'# numeric covariate 'cov'
#'tmpData2 <- dataEP05A2_3[1:10,] 
#'tmpData2$cov <- 10+rnorm(10,,3)
#'model.matrix(~day*cov, tmpData2)
#'}
#'

getMM <- function(form, Data, keep.order=TRUE)
{
	stopifnot(class(form) == "formula")
	stopifnot(identical(class(Data),"data.frame"))
	tform <- terms(form, simplify=TRUE, data=Data, keep.order=keep.order)
	int   <- attr(tform, "intercept") == 1
	form  <- as.character(tform)
	fmat  <- attr(tform, "factors")
	if(length(fmat) == 0)											# case y~1 (y is response here)
		return(Matrix(1, ncol=1, nrow=nrow(Data)))
	tlabs <- rownames(fmat)[which(apply(fmat, 1, sum) > 0)]			# excludes response if available
	tlab.cls <- sapply(Data[,tlabs,drop=F], class)
	fac.obs  <- tlab.cls[tlab.cls == "factor"]						# factor variables in terms of form
	num.obs  <- tlab.cls[tlab.cls != "factor"]						# non-factor variables
	N <- nrow(Data)
	lvls <- strsplit(form[ifelse(length(form) == 3, 3, 2)], "\\+")  # obtain single factors
	lvls <- gsub(" ", "", unlist(lvls))  
	
	if(int)
	{
		mm <- matrix(1, nrow=N, ncol=1)                             # include intercept
		colnames(mm) <- "int"
		assign <- 0
	}
	else
	{
		mm <- matrix(nrow=N, ncol=0)
		assign <- NULL
	}
	
	for(i in 1:length(lvls))                                        # over terms in the formula
	{   
		if(grepl("-.?1", lvls[i]))
			lvls[i] <- sub("-.?1", "", lvls[i])
		
		if(grepl(":", lvls[i]))                                     # crossed terms
		{
			tmpVar  <- unlist(strsplit(lvls[i], ":"))
			facVar  <- tmpVar[tmpVar %in% names(fac.obs)]
			numVar  <- tmpVar[tmpVar %in% names(num.obs)]
			VarNum  <- length(numVar) > 0
			if(length(facVar) > 0)
			{
				tmpData <- Data[,facVar,drop=F]                                         # cols of factor variables in lvls[i]
				VarFac  <- TRUE
				if(ncol(tmpData) > 1)
				{
					tmpName <- t(apply(tmpData, 1, function(x) paste(facVar, x, sep="")))   # terms = factor levels concatenated to variable names 
					tmpName <- apply(tmpName, 1, paste, collapse=":")                       # terms connected by ":"
					tmpName <- unique(tmpName)                              
				}
				else
				{
					tmpName <- paste(facVar, unique(tmpData[,1]), sep="")#":")
				}
				
				if(VarNum)
				{
					nam0 <- rep(lvls[i], length(tmpName))
					
					for(j in 1:length(facVar))
					{
						mat  <- tmpName
						re   <- regexpr(paste0(facVar[j],".*:"), tmpName)				# string in the middle
						sub  <- 1
						if(re[1] == - 1)
						{																# string at the end
							re <- regexpr(paste0(facVar[j],".*"), tmpName)
							sub <- 0
						}
						mat  <- rbind(mat, re, attr(re, "match.length"))
						
						ss <- apply(mat, 2, function(x)substr(x[1], as.numeric(x[2]), as.numeric(x[3])-sub) )
						
						for(k in 1:length(ss))
							nam0[k] <- sub(facVar[j], ss[k], nam0[k])
					}
					tmpName <- nam0
				}
				
				fac  <- apply(tmpData, 1, paste, collapse="")	
				Data <- cbind(Data, as.factor(fac))                                     # add new factor-variable to Data
				colnames(Data)[ncol(Data)] <- lvls[i]
			}
			else																		# interaction of numeric variables
			{
				VarFac <- FALSE
				tmpName <- lvls[i]
			}
		}    
		else														# main factors (no interaction)
		{
			if(lvls[i] %in% names(fac.obs))
			{
				tmpName <- paste(lvls[i], unique(Data[,lvls[i]]), sep="")   
				VarNum  <- FALSE
				VarFac  <- TRUE
			}
			else
			{
				tmpName <- numVar <- lvls[i]
				VarNum  <- TRUE
				VarFac  <- FALSE
			}
		}
		
		if(VarFac)
		{
			eff <- unique(Data[,lvls[i]])
			Ni  <- length(eff)                                          # number of effects for the i-th factor
			tmp <- matrix(0, N, Ni)
			
			colnames(tmp) <- unique(tmpName)
			
			for(j in 1:Ni)                                              # over effects of the i-th factor 
			{
				tmp[which(Data[,lvls[i]] == eff[j]),j] <- 1
			}
		}
		else
		{
			Ni  <- 1
			tmp <- matrix(1,nrow=N)
			colnames(tmp) <- tmpName
		}        
		
		if(VarNum)										
		{
			for(j in 1:ncol(tmp))
			{
				for(k in 1:length(numVar))
					tmp[,j] <- tmp[,j] * Data[,numVar[k]]	
			}
		}
		mm <- cbind(mm, tmp)
		assign <- c(assign, rep(i, Ni))
	}
	mm <- Matrix(mm)
	attr(mm, "assign") <- assign
	
	return(mm)
}

#'Moore-Penrose Generalized Inverse of a Matrix
#'
#'This function is originally impelemented in package 'MASS' as function \code{ginv}. 
#'It was adapted to be able to deal with matrices from the 'Matrix' package,
#'e.g. sparse matrices.
#'
#'@param X			(object) two-dimensional, for which a Moore-Penrose inverse
#'					has to be computed
#'@param tol		(numeric) tolerance value to be used in comparisons
#'
#'@return (object) A Moore-Penrose inverse of X.
#'
#'@author Authors of the 'MASS' package.

MPinv <- function (X, tol = sqrt(.Machine$double.eps)) 
{
	if (length(dim(X)) > 2L) 
		stop("'X' must be two-dimensional")
	Xsvd <- svd(X)
	Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
	if (all(Positive)) 
		Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
	else if (!any(Positive)) 
		array(0, dim(X)[2L:1L])
	else 
		Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop = FALSE]))
}



#'Least Squares Means of Fixed Effects
#'
#'Computes Least Squares Means (LS Means) of fixed effects for fitted mixed models of class 'VCA'.
#'
#'Function computes LS Means of fixed effects and their corresponding 
#'standard errors. In case of setting argument 'type' equal to "complex" (or any
#'abbreviation) a \eqn{t}-test is performed on each LS Mean, returning degrees
#'of freedom, t-statistic and corresponding p-values. One can choose from one of three
#'denominator degrees of freedom ('ddfm')-methods.
#'
#'Actually, function \code{\link{test.fixef}} is called with the "no intercept" 
#'version of the fitted model. The "complex" option is significantly slower for unbalanced
#'designs (see \code{\link{test.fixef}} for details). In case that the 'VarCov' element of
#'the 'VCA' object already exists (calling \code{\link{vcovVC}}), which is the most time 
#'consuming part, results can be obtained in less amount of time.
#'
#'Standard Errors of LS Means are computed as \eqn{TPT^{T}}{T * P * T'}, where \eqn{T} is the
#'LS Means generating contrast matrix and \eqn{P} is the variance-covariance matrix of
#'fixed effects.
#'
#'Argument \code{at} can be used to modify the values of covariables when computing LS Means and/or
#'to apply different weighting schemes for (fixed) factor varialbes in the model, e.g. when the prevelance
#'of factor-levels differs from a uniform distribution. Usually, if the weighting scheme is not modified,
#'each factor-level will contribute \eqn{1/N} to the LS Mean, where \eqn{N} corresponds to the number of factor-levels. 
#'
#'Covariables have to be specified as 'name=value', where value can be a vector of length > 1. 
#'Each value will be evaluated for each row of the original LS Means contrast matrix. 
#'If multiple covariables are specified, the i-th element of covariable 1 will be matched with
#'the i-th element of covariable(s) 2...M, where \eqn{M} is the number of covariables in the model.
#'
#'To apply a different weighting scheme for factor-variables one has to specify 'factor-name=c(level-name_1=value_1,
#'level-name_2=value_2, ..., level-name_N=value_N)'. The sum of all 'value_i' elements must be equal to 1, otherwise,
#'this factor-variable will be skipped issuing a warning. If any levels 'level-name_i' cannot be found for 
#'factor-variable 'factor-name', this variable will also be skipped and a warning will be issued.
#'See the examples section to get an impression of how this works.
#'
#'@param obj			(VCA) object having at least one fixed effect
#'@param var			(character) string specifying a fixed effects variable for which
#'LS Means should be computed, defaults to all fixed effects, i.e. for
#'each level of a fixed effects variable ls means will be computed
#'@param type			(character) "simple" = fast version of computing LS means
#'@param ddfm			(character) string specifying the method used for computing the 
#'degrees of freedom of the t-statistic. Only used when type="complex".
#'Available methods are "contain", "residual", and "satterthwaite".
#'@param at			(list) where each element corresponds either to a (numeric) covariable or
#'to a factor-variable for which the weighting scheme should be adjusted.
#'See details section for a thorough description of how argument 'at' works
#'and also see the examples.
#'@param contr.mat		(logical) TRUE = the LS Means generating contrast-matrix will be added to the
#'result as attribute \code{contrasts}
#'@param quiet			(logical) TRUE = suppress warning messages, e.g. for non-estimable contrasts
#'
#'@return (matrix) with LS Means of fixed effects and respective standard errors,
#'in case of 'type="complex"'
#'
#'@author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#'
#'@examples #
#'\dontrun{
#'data(dataEP05A2_2)
#'fit1 <- anovaMM(y~day/(run), dataEP05A2_2)
#'lsmeans(fit1)
#'lsmeans(fit1,, "complex")
#'
#'# a more complex model
#'data(VCAdata1)
#'fit2 <- anovaMM(y~(lot+device)/(day)/(run), VCAdata1[VCAdata1$sample==2,])
#'lsmeans(fit2, "lot")
#'lsmeans(fit2, "device", "complex")
#'
#'# pre-computed 'VarCov' element saves time
#'system.time(lsm1 <- lsmeans(fit2, "device", "complex"))
#'fit2$VarCov <- vcovVC(fit2)
#'system.time(lsm2 <- lsmeans(fit2, "device", "complex"))
#'lsm1
#'lsm2 
#'
#'# simulate some random data 
#'set.seed(212)
#'id <- rep(1:10,10)
#'x <- rnorm(200)
#'time <- sample(1:5,200,replace=T)
#'y <- rnorm(200)+time
#'snp <- sample(0:1,200,replace=T)
#'dat <- data.frame(id=id,x=x,y=y,time=time,snp=snp)
#'dat$snp <- as.factor(dat$snp)
#'dat$id <- as.factor(dat$id)
#'dat$time <- as.numeric(dat$time)
#'dat$sex <- gl(2, 100, labels=c("Male", "Female"))
#'dat$y <- dat$y + rep(rnorm(2, 5, 1), c(100, 100))
#'
#'fit3 <- remlMM(y~snp+time+snp:time+sex+(id)+(id):time, dat)
#'
#'# compute standard LS Means for variable "snp"
#'lsmeans(fit3, var="snp")
#'lsmeans(fit3, var="snp", type="c")    # comprehensive output
#'
#'# compute LS Means at timepoints 1, 2, 3, 4
#'# Note: original LS Means are always part of the output
#'lsmeans(fit3, var="snp", at=list(time=1:4))
#'
#'# compute LS Means with different weighting scheme
#'# for factor-variable 'sex'
#'lsmeans(fit3, var="snp", at=list(sex=c(Male=.3, Female=.7)))
#'
#'# combine covariables at some value and altering the
#'# weighting scheme
#'lsmeans(fit3, var="snp", at=list(time=1:4, sex=c(Male=.3, Female=.7)))
#'
#'# now with comprehensive output and requesting the
#'# LS Means generating contrast matrix
#'lsmeans(fit3, var="snp", type="complex", contr.mat=TRUE,
#'at=list(time=1:4, sex=c(Male=.3, Female=.7)))
#'}

lsmeans <- function(obj, var=NULL, type=c("simple", "complex"), ddfm=c("contain", "residual", "satterthwaite"), 
		at=NULL, contr.mat=FALSE, quiet=FALSE)
{
	Call <- match.call()
	
	if(is.list(obj) && !is(obj, "VCA"))
	{
		if(!all(sapply(obj, class) == "VCA"))
			stop("Only lists of 'VCA' object are accepted!")
		
		obj.len <- length(obj)
		
		if(is.null(var))
		{
			res <- mapply(	FUN=lsmeans, obj=obj, 
					type=type[1], ddfm=ddfm[1], quiet=quiet,
					SIMPLIFY=FALSE)
		}
		else
		{
			res <- mapply(	FUN=lsmeans, obj=obj, var=var, 
					type=type, ddfm=ddfm, quiet=quiet,
					SIMPLIFY=FALSE)
		}
		
		names(res) <- names(obj)
		
		if(obj.len == 1)			# mapply returns a list of length 2 in case that length(obj) was equal to 1
			res <- res[1]
		
		return(res)
	}	
	
	stopifnot(class(obj) == "VCA")
	stopifnot(obj$Type %in% c("Linear Model", "Mixed Model"))		# won't work for random models
	
	type <- match.arg(type)
	
	if(!is.null(var))												# does this fixed effects variable exist
		stopifnot(var %in% obj$fixed)
	
	if(length(ddfm) > 1 && type == "complex")
	{		
		ddfm <- "satterthwaite"
		if(!quiet)
			message("Note: 'ddfm' not specified, option \"satterthwaite\" was used!")
	}
	ddfm <- match.arg(ddfm)
	
	if(obj$EstMethod == "REML" && ddfm == "contain" && type=="complex")
	{
		ddfm <- "satterthwaite"
		if(!quiet)
			message("Note: 'ddfm' set to \"satterthwaite\", currently option \"contain\" does not work for REML-estimation!")
	}
	
	if(obj$Type != "Mixed Model" && !obj$intercept)
	{
		if(!quiet)
			warning("There are no fixed effects for which LS Means could be calculated!")
		return(NA)
	}
	
	T <- lsmMat(obj, var=var, quiet=quiet)					# LS Means generating contrast matrix
	
	id.mat <- NULL											# being non-NULL means that something was done via 'at'
	
	if(!is.null(at))										# LS Means at some fixed values of covariates and/or with different weighting scheme for fixed factor-variables
	{
		Means <- attr(T, "means")
		
		dat.class <- attr(T, "var.class")
		
		nam <- names(at)
		
		T2  <- T
		cnT <- colnames(T)
		fe.terms <- attr(obj$fe.assign, "terms")
		
		num.at <- nam[dat.class[nam] == "numeric"]
		fac.at <- nam[dat.class[nam] == "factor"]
		
		if(all(is.na(c(num.at, fac.at))) && !quiet)
			warning("Argument 'at' was not correctly specified! Neither covariables nor factor variables could be matched!")
		
		if(length(num.at) == 1 && is.na(num.at))								
			num.at <- character()							# ensure feasibility tests below
		if(length(fac.at) == 1 && is.na(fac.at))
			fac.at <- character()
		
		if(length(num.at) > 0)								# if there are any covariables
		{
			at.mat <- matrix(nrow=max(sapply(at[num.at], length)), ncol=length(num.at))
			colnames(at.mat) <- num.at
			
			for(i in 1:length(num.at))						# fill matrix, could also be user-defined
			{
				if(!num.at[i] %in% cnT)
				{
					if(!quiet)
						warning("There is no fixed term ", paste0("'", nam[i],"'!"),"Possible terms are:", paste(cnT, collapse=", "))
					next
				}
				
				if(length(at[[num.at[i]]]) < nrow(at.mat))
				{
					if(!quiet)
						warning("at[[",num.at[i],"]] does not match the max-length element of 'at', it will be replicated as necessary!")
					
					at.mat[,i] <- rep(at[[num.at[i]]], ceiling(nrow(at.mat)/length(at[[num.at[i]]])))[1:nrow(at.mat)]	# replicate
				}
				else
					at.mat[,i] <- at[[num.at[i]]]
			}
		}
		else
			at.mat <- matrix(nrow=length(fac.at), ncol=0)					# empty matrix
		
		if(length(fac.at) > 0)										# treat factor-variables, i.e. different weighting scheme
		{
			fac.lst <- vector("list", length(fac.at))
			names(fac.lst) <- fac.at
			
			for(i in 1:length(fac.at))								# over all specified factor variables
			{	
				if(sum(at[[fac.at[i]]]) != 1)
				{
					if(!quiet)
						warning("Sum of all coefficients of factor-variable ", paste0("'", fac.at[i],"'"), " is not equal to 1! It will be skipped!")
					next
				}
				
				skip <- FALSE
				
				tmp.mat <- matrix(nrow=nrow(at.mat), ncol=0)
				
				for(j in 1:length(at[[fac.at[i]]]))					# over all levels of the i-th factor variable
				{
					tmp <- matrix(at[[fac.at[i]]][j], ncol=1, nrow=nrow(at.mat))
					colnames(tmp) <- paste0(fac.at[i], names(at[[fac.at[i]]])[j])
					
					if(!colnames(tmp) %in% colnames(T))
					{
						if(!quiet)
							warning("Factor-level ", paste0("'",colnames(tmp), "'"),
									" of variable ",paste0("'", fac.at[i], "'"),
									" does not correspond to a fixed effect!\n  This element of 'at' will be skipped!")
						skip <- TRUE
					}
					tmp.mat <- cbind(tmp.mat, tmp)
				}
				
				if(skip)
					next
				else
				{
					at.mat <- cbind(at.mat, tmp.mat)
					fac.lst[[i]] <- tmp.mat
				}
			}
		}
		at.mat <- na.omit(at.mat)
		
		id.mat <- matrix(nrow=nrow(T), ncol=ncol(at.mat))
		cn.id  <- colnames(at.mat)
		colnames(id.mat) <- cn.id
		
		if(length(num.at) > 0)
		{
			for(i in 1:length(num.at))								# identifies combination of covar-levels
				id.mat[,num.at[i]] <- Means[[num.at[i]]]["mean"]
			
			tmp.ind <- cn.id[!cn.id %in% num.at]					# those columns not corresponding to covariables
		}
		else
			tmp.ind <- cn.id										# only columns corresponding to factor variable weights
		
		if(length(fac.at) > 0)
			id.mat[,tmp.ind] <- T[,tmp.ind]
		
		if(ncol(at.mat) > 0)								# only if there is anything to evaluate
		{
			for(i in 1:nrow(at.mat))						# for each combination specified
			{
				T3  <- T
				
				if(length(num.at) > 0)						# if there any covariables
				{
					for(j in 1:length(num.at))				# over covariables
					{
						cnTi <- cnT[grepl(num.at[j], cnT)]					# all relevant columns in matrix T					
						
						for(k in 1:length(cnTi))							# over all columns in which the current covariable appears (main-effect or interaction)
						{
							if(grepl(":", cnTi[k]))							# interaction
							{
								splt <- unlist(strsplit(cnTi[k], ":"))		# split interaction into atomic terms
								
								tmp.val <- 1
								
								for(l in 1:length(splt))							# over all terms in the interaction
								{
									if(splt[l] == num.at[j])						# user-specified level	
										tmp.val <- tmp.val * at.mat[i,num.at[j]]
									else
									{
										if(splt[l] %in% names(Means))						# another numeric covariable -> use mean of this value
											tmp.val <- tmp.val * Means[[splt[l]]]["mean"]
										else												# check with rowname
										{
											if(splt[l] %in% rownames(T3))
											{
												tmp.fac <- rep(0, nrow(T3))
												tmp.fac[grepl(splt[l], rownames(T3))] <- 1		# only that LS Mean affected, which is part of the interaction in column cnTi[k]
												tmp.val <- tmp.val * tmp.fac
											}
											else										# not part of the current interaction
												tmp.val <- tmp.val * rep(0, nrow(T3))
										}
									}
								}
								T3[,cnTi[k]] <- tmp.val
							}
							else
								T3[,cnTi[k]] <- at.mat[i,num.at[j]]					# main effects		
						}					
					}
				}
				
				if(length(fac.at) > 0)
				{
					for(j in 1:length(fac.at))									# over factor variables
					{
						tmp.row <- fac.lst[[fac.at[j]]][i,]
						T3[, tmp.ind] <- rep(tmp.row, rep(nrow(T3), length(tmp.row)))  
					}
				}
				
				T2 <- rbind(T2, T3)
				
				tmp.id.mat <- matrix(rep(at.mat[i,], nrow(T)), ncol=ncol(at.mat), byrow=TRUE)
				colnames(tmp.id.mat) <- colnames(at.mat)
				
				id.mat <- rbind(id.mat, tmp.id.mat)
			}
		}		
		T <- T2										# T2 might be identical to T
	}
	
	attr(T, "var.class") <- NULL
	attr(T, "means") 	 <- NULL
	
	x   <- obj
	if(x$intercept)
	{
		x$intercept <- FALSE		
	}
	
	if(type == "complex")						# complex output --> slower
	{	
		if(is.null(obj$VarCov))
			obj$VarCov <- vcovVC(obj)			# determine vcov of VCs if missing
		lsm <- test.fixef(obj, T, ddfm=ddfm, quiet=TRUE, lsmeans=TRUE)
		rownames(lsm) <- rownames(T)
	}
	else										# simple output --> faster
	{
		vc  <- vcov(obj)	
		vc  <- T %*% vc %*% t(T)
		se  <- sqrt(diag(vc))
		lsm <- T %*% obj$FixedEffects
		lsm <- cbind(as.matrix(lsm), SE=se)
	}
	
	if(!is.null(id.mat) && ncol(id.mat) > 0)
		lsm <- cbind(id.mat, lsm)
	
	if(contr.mat)
		attr(lsm, "contrasts") <- T
	
	return(lsm)
}


#'Determine V-Matrix for a 'VCA' Object
#'
#'Determine the estimated variance-covariance matrix of observations \eqn{y}.
#'
#'A linear mixed model can be written as \eqn{y = Xb + Zg + e}, where \eqn{y} is the column
#'vector of observations, \eqn{X} and \eqn{Z} are design matrices assigning fixed (\eqn{b}),
#'respectively, random (\eqn{g}) effects to observations, and \eqn{e} is the column vector of 
#'residual errors.
#'The variance-covariance matrix of \eqn{y} is equal to \eqn{Var(y) = ZGZ^{-T} + R}{Var(y) = ZGZ' + R}, where \eqn{R}
#'is the variance-covariance matrix of \eqn{e} and \eqn{G} is the variance-covariance matrix of \eqn{g}.
#'Here, \eqn{G} is assumed to be a diagonal matrix, i.e. all random effects \eqn{g} are mutually independent
#'(uncorrelated).
#'
#'@param obj			(VCA) object
#'
#'@return (VCA) object with additional elements in the 'Matrices' element, including matrix \eqn{V}.
#'
#'@author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}

getV <- function(obj)
{
	stopifnot(class(obj) == "VCA")
	mats <- obj$Matrices
	R <- Diagonal(obj$Nobs) * obj$aov.tab["error", "VC"]
	
	if(as.character(obj$terms)[3] == "1")
	{
		mats$V <- mats$R <- R
		mats$G <- mats$Zre <- NULL
		obj$Matrices <- mats
		return(obj)
	}
	Z <- Vs <- nam <- VCnam <- assign <- NULL
	Zi <- mats$Z
	VC <- mats$VCall									# all VCs, i.e. as if the model were a random model
	rf.ind <- mats$rf.ind
	
	rf.ind <- rf.ind[-length(rf.ind)]					# remove index to error VC
	ind <- 1
	for(i in rf.ind)									# constructing complete Z-matrix from VC-wise Z-matrices
	{
		if(ind == 1)
			Z <- cbind(Z, as.matrix(Zi[[i]]))
		else
			Z <- cbind(as.matrix(Z), as.matrix(Zi[[i]]))
		Vs <- c(Vs, rep(VC[i], ncol(Zi[[i]])))
		nam <- c(nam, colnames(Zi[[i]]))
		VCnam <- c(VCnam , attr(Zi[[i]], "term"))
		assign <- c(assign, rep(ind, ncol(Zi[[i]])))
		ind <- ind + 1									# ind vector specifies which REs belong to which terms in the formula
	}
	
	if(is.null(Z))
	{
		V <- R
		G <- NULL
	}
	else
	{
		Z <- Matrix(Z)
		assign <- list(ind=assign, terms=VCnam)
		G <- Diagonal(ncol(Z)) * Vs						# estimated variance-covariance matrix of random effects
		rownames(G) <- colnames(G) <- nam
		V <- Z %*% G %*% t(Z) + R						# compute V
	}
	colnames(Z) 	<- nam
	mats$Zre 		<- Z								# complete Z-matrix
	mats$G	 		<- G
	mats$V	 		<- V
	mats$R   		<- R
	obj$Matrices 	<- mats
	obj$re.assign 	<- assign
	
	return(obj)
}


#'Extract a Specific Matrix from a 'VCA' Object
#'
#'For convenience only, extracting a specific matrix from the 
#'"Matrices" element of a 'VCA' object if this matrix exists.
#'
#'When 'mat="Z"' the design matrix of random effects will be returned.
#'If one is interested in the design matrix of random effects for a specific
#'variance component use a name like "Z" + NAME, where NAME has to be equal to
#'the name of the VC in the 'VCA' object (see examples). The same applies to 
#'the A-matrices in the quadratic forms, use "A" + NAME for extracting a specific 
#'A-matrix.
#'
#'@param obj			... (VCA) object
#'@param mat			... (character) string specifying the matrix to be extracted
#'
#'@return (matrix) as requested by the user
#'
#'@examples
#'\dontrun{
#'data(dataEP05A2_1)
#'fit <- anovaVCA(y~day/run, dataEP05A2_1)
#'getMat(fit, "Z")
#'getMat(fit, "Zday")
#'getMat(fit, "Zday:run")
#'getMat(fit, "Zerror")
#'fit2 <- anovaMM(y~day/(run), dataEP05A2_1)
#'getMat(fit2, "V")			 	# Var(y)
#'getMat(fit2, "G")				# Var(re)
#'}

getMat <- function(obj, mat)
{
	stopifnot(class(obj) == "VCA")
	mats <- obj$Matrices
	if(is.null(mats))
		return(NULL)
	else
	{		
		if(grepl("^A", mat) || grepl("^Z", mat))
		{
			if(mat == "Z")
			{
				return(mats$Zre)
			}
			else
			{
				for(i in 1:length(mats$Z))					# same length as mats$A
				{
					Znam <- paste("Z", attr(mats$Z[[i]], "term"), sep="")
					Anam <- paste("A", attr(mats$A[[i]], "term"), sep="")
					
					if(mat == Znam)
					{
						return(mats$Z[[i]])
					}
					if(mat == Anam)
					{
						return(mats$A[[i]])
					}
				}
				return(NULL)
			}
		}
		else
		{
			return(mats[[mat]])
		}
	}
}


#'Extract Degrees of Freedom from Linear Hypotheses of Fixed Effects or LS Means
#'
#'Determine degrees of freedom for custom linear hypotheses of fixed effects or LS Means
#'using one of three possible approximation methods. 
#'
#'This is a convenience function to determine the DFs for linear hypotheses in the same way
#'as function \code{\link{test.fixef}}. Only the "DF" part is returned here which can be passed
#'to other functions expecting DFs as input.
#'
#'@param obj			(VCA) object
#'@param L				(matrix) specifying one or multiple linear hypothese, as returned by function
#'\code{\link{getL}} 
#'@param method		(character) the method to be used to determine the degrees of freedom for a 
#'linear hypothesis
#'@param ...			additional parameters
#'
#'@return (numeric) vector with the DFs for each row of 'L'
#'
#'@author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#'
#'@examples 
#'\dontrun{
#'data(VCAdata1)
#'tmpDat <- VCAdata1[VCAdata1$sample==1,]
#'tmpDat <- tmpDat[-c(11,51,73:76),]
#'fitMM <- anovaMM(y~(lot+device)/(day)/(run), tmpDat)
#'fitMM
#'L <- getL(fitMM, c("lot1-lot2", "device1-device2"))
#'getDF(fitMM, L)						# method="contain" is Default
#'getDF(fitMM, L, method="res")
#'
#'getDF(fitMM, L, method="satt")		# takes quite long for this model
#'}

getDF <- function(obj, L, method=c("contain", "residual", "satterthwaite"), ...)
{
	stopifnot(class(obj) == "VCA")
	method <- match.arg(method)
	return(test.fixef(obj, L=L, ddfm=method, onlyDF=TRUE))
}


#'Calculate Variance-Covariance Matrix and Standard Errors of Fixed Effects for an 'VCA' Object
#'
#'The variance-covariance matrix of fixed effects for the linear mixed model in 'obj' is calculated.
#'
#'The variance-covariance matrix of fixed effects for a linear mixed model corresponds to matrix
#'\eqn{(X^{T}V^{-1}X)^{-}}{(X'V"X)`}, where >\eqn{^{T}}{'}< denotes the transpose operator, >\eqn{^{-1}}{"}< 
#'the regular matrix inverse, and >\eqn{^{-}}{`}< the generalized (Moore-Penrose) inverse of a matrix.
#'
#'@param obj			(VCA) object for which the variance-covariance matrix of
#'fixed effects shall be calculated
#'@param quiet			(logical) TRUE = will suppress any warning, which will be issued otherwise 
#'
#'@return (matrix) corresponding to the variance-covariance matrix of fixed effects
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

vcovFixed <- function(obj, quiet=FALSE)
{
	Call <- match.call()
	if(!obj$intercept && length(obj$fixed) == 0)
	{
		if(!quiet)
			warning("There is no variance-convariance matrix of fixed effects for this object!")
		return(NA)
	}
	
	X  <- getMat(obj, "X")
	if(is.null(obj$Matrices$Vi))			# MME not yet been solved
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
			if(!quiet)
				message("Some required information missing! Usually solving mixed model equations has to be done as a prerequisite!")
		}
	}
	Vi <- getMat(obj, "Vi")
	VCov <- obj$VarFixed
	if(is.null(VCov))
		VCov <- Solve(t(X) %*% Vi %*% X, quiet=TRUE)
	#VCov <- MPinv(t(X) %*% Vi %*% X)
	rownames(VCov) <- colnames(VCov) <- rownames(obj$FixedEffects)
	
	return(VCov)
}


#'Variance-Covariance Matrix of Fixed Effects as Function of Covariance Parameter Estimates
#'
#'This is a helper function for function \code{\link{test.fixef}} approximating degrees of freedom for 
#'linear contrasts of fixed effect parameter estimates.
#'
#'@param obj			(VCA) object
#'@param x				(numeric) vector of covariance parameter estimates
#'
#'@return (matrix) corresponding to the variance-covariance matrix of fixed effects

DfSattHelper <- function(obj, x)
{
	stopifnot(class(obj) == "VCA")
	
	rf.ind <- obj$Matrices$rf.ind
	Zi <- obj$Matrices$Zre
	Z <- Vs <- nam <- NULL 
	
	ura <- unique(obj$re.assign[[1]])
	
	for(i in 1:length(x))
	{
		if(i != length(x))
		{
			ind <- which(obj$re.assign[[1]] == i)
			Z <- cbind(Z, as.matrix(Zi[,ind]))
			Vs <- c(Vs, rep(x[i], length(ind)))
			nam <- c(nam, colnames(Zi)[ind])
		}
		else
		{
			Z  <- cbind(Z, diag(obj$Nobs))
			Vs <- c(Vs, rep(x[i], obj$Nobs)) 
		}
	}
	
	G  <- diag(Vs)												# estimated variance-covariance matrix of random effects
	Z  <- Matrix(Z)
	R  <- getMat(obj, "R")
	X  <- getMat(obj, "X")
	V  <- Z %*% G %*% t(Z) + R
	Vi <- Solve(V)
	P  <- Solve(t(X) %*% Vi %*% X, quiet=TRUE)
	P <- as.matrix(P)	
	return(P)
}


#'Perform t-Tests for Linear Contrasts on LS Means
#'
#'Perform custom hypothesis tests on Least Squares Means (LS Means) of fixed effect.
#'
#'This function is similar to function \code{\link{test.fixef}} and represents a convenient way of specifying
#'linear contrasts of LS Means. 
#'
#'@param obj			(VCA) object
#'@param L				(matrix) specifying one or multiple custom hypothesis tests as linear contrasts of LS Means.
#'Which LS Means have to be used is inferred from the column names of matrix \eqn{L}, which need to 
#'be in line with the naming of LS Means in function \code{\link{lsmeans}}.
#'@param ddfm			(character) string specifying the method used for computing the denominator
#'degrees of freedom of t-tests of LS Means. Available methods are "contain", 
#'"residual", and "satterthwaite".
#'@param quiet			(logical) TRUE = will suppress any warning, which will be issued otherwise 
#'
#'@author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#'
#'@seealso \code{\link{test.fixef}}, \code{\link{lsmeans}}
#'
#'@examples
#'\dontrun{
#'data(dataEP05A2_2)
#'ub.dat <- dataEP05A2_2[-c(11,12,23,32,40,41,42),]
#'fit1 <- anovaMM(y~day/(run), ub.dat)
#'fit2 <- remlMM(y~day/(run), ub.dat)
#'lsm1 <- lsmeans(fit1)
#'lsm2 <- lsmeans(fit2)
#'lsm1
#'lsm2
#'
#'lc.mat <- getL(fit1, c("day1-day2", "day3-day6"), "lsm")
#'lc.mat[1,c(1,2)] <- c(1,-1)
#'lc.mat[2,c(3,6)] <- c(1,-1)
#'lc.mat
#'test.lsmeans(fit1, lc.mat) 
#'test.lsmeans(fit2, lc.mat)
#'
#'# fit mixed model from the 'nlme' package
#'
#'library(nlme)
#'data(Orthodont)
#'fit.lme <- lme(distance~Sex*I(age-11), random=~I(age-11)|Subject, Orthodont) 
#'
#'# re-organize data for using 'anovaMM'
#'Ortho <- Orthodont
#'Ortho$age2 <- Ortho$age - 11
#'Ortho$Subject <- factor(as.character(Ortho$Subject))
#'
#'# model without intercept
#'fit.anovaMM <- anovaMM(distance~Sex+Sex:age2+(Subject)+(Subject):age2-1, Ortho)
#'fit.remlMM1 <- remlMM( distance~Sex+Sex:age2+(Subject)+(Subject):age2-1, Ortho)
#'fit.remlMM2 <- remlMM( distance~Sex+Sex:age2+(Subject)+(Subject):age2-1, Ortho, cov=FALSE)
#'lsm0 <- lsmeans(fit.anovaMM)
#'lsm1 <- lsmeans(fit.remlMM1)
#'lsm2 <- lsmeans(fit.remlMM2)
#'lsm0
#'lsm1
#'lsm2
#'
#'lc.mat <- matrix(c(1,-1), nrow=1, dimnames=list("int.Male-int.Female", c("SexMale", "SexFemale")))
#'lc.mat
#'test.lsmeans(fit.anovaMM, lc.mat)	
#'test.lsmeans(fit.remlMM1, lc.mat)
#'test.lsmeans(fit.remlMM2, lc.mat)
#'}	

test.lsmeans <- function(obj, L, ddfm=c("contain", "residual", "satterthwaite"), quiet=FALSE)
{
	stopifnot(class(obj) == "VCA")
	lsmm <- lsmMat(obj, NULL, quiet=quiet)
	stopifnot(all(colnames(L) %in% rownames(lsmm)))			# only existing LS Means can be used
	
	lsmm <- lsmm[which(colnames(L) %in% rownames(lsmm)),]
	lsmm <- as.matrix(lsmm)
	
	U <- matrix(nrow=0, ncol=ncol(lsmm))
	colnames(U) <- colnames(lsmm)
	
	for(i in 1:nrow(L))
	{
		U <- rbind(U, matrix(L[i,], nrow=1) %*% lsmm)
	}
	
	if(is.null(obj$VarCov))
		obj$VarCov 	<- vcovVC(obj)
	
	res <- test.fixef(obj, U, ddfm=ddfm, quiet=quiet)
	rownames(res) <- rownames(L)
	return(res)
}


#'Perform t-Tests for Linear Contrasts on Fixed Effects
#'
#'This function performs t-Tests for one or multiple linear combinations (contrasts) of estimated 
#'fixed effects.
#'
#'Here, the same procedure as in \code{SAS PROC MIXED ddfm=satterthwaite} (sat) is implemented. 
#'This implementation was inspired by the code of function 'calcSatterth' of R-package 'lmerTest'. 
#'Thanks to the authors for this nice implementation. \cr
#'Note, that approximated Satterthwaite degrees of freedom might differ from 'lmerTest' and SAS PROC MIXED.
#'Both use the inverse Fisher-information matrix as approximation of the variance-covariance matrix
#'of variance components (covariance parameters). Here, either the exact algorithm for ANOVA-estimators of
#'variance components, described in Searle et. al (1992) p. 176, or the approximation presented in Giesbrecht and 
#'Burns (19985) are used. For balanced designs their will be no differences, usually. 
#'In case of balanced designs, the Satterthwaite approximation is equal to the degrees of freedom of the highest
#'order random term in the model (see examples).
#'
#'@param obj			(VCA) object
#'@param L				(numeric) vector or matrix, specifying linear combinations of the fixed effects, in the latter case,
#'each line represents a disctinct linear contrast
#'@param ddfm			(character) string specifying the method used for computing the denominator
#'degrees of freedom for tests of fixed effects or LS Means. Available methods are
#'"contain", "residual", and "satterthwaite".
#'@param method.grad	(character) string specifying the method to be used for approximating the gradient
#'of the variance-covariance matrix of fixed effects at the estimated covariance parameter
#'estimates (see function 'grad' (numDeriv) for details)
#'@param tol			(numeric) value specifying the numeric tolerance for testing equality to zero
#'@param quiet			(logical) TRUE = suppress warning messages, e.g. for non-estimable contrasts
#'@param opt			(logical) TRUE = tries to optimize computation time by avoiding unnecessary computations
#'for balanced datasets (see details). 
#'@param onlyDF		(logical) TRUE = only the specified type of degrees of freedom are determined without carrying out
#'the actual hypothesis test(s)
#'@param ...			further parameters (for internal use actually)
#'
#'@return (numeric) vector or matrix with 4 elements/columns corresponding to "Estimate", "t Value", "DF", and
#'"Pr > |t|".
#'
#'@author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com} inspired by authors of R-package 'lmerTest'
#'
#'@references 
#'
#'Searle, S.R, Casella, G., McCulloch, C.E. (1992), Variance Components, Wiley New York	
#'
#'Giesbrecht, F.G. and Burns, J.C. (1985), Two-Stage Analysis Based on a Mixed Model: Large-Sample
#'Asymptotic Theory and Small-Sample Simulation Results, Biometrics 41, p. 477-486
#'
#'SAS Help and Documentation PROC MIXED (MODEL-statement, Option 'ddfm'), SAS Institute Inc., Cary, NC, USA
#'
#'@seealso \code{\link{test.lsmeans}}, \code{\link{getL}}
#'
#'@examples
#'\dontrun{
#'data(dataEP05A2_2)
#'ub.dat <- dataEP05A2_2[-c(11,12,23,32,40,41,42),]
#'fit1 <- anovaMM(y~day/(run), ub.dat)
#'fit2 <- remlMM(y~day/(run), ub.dat)
#'fe1 <- fixef(fit1)
#'fe1
#'fe2 <- fixef(fit2)
#'fe2
#'lc.mat <- getL( fit1, c("day1-day2", "day3-day6"))
#'lc.mat
#'test.fixef(fit1, lc.mat, ddfm="satt") 
#'test.fixef(fit2, lc.mat, ddfm="satt")
#'
#'# some inferential statistics about fixed effects estimates
#'L <- diag(nrow(fe1))
#'rownames(L) <- colnames(L) <- rownames(fe1)
#'test.fixef(fit1, L)
#'test.fixef(fit2, L)
#'
#'# using different "residual" method determining DFs
#'test.fixef(fit1, L, ddfm="res")
#'test.fixef(fit2, L, ddfm="res")  
#'
#'# having 'opt=TRUE' is a good idea to save time 
#'# (in case of balanced designs)
#'data(VCAdata1)
#'datS3 <- VCAdata1[VCAdata1$sample==3,]
#'fit3 <- anovaMM(y~(lot+device)/(day)/run, datS3)
#'fit4 <- remlMM(y~(lot+device)/(day)/run, datS3)  
#'fit3$VarCov <- vcovVC(fit3)
#'fe3 <- fixef(fit3)
#'fe4 <- fixef(fit4)
#'L <- diag(nrow(fe3))
#'rownames(L) <- colnames(L) <- rownames(fe3)
#'system.time(tst1 <- test.fixef(fit3, L))
#'system.time(tst2 <- test.fixef(fit3, L, opt=FALSE))
#'system.time(tst3 <- test.fixef(fit4, L, opt=FALSE))
#'tst1
#'tst2
#'tst3
#'}

test.fixef <- function(	obj, L, ddfm=c("contain", "residual", "satterthwaite"),
		method.grad="simple", tol=1e-12, quiet=FALSE, opt=TRUE,
		onlyDF=FALSE, ...)
{
	stopifnot(class(obj) == "VCA")
	ddfm <- match.arg(ddfm)
	
	if(obj$EstMethod == "REML" && ddfm == "contain")
	{
		ddfm <- "satterthwaite"
		if(!quiet)
			message("Note: 'ddfm' set to \"satterthwaite\", currently option \"contain\" does not work for REML-estimation!")
	}
	
	args <- list(...)
	lsmeans <- args$lsmeans						# this function is called when LS Means with complex output is requested
	if(is.null(lsmeans))
		lsmeans <- FALSE
	
	b <- fixef(obj)[,"Estimate", drop=F]
	
	if(is.null(dim(L)))
		L <- matrix(L, nrow=1)
	
	stopifnot(ncol(L) == nrow(b))
	
	if(nrow(L) > 1)
	{
		res <- matrix(ncol=5, nrow=nrow(L))
		colnames(res) <- c("Estimate", "DF", "SE", "t Value", "Pr > |t|")
		
		for(i in 1:nrow(L))
			res[i,] <- test.fixef(obj=obj, L=L[i,,drop=F], ddfm=ddfm, method.grad=method.grad, quiet=quiet, opt=opt, onlyDF=onlyDF, lsmeans=lsmeans)
		
		rownames(res) <- rownames(L)
		if(onlyDF)
			return(res[,"DF"])
		attr(res, "ddfm") <- ddfm
		return(res)
	}	
	else
	{
		ind <- which(L != 0)						# test whether fixef(obj, "complex") was called
		
		if(length(ind) == 1 && any(abs(b[ind,]) < tol) && !lsmeans)
		{
			if(!quiet)
				warning("Linear contrast'", L ,"' not estimable!")
			res <- rep(NA, 5)
			names(res) <- c("Estimate", "DF", "SE", "t Value", "Pr > |t|")
			return(res)
		}		
		est		<- as.numeric(L %*% b)
		sgn		<- sign(est)
		X		<- getMat(obj, "X")
		P    	<- vcov(obj)
		lPl  	<- L %*% P %*% t(L)
		lPli 	<- solve(lPl)
		r       <- as.numeric(rankMatrix(lPli))
		SVD 	<- eigen(lPl)		
		items   <- list(SVD=SVD, r=r, b=b)	
		
		se <- L %*% P %*% t(L)						# compute standard error of linear contrast
		se <- sqrt(diag(se))
		
		DF		<- getDDFM(obj, L, ddfm, tol=tol, method.grad=method.grad, opt=opt, items=items)
		if(onlyDF)
		{
			return(c(NA, DF, NA, NA, NA))
		}
		t.stat 	<- as.numeric(sqrt((t(L %*% b) %*% lPli %*% (L %*% b))/r))
		
		res <- matrix(c(est, DF, se, sgn*t.stat, 2*pt(abs(t.stat), df=DF, lower.tail=FALSE)), nrow=1,
				dimnames=list(NULL, c("Estimate", "DF", "SE", "t Value", "Pr > |t|"))) 
		
		attr(res, "ddfm") <- ddfm
		
		return( res )
	}	
}


#'Construct Linear Contrast Matrix for Hypothesis Tests
#'
#'Function constructs coefficient/contrast matrices from a string-representation of linear hypotheses.
#'
#'Function constructs matrices expressing custom linear hypotheses of fixed effects or
#'LS Means. The user has to specify a string denoting this contrast which is then 
#'transformed into a coefficient/contrast matrix. This string may contain names of fixed effects
#'belonging to same same fixed term, numeric coefficients and mathematical operators "+"
#'and "-" (see examples).
#'
#'@param obj			(VCA) object
#'@param s				(character) string or vector of strings, denoting one or multiple
#'linear contrasts
#'@param what			(character) string specifying whether to construct contrast matrices
#'of fixed effects ("fixed") or LS Means ("lsmeans"), abbreviations are allowed.
#'
#'@return (matrix) representing one linear hypothesis of fixed effects or LS Means per row 
#'
#'@author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#'
#'@examples
#'\dontrun{
#'data(dataEP05A2_2)
#'fit <- anovaMM(y~day/(run), dataEP05A2_2)
#'L <- getL(fit, c("day1-day2", "day5-day10"), what="fixef")
#'L
#'test.fixef(fit, L=L)
#'
#'# another custom hypothesis
#'L2 <- getL(fit, "0.25*day1+0.25*day2+0.5*day3-0.5*day4-0.5*day5")
#'L2
#'
#'# more complex model
#'data(VCAdata1)
#'dataS2 <- VCAdata1[VCAdata1$sample==2,]
#'fit.S2 <- anovaMM(y~(lot+device)/day/(run), dataS2)
#'L3 <- getL(fit.S2, c("lot1-lot2", "lot1:device3:day19-lot1:device3:day20", 
#'"lot1:device1:day1-lot1:device1:day2"))
#'L3
#'test.fixef(fit.S2, L3)
#'}

getL <- function(obj, s, what=c("fixef", "lsmeans"))
{
	what <- match.arg(what)
	nwhat <- c(fixef="fixed effects", lsmeans="LS Means")
	
	if(what == "fixef")
		b <- fixef(obj, quiet=TRUE)
	else
		b <- lsmeans(obj, quiet=TRUE)
	n <- rownames(b)
	
	if(length(s) > 1)
	{
		L <- try(t(sapply(s, function(x) getL(obj, x, what=what))), silent=TRUE)
		if(is(L[1], "try-error"))
			stop(L[1])
		colnames(L) <- n
		return(L)
	}	
	contr <- rep(0, nrow(b))
	cname <- names(s)
	s <- gsub("\\*", "", s)
	
	splt <- sapply(unlist(strsplit(s, "\\+")), strsplit, "\\-")
	
	fac <- NULL
	sgn <- NULL
	
	for(i in 1:length(splt))
	{
		if(i == 1)
		{
			if(splt[[1]][1] == "")
				sgn <- -1
			else
				sgn <- 1 
		}
		else
			sgn <- 1									# 1st element always "+" since it was firstly used as split-char
		
		sgn <- c(sgn, rep(-1, length(splt[[i]])-1))
		
		for(j in 1:length(splt[[i]]))
		{
			tmp <- regexpr("^[[:digit:]]*\\.?[[:digit:]]*", splt[[i]][j])
			
			if(tmp > -1 && attr(tmp, "match.length") > 0)				# there is a factor at the beginning of the string
			{
				tmp.fac <- substr(splt[[i]][j], 1, attr(tmp, "match.length"))
				fac <- c(fac, as.numeric(tmp.fac) * sgn[j])	
				splt[[i]][j] <- sub(tmp.fac, "", splt[[i]][j])			# remove factor sub-string
			}
			else
				fac <- c(fac, 1 * sgn[j])
		}
	}
	
	splt <- unlist(splt)
	names(splt) <- NULL
	
	if(!all(splt %in% n))
		stop("\nError: There are terms which do not belong to the '",nwhat[what],"'!")
	
	res <- sapply(n, regexpr, splt)
	rownames(res) <- splt
	cnl <- nchar(colnames(res))	
	
	tms <- apply(res, 1, function(x){
				ind <- which(x > 0)
				len <- cnl[ind]
				return(ind[which(len == max(len))])
			})
	contr[tms] <- fac	
	return(matrix(contr, nrow=1, dimnames=list(cname, n)))
}


