###############################################################################
#
# TODO: Add comment
# 
# Author: schueta6
#
# 'Fehler in solve.default(VCvar[nze, nze])' due to missing variation within a
# factor-level --> check if all variances are zero for a term
#
###############################################################################



#' Check Tow Formula Terms for Potential Problems.
#' 
#' Is is checked whether 'term2' is different from 'term1'
#' in adding information to the model. If both are main
#' factors, i.e. no interactions terms, it is checked whether
#' levels of 'term2' differ from those of 'term1'. Otherwise,
#' 'term2' is an interaction with a part being different from
#' 'term1'.
#' 
#' @param Data			(data.frame) containing all variables of
#' 						'term1' and 'term2'
#' @param term1			(character) term of a model formula as
#' 						returned by 'attr(terms(form), \"term.labels\")')
#' @param term2			(character) 2nd term of a model formula as
#' 						returned by 'attr(terms(form), \"term.labels\")')
#' 						to check against 'term1'
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @return (list) with components \"Diff\"=part of 'term2' distinguishing
#' 			it from 'term1', \"AddInfo\"=message informing about potential
#' 			problems with 'term2'

checkVars <- function(Data, term1, term2) {
	ia1  <- grepl(":", term1)				# term 1 interaction?
	ia2  <- grepl(":", term2)				# term 2?
	Diff <- term2							# differentiator probably term2						
	if(ia2) {								# term2 is interaction, term1 might be interaction
		if(grepl(term1, term2)) {			# term1 part of term2
			Diff <- sub(term1, "", term2)
			Diff <- sub("^:|:$", "", Diff)	# remove leading or trailing ':'
		} 
	}
	unq1 <- unique(Data[,unlist(strsplit(term1, ":")), drop=FALSE])
	unq2 <- unique(Data[,unlist(strsplit(term2, ":")), drop=FALSE])
	# same number of rows --> no additional information, yes, otherwise
	return(list(Diff=Diff, AddInfo=!(nrow(unq1) == nrow(unq2))))
}

#' Check Random Model for Given Dataset.
#' 
#' This function is intended to check a variance component analysis
#' either before or after performing it. This is particularily important
#' for less experienced users who my not exactly know where error messages 
#' come from. External software using functions \code{\link{anovaVCA}}
#' or \code{\link{remlVCA}} also via function \code{\link{fitVCA}} may
#' also benefit from more user-friendly error messages.
#' 
#' @param form			(formula) describing the model to be analyzed
#' @param Data			(data.frame) with all variables used in 'form'
#' 
#' @return  (list) of length equal to the number of terms in 'form' 
#' 			each element being a list of messages with identified, 
#' 			potential problems.
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @examples 
#' \dontrun{
#' data(dataEP05A2_1)
#' dat0 <- dataEP05A2_1[1:16,]
#' # everything should be ok
#' checkData(y~day/run, dat0)
#' # data identical response for all obs 
#' dat1 <- dat0
#' dat1$y <- dat1[1,"y"]
#' remlVCA(y~day/run, dat1)
#' checkData(y~day/run, dat1)
#' # now factor-levels have identical values
#' dat2 <- dat0
#' dat2$y <- dat2$y[rep(seq(1,7,2), rep(2,4))] 
#' checkData(y~day/run, dat2)
#' remlVCA(y~day/run, dat2, quiet=TRUE)
#' # indistinguishable terms are also problematic
#' dat3 <- data.frame(	y=rnorm(8,10),
#' 						day=paste("day",rep(c(1,2),c(4,4))), 
#' 						run=rep(c(2,1), c(4,4)))
#' checkData(y~day/run, dat3)
#' anovaVCA(y~day/run, dat3)
#' # no replicates, thus, no error variability
#' dat4 <- dat0[seq(1,15,2),]
#' dat4$day <- factor(dat4$day)
#' dat4$run <- factor(dat4$run)
#' checkData(y~day/run, dat4)
#' remlVCA(y~day/run, dat4)
#' }

checkData <- function(form, Data) {
	if(!is(form, "formula"))
		return(NULL)
	trm <- terms(form)
	res <- as.character(trm)[2]
	lab <- attr(trm, "term.labels")
	Nlv <- rep(0, length(lab))				# record number of levels for each term, interaction terms with the same number as main terms being part of the interaction do not add information
	ret <- vector("list", length(lab)+1)
	names(Nlv) <- lab
	names(ret) <- c(res, lab)
	# all equal values? if so, report back and do nothing else
	if(var(Data[,res]) == 0) {
		ret[[res]] <- "No variability among values of the response variable!"
		return(ret)
	}
	
	for(i in 1:length(lab)) {
		vrs <- unlist(strsplit(lab[i], ":"))
		tmp1 <- apply(Data[,vrs,drop=FALSE], 1, paste, collapse="_")
		# check if there are more observations than group-levels, will make problems in REML-estimation
		if(length(unique(tmp1)) == nrow(Data)) {
			ret[[lab[i]]] <- "Number of levels of each grouping factor must be < number of observations!"
			return(ret)
		}
		# check whether two variables are indistinguishable, will make problems in REML-estimation
		Nlv[lab[i]] <- length(unique(tmp1))
		if(i > 1) {
			for(j in 1:(i-1)) {
				chk <- checkVars(Data=Data, term1=lab[j], term2=lab[i])
				if(!chk$AddInfo)
					ret[[lab[i]]] <- c(ret[[lab[i]]], paste0("Levels of '", chk$Diff, "' are indistinguishable from levels of '", lab[j], "'!"))
			}
		}
		Var <- tapply(Data[,res], tmp1, var)
		# check potential problems with matrix inversion of VCvar in getGB
		if(all(Var < 2*.Machine$double.eps))
			ret[[lab[i]]] <- c(ret[[lab[i]]], paste0("REML-estimation error! Zero variation within factor-levels of '",lab[i],"'!"))
	}
	ret
}


#' Convert Objects to Detailed Error Message.
#' 
#' Function takes one or multiple objects and converts them to a single
#' error-message. Objects can be output of functions \code{\link{try}} or
#' \code{\link{checkData}}.
#' 
#' @param ...			one or multiple objects separated by comma
#' 
#' @return (characer) string combining information from input-objects
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' @examples 
#' \dontrun{
#' data(dataEP05A2_1)
#' dat2 <- dataEP05A2_1[1:16,]
#' dat2$y <- dat2$y[rep(seq(1,7,2), rep(2,4))] 
#' errorMessage(try(1/"a"), checkData(y~day/run, dat2))
#' }

errorMessage <- function(...) {
	args <- list(...)
	out  <- NULL
	for(i in 1:length(args)) {
		if(is(args[[i]], "try-error")) {
			out <- c(out, paste0("\n try-error:", as.character(args[[i]])))
		} else if(is(args[[i]], "list")) {
			frst <- TRUE
			for(j in 1:length(args[[i]])){
				if(!is.null(args[[i]][[j]])){
					if(frst)
						out <- c(out, " Additional error checks:\n")
					out <- c(out, paste0("  ", args[[i]][[j]], "\n"))	
				}
			}
		}
	}
	out
}


#' Wrap Function-Calls to Execute Additional Checks.
#' 
#' Function can be used to wrap function-calls, here, intended for model fitting
#' functions \code{\link{anovaVCA}}, \code{\link{anovaMM}}, \code{\link{remlVCA}}, \code{\link{remlMM}},
#' \code{\link{fitVCA}}, and \code{\link{fitLMM}}. When wrapped, there is the option to 
#' perform additional checks and reporting back identified problems by setting 'ErrorType="Detailed"'.
#' There is no error-handling provided by this function, i.e. any error issued will remain an error.
#' It would need to be handled by \code{\link{try}}, \code{\link{tryCatch}} or similar.
#' Note, that inline definition of datasets within 'expr' is not supported and will issue an error. 
#' 
#' @param expr			(expression) to be protected, typically, a call to a model-fitting 
#' 						function from this package (see details)
#' @param ErrorType		(ErrorType) "Simple"=default error-messages, "Detailed"= additional
#' 						data consistency checks will be performed
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @examples 
#' \dontrun{
#' # nothing happens if no error occurs
#' data(dataEP05A2_1)
#' res <- protectedCall(anovaVCA(form=y~day/run, Data=dataEP05A2_1))
#' res
#' # error message without additional consistency checks (default)
#' dat3 <- data.frame(	y=rnorm(8,10),
#' 						day=rep(c(1,2),c(4,4)), 
#' 						run=rep(c(2,1), c(4,4)))
#' protectedCall(anovaVCA(form=y~day/run, Data=dat3), ErrorType="Simple")
#' # error message with additional consistency checks hopefully helpful for the user
#' protectedCall(anovaVCA(form=y~day/run, Data=dat3), ErrorType="Detailed")
#' 
#' # handle error
#' res <- try(protectedCall(anovaVCA(form=y~day/run, Data=dat3), ErrorType="Detailed"), silent=TRUE)
#' if(is(res, "try-error"))
#' 	cat(sub(", ErrorType .*\\)", "", sub("protectedCall\\(", "", res)))
#' 
#' # inline-definition of data.frames issues an error
#' protectedCall(anovaVCA(	form=y~day/run, 
#' 							Data=data.frame(y=rnorm(8,10),
#' 								day=rep(c(1,2),c(4,4)), 
#' 								run=rep(c(2,1), c(4,4)))))
#' }

protectedCall <- function(expr, ErrorType=c("Simple", "Detailed")){
	call 		<- as.list(match.call())
	ErrorType	<- match.arg(ErrorType)
	callS		<- deparse(call$expr)
	if(length(callS) > 1)
		callS <- paste(callS, collapse="")
	callS		<- gsub(" ", "", callS)
	fun  		<- sub("\\(.*$", "", callS)[1]
	fun			<- sub("^.*=", "", sub("^.*<-", "", fun))
	fmls 		<- formals(fun) 
	args		<- sub(paste0("^.*", fun), "", callS)
	if(grepl("data.frame", args))
		stop("Inline-definition of datasets is not supported in 'protectedCall'!")
	args		<- sub("\\){1}$", "", sub("^\\({1}", "", args))
	args		<- unlist(strsplit(args, ","))
	args		<- gsub(" ", "", args)
	form		<- args[which(sapply(args, grepl, pattern="form"))]	
	if(length(form) == 0)
		form <- args[1]
	else
		form <- sub("form=", "", form)

	form <- as.formula(form)
	
	Data <- args[which(sapply(args, grepl, pattern="Data"))]	
	if(length(Data) == 0)
		Data <- args[2]
	else
		Data <- sub("Data=", "", Data)
	
	Data <- eval(parse(text=Data))
	
	### catch any errors and call consistency-check function on formula and data
	tryObj <- try(eval(parse(text=callS)), silent=TRUE)
	if(inherits(tryObj, "try-error")) {
		chk <- checkData(form, Data)
		if(ErrorType == "Simple") {
			stop(tryObj)
		} else {
			stop(errorMessage(tryObj, chk))
		}
	} else {
		return(tryObj)
	}
}




