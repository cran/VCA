# TODO: R-source file for graphical methods of the R-package VCA.
# 
# Author: schueta6
###############################################################################



#' Plot Random Variates of a Mixed Model ('VCA' Object).
#' 
#' Plots, possibly transformed, random variates of a linear mixed model (random effects, contitional or marginal
#' residuals).
#' 
#' Function plots either random effects of a 'VCA' object or residuals. Parameter 'term' is used to specify either
#' one. If 'term' is one of c("conditional", "marginal") corresponding residuals will be plotted 
#' (see \code{\link{resid}} for details). If 'term' is either the name of a random term in the formula of the 'VCA'
#' object or an integer specifying the i-th random term, corresponding random effects will be plotted. Both types
#' of random variates (random effects, residuals) can be plotted untransformed ("raw"), "studentized" or "standardized".
#' In case of residuals, one can also use the "Pearson"-type transformation.
#' 
#' @param obj			(VCA) object
#' @param term			(character, integer) specifying a type of residuals if one of c("conditional",
#' 						"marginal"), or, the name of a random term (one of obj$re.assign$terms). If 'term'
#' 						is a integer, it is interpreted as the i-th random term in 'obj$re.assign$terms'. 
#' @param mode			(character) string specifying a possible transformation of random effects or 
#'                      residuals (see \code{\link{residuals.VCA}} and \code{\link{ranef.VCA}}for details)
#' @param main			(character) string used as main title of the plot, if NULL, it will be automatically
#'                      generated
#' @param Xlabels		(list) passed to function \code{\link{text}} adding labels to the bottom margin at
#'                      x-coordinates 1:N, where N is the number of random variates. Useful for customization.
#' @param Points		(list) passed to function \code{\link{points}} for customization of plotting symbols
#' @param Vlines		(list) passed to function (abline) adding vertical lines, separating random variates
#'                      for better visual separation, set to NULL for omitting vertical lines.
#' @param pick			(logical) TRUE = lets the user identify single points using the mouse, useful, when many,
#'                      points were drawn where the X-labels are not readable.
#' @param ...			additional arguments to be passed to methods, such as graphical parameters (see \code{\link{par}})
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @examples 
#' 
#'  \dontrun{
#' data(dataEP05A2_1)
#' fit <- anovaVCA(y~day/run, dataEP05A2_1)
#' # solve mixed model equations including random effects
#' fit <- solveMME(fit)
#' plotRandVar(fit, "cond", "stand")	
#' plotRandVar(fit, 1, "stud")						# 1st random term 'day'
#' plotRandVar(fit, "day", "stud")					# equivalent to the above
#' 
#' # for larger datasets residuals can hardly be identified
#' # pick out interesting points with the mouse
#'
#' plotRandVar(fit, "marg", "stud", pick=TRUE)
#' 
#' # customize the appearance
#' plotRandVar( fit, 1, "stud", Vlines=list(col=c("red", "darkgreen")), 
#' 	Xlabels=list(offset=.5, srt=60, cex=1, col="blue"),
#' 	Points=list(col=c("black", "red", rep("black", 18)),
#'	pch=c(3,17,rep(3,18)), cex=c(1,2,rep(1,18))))	
#' } 

plotRandVar <- function(obj, term=NULL, mode=c("raw", "student", "standard", "pearson"), main=NULL, 
						Xlabels=list(), Points=list(), Vlines=list(), pick=FALSE, ...)
{
	Call <- match.call()
	
	stopifnot(class(obj) == "VCA")
	stopifnot(!is.null(term))
	
	if(length(mode) > 1)
		mode <- mode[1]
	
	ObjNam  <- as.character(as.list(Call)$obj)	
	
	if(is.null(obj$RandomEffects))					# compute required matrices
	{
		obj  <- solveMME(obj)
		
		if(!grepl("\\(", ObjNam))
		{
			expr <- paste(ObjNam, "<<- obj")		# update object missing MME results
			eval(parse(text=expr))
		}
		else
			warning("Some required information missing! Usually solving mixed model equations has to be done as a prerequisite!")
	}
	
	if(grepl(mode, "student") && is.null(obj$Matrices$Q))		# might be required when solveMME has been called yet
	{
		mats <- obj$Matrices
		X  <- mats$X
		T  <- mats$T
		Vi <- mats$Vi
		mats$H <- H  <- X %*% T
		mats$Q <- Q  <- Vi %*% (diag(nrow(H))-H)
		obj$Matrices <- mats
		
		if(!grepl("\\(", ObjNam))
		{
			expr <- paste(ObjNam, "<<- obj")		# update object missing MME results
			eval(parse(text=expr))
		}
		else
			warning("Some required information missing! Usually solving mixed model equations has to be done as a prerequisite!")
	}
	
	if(!is.character(term))
	{
		term <- obj$re.assign$terms[as.integer(term)]
		RE <- TRUE
		mode.opts <- c("raw", "student", "standard")		
	}
	else
	{
		term <- match.arg(term, choices=c("conditional", "marginal", obj$re.assign$terms))
		
		if(term %in% c("conditional", "marginal"))
		{
			RE <- FALSE
			mode.opts <- c("raw", "student", "standard", "pearson")
		}
		else
		{
			RE <- TRUE
			mode.opts <- c("raw", "student", "standard")			
		}
	}
	
	mode <- match.arg(mode, choices=mode.opts)
	m.names <- c(raw="raw", student="studentized", standard="standardized", pearson="Pearson-type")
	
	if(is.null(main))
		main <- paste(m.names[mode], if(RE)
							paste(" Random Effects: ",  "'", term, "'", sep="")
						else
							paste(term, "Residuals") )
	
	if(RE)
	{
		vec  <- ranef(obj, term=term, mode=mode)
		Ylab <- paste(m.names[attr(vec, "mode")], "Random Effects")
		nam  <- rownames(vec)
		vec  <- vec[,1]
		names(vec) <- nam
	}
	else
	{
		vec  <- resid(obj, type=term, mode=mode)
		Ylab <- paste(attr(vec, "mode"), "Residuals")
	}
	Nvec <- length(vec)
	
	nam <- names(vec)
	
	plot(1:Nvec, vec, main=main, axes=FALSE, xlab=NA, ylab=Ylab, type="n", ...)
	
	axis(2, las=1)
	axis(1, at=1:Nvec, labels=NA)
	
	if(!is.null(Vlines))
	{
		Vlines.def <- list(v=seq(1.5, Nvec-.5), col="lightgray")
		Vlines.def[names(Vlines)] <- Vlines
		Vlines <- Vlines.def
		do.call("abline", Vlines)
	}
	
	box()
	
	Points.def <- list(x=1:Nvec, y=vec, pch=3, cex=1, col="black")
	Points.def[names(Points)] <- Points
	Points <- Points.def
	do.call("points", Points)
	
	par(xpd=TRUE)
	
	usr <- par("usr")
	Xlabels.def <- list(x=1:Nvec, y=usr[3]-0.03*abs(usr[4]-usr[3]), srt=45, cex=.75, labels=nam, pos=2, offset=0)
	Xlabels.def[names(Xlabels)] <- Xlabels
	Xlabels <- Xlabels.def
	do.call("text", Xlabels)
	
	par(xpd=FALSE)
	
	if(pick)
		identify(x=1:Nvec, vec, nam)
}




#' Variability Chart for Hierarchical Models.
#' 
#' Function \code{varPlot} determines the sequence of variables in the model formula and uses this information to construct
#' the variability chart. 
#' 
#' This function implements a variability-chart, known from, e.g. JMP (JMP, SAS Institute Inc., Cary, NC). 
#' Arbitrary models can be specified via parameter 'form'. Formulas will be reduced to a simple hierarchical structure 
#' ordering factor-variables according to the order of appearance in 'form'. This is done to make function \code{varPlot} 
#' applicable to any random model considered in this package. 
#' Even if there are main factors, neither one being above or below another main factor, these are forced into a hierachy.
#' Besides the classic scatterplot, where observations are plotted in sub-classes emerging from the model formula, a plot of
#' standard deviations (SD) or coefficients of variation (CV) is provided (type=2) or both types of plots together (type=3).
#' 
#' @param Data          (data.frame) with the data
#' @param form          (formula) object specifying the model, NOTE: any crossed factors are reduced to last term of the crossing structure, i.e.
#'                      "a:b" is reduced to "b", "a:b:c" is reduced to "c". 
#' @param Title			(list) specifying all parameters applicable in function \code{\link{title}} for printing main- or sub-titles to plots. If 'type==3',
#'                      these settings will apply to each plot. For individual settings specify a list with two elements, where each element is a
#' 						list itself specifying all parameters of function 'title'. The first one is used for the variability chart, 
#' 						the second one for the SD or CV plot. Set to NULL to omit any titles.
#' @param keep.order    (logical) TRUE = the ordering of factor-levels is kept as provided by 'Data', FALSE = factor-levels are sorted on 
#'                      and within each level of nesting.
#' @param type          (integer) specifying the type of plot to be used, options are 1 = regular scatterplot, 2 = plot of the standard deviation,
#'                      3 = both type of plots.
#' @param VSpace		(numeric) vector of the same length as there are variance components, specifying the proportion of vertical space assigned
#'                      to each variance component in the tabular indicating the model structure. These elements have to sum to 1, otherwise equal
#' 						sizes will be used for each VC.
#' @param VARtype       (character) either "SD" (standard deviation) or "CV" (coefficient of variation), controls which type of measures is used to
#'                      report variability in plots when 'type' is set to either 2 or  (see 'type' above). Note that all parameters which apply to
#'                      the SD-plot will be used for the CV-plot in case 'VARtype="CV"'.
#' @param htab          (numeric) value 0 < htab < 1 specifying the height of the table representing the experimental design. This value represents
#'                      the proportion in relation to the actual plotting area, i.e. htab=1 mean 50\% of the vertical space is reserved for the table.
#' @param VarLab        (list) specifying all parameters applicable in function \code{\link{text}}, used to add labels within the table environment refering
#'                      to the nesting structure. This can be a list of lists, where the i-th list corresponds to the i-th variance component, counted
#'                      in bottom-up direction, i.e. starting from the most general variance component ('day' in the 1st example).
#' @param YLabel        (list) specifying all parameters applicable in function \code{\link{mtext}}, used for labelling the Y-axis.
#' @param SDYLabel      (list) specifying all parameters applicable in function \code{\link{mtext}}, used for labelling the Y-axis.
#' @param Points        (list) specifying all parameters applicable in function \code{\link{points}}, used to specify scatterplots per lower-end factor-level
#'                      (e.g. 'run' in formula run/day). If list-elements "col", "pch", "bg" and "cex" are lists themselves with elements "var" and "col"/"pch"/"bg"/"cex", 
#' 						where the former specifies a variable used for assigning colors/symbols/backgrounds/sizes according to the class-level of
#' 						variable "var", point-colors/plotting-symbols/plotting-symbol backgrounds/plotting-symbol sizes can be used for indicating
#' 						specific sub-classes not addressed by the model/design or indicate any sort of information (see examples).
#' 						Note the i-th element of 'col'/'pch' refers of the i-th element of unique(Data$var), even if 'var' is an integer variable.
#' @param SDs           (list) specifying all parameters applicable in function \code{\link{points}}, used to specify the appearance of SD-plots.
#' @param SDline        (list) specifying all parameters applicable in function \code{\link{lines}}, used to specify the (optional) line joining individual SDs,
#'                      Set to NULL to omit.
#' @param BG            (list) specifying the background for factor-levels of a nested factor. This list is passed on to function \code{\link{rect}} after element
#'                      'var', which identifies the factor to be used for coloring, has been removed. If not set to NULL and no factor has been specified by the
#'                      user, the top-level factor is selected by default. If this list contains element 'col.table=TRUE', the same coloring schema is used
#' 						in the table below at the corresponding row/factor (see examples). Addionally, list-elment 'col.bg=FALSE' can be used to turn off
#' 						BG-coloring, e.g. if only the the respective row in the table below should be color-coded (defaults to 'col.bg=TRUE').
#' 						When specifying as many colors as there are factor-levels, the same color will be applied to a factor-level automatically. 
#' 						This is relevant for factors, which are not top-level (bottom in the table).
#'                      Example: BG=list(var="run", col=c("white", "lightgray"), border=NA) draws the background for alternating levels of factor "run" 
#'                      white and gray for better visual differentiation. Set to NULL to omit. Use list( ..., col="white", border="gray") for using gray 
#'                      vertical lines for separation. See argument 'VLine' for additional highlighting options of factor-levels.
#' @param VLine			(list) specifying all parameters applicable in \code{\link{lines}} optionally separating levels of one or multiple variables
#' 						as vertical lines. This is useful in addition to 'BG' (see examples), where automatically 'border=NA' will be set that 'VLine' will
#' 						take full effect. If this list contains element 'col.table=TRUE', vertical lines will be extended to the table below the plot.
#' @param HLine			(list) specifying all parameters applicable in function \code{\link{abline}} to add horizontal lines. Only horizontal lines can be
#' 						set specifying the 'h' parameter. 'HLine=list()' will use default settings. 'HLine=NULL' will omit horizontal lines.
#' 						In case 'type=3', two separate lists can be specified where the first list applies to the variability chart and the second list
#' 						to the SD-/CV-chart.
#' @param ylim          (numeric) vector of length two, specifying the limits in Y-direction, if not set these values will be determined automatically.
#' 						In case of plot 'type=3' this can also be a list of two ylim-vectors, first corresponding to the variability chart, second to the
#' 						plot of error variability per replicate group
#' @param Join          (list) specifying all parameter applicable in function \code{\link{lines}} controlling how observed values within lower-level factor-levels,
#' 						are joined. Set to NULL to omit.
#' @param JoinLevels	(list) specifying all arguments applicable in function \code{\link{lines}}, joining factor-levels nested within higher order factor levels,
#' 						list-element "var" specifies this variable
#' @param Mean          (list) passed to function \code{\link{points}} specifying plotting symbols used to indicate mean values per lower-level factor-level, 
#' 						set equal to NULL to omit. 
#' @param MeanLine		(list) passed to function \code{\link{lines}} specifying the appearance of horizontal lines indicating mean values of factor levels. 
#' 						The factor variable for which mean-values of factor-levels are plotted can be specified via list-element "var" accepting any factor 
#' 						variable specified in 'form'. List element "mar" takes values in [0;.5] setting the left and right margin size of mean-lines.
#' 						Set equal to NULL to omit. Use 'var="int"' for specifying the overall mean (grand mean, intercept).
#' 						If this list contains logical 'join' which is set to TRUE, these mean lines will be joined. If list-element "top" is set to TRUE, 
#' 						these lines will be plotted on top, which is particularily useful for very large datasets.
#' @param Boxplot		(list) if not NULL, a boxplot of all values within the smallest possible subgroup (replicates) will be added to the plot,
#' 						On can set list-elements 'col.box="gray65"', 'col.median="white"', 'col.whiskers="gray65"' specifying different colors and 'lwd=3'
#' 						for the line width of the median-line and whiskers-lines as well as 'jitter=1e3' controlling the jittering of points around the
#' 						center of the box in horizontal direction, smallest possible value is 5 meaning the largest amount of jittering (1/5 in both directions)
#' 						value is)
#' @param VCnam         (list) specifying the text-labels (names of variance components) appearing as axis-labels. These parameters are passed to function
#'                      \code{\link{mtext}}. Parameter 'side' can only be set to 2 (left) or 4 (right) controlling where names of variance components appear. 
#' 						Set to NULL to omit VC-names. 
#' @param useVarNam     (logical) TRUE = each factor-level specifier is pasted to the variable name of the current variable and used as list-element name, 
#'                               FALSE = factor-level specifiers are used as names of list-elements; the former is useful when factor levels are indicated
#'                               as integers, e.g. days as 1,2,..., the latter is useful when factor levels are already unique, e.g. day1, day2, ... .
#' @param max.level		(integer) specifying the max. number of levels of a nested factor in order to draw vertical lines. If there are too many levels a black
#' 						area will be generated by many vertical lines. Level names will also be omitted.
#' @param ...			further graphical parameters passed on to function 'par', e.g. use 'mar' for specification of margin widths. Note, that not all of them
#' 						will have an effect, because some are fixed ensuring that a variability chart is drawn.
#' 
#' @return (invisibly) returns 'Data' with additional variable 'Xcoord' giving X-coordinates of each observation
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @examples
#' \dontrun{
#' 
#' # load data (CLSI EP05-A2 Within-Lab Precision Experiment)
#' data(dataEP05A2_3)
#' 
#' # two additional classification variables (without real interpretation)
#' dataEP05A2_3$user <- sample(rep(c(1,2), 40))
#' dataEP05A2_3$cls2 <- sample(rep(c(1,2), 40))
#' 
#' # plot data as variability-chart, using automatically determined parameter 
#' # settings (see 'dynParmSet')
#' varPlot(y~day/run, dataEP05A2_3)
#' 
#' # display intercept (total mean)
#' varPlot(y~day/run, dataEP05A2_3, MeanLine=list(var="int"))
#' 
#' # use custom VC-names
#' varPlot(y~day/run, dataEP05A2_3, VCnam=list(text=c("_Day", "_Run")))
#' 
#' # re-plot now also indicating dayly means as blue horizontal lines
#' varPlot(y~day/run, dataEP05A2_3, MeanLine=list(var=c("day", "int"), col="blue"))
#' 
#' # now use variable-names in names of individual factor-levels and use a different 
#' # notation of the nesting structure
#' varPlot(y~day+day:run, dataEP05A2_3, useVarNam=TRUE)
#' 
#' # rotate names of 2nd VC
#' varPlot(y~day+day:run, dataEP05A2_3, useVarNam=TRUE, VarLab=list(list(font=2), list(srt=60)))
#' 
#' # use alternating backgrounds for each level of factor "day" 
#' # (top-level factor is default) 
#' # use a simplified model formula (NOTE: only valid for function 'varPlot')
#' varPlot(y~day+run, dataEP05A2_3, BG=list(col=c("gray70", "gray90"), border=NA))
#' 
#' # now also color the corresponding row in the table accordingly
#' varPlot( y~day+run, dataEP05A2_3, 
#'          BG=list(col=c("gray70", "gray90"), border=NA, col.table=TRUE))
#' 
#' # assign different point-colors according to a classification variable
#' # not part of the model (artificial example in this case)
#' varPlot( y~day+day:run, dataEP05A2_3, mar=c(1,5,1,7), VCnam=list(side=4),
#'          Points=list(col=list(var="user", col=c("red", "green"))) )
#' 
#' # always check order of factor levels before annotating
#' order(unique(dataEP05A2_3$user))
#' 
#' # add legend to right margin
#' legend.m(fill=c("green", "red"), legend=c("User 1", "User 2"))
#' 
#' # assign different plotting symbols according to a classification
#' # variable not part of the model
#' varPlot( y~day+day:run, dataEP05A2_3, mar=c(1,5,1,7), VCnam=list(side=4), 
#'          Points=list(pch=list(var="user", pch=c(2, 8))) )
#' 
#' # add legend to right margin
#' legend.m(pch=c(8,2), legend=c("User 1", "User 2"))
#' 
#' # assign custom plotting symbols by combining 'pch' and 'bg'
#' varPlot( y~day+day:run, dataEP05A2_3, 
#'          Points=list(pch=list(var="user", pch=c(21, 24)),
#'                      bg=list( var="user", bg=c("lightblue", "yellow"))) )
#' 
#' # assign custom plotting symbols by combining 'pch', 'bg', and 'cex'      
#' varPlot( y~day+day:run, dataEP05A2_3,                                    
#'          Points=list(pch=list(var="user", pch=c(21, 24)),                
#'                      bg =list(var="user", bg=c("lightblue", "yellow")),
#'                      cex=list(var="user",  cex=c(2,1))) )
#' 
#' # now combine point-coloring and plotting symbols
#' # to indicate two additional classification variables
#' varPlot( y~day+day:run, dataEP05A2_3, mar=c(1,5,1,10), 
#'          VCnam=list(side=4, cex=1.5),
#'          Points=list(col=list(var="user", col=c("red", "darkgreen")), 
#'                      pch=list(var="cls2", pch=c(21, 22)),
#'                      bg =list(var="user", bg =c("orange", "green"))) )
#' 
#' # add legend to (right) margin
#'  legend.m( margin="right", pch=c(21, 22, 22, 22), 
#'            pt.bg=c("white", "white", "orange", "green"), 
#'            col=c("black", "black", "white", "white"), 
#'            pt.cex=c(1.75, 1.75, 2, 2), 
#'            legend=c("Cls2=1", "Cls2=2", "User=2", "User=1"),
#'            cex=1.5)
#' 
#' # use blue lines between each level of factor "run" 
#' varPlot(y~day/run, dataEP05A2_3, BG=list(var="run", border="blue"))
#' 
#' # plot SDs for each run
#' varPlot(y~day+day:run, dataEP05A2_3, type=2)
#' 
#' # use CV instead of SD
#' varPlot(y~day/run, dataEP05A2_3, type=2, VARtype="CV")
#' 
#' # now plot variability-chart and SD-plot in one window
#' varPlot(y~day/run, dataEP05A2_3, type=3, useVarNam=TRUE)
#' 
#' # now further customize the plot
#' varPlot( y~day/run, dataEP05A2_3, BG=list(col=c("lightgray", "gray")),
#'          YLabel=list(font=2, col="blue", cex=1.75, text="Custom Y-Axis Label"),
#'          VCnam=list(col="red", font=4, cex=2),
#'          VarLab=list(list(col="blue", font=3, cex=2), list(cex=1.25, srt=-15)))
#' 
#' # create variability-chart of the example dataset in the CLSI EP05-A2 
#' # guideline (listed on p.25)
#' data(Glucose,package="VCA")
#' varPlot(result~day/run, Glucose, type=3)
#' 
#' # use individual settings of 'VarLab' and 'VSpace' for each variance component
#' varPlot(result~day/run, Glucose, type=3, 
#'         VarLab=list(list(srt=45, col="red", font=2), 
#'         list(srt=90, col="blue", font=3)), VSpace=c(.25, .75))
#' 
#' # set individual titles for both plot when 'type=3'
#' # and individual 'ylim' specifications
#' varPlot(	result~day/run, Glucose, type=3, 
#'          Title=list(list(main="Variability Chart"), 
#'                     list(main="Plot of SD-Values")),
#' 			ylim=list(c(230, 260), c(0, 10)))
#' 
#' # more complex experimental design
#' data(realData)
#' Data <- realData[realData$PID == 1,]
#' varPlot(y~lot/calibration/day/run, Data, type=3)
#' 
#' # improve visual appearance of the plot
#' varPlot(y~lot/calibration/day/run, Data, type=3, keep.order=FALSE,
#'         BG=list(var="calibration", col=c("white", "lightgray")))
#' 
#' # add horizontal lines indicating mean-value for each factor-level of all variables
#' varPlot(y~lot/calibration/day/run, Data, type=3, keep.order=FALSE,
#'         BG=list(var="calibration", 
#'                 col=c("aquamarine","antiquewhite2","antiquewhite4",
#'                       "antiquewhite1","aliceblue","antiquewhite3",
#'                       "white","antiquewhite","wheat" ), 
#'                 col.table=TRUE),
#'         MeanLine=list(var=c("lot", "calibration", "day", "int"), 
#'                       col=c("orange", "blue", "darkgreen", "yellow"), 
#'                       lwd=c(2,2,2,2)))
#' 
#' # now also highlight bounds between factor levels of "lot" and "day" 
#' # as vertical lines and extend them into the table (note that each 
#' # variable needs its specific value for 'col.table')
#' varPlot(y~lot/calibration/day/run, Data, type=3, keep.order=FALSE,
#'         BG=list(var="calibration", 
#'                 col=c("aquamarine","antiquewhite2","antiquewhite4",
#'                       "antiquewhite1","aliceblue","antiquewhite3",
#'                       "white","antiquewhite","wheat" ), 
#'                  col.table=TRUE),
#'         MeanLine=list(var=c("lot", "calibration", "day", "int"), 
#'                       col=c("orange", "blue", "darkgreen", "yellow"), 
#'                       lwd=c(2,2,2,2)),
#'         VLine=list(var=c("lot", "day"), col=c("black", "skyblue1"),
#'                    lwd=c(2, 1), col.table=c(TRUE, TRUE)))
#' 
#' # one can use argument 'JoinLevels' to join factor-levels or a variable
#' # nested within a higher-level factor, 'VLine' is used to separate levels
#' # of variables "calibration" and "lot" with different colors
#'  varPlot(y~calibration/lot/day/run, Data, 
#'          BG=list(var="calibration", 
#'                  col=c("#f7fcfd","#e5f5f9","#ccece6","#99d8c9",
#'                        "#66c2a4","#41ae76","#238b45","#006d2c","#00441b"), 
#'                  col.table=TRUE), 
#'          VLine=list(var=c("calibration", "lot"), 
#'                     col=c("black", "darkgray"), lwd=c(2,1), col.table=TRUE), 
#'          JoinLevels=list(var="lot", col=c("#ffffb2","orangered","#feb24c"), 
#'                          lwd=c(2,2,2)), 
#'          MeanLine=list(var="lot", col="blue", lwd=2))
#' 
#' # same plot demonstrating additional features applicable via 'Points' 
#'  varPlot(y~calibration/lot/day/run, Data, 
#'          BG=list(var="calibration", 
#'                  col=c("#f7fcfd","#e5f5f9","#ccece6","#99d8c9",
#'                        "#66c2a4","#41ae76","#238b45","#006d2c","#00441b"), 
#'                  col.table=TRUE), 
#'          VLine=list(var=c("calibration", "lot"), 
#'                     col=c("black", "mediumseagreen"), lwd=c(2,1), 
#'                     col.table=c(TRUE,TRUE)), 
#'          JoinLevels=list(var="lot", col=c("lightblue", "cyan", "yellow"), 
#'                          lwd=c(2,2,2)), 
#'          MeanLine=list(var="lot", col="blue", lwd=2),
#'          Points=list(pch=list(var="lot", pch=c(21, 22, 24)), 
#'                      bg =list(var="lot", bg=c("lightblue", "cyan", "yellow")), 
#'                      cex=1.25))
#' 
#' # depict measurements as boxplots
#' datS5 <- subset(VCAdata1, sample==5)
#' varPlot(y~device/day, datS5, Boxplot=list()) 
#' 
#' # present points as jitter-plot around box-center
#' varPlot(	y~device/day, datS5, 
#'          Boxplot=list(jitter=1, col.box="darkgreen"),
#' 			BG=list(var="device", col=paste0("gray", c(60, 70, 80)),
#' 					col.table=TRUE),
#'          Points=list(pch=16, 
#' 						col=list(var="run", col=c("blue", "red"))), 
#' 			Mean=list(col="black", cex=1, lwd=2)) 
#' legend( "topright", legend=c("run 1", "run 2"),
#'         fill=c("blue", "red"), box.lty=0, border="white")
#' }

varPlot <- function(form, Data, keep.order=TRUE, 
					type=c(1L, 2L, 3L)[1], VARtype="SD", htab=.5,
					Title=NULL, VSpace=NULL,
					VarLab=list(cex=.75, adj=c(.5, .5)),
					YLabel=list(text="Value", side=2, line=3.5, cex=1.5),
					SDYLabel=list(side=2, line=2.5),
					Points=list(pch=16, cex=.5, col="black"),
					SDs=list(pch=16, col="blue", cex=.75),
					SDline=list(lwd=1, lty=1, col="blue"),
					BG=list(border="lightgray", col.table=FALSE),
					VLine=list(lty=1, lwd=1, col="gray90"),
					HLine=NULL,
					Join=list(lty=1, lwd=1, col="gray"),
					JoinLevels=NULL,
					Mean=list(pch=3, col="red", cex=.5),
					MeanLine=NULL,
					Boxplot=NULL,
					VCnam=list(cex=.75, col="black", line=0.25),
					useVarNam=FALSE, 
					ylim=NULL, max.level=25, ...)
{        
	stopifnot(is.data.frame(Data))
	stopifnot(class(form) == "formula")
	stopifnot(nrow(Data) > 2)                                               # at least 2 observations for estimating a variance
	
	args <- list(...)
	
	VARtype <- match.arg(toupper(VARtype), c("SD", "CV"))
	
	form <- terms(form, simplify=TRUE)
	stopifnot(attr(form, "response") == 1)
	resp <- rownames(attr(form, "factors"))[1]
	
	IntVal <- mean(Data[,resp], na.rm=TRUE)									# grand mean, possibly needed for MeanLine-functionality
	
	nest <- attr(form, "term.labels")
	if(length(nest) == 0)
		stop("Models, where the intercept is the only coefficient, are not supported!")
	
	nest <- sapply(nest, function(x){
				x <- unlist(strsplit(x, ":"))
				return(x[length(x)])
			})
	
	tab <- table(nest)
	if(any(tab > 1))
	{
		ind <- which(tab > 1)
		warning("Term(s) ", paste("'", nest[ind], "'", sep=""), " occurring multiple times as last sub-term in: ",paste(paste("'",names(nest), "'", sep=""), collapse=", "),". It will be kept only once!")
		nest <- unique(nest)
	}
	
	Nvc  <- length(nest)
	
	if(!is.null(Title))
	{
		stopifnot(class(Title) == "list")
		
		if(class(Title[[1]]) == "list")
			stopifnot( class(Title[[2]]) == "list" )
		else																# replicate title settings for both possible plots
			Title <- list(Title, Title)
	}
	
	if( is.null(VSpace) || !is.numeric(VSpace) || sum(VSpace) != 1 )
		VSpace <- rep(1/Nvc, Nvc)                					 		# fraction of the table reserved for vertical names is automatically set (e.g. run/part)
	
	
	
	if( !is.null(BG) && is.null(BG$var) )
		BG$var <- nest[1]                                   				# use top-level factor for separating factor-levels if not otherwise specified 
	
	if(nest[1] != "1")                                         				# travers nesting structure and build nested list
	{
		lst <- buildList(Data=Data, Nesting=nest, Current=nest[1], resp=resp,
				keep.order=keep.order, useVarNam=useVarNam, Points=Points)
	}
	else
	{
		lst <- as.list(Data[,resp])
		
		names(lst) <- paste("rep", 1:nrow(Data), sep="")
		attr(lst, "Nelem") <- rep(1, length(lst))
		
		# point colors used to indicated class levels?
		
		if( (is.list(Points$col)) && ("var" %in% names(Points$col)) && (Points$col$var %in% colnames(Data)) )
		{
			if(length(unique(Data[,Points$col$var])) > length(unique(Points$col$col)))          # not enough colors, need to be recycled
			{
				Points$col$col <- rep(Points$col$col, ceiling(length(unique(Data[,Points$col$var])) / length(unique(Points$col$col))))
			}
			tmp <- 1:length(unique(Data[,Points$col$var]))
			names(tmp) <- unique(Data[,Points$col$var])
			attr(lst, "Color") <- Points$col$col[tmp[Data[,Points$col$var]]]                   # assign colors to class-levels of Points$col$var variable
		}
		else
		{
			attr(lst, "Color") <- rep(Points$col, nrow(Data))
		}
		
		# plotting symbols used to indicated class levels?
		
		if( (is.list(Points$pch)) && ("var" %in% names(Points$pch)) && (Points$pch$var %in% colnames(Data)) )
		{
			if(length(unique(Data[,Points$pch$var])) > length(unique(Points$pch$pch)))          # not enough colors, need to be recycled
			{
				Points$pch$pch <- rep(Points$pch$pch, ceiling(length(unique(Data[,Points$pch$var])) / length(unique(Points$pch$pch))))
			}
			tmp <- 1:length(unique(Data[,Points$pch$var]))
			names(tmp) <- unique(Data[,Points$pch$var])
			attr(lst, "Symbol") <- Points$pch$pch[tmp[Data[,Points$pch$var]]]                   # assign colors to class-levels of Points$col$var variable
		}
		else
		{
			attr(lst, "Symbol") <- rep(Points$pch, nrow(Data))
		}
		
		# plotting symbol backgrounds used to indicated class levels?
		
		if( (is.list(Points$bg)) && ("var" %in% names(Points$bg)) && (Points$bg$var %in% colnames(Data)) )
		{
			if(length(unique(Data[,Points$bg$var])) > length(unique(Points$bg$bg)))          # not enough colors, need to be recycled
			{
				Points$bg$bg <- rep(Points$bg$bg, ceiling(length(unique(Data[,Points$bg$var])) / length(unique(Points$bg$bg))))
			}
			tmp <- 1:length(unique(Data[,Points$bg$var]))
			names(tmp) <- unique(Data[,Points$bg$var])
			attr(lst, "BG") <- Points$bg$bg[tmp[Data[,Points$bg$var]]]                   # assign colors to class-levels of Points$col$var variable
		}
		else
		{
			attr(lst, "BG") <- rep(Points$bg, nrow(Data))
		}
		
		# plotting symbol size used to indicated class levels?
		
		if( (is.list(Points$cex)) && ("var" %in% names(Points$cex)) && (Points$cex$var %in% colnames(Data)) )
		{
			if(length(unique(Data[,Points$cex$var])) > length(unique(Points$cex$cex)))          # not enough colors, need to be recycled
			{
				Points$cex$cex <- rep(Points$cex$cex, ceiling(length(unique(Data[,Points$cex$var])) / length(unique(Points$cex$cex))))
			}
			tmp <- 1:length(unique(Data[,Points$cex$var]))
			names(tmp) <- unique(Data[,Points$cex$var])
			attr(lst, "CEX") <- Points$cex$cex[tmp[Data[,Points$cex$var]]]                   # assign colors to class-levels of Points$col$var variable
		}
		else
		{
			attr(lst, "CEX") <- rep(Points$cex, nrow(Data))
		}
	}
	
	Nelem <- attr(lst, "Nelem")                                             # number of bottom-level sub-classes per top-level factor-level  
	Nbin <- attr(lst, "Nbin")                                               # number of bins, i.e. scatterplots
	Xdiff <- length(lst)/sum(Nelem)                                         # basic cell-width of the tabular        
	
	Range 		<- range(Data[, resp], na.rm=TRUE)                          # plotting limits (Y-direction) for the variability-plot    
	Range[1] 	<- Range[1] - diff(Range) * .1
	Range[2] 	<- Range[2] + diff(Range) * .1
	Range 		<- range(pretty(Range))
	Range2		<- NULL
		
	if( !is.null(ylim) )                                                    # user-defined Y-limits
	{
		if (type == 1) {
			if( is.numeric(ylim) && length(ylim) > 1)
				Range <- ylim[1:2]
		} else if(type == 2) {
			Range2 <- ylim[1:2]
		} else {
			if(is.list(ylim))
			{
				if( is.numeric(ylim[[1]]) && length(ylim[[1]]) > 1)
					Range <- ylim[[1]][1:2]
				
				if( is.numeric(ylim[[2]]) && length(ylim[[2]]) > 1)
					Range2 <- ylim[[2]][1:2]
				
				ylim <- ylim[[1]]
			} else {						# assume that only variability chart is meant
				Range <- ylim[1:2]
			}
		}
	}
	
	YLim <- Range
	
	YLim[1] <- YLim[1] - diff(YLim) * htab
	
	if(VARtype == "SD")
		SDrange <- attr(lst, "SDrange")                                     # plotting limits (Y-direction) for the SD-plot
	if(VARtype == "CV")
		SDrange <- attr(lst, "CVrange")                                     # use CV instead of SD
	
	SDrange[1] 	<- SDrange[1] - diff(SDrange) * .05
	SDrange[2] 	<- SDrange[2] + diff(SDrange) * .1
	SDrange 	<- range(pretty(SDrange))
	SDYLim 		<- SDrange
	if(!is.null(Range2))
	{
		SDrange <- Range2
		SDYLim  <- Range2
	} else {
		Range2 <- SDrange
	}
	
	SDYLim[1] 	<- SDYLim[1] - diff(SDYLim) * htab

	if(type == 3L)
		MFROW <- c(2,1)
	#else
	#	MFROW <- c(1,1)

	if(!is.null(VCnam))
	{
		VCnam.default <- list(cex=.75, col="black", line=0.35, side=2)
		VCnam.default[names(VCnam)] <- VCnam
		VCnam <- VCnam.default
		if(!VCnam$side %in% c(2, 4))
			VCnam$side <- 2
		VCnam$adj <- ifelse(VCnam$side == 2, 1, 0)
	}

	mar <- c(	1,
				5.1, 
				ifelse(is.null(Title), 1, 3.5),
				ifelse(VCnam$side == 2, 1, 5.1))
	
#	if(is.null(Title))
#		mar <- c(1,4.1,1,1) 
#	else
#		mar <- c(2,4.1,3.5,1)
	
	if(!is.null(args))
	{
		if(!"mar" %in% names(args))
			args[["mar"]] <- mar
		args[["xaxs"]]  <- "i"
		args[["yaxs"]]  <- "i"
		
		if(type == 3L)						# only if two plots are requested at once, otherwise it may interfer with 'layout'
			args[["mfrow"]] <- MFROW
	}
	else
	{
		if(type == 3L)
			args <- list(mar=mar, xaxs="i", yaxs="i", mfrow=MFROW)
		else
			args <- list(mar=mar, xaxs="i", yaxs="i")	
	}
	
	old.par <- do.call("par", args)											# set graphical parameters
	
	Ybound   <- YLim[1]														# generate Y-limits for the tabular environment(s)
	SDYbound <- SDYLim[1]
	
	for(i in 1:Nvc)
	{
		Ybound   <- c(Ybound,     Ybound[length(  Ybound)] + VSpace[i] * abs(diff(c(Range[1], YLim[1]))) )
		SDYbound <- c(SDYbound, SDYbound[length(SDYbound)] + VSpace[i] * abs(diff(c(SDrange[1], SDYLim[1]))) )
	}	
	
	Xbound <- c(0, cumsum(Xdiff * Nelem))
	
	VarLab.default <- list(cex=.75, adj=c(.5, .5))
	tmpList <- list()
	for(i in 1:Nvc)																			# replicate default-settings list as many times as there are VCs
		tmpList[[i]] <- VarLab.default
	VarLab.default <- tmpList
	
	if(any(sapply(VarLab, class) == "list"))												# VarLab consists of element(s) which is/are list(s)
	{
		for(i in 1:length(VarLab.default))
		{
			if(class(VarLab[[i]]) == "list")
				VarLab.default[[i]][names(VarLab[[i]])] <- VarLab[[i]]
		}
	}
	else																					# default settings used for each VC (no user-specification as list of lists)
	{
		for(i in 1:length(VarLab.default))
		{
			VarLab.default[[i]][names(VarLab)] <- VarLab
		}
	}
	VarLab <- VarLab.default
	
	YLabel.default <- list(text="Value", side=2, line=3.5, cex=1.5, at=mean(Range))
	YLabel.default[names(YLabel)] <- YLabel
	YLabel <- YLabel.default
	
	SDYLabel.default <- list(text=ifelse(VARtype=="SD", "SD", "CV %"), side=2, line=2, at=mean(SDrange))
	SDYLabel.default[names(SDYLabel)] <- SDYLabel
	SDYLabel <- SDYLabel.default
	
	Points.default <- list(pch=16, cex=.5, col="black")
	Points.default[names(Points)] <- Points
	Points <- Points.default
	
	if(!is.null(Boxplot))
	{
		Boxplot.default <- list( col.box="gray65", col.median="white",
								 col.whiskers="gray65", lwd=2, jitter=1e3)
		Boxplot.default[names(Boxplot)] <- Boxplot
		Boxplot <- Boxplot.default
		if(!is.null(Boxplot$jitter) && Boxplot$jitter < 5)
			Boxplot$jitter <- 5
	}
	
	SDs.default <- list(pch=16, col="blue", cex=.75)
	SDs.default[names(SDs)] <- SDs
	SDs <- SDs.default
	
	if(!is.null(BG))
	{
		BG.default <- list(border=NA, col.table=FALSE, col.bg=TRUE)
		BG.default[names(BG)] <- BG
		BG <- BG.default
		
		if(!is.null(BG$col) && !is.null(BG$var))			# color(s) and variable specified
		{
			Lvls <- unique(Data[,BG$var])
			
			if(length(BG$col) == length(Lvls))
			{
				names(BG$col) <- Lvls						# same color will be used for same level
			}
		}
	}
	
	if(!is.null(SDline))
	{
		SDline.default <- list(lwd=1, lty=1, col="blue")
		SDline.default[names(SDline)] <- SDline
		SDline <- SDline.default
	}
	
#	if(!is.null(VCnam))
#	{
#		VCnam.default <- list(cex=.75, col="black", line=0.25, side=2)
#		VCnam.default[names(VCnam)] <- VCnam
#		VCnam <- VCnam.default
#	}
	
	if(!is.null(Join))
	{
		Join.default <- list(lty=1, lwd=1, col="gray")
		Join.default[names(Join)] <- Join   
		Join <- Join.default
	}
	
	if (!is.null(Mean)) {
		Mean.default <- list(pch = 3, col = "yellow", cex = 0.35)
		Mean.default[names(Mean)] <- Mean
		Mean <- Mean.default
	}
	
	if(!is.null(JoinLevels) && is.list(JoinLevels) && JoinLevels$var %in% nest[2:length(nest)])
	{
		JoinLevelsLevels <- unique(as.character(Data[, JoinLevels$var]))
		if(is.null(JoinLevels$col) || length(JoinLevels$col) < length(JoinLevelsLevels))
		{
			JoinLevels$col <- sample(colors()[1:length(JoinLevelsLevels)])
			warning("Random selection of colors used because too few user-specified colors detected!")
		}
		JLcolorMap <- JoinLevels$col
		names(JLcolorMap) <- JoinLevelsLevels
		JoinLevels$col <- NULL									# information not needed any more
	}
	
	# environment for recording plotting coordinates and used as interims storage for plotting information
	rec.env <- environment()
	rec.env$dat <- list()
	rec.env$VLineCollection <- list()							# also use it for collecting vertical line information
	rec.env$JoinLevelsCollection <- list()						# and for joining factor level means
	rec.env$JoinCollection <- list()
	rec.env$MeanLineCollection <- list()
	rec.env$MeansCollection <- list()
	rec.env$PointCollection <- list()
	
	#abline(h=Ybound[Nvc+1])
	
	tlNames <- names(lst)
	
	global.Points <- list()
	
	# adding tabluar environment at the bottom and scatterplot of observed values
	#
	# lst       	... (list) to be processed 
	# xlim      	... (numeric) vector specifying the X-limits for the current list
	# yval      	... (numeric) vector specifying the Y-coordinates of horizontal lines
	# Xdiff     	... (numeric) value specifying the width of a base cell (cell corresponding to a node)
	# type      	... (integer) see 'type' above
	# VARtype   	... (character) see 'VARtype' above
	# StatsList 	... (list) for descriptive statistics
	# VarLab		... (list of lists) specifying the labelling style in the tabular-environment
	# MeanLine		... (list) of function arguments specifying factor-levels means
	# VLine			... (list) of function arguments for vertically separating factor-levels
	# Intercept 	... (list) of function arguments specifying the horizontal intercept-line
	# JoinLevels	... (list) of functin arguments specifying lines joining factor-level means, nested within higher order factor
	# draw.vertical ... (logical) indicating whether vertical lines in the table should be drawn or not, depends on 'max.level'
	
	processList <- function(lst, xlim, yval, Xdiff, type, VARtype, StatsList, index=1, 
			VarLab, MeanLine, VLine, Intercept, JoinLevels, draw.vertical=TRUE)
	{
		Xlower <- xlim[1]    
		
		Stats  <- attr(lst, "Stats")
		Nelem  <- attr(lst, "Nelem")
		Nlevel <- attr(lst, "Nlevel")
		
		Names <- names(lst)
		
		col.table <- FALSE
		
		if( !is.null(BG))                                                           # draw vertical lines between or color background for factor-levels
		{
			if( BG$var == attr(lst, "factor") )                                  	# apply BG-list to the current factor
			{
				tmpBG <- BG
				ColNames <- names(BG$col)											# might be NULL
				
				col.table <- tmpBG$col.table										# color respective row in table below?
				col.bg    <- tmpBG$col.bg											# color the background?
				tmpBG$col.table <- tmpBG$col.bg <- NULL
			}               
			else
				tmpBG <- NULL
		}
		else																		# background will not be colored
		{
			tmpBG <- NULL
			
			if(!is.null(Intercept))													# horizontal intercept line will not be covered by 'BG'
				do.call("lines", Intercept)
		}
		
		
		Factor <- attr(lst, "factor")												# current factor-variable to be processed
		
		for(i in 1:length(lst))
		{  
			tmp <- lst[[i]]    
			tmp.draw.vertical <- Nlevel[i] <= max.level								# no further sub-division if too many levels
			
			Xupper <- Xlower + Xdiff * Nelem[i]
			
			if( !is.null(tmpBG) )                                                   # apply BG-settings
			{
				tmpBG$factor    <- NULL
				tmpBG$xleft     <- Xlower
				tmpBG$xright    <- Xupper
				
				if(type %in% c(1,3))
				{
					tmpBG$ybottom   <- Range[1]										# outside of 'processList' defined
					tmpBG$ytop      <- Range[2]
				}
				else if(type %in% c(2,3))
				{
					tmpBG$ybottom   <- SDYbound[length(SDYbound)] 
					tmpBG$ytop      <- SDYLim[2] * 2
				}
				
				if( length(BG$col) > 1 )
				{
					if(!is.null(ColNames))
						tmpBG$col <- BG$col[Names[i]]
					else
						tmpBG$col <- BG$col[index]
					last.index <- index
				}
				
				tmpBG$var <- NULL													# "var" is no graphical parameter
				if(col.bg)
					do.call("rect", tmpBG)                                          # draw the rectangle as background for the specified group
				
				if(col.table)														# color levels for the factor determining BG-colors
				{
					tmpBG$ybottom <- yval[1]
					tmpBG$ytop    <- yval[2]
					
					if(is.na(tmpBG$border))											# ensure that border is always drawn in the tabular part
					{						
						tmpBG$border <- "black"
						do.call("rect", tmpBG)
						tmpBG$border <- NA
					}
					else
						do.call("rect", tmpBG)
				}
				
				if( length(BG$col) > 1 )
				{
					if( index + 1 > length(BG$col) )
						index <- 1
					else
						index <- index + 1                                          # use color (i+1) or go back to 1st color
				}
			}
			
			if(draw.vertical || i == length(lst))
				lines(x=rep(Xupper, 2), y=yval[1:2])                                # draw vertical lines of the table       
			
			if(is.list(tmp))                                                        # recursive descent
			{	
				VarLabel <- VarLab[[1]]												# use different variable name because VarLab is passed down and is not allowed to be manipulated
				
				VarLabel$x=mean(c(Xlower, Xupper))
				VarLabel$y=mean(yval[1:2])
				VarLabel$labels=Names[i]
				
				if(draw.vertical)
					do.call("text", args=VarLabel)
				
				if(!tmp.draw.vertical)
				{
					text(mean(c(Xlower, Xupper)),  mean(yval[2:3]), 
							paste("N=", Nlevel[i], sep=""), cex=.75)
				}
				
				StatsList <- processList(	lst=tmp, xlim=c(Xlower, Xupper), yval=yval[-1], Xdiff=Xdiff, type=type, 
											VARtype=VARtype, StatsList=StatsList, index=index, VarLab=VarLab[-1],
											MeanLine=MeanLine, VLine=VLine, Intercept=Intercept, JoinLevels=JoinLevels,
											draw.vertical=tmp.draw.vertical)
				index <- attr(StatsList, "index")
			}
			else                                                                    # leaf-node reached (numeric vector)
			{
				VarLabel <- VarLab[[1]]
				VarLabel$x=mean(c(Xlower, Xupper))
				VarLabel$y=mean(yval[1:2])
				VarLabel$labels=Names[i]
				
				if(draw.vertical)
					do.call("text", args=VarLabel)
				
				if(type == 1)
				{
					Points$x=rep(mean(mean(c(Xlower, Xupper))), length(tmp))        # add points to plot (X-coordinates)
					
					stylim <- tmp < ylim[1]											# smaller than ylim
					gtylim <- tmp > ylim[2]
					
					tmp.points <- tmp
					
					if(any( stylim | gtylim ))
					{
						tmp.points[stylim | gtylim] <- NA
					}
					
					Points$y <- tmp.points						
					
					rec.env$dat$x <- c(rec.env$dat$x, Points$x)						# record coordinates for later use
				}
				else
				{
					SDs$x=mean(c(Xlower, Xupper))
					SDs$y=ifelse(VARtype=="SD", Stats$SD[i], Stats$CV[i])
					StatsList$SDvec <- c(StatsList$SDvec, SDs$y)
					names(StatsList$SDvec)[length(StatsList$SDvec)] <- SDs$x
				}
				
				StatsList$Nobs <- c(StatsList$Nobs, length(na.omit(tmp)))           # record number of non-NA observation
				
				if(!is.null(Join) && type == 1)
				{
					Join$x=Points$x
					
					if(any( stylim | gtylim ))										# set extreme values to ylim
					{
						tmp.points <- tmp
						
						if(any(stylim))
							tmp.points[stylim] <- ylim[1]
						if(any(gtylim))
							tmp.points[gtylim] <- ylim[2]
						
						Join$y <- tmp.points
					}
					else
						Join$y=Points$y
					
					#do.call("lines", args=Join)
					rec.env$JoinCollection[[length(rec.env$JoinCollection)+1]] <- Join
				}
				
				if(!is.null(Mean) && type == 1)
				{
					Mean$x=Points$x[1]
					Mean$y=Stats$Mean[i]
					mstylim <- Mean$y < ylim[1]
					mgtylim <- Mean$y > ylim[2]
					if(any(mstylim))												# do plot mean-value if outside the plotting region
						Mean$y[mstylim] <- NA
					if(any(mgtylim))
						Mean$y[mgtylim] <- NA
					
					#do.call("points", args=Mean)
					rec.env$MeanCollection[[length(rec.env$MeanCollection) + 1]] <- Mean
				}
				
				Points$col <- attr(tmp, "Color")
				Points$pch <- attr(tmp, "Symbol")
				Points$bg  <- attr(tmp, "BG")
				Points$cex <- attr(tmp, "CEX")
				
				if(type == 1)
				{
					#do.call("points", args=Points)
					rec.env$PointsCollection[[length(rec.env$PointsCollection) + 1]] <- Points
				}
				else
				{
					if(!is.null(Range2))					# SD-points outside ylim will not be plotted 
					{
						SD.idx <- which(SDs$y < Range2[1])
						if(length(SD.idx) > 0)
							{
								SDs$x <- SDs$x[-SD.idx]
								SDs$y <- SDs$y[-SD.idx]
							}
					}
		
					do.call("points", args=SDs)
				}
				
#################
				if(!is.null(Boxplot))
					rec.env$BoxCollection <- Boxplot
			}
			
			if( !is.null(MeanLine) && !is.null(Factor) && Factor %in% MeanLine$var)
			{
				var.ind <- which(MeanLine$var == Factor)
				MeanLineClone <- MeanLine
				for(e in 1:length(MeanLineClone))
				{
					MeanLineClone[[e]] <- MeanLineClone[[e]][var.ind]				# only keep parameter values for current factor-variable
				}
				MeanLine.default <- list(lty=1, lwd=1, col="red")
				MeanLine.default[names(MeanLineClone)] <- MeanLineClone
				MeanLineClone <- MeanLine.default
				
				tmp.Mean     <- MeanLineClone
				if("mar" %in% names(MeanLineClone) && is.numeric(MeanLineClone$mar) && MeanLineClone$mar >= 0 && MeanLineClone$mar < .5 )
				{
					tmpXdiff   <- diff(c(Xlower, Xupper))	
					tmp.Mean$x <- c(Xlower+MeanLineClone$mar*tmpXdiff, Xupper-MeanLineClone$mar*tmpXdiff)
				}
				else
					tmp.Mean$x <- c(Xlower, Xupper)
				tmp.Mean$y   <- rep(Stats$Mean[i], 2)
				tmp.Mean$var <- NULL
				
				Mstylim <- tmp.Mean$y < ylim[1]
				Mgtylim <- tmp.Mean$y > ylim[2]
				
				if(any(Mstylim))
					tmp.Mean$y[Mstylim] <- NA
				if(any(Mgtylim))
					tmp.Mean$y[Mgtylim] <- NA
				
				#do.call("lines", tmp.Mean)
				rec.env$MeanLineCollection[[length(rec.env$MeanLineCollection) + 1]] <- tmp.Mean
			}
			
			nestInd <- which(nest == Factor)
			
			if( !is.null(VLine) && !is.null(Factor) && Factor %in% VLine$var)		# skip last vertical line
			{
				drawVLine <- TRUE
				
				if(nestInd > 1)														# not upper-most factor
				{
					if(nest[nestInd-1] %in% VLine$var)								# preceding factor in nesting structure also specified
					{
						drawVLine <- i < length(lst)
					}
				}
				else
				{
					if(i == length(lst))
						drawVLine <- FALSE											# omit last vertical line for upper-most factor
				}
				
				if(drawVLine)														# do not draw last vertical line, this separates levels of the factor one level above
				{
					var.ind <- which(VLine$var == Factor)
					VLineClone <- VLine
					for(e in 1:length(VLineClone))									# only keep parameter values for current factor-variable
					{
						VLineClone[e] <- VLineClone[[e]][var.ind]
						if(is.na(VLineClone[e]))									# not enough arguments provided, use last value
							VLineClone[e] <- tail(VLineClone[[e]])
					}
					
					VLine.default <- list(var=nest[1], lty=1, lwd=1, col="gray90")
					VLine.default[names(VLineClone)] <- VLineClone
					VLineClone <- VLine.default
					
					tmpVLine <- VLineClone
					tmpVLine$x <- c(Xupper, Xupper)
					if("col.table" %in% names(tmpVLine) && is.logical(tmpVLine$col.table) && isTRUE(tmpVLine$col.table))
					{
						if(type %in% c(1, 3))
							tmpVLine$y <- c(yval[1], Range[2]+abs(diff(Range)))
						if(type %in% c(2, 3))
							tmpVLine$y <- c(yval[1], SDYLim[2]*2)
					}
					else
					{
						if(type %in% c(1, 3))
							tmpVLine$y <- c(Range[1], Range[2]+abs(diff(Range)))
						if(type %in% c(2, 3))
							tmpVLine$y <- c(SDrange[1], SDYLim[2]*2)
					}
					tmpVLine$col.table <- NULL
					tmpVLine$var <- NULL
					
					#do.call("lines", tmpVLine)
					rec.env$VLineCollection[[length(rec.env$VLineCollection) + 1]] <- tmpVLine
				}
			}
			
			
			precFactor <- NULL
			if(nestInd > 1)
				precFactor <- nest[nestInd-1]											# preceding factor in nesting hierarchy
			
			if( !is.null(JoinLevels) && !is.null(Factor) && !is.null(precFactor) && Factor %in% JoinLevels$var)	# joining mean-values of current factor-variable
			{
				if(length(rec.env$JoinLevelsCollection[[Names[i]]]) == 0)
				{
					JoinLevelsClone <- JoinLevels
					
					JoinLevels.default <- list(lty=1, lwd=1)
					JoinLevels.default[names(JoinLevelsClone)] <- JoinLevelsClone
					JoinLevelsClone <- JoinLevels.default
					JoinLevelsClone$col <- JLcolorMap[Names[i]]							# use correct color for current factor-level
					
					JoinLevelsClone$x <- mean(c(Xlower, Xupper))
					JoinLevelsClone$y <- Stats$Mean[i]
					JoinLevelsClone$var <- NULL
					
					rec.env$JoinLevelsCollection[[Names[i]]] <- JoinLevelsClone
				}
				else																	# only append X- and Y-coordinates
				{
					rec.env$JoinLevelsCollection[[Names[i]]]$x <- c(rec.env$JoinLevelsCollection[[Names[i]]]$x, mean(c(Xlower, Xupper)))
					rec.env$JoinLevelsCollection[[Names[i]]]$y <- c(rec.env$JoinLevelsCollection[[Names[i]]]$y, Stats$Mean[i])
					rec.env$JoinLevelsCollection[[Names[i]]]$col <- JLcolorMap[Names[i]]
				}
			}
			
			if(!is.null(tmpBG))														# all BG-drawing finished
			{
				if(!is.null(Intercept))												# horizontal intercept line will not be covered by 'BG'
					do.call("lines", Intercept)
			}
			
			Xlower <- Xupper
		}
		abline(h=yval[2])
		
		attr(StatsList, "index") <- index                                           # remember last values of index
		invisible(StatsList)
	}
	
	MeanLineOnTop <- FALSE															# default
	
	if(!is.null(MeanLine$top))
	{
		MeanLineOnTop <- MeanLine$top
		MeanLine$top  <- NULL 
	}
	
	InterceptMeanLine <- NULL
	
	if(is.list(MeanLine) && "int" %in% MeanLine$var)								# draw horizontal line indicating the intercept if specified
	{		
		var.ind <- which(MeanLine$var == "int")
		
		MeanLineClone <- MeanLine
		for(e in 1:length(MeanLineClone))
		{
			MeanLineClone[[e]] <- MeanLineClone[[e]][var.ind]						# only keep parameter values for current factor-variable
			
			if(is.na(MeanLineClone[[e]]))											# remove if not set by the user (defaults will be applied)
				MeanLineClone[[e]] <- NULL
		}
		Intercept <- list(lty=3, lwd=2, col="gray")									# default settings
		Intercept[names(MeanLineClone)] <- MeanLineClone
		Intercept$var <- NULL
		
		if("mar" %in% names(Intercept) && is.numeric(Intercept$mar) && Intercept$mar >= 0 && Intercept$mar < .5 )
		{
			Xlower <- Xbound[1]
			Xupper <- Xbound[2]
			tmpXdiff   <- diff(c(Xlower, Xupper))	
			Intercept$x <- c(Xlower+MeanLineClone$mar*tmpXdiff, Xupper-Intercept$mar*tmpXdiff)
		}
		else
			Intercept$x <- range(Xbound)
		
		Intercept$y <- rep(IntVal, 2)	
		
		if(MeanLineOnTop)															# plot it on top of measurements
		{
			InterceptMeanLine <- Intercept
			Intercept <- NULL
		}
	}
	else
		Intercept <- NULL
	
	if(type %in% c(1,3))
	{
		plot(1, axes=FALSE, xlab="", ylab="", ylim=YLim, xlim=c(0, length(lst)), type="n")
		
#		if(!is.null(Grid))
#		{			
#			USR	<- par("usr")
#			Grid.def <- list(v=pretty(USR[1:2]), h=pretty(USR[3:4]), lty=3, col="gray90")
#			Grid.def[names(Grid)] <- Grid
#			Grid <- Grid.def
#			do.call("abline", Grid)			
#			rect(USR[1], USR[3], USR[2], Range[1], col="white", border="white")
#		}
#		
		tmp.at <- pretty(Range)
		if(any(tmp.at < min(Range)))
			tmp.at <- tmp.at[-which(tmp.at < min(Range))]
		axis(2, at=tmp.at, las=1) 
		#axis(2, at=pretty(Range))
		
		do.call("mtext", args=YLabel)
		
		if(!is.null(Title))
			do.call("title", args=Title[[1]])
		
		StatsList <- list(SDvec=numeric(), Nobs=integer())
		
		StatsList <- processList(lst=lst, xlim=range(Xbound), yval=Ybound, Xdiff=Xdiff, type=1, 
				StatsList=StatsList, VARtype=VARtype, VarLab=VarLab, MeanLine=MeanLine, 
				Intercept=Intercept, VLine=VLine, JoinLevels=JoinLevels)   		  		
		
		if(!is.null(HLine))
		{
			HLine.def1 <- list(h=tmp.at, col="gray90", lty=3)
			if(type == 3 && class(HLine[[1]]) == "list")
				HLine.def1[names(HLine[[1]])] <- HLine[[1]]
			else
				HLine.def1[names(HLine)] <- HLine
			HLine1 		<- HLine.def1
			HLine1$v 	<- NULL
			oor   		<- which(HLine1$h <= min(Range))
			if(length(oor) > 0)
				HLine1$h <- HLine1$h[-oor]
			do.call("abline", HLine1)
		}
		
		box()
		
		if(length(rec.env$VLineCollection) > 0)
		{
			for(i in 1:length(rec.env$VLineCollection))								# now actually draw vertical lines, avoiding BG-coloring to partially cover them
				do.call("lines", rec.env$VLineCollection[[i]])
			
			rec.env$VLineCollection <- NULL											# should not be returned
		}
		
		if(length(rec.env$JoinLevelsCollection) > 0)
		{
			for(i in 1:length(rec.env$JoinLevelsCollection))						# now actually draw lines joining factor-level means
				do.call("lines", rec.env$JoinLevelsCollection[[i]])
			
			rec.env$JoinLevelsCollection <- NULL
		}
		
		if(length(rec.env$JoinCollection) > 0)										# vertical lines joining observations
		{
			for(i in 1:length(rec.env$JoinCollection))
				do.call("lines", rec.env$JoinCollection[[i]])
			
			rec.env$JoinCollection <- NULL
		}
		
		if(!MeanLineOnTop)															# mean-lines in front
		{
			if(length(rec.env$MeanLineCollection) > 0)									# mean line per factor-level
			{
				for(i in 1:length(rec.env$MeanLineCollection))
					do.call("lines", rec.env$MeanLineCollection[[i]])
				
				rec.env$MeanLineCollection <- NULL
			}
		}
		
		if(length(rec.env$BoxCollection) > 0)
		{
			drawBox <- function(y, x, xdiff, col.box="gray65", col.median="white",
								col.whiskers="gray65", lwd=3, jitter=FALSE)
			{
				bp 		<- boxplot(y, plot=F)									# get stats
				Xlim 	<- x + c(-1, 1) * xdiff/3
				rect(	Xlim[1], bp$stats[2,1], Xlim[2], bp$stats[4,1], 		# box
						col=col.box, border=NA)
				lines(	Xlim, rep(bp$stats[3,1], 2), col=col.median, lwd=lwd)	# median
				lines(	rep(x, 2), bp$stats[1:2,1], col=col.whiskers, lwd=lwd)	# whiskers
				lines(	rep(x, 2), bp$stats[4:5,1], col=col.whiskers, lwd=lwd)
				lines(	x + c(-1, 1) * xdiff/5, rep(bp$stats[1,1], 2), 
						col=col.whiskers, lwd=lwd)
				lines(	x + c(-1, 1) * xdiff/5, rep(bp$stats[5,1], 2), 
						col=col.whiskers, lwd=lwd)
			}
			for(i in 1:length(rec.env$PointsCollection))
			{
				tmp <- c(list(	y=rec.env$PointsCollection[[i]]$y, 
						 		x=rec.env$PointsCollection[[i]]$x[1],
								xdiff=Xdiff),
						 rec.env$BoxCollection )

				 do.call("drawBox", tmp)
			}
			rec.env$BoxCollection <- NULL
		}
		
		if(length(rec.env$PointsCollection) > 0)									# observations
		{
			for(i in 1:length(rec.env$PointsCollection))
			{
				if(!is.null(Boxplot))
				{
					if(!is.null(Boxplot$jitter) && Boxplot$jitter)
					{
						jitter <- runif(length(rec.env$PointsCollection[[i]]$y),
										-Xdiff/Boxplot$jitter, Xdiff/Boxplot$jitter)
						jitter <- sample(jitter)
						rec.env$PointsCollection[[i]]$x <-  rec.env$PointsCollection[[i]]$x + jitter
					}
				}
				do.call("points", rec.env$PointsCollection[[i]])
			}
			rec.env$PointsCollection <- NULL
		}
		
		if(length(rec.env$MeanCollection) > 0)										# Mean values in factor-levels of variable one above error
		{
			for(i in 1:length(rec.env$MeanCollection))
				do.call("points", rec.env$MeanCollection[[i]])
			
			rec.env$MeanCollection <- NULL
		}
		
		if(MeanLineOnTop)
		{
			if(!is.null(InterceptMeanLine))
				do.call("lines", InterceptMeanLine)
			
			if(length(rec.env$MeanLineCollection) > 0)									# mean line per factor-level
			{
				for(i in 1:length(rec.env$MeanLineCollection))
					do.call("lines", rec.env$MeanLineCollection[[i]])
				
				rec.env$MeanLineCollection <- NULL
			}
		}
		
		if(!is.null(VCnam))
		{
			VCnam$at <- Ybound[-length(Ybound)]+diff(Ybound[1:2])/2
			if(is.null(VCnam$text))																			# do not overwrite user-specified text
				VCnam$text <- nest
			VCnam$las <- 1

			do.call("mtext", args=VCnam)
		}
	}
	
	if(type %in% c(2,3))                                                                        			# add a 2nd plot
	{
		plot(1, axes=FALSE, xlab="", ylab="", ylim=SDYLim, xlim=c(0, length(lst)), type="n")
				
		if(SDrange[1] == 0)
			abline(h=0, lty=2, col="gray")
		
		tmp.at <- pretty(Range2)
		if(any(tmp.at < min(Range2)))
			tmp.at <- tmp.at[-which(tmp.at < min(Range2))]
		
		axis(2, at=tmp.at, las=1) 
		
		abline(h=0, lty=2, col="gray")
		
		if(!is.null(HLine))
		{
			HLine.def2 <- list(h=tmp.at, col="gray90", lty=3)
			if(type == 3 && class(HLine[[2]]) == "list")
				HLine.def2[names(HLine[[2]])] <- HLine[[2]]
			else
				HLine.def2[names(HLine)] <- HLine
			HLine2 		<- HLine.def2
			HLine2$v 	<- NULL
			oor   		<- which(HLine2$h <= min(SDrange))
			if(length(oor) > 0)
				HLine2$h <- HLine2$h[-oor]

			do.call("abline", HLine2)
		}
		
		#axis(2, at=pretty(SDrange)[-1])   
		
		do.call("mtext", args=SDYLabel)
		
		if(!is.null(Title))
			do.call("title", args=Title[[2]])
		
		StatsList <- list(SDvec=numeric(), Nobs=integer())
		#SDvec <- numeric()
		StatsList <- processList(	lst=lst, xlim=range(Xbound), yval=SDYbound, Xdiff=Xdiff, type=2, StatsList=StatsList, 
									VARtype=VARtype, VarLab=VarLab, MeanLine=MeanLine, Intercept=Intercept, VLine=VLine,
									JoinLevels=JoinLevels)   
		
		if(length(rec.env$VLineCollection) > 0)
		{
			for(i in 1:length(rec.env$VLineCollection))								# now actually draw vertical lines, avoiding BG-coloring to partially cover them
				do.call("lines", rec.env$VLineCollection[[i]])
			
			rec.env$VLineCollection <- NULL											# should not be returned
		}
		
		box()
		


		if(!is.null(Range2))
		{
			# adjust plotting region to user-specification by restricting the plotting region
			USR 	<- par("usr")
			PLT 	<- par("plt")
			PLT2 	<- PLT
			PLT2[3] <- PLT[4] - diff(Range2) * diff(PLT[3:4]) / (diff(USR[3:4]))
			par(plt=PLT2)
		}	
		if(!is.null(SDline))
		{
			SDline$x <- as.numeric(names(StatsList$SDvec))
			SDline$y <- StatsList$SDvec
			do.call("lines", args=SDline)
		}
		
		if(!is.null(Range2))			# 
			par(plt=PLT)
		
		if(!is.null(VCnam))
		{
#			VCnam$side=2
			VCnam$at <- SDYbound[-length(SDYbound)]+diff(SDYbound[1:2])/2
			if(is.null(VCnam$text))																			# do not overwrite user-specified text
				VCnam$text <- nest
			VCnam$las <- 1
			VCnam$adj <- 1
			do.call("mtext", args=VCnam)
		}
	}  
	
	#if(type == 3)                                          # adds number of observation per bin to the plot (between variability chart and SD/CV-plot)
	#{
	#    tmp <- StatsList$Nobs
	#    len <- length(tmp)
	#    rxb <- range(Xbound)
	#    dx <- diff(rxb)/(len+1)
	#    mtext(side=3, line=1.25, at=0-diff(rxb)*.1, text="Obs:", cex=.75)
	#    mtext(side=3, at=seq(rxb[1]+dx/2, rxb[2]-dx/2, length.out=len), text=tmp, cex=.75, line=1.25)
	#}
	
	old.par$yaxs <- old.par$xaxs <- "i"						# unfortunatly need to be maintained in order to further add elements to the plot
	old.par$mar <- par("mar")
	par(old.par)
	
	Data$Xcoord <- rec.env$dat$x
	
	invisible(Data)
}








#' Build a Nested List.
#'
#' Function \code{buildList} creates a nested-list reflecting the hierarchical structure of a fully-nested model, respectively, the imputed
#' hierarchical structure of the data (see details).
#' 
#' This function is not intended to be used directly and serves as helper function for \code{\link{varPlot}}.
#' Each factor-level, on each level of nesting is accompanied with a set of descriptive statistics, such as mean, median, var, sd, ... which can be evaluated
#' later on. These information are used in function \code{varPlot}, which implements a variability chart.
#' Note, that this function is also used if data does not correspond to a fully-nested design, i.e. the hierarchical structure is
#' inferred from the model formula. The order of main factors (not nested within other factors) appearing in the model formula determines
#' the nesting structure imputed in order to plot the data as variability chart.
#' 
#' @param Data          (data.frame) with the data
#' @param Nesting       (character) vector specifying the nesting structure with the top-level variable name
#'                      as 1st element and the variance component one above the residual error as last element
#' @param Current       (character) string specifying the current level which has to be processed
#' @param resp          (character) string specifying the name of the response variable (column in 'Data')
#' @param keep.order    (logical) TRUE = the ordering of factor-levels is kept as provided by 'Data', FALSE = factor-levels are sorted on 
#'                      and within each level of nesting 
#' @param useVarNam     (logical) TRUE = each factor-level specifier is pasted to the variable name of the current variable and used as list-element name, 
#'                               FALSE = factor-level specifiers are used as names of list-elements; the former is useful when factor levels are indicated
#'                               as integers, e.g. days as 1,2,..., the latter is useful when factor levels are already unique, e.g. day1, day2, ...
#' @param sep           (character) string specifying the separator-string in case useVarNam=TRUE
#' @param na.rm         (logical) TRUE = NAs will be removed before computing the descriptive statistics AND NAs will be omitted when counting number of elements, 
#'                               FALSE = if there are NAs, this will result in NAs for the descriptive statistics  
#' @param Points        (list) specifying all parameters applicable to function 'points', used to specify scatterplots per lower-end factor-level
#'                      (e.g. run/part in EP05-A2 experiments). If list-element "col" is itself a list with elements "var" and "col", where the former
#'                      specifies a variable used for assigning colors "col" according to the class-level of "var", point-colors can be used for indicating
#'                      specific sub-classes not addressed by the model/design (see examples).
#' 
#' @return (list) which was recursively built, representing the data of the fully-nested as hierarchy
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @examples
#' 
#' \dontrun{
#' # load data (CLSI EP05-A2 Within-Lab Precision Experiment)
#' data(dataEP05A2_3)
#' 
#' # build a list representing the hierarichal structure of a fully-nested model
#' # there needs to be a distinct hierarchy for being able to plot the data
#' # as variability chart (this function is not exported)
#' lst <- VCA:::buildList(Data=dataEP05A2_3, Nesting=c("day", "run"), Current="day", resp="y")
#' }

buildList <- function(Data, Nesting, Current, resp, keep.order=TRUE, useVarNam=TRUE, sep="", na.rm=TRUE, Points=list(pch=16, cex=.5, col="black"))
{
	lev <- unique(as.character(Data[,Current]))
	
	if( is.list(Points))                                            # initial call
	{
		# check whether point-colors have to be used for indicating class-levels
		
		if( (is.list(Points$col)) && ("var" %in% names(Points$col)) && (Points$col$var %in% colnames(Data)) )
		{
			if(length(unique(Data[,Points$col$var])) > length(unique(Points$col$col)))          # not enough colors, need to be recycled
			{
				Points$col$col <- rep(Points$col$col, ceiling(length(unique(Data[,Points$col$var])) / length(unique(Points$col$col))))
			}
			tmp <- 1:length(unique(Data[,Points$col$var]))
			names(tmp) <- unique(Data[,Points$col$var])
			Data$PointsColor <- Points$col$col[tmp[as.character(Data[,Points$col$var])]]                         # assign colors to class-levels of Points$col$var variable
		}
		else
		{
			Data$PointsColor <- rep(Points$col, nrow(Data))
		}
		
		# check whether plotting-symbols have to be used for indicating class-levels
		
		if( (is.list(Points$pch)) && ("var" %in% names(Points$pch)) && (Points$pch$var %in% colnames(Data)) )
		{
			if(length(unique(Data[,Points$pch$var])) > length(unique(Points$pch$pch)))          # not enough plotting symbols, need to be recycled
			{
				Points$pch$pch <- rep(Points$pch$pch, ceiling(length(unique(Data[,Points$pch$var])) / length(unique(Points$pch$pch))))
			}
			tmp <- 1:length(unique(Data[,Points$pch$var]))
			names(tmp) <- unique(Data[,Points$pch$var])
			Data$PointsPCH <- Points$pch$pch[tmp[as.character(Data[,Points$pch$var])]]                         # assign colors to class-levels of Points$col$var variable
		}
		else
		{
			Data$PointsPCH <- rep(Points$pch, nrow(Data))
		}
		
		# check whether plotting-symbol backgrounds have to be used for indicating some sort of information
		
		if( (is.list(Points$bg)) && ("var" %in% names(Points$bg)) && (Points$bg$var %in% colnames(Data)) )
		{
			if(length(unique(Data[,Points$bg$var])) > length(unique(Points$bg$bg)))          # not enough background colors, need to be recycled
			{
				Points$bg$bg <- rep(Points$bg$bg, ceiling(length(unique(Data[,Points$bg$var])) / length(unique(Points$bg$bg))))
			}
			tmp <- 1:length(unique(Data[,Points$bg$var]))
			names(tmp) <- unique(Data[,Points$bg$var])
			Data$PointsBG <- Points$bg$bg[tmp[as.character(Data[,Points$bg$var])]]                         # assign colors to class-levels of Points$col$var variable
		}
		else
		{
			Data$PointsBG <- rep(Points$bg, nrow(Data))
		}
		
		# check whether plotting-symbol sizes have to be used for indicating some sort of information
		
		if( (is.list(Points$cex)) && ("var" %in% names(Points$cex)) && (Points$cex$var %in% colnames(Data)) )
		{
			if(length(unique(Data[,Points$cex$var])) > length(unique(Points$cex$cex)))          # not enough background colors, need to be recycled
			{
				Points$cex$cex <- rep(Points$cex$cex, ceiling(length(unique(Data[,Points$cex$var])) / length(unique(Points$cex$cex))))
			}
			tmp <- 1:length(unique(Data[,Points$cex$var]))
			names(tmp) <- unique(Data[,Points$cex$var])
			Data$PointsCEX <- Points$cex$cex[tmp[as.character(Data[,Points$cex$var])]]                         # assign colors to class-levels of Points$col$var variable
		}
		else
		{
			Data$PointsCEX <- rep(Points$cex, nrow(Data))
		}
		
		Points <- NULL
	}
	
	if(!keep.order)
	{
		suppressWarnings(levInt <- as.integer(lev))
		if(!any(is.na(levInt)))
			lev <- lev[sort(levInt, index.return=TRUE)$ix]
		else
			lev <- sort(lev)
	}
	
	lst <- vector("list", length=length(lev))
	attr(lst, "factor") <- Current
	
	Mean <- Median <- SD <- CV <- Nelem <- Nlevel <- numeric(length(lev))
	Nbin <- 0
	SDrange <- c(Inf, -Inf)
	CVrange <- c(Inf, -Inf)
	
	if(useVarNam)
	{
		names(lst) <- paste(Current, lev, sep=sep)
	}
	else
		names(lst) <- lev
	
	for(i in 1:length(lev))
	{
		tmpData <- Data[which(Data[,Current] == lev[i]),]
		
		Mean[i] <- mean(tmpData[,resp], na.rm=TRUE)
		Median[i] <- median(tmpData[,resp], na.rm=TRUE)
		SD[i] <- sd(tmpData[,resp], na.rm=TRUE)        
		CV[i] <- SD[i]*100/Mean[i]
		
		if(Current == Nesting[length(Nesting)])                     # last level above residual error
		{
			lst[[i]] <- tmpData[, resp]
			attr(lst[[i]], "Color")  <- tmpData$PointsColor
			attr(lst[[i]], "Symbol") <- tmpData$PointsPCH 
			attr(lst[[i]], "BG")	 <- tmpData$PointsBG
			attr(lst[[i]], "CEX")	 <- tmpData$PointsCEX
			
			if(na.rm)                                                # na.rm=TRUE?
				lst[[i]] <- na.omit(lst[[i]])
			
			Nelem[i] <- length(unique(tmpData[,Current]))            # number of bottom-level classes
			Nbin <- Nbin + 1
			
#            Mean[i] <- mean(lst[[i]])
#            Median[i] <- median(lst[[i]])
#            SD[i] <- sd(lst[[i]])            
#            CV[i] <- SD[i]*100/Mean[i]
			
			if(!is.na(SD[i]) && SD[i] < SDrange[1])                                   # range of all SD-values   
				SDrange[1] <- SD[i] 
			if(!is.na(SD[i]) && SD[i] > SDrange[2])
				SDrange[2] <- SD[i]
			
			if(!is.na(CV[i]) && CV[i] < CVrange[1])                                   # range of all CV-values   
				CVrange[1] <- CV[i] 
			if(!is.na(CV[i]) && CV[i] > CVrange[2])
				CVrange[2] <- CV[i]
		}   
		else
		{   
			lst[[i]] <- buildList(Data=tmpData, Nesting=Nesting, Current=Nesting[which(Nesting == Current)+1], 
					resp=resp, keep.order=keep.order, useVarNam, Points=Points)
			Nelem[i] <- sum(attr(lst[[i]], "Nelem"))
			Nbin <- Nbin + attr(lst[[i]], "Nbin")
			
			tmpSD <- attr(lst[[i]], "SDrange")
			
			if(tmpSD[1] < SDrange[1])                               # refresh SD-range values
				SDrange[1] <- tmpSD[1]
			if(tmpSD[2] > SDrange[2])
				SDrange[2] <- tmpSD[2]
			
			tmpCV <- attr(lst[[i]], "CVrange")
			
			if(tmpCV[1] < CVrange[1])                                  # range of all CV-values   
				CVrange[1] <- tmpCV[1] 
			if(tmpCV[2] > CVrange[2])
				CVrange[2] <- tmpCV[2]
		} 
		Nlevel[i] <- length(lst[[i]])
	}
	attr(lst, "Nelem")   <- Nelem
	attr(lst, "Nlevel")  <- Nlevel
	attr(lst, "Nbin")    <- Nbin
	attr(lst, "SDrange") <- SDrange
	attr(lst, "CVrange") <- CVrange
	
	attr(lst, "Stats") <- list(Mean=Mean, Median=Median, SD=SD, CV=CV)
	
	return(lst)
}


#' Standard 'plot' Method for 'VCA' S3-Objects.
#' 
#' Create a variability chart from a 'VCA'-object, i.e. from a fitted model.
#' 
#' This function extracts the data and the model-formula from a fitted 'VCA'-object and calls function \code{\link{varPlot}}
#' accepting all its arguments. Please see the documention of function \code{\link{varPlot}} for a detailed description.
#' 
#' It will not be differentiated between fixed and random effects when calling this function on a fitted linear mixed model.
#' 
#' @param x         (VCA) object 
#' @param ...       additional arguments to be passed to or from methods.
#' 
#' @return 	nothing, instead a plot is generated
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @seealso \code{\link{varPlot}}, \code{\link{anovaVCA}},\code{\link{remlVCA}}, \code{\link{anovaMM}},\code{\link{remlMM}}
#' 
#' @method plot VCA
#' 
#' @examples 
#' \dontrun{
#' data(dataEP05A2_1)
#' fit <- anovaVCA(y~day/run, dataEP05A2_1)
#' 
#' # standard plot without any extras
#' plot(fit)
#' 
#' # plot with some additional features
#' plot(fit, MeanLine=list(var=c("int", "day"), col=c("cyan", "blue"), lwd=c(2,2)))
#' 
#' # more complex model
#' data(realData)
#' Data <- realData[realData$PID == 1,]
#' fit2 <- anovaVCA(y~(calibration+lot)/day/run, Data)
#' plot(fit2, 
#' 		BG=list(var="calibration",
#' 				col=c("#f7fcfd","#e5f5f9","#ccece6","#99d8c9",
#' 				      "#66c2a4","#41ae76","#238b45","#006d2c","#00441b"),
#' 				col.table=TRUE),
#'		VLine=list(var=c("calibration", "lot"),
#' 				   col=c("black", "darkgray"), lwd=c(2,1), col.table=TRUE),
#' 		JoinLevels=list(var="lot", col=c("#ffffb2","orangered","#feb24c"),
#' 				        lwd=c(2,2,2)),
#' 		MeanLine=list(var="lot", col="blue", lwd=2))
#'
#' }

plot.VCA <- function(x, ...)
{
	Data <- x$data
	form <- x$formula
	varPlot(form=form, Data=Data, ...)
}



#' Add Legend to Margin.
#' 
#' This function accepts all parameters applicable in and forwards them to function \code{\link{legend}}.
#' There will be only made some modifications to the X-coordinate ensuring that the legend is plotted in
#' the right margin of the graphic device. Make sure that you have reserved sufficient space in the right
#' margin, e.g. 'plot.VFP(....., mar=c(4,5,4,10))'.
#' 
#' @param x			(character, numeric) either one of the character strings "center","bottomright", "bottom", "bottomleft", 
#' 					"left", "topleft", "top", "topright", "right" or a numeric values specifying the X-coordinate in user
#' 					coordinates
#' @param y			(numeric) value specifying the Y-coordiante in user coordinates, only used in case 'x' is numeric
#' @param margin	(character) string specifying in which part of the margin the legend shall be added, choices are
#' 					"right", "bottomright", "bottom", "bottomleft", 
#' 					"left", "topleft", "top", "topright" with "right" being the default
#' @param offset	(numeric) value in [0, 0.5] specifying the offset as fraction in regard to width of the right margin
#' @param ...		all parameters applicable in function \code{\link{legend}}
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @examples 
#' \dontrun{
#' 
#' par( mar=c(10,10,10,10) )
#' plot(1, type="n", axes=FALSE, xlab="", ylab="")
#' box()
#' # add legend to different regions within the 'margin'
#' legend.m(margin="topleft", 		fill="black",	legend=c("topleft"))
#' legend.m(margin="top", 			fill="red", 	legend=c("top"))
#' legend.m(margin="topright", 		fill="blue",	legend=c("topright"))
#' legend.m(margin="right", 		fill="green",	legend=c("right"))
#' legend.m(margin="bottomright", 	fill="yellow",	legend=c("bottomright"))
#' legend.m(margin="bottom", 		fill="orange",	legend=c("bottom"))
#' legend.m(margin="bottomleft", 	fill="cyan",	legend=c("bottomleft"))
#' legend.m(margin="left", 			fill="magenta", legend=c("left"))
#' 
#' data(dataEP05A2_3)
#' dataEP05A2_3$user <- sample(rep(c(1,2), 40))
#'  
#' varPlot( y~day+day:run, dataEP05A2_3, mar=c(1,5,1,7), VCnam=list(side=4), 
#' 	        Points=list(pch=list(var="user", pch=c(2, 8))) )
#' # always check order of factor levels before annotating
#' order(unique(dataEP05A2_3$user))
#' legend.m(pch=c(8,2), legend=c("User 1", "User 2"))
#' 
#' # using different colors 
#' varPlot( y~day+day:run, dataEP05A2_3, mar=c(1,5,1,7), VCnam=list(side=4),
#'          Points=list(col=list(var="user", col=c("red", "green"))) )
#' legend.m(fill=c("green", "red"), legend=c("User 1", "User 2"))
#'
#' # two additional classification variables
#' dataEP05A2_3$user <- sample(rep(c(1,2), 40))
#' dataEP05A2_3$cls2 <- sample(rep(c(1,2), 40))
#' 
#' # now combine point-coloring and plotting symbols
#' # to indicate two additional classification variables
#' varPlot( y~day+day:run, dataEP05A2_3, mar=c(1,5,1,10),
#'          VCnam=list(side=4, cex=1.5),
#'          Points=list(col=list(var="user", col=c("red", "darkgreen")),
#'                      pch=list(var="cls2", pch=c(21, 22)),
#'                      bg =list(var="user", bg =c("orange", "green"))) )
#'
#' # add legend to (right) margin
#' legend.m(margin="right", pch=c(21, 22, 22, 22), 
#'          pt.bg=c("white", "white", "orange", "green"),  
#'          col=c("black", "black", "white", "white"), 
#'          pt.cex=c(1.75, 1.75, 2, 2),
#'          legend=c("Cls2=1", "Cls2=2", "User=2", "User=1"), 
#'          cex=1.5)
#' }

legend.m <- function(	x=c("center","bottomright", "bottom", "bottomleft", 
							"left", "topleft", "top", "topright", "right"), 
						y=NULL, 
						margin=c(	"right", "bottomright", "bottom", "bottomleft", 
									"left", "topleft", "top", "topright"),
						offset=.05, ...)
{
	stopifnot(	is.numeric(x) || is.character(x) )
	if(is.character(x))
	{
		x <- match.arg(x[1], choices=c("center","bottomright", "bottom", "bottomleft", 
						"left", "topleft", "top", "topright", "right"))
	}else
	{
		stopifnot(is.numeric(y))
	}
	margin <- match.arg(margin[1], 
						choices=c(	"right", "bottomright", "bottom", "bottomleft", 
									"left", "topleft", "top", "topright"))
	
	par(xpd=TRUE)
	args <- list(...)
	
	USR  <- par("usr")
	PLT  <- par("plt")
	FIG  <- par("fig")
	
	if(margin == "right")
	{
		wm 		<- FIG[2] - PLT[2]			# width of the margin
		hm 		<- PLT[4] - PLT[3]			# height of the margin
		Left	<- PLT[2]
		Bottom	<- PLT[3]
	}
	else if(margin == "topright")
	{
		wm 		<- FIG[2] - PLT[2]						
		hm 		<- FIG[4] - PLT[4]
		Left	<- PLT[2]
		Bottom	<- PLT[4]
	}
	else if(margin == "bottomright")
	{
		wm 		<- FIG[2] - PLT[2]						
		hm 		<- PLT[3] - FIG[3]
		Left	<- PLT[2]
		Bottom	<- FIG[3]
	}
	else if(margin == "bottomleft")
	{
		wm 		<- PLT[1] - FIG[1]						
		hm 		<- PLT[3] - FIG[3]
		Left	<- FIG[1]
		Bottom	<- FIG[3]
	}
	else if(margin == "bottom")
	{
		wm 		<- PLT[2] - PLT[1]						
		hm 		<- PLT[3] - FIG[3]
		Left	<- PLT[1] 
		Bottom	<- FIG[3]
	}
	else if(margin == "top")
	{
		wm 		<- PLT[2] - PLT[1]
		hm 		<- FIG[4] - PLT[4]
		Left	<- PLT[1]
		Bottom	<- PLT[4]
	}
	else if(margin == "left")
	{
		wm 		<- PLT[1] - FIG[1]	
		hm 		<- PLT[3] - FIG[3]
		Left	<- FIG[1]
		Bottom	<- PLT[3]
	}
	else if(margin == "topleft")
	{
		wm 		<- PLT[1] - FIG[1]
		hm 		<- FIG[4] - PLT[4]
		Left	<- FIG[1]
		Bottom	<- PLT[4]
	}	
	
	if(is.character(x))
	{
		X.orig	<- x
		xjust 	<- 0.5							# defaults to center
		x 		<- Left + 0.5 * wm
		yjust   <- 0.5
		y		<- Bottom + 0.5 * hm
		
		if(grepl("left", X.orig))
		{
			xjust 	<- 0
			x 		<- Left + offset * wm
		}		
		if(grepl("right", X.orig))
		{
			xjust 	<- 1
			x 		<- Left + (1-offset) * wm
		}
		if(grepl("top", X.orig))
		{
			yjust 	<- 1
			y 		<- Bottom + (1-offset) * hm
		}
		if(grepl("bottom", X.orig))
		{
			yjust 	<- 0
			y 		<- Bottom + offset * hm
		}
	}
	
	x <- grconvertX(x, from="nic", to="user")
	y <- grconvertY(y, from="nic", to="user")
	
	args$x 		<- x
	args$y 		<- y
	
	if(!"xjust" %in% names(args))
		args$xjust 	<- xjust
	if(!"yjust" %in% names(args))
		args$yjust 	<- yjust
	
	do.call(legend, args)
	par(xpd=FALSE)
}



