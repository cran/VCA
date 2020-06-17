## ----global_options, echo=FALSE, eval=TRUE------------------------------------
knitr::opts_chunk$set(fig.width=7, fig.height=5, fig.align='center', fig.path='figures/',
                      echo=TRUE, eval=TRUE, warning=FALSE, message=FALSE)

## ----create_data, echo=FALSE--------------------------------------------------
library(VCA)
library(STB)
data(VCAdata1)
datS5 <- subset(VCAdata1, sample==5)
# limit number of decimal places
options(digits=4)

## ----str_VCAdata1, echo=FALSE-------------------------------------------------
str(VCAdata1)

## ----create_data_fake, eval=FALSE---------------------------------------------
#  library(VCA)
#  data(VCAdata1)
#  datS5 <- subset(VCAdata1, sample==5)

## ----varPlot_plain, fig.width=10, fig.height=7--------------------------------
varPlot(form=y~(device+lot)/day/run, Data=datS5)

## ----varPlot_full, fig.width=10, fig.height=7---------------------------------
varPlot(y~(device+lot)/day/run, datS5, 
		MeanLine=list(var=c("int", "device", "lot"), 
					  col=c("white", "blue", "magenta"), lwd=c(2,2,2)), 
		BG=list(var="lot", col=paste0("gray", c(70,80,90))))

## ----varPlot_Outlier, fig.width=10, fig.height=7------------------------------
#indicate outliers by new variable
datS5$out <- 0
datS5$out[which(datS5$device==1 & datS5$lot==2 & datS5$day==4 &datS5$run==2)] <- 1

varPlot(y~(device+lot)/day/run, datS5,
		# plots horizontal lines for sub-sets
		MeanLine=list(	var=c("int", "device", "lot"),
						col=c("white", "blue", "magenta"),
						lwd=c(3,3,3)),
		# colors the background according to a variable's levels
		# 'col.table=TRUE' indicates the variable in the table 
		BG=list(var="lot", col=paste0("gray", c(70,80,90)),
				col.table=TRUE),
		# tailoring the appearance of measurements (points)
		Points=list(pch=list(var="out", pch=c(21, 24)),
					col=list(var="out", col=c("black", "red")),
					bg= list(var="out", bg=c("white", "yellow")),
					cex=list(var="out", cex=c(1, 1.5))),
		# specification of the text within the table below
		VarLab=list(list(cex=1.75), list(cex=1.5),
					list(cex=1, srt=90), list()),
		# variable names right of the table (since 'side=4')
		VCnam=list(cex=1.25, side=4),
		# use 33% of the height of the upper part for the table
		htab=.33)  

## ----plotRandVar_fake, eval=FALSE---------------------------------------------
#  fitS5 <- anovaVCA(y~(device+lot)/day/run, datS5)
#  # if varPlot was called before, better shut down the old graphics device
#  # graphical parameters were not reset to allow the user to add further
#  # information to the variability chart, which would not be possible otherwise
#  dev.off()
#  plotRandVar(fitS5, term="cond", mode="student")

## ----plotRandVar, echo=FALSE, fig.width=10, fig.height=7----------------------
fitS5 <- anovaVCA(y~(device+lot)/day/run, datS5)
plotRandVar(fitS5, term="cond", mode="student")
abline(h=c(-3, 3), lty=2, col="red")
mtext(side=4, at=c(-3, 3), col="red", line=.25, las=1, text=c(-3, 3))

## ----plotRandVar_pick_fake, eval=FALSE----------------------------------------
#  plotRandVar(fitS5, term="cond", mode="student", pick=TRUE)
#  abline(h=c(-3, 3), lty=2, col="red", lwd=2)
#  mtext(side=4, at=c(-3, 3), col="red", line=.25, las=1, text=c(-3, 3))

## ---- plotRandVar_picked, echo=FALSE------------------------------------------
knitr::include_graphics("figures/plotRandVar_pick.png")

## ----picked_observations------------------------------------------------------
datS5[c("1191","1192"),]

## ----define_MD68--------------------------------------------------------------
# Function computing the MD68 using SAS PCTLDEF5 quantile definition.
# x (numeric) values from which the MD68 shall be computed
# na.rm (logical) TRUE = missing values will be excluded automatically
md68 <- function(x, na.rm=FALSE)
{
	stopifnot(is.numeric(x))
	Med 	<- median(x, na.rm=na.rm)
	Diff 	<- abs(x-Med)
	MD68 	<- quantile(Diff, probs=.68, type=2)
	MD68
}

## ----Outlier_Detection_Algorithm----------------------------------------------
# identify groups of intermediate precision (IP) measuring conditions
datS5$IPgroup <- paste(datS5$device, datS5$lot, sep="_")
# uniquely identify replicate groups within IP-groups
datS5$RepGroup <- paste(datS5$day, datS5$run, sep="_")

# Define a function performing the MD68-based outlier-algorithm.
# The data.frame will be returned with two additional variables,
# "diff" (absolute differences from the replicate-group median) and
# "outlier" (TRUE = outlier, FALSE = no outlier).
#  
# obj 		(data.frame) with at least two variables
# resp 		(character) name of the numeric response variable
# RepGroup 	(character) name of the replicate-grouping variable
  
md68OutlierDetection <- function(obj=NULL, resp=NULL, RepGroup=NULL)
{
	stopifnot(is.data.frame(obj))
	cn 		<- colnames(obj)
	stopifnot(is.character(resp) && resp %in% cn && is.numeric(obj[,resp]))
	stopifnot(is.character(RepGroup) && RepGroup %in% cn)
	MD68 	<- md68(obj[, resp])
	Crit 	<- MD68*2.75
	# tapply returns a list with as many elements as there are unique
	# elements in obj[,RepGroup], which must be converted to a vector
	obj$diff <- unlist(
					tapply(obj[,resp], obj[,RepGroup],
					function(x)
					{
						m <- median(x)
						d <- abs(x - m)
						d
					}))
	obj$threshold 	<- Crit
	obj$outlier 	<- obj$diff > Crit
	obj
}

## ----Outlier_Detection_Example------------------------------------------------
IPgroups <- unique(datS5$IPgroup)
for(i in 1:length(IPgroups))
{
	tmpData <- subset(datS5, IPgroup == IPgroups[i])
	if(i == 1)
		out <- md68OutlierDetection(tmpData, resp="y", RepGroup="RepGroup")
	else
		out <- rbind(out, md68OutlierDetection(tmpData, resp="y", RepGroup="RepGroup"))
}
head(out)


# Were there any outliers detected?
any(out$outlier)

## ----Show_VCA_Result----------------------------------------------------------
# use non-default number of decimal places 
print(fitS5, digits=4)

## ----RandomEffects_Day, fig.height=7, fig.width=10----------------------------
plotRandVar(fitS5, term="device:lot:day", mode="student")

## ----RandomEffects_Run, fig.height=7, fig.width=10----------------------------
plotRandVar(fitS5, term="device:lot:day:run", mode="student")

## ----ANOVAtable_fake, eval=FALSE----------------------------------------------
#  fitS5$aov.tab

## ----ANOVAtable, echo=FALSE---------------------------------------------------
as.data.frame(fitS5$aov.tab)

## ----STB_construction, eval=FALSE---------------------------------------------
#  set.seed(23)
#  fitS5.LMM <- anovaMM(y~((device)+(lot))/(day)/(run), datS5)
#  STB.res <- stb(fitS5.LMM, term="cond", mode="student", N=5000)

## ---- STB_picked1, echo=FALSE, fig.align="center"-----------------------------
knitr::include_graphics("figures/STB1.png")

## ----STB_picked_obs_fake, eval=FALSE------------------------------------------
#  plot(STB.res, pick=TRUE)

## ---- STB_picked2, echo=FALSE, fig.align="center"-----------------------------
knitr::include_graphics("figures/STB1_pick.png")

## ----delete_obs1_1, eval=TRUE-------------------------------------------------
datS5.reduced <- datS5[!rownames(datS5) == "1192",]
fitS5.reduced <- anovaMM(y~((device)+(lot))/(day)/(run), datS5.reduced)

## ----delete_obs1_2, eval=FALSE------------------------------------------------
#  STB.res2 <- stb(fitS5.reduced, term="cond", mode="student", N=5000)

## ---- STB_obs_removed, echo=FALSE, fig.align="center"-------------------------
knitr::include_graphics("figures/STB2.png")

## ----plotRandVar2, fig.height=7, fig.width=10---------------------------------
plotRandVar(fitS5.reduced, term="cond", mode="student")

## ----ANOVAtable_full, echo=FALSE----------------------------------------------
as.data.frame(fitS5$aov.tab)

## ----ANOVAtable_reduced, echo=FALSE-------------------------------------------
as.data.frame(fitS5.reduced$aov.tab)

## ----WithinLab_Example_plot, echo=TRUE, fig.height=7, fig.width=10------------
# Function converts a color-string into RGB-code
# col (character) string specifying an R-color
# alpha (numeric) degree of transparency in [0, 1], 0=fully transparency, 1=opaque
asRGB <- function(col, alpha)
			rgb(t(col2rgb(col))/255, alpha=alpha)
			
data(dataEP05A2_3)

varPlot(y~day/run, dataEP05A2_3,
	# controls horizontal mean lines
	MeanLine=list(var=c("int", "day"), col=c("gray75", "blue"), lwd=c(2,2)),
	# controls how points (concentrations) are plotted, here using semi-transparency
	# to see overlayed points
	Points=list(pch=16, col=asRGB("black", .5), cex=1.25),
	# controls how replicate-means are plotted
	Mean=list(col="magenta", cex=1.25, lwd=2),
	# controls how the title is shown
	Title=list(main="20 x 2 x 2 Single-Site Evaluation", cex.main=1.75),
	# controls plotting of levels per VC, if as many lists as there are VCs are
	# specified, each VC can be specified individually
	VarLab=list(list(cex=1.5), list(cex=1.25)),
	# controls how names of VCs are plotted
	VCnam=list(font=2, cex=1.5),
	# controls appearance of the Y-axis label
	YLabel=list(text="Concentation [mg/dL]", las=0, line=3, font=2, cex=1.25),
	# Y-axis labels rotated
	las=1)

## ----WithinLab_Example_fit, echo=TRUE-----------------------------------------
# fit 20 x 2 x 2 model to data
fit.SS3 <- fitVCA(y~day/run, dataEP05A2_3)
fit.SS3

# estimate 95% confidence intervals, request CI for
# all variance components via 'VarVC=TRUE'
inf.SS3 <- VCAinference(fit.SS3, VarVC=TRUE)
inf.SS3

## ----Multi_Site_Design_1_Plot, echo=TRUE, fig.height=7, fig.width=10----------
data(dataEP05A3_MS_1)
varPlot(y~site/day, dataEP05A3_MS_1,
		BG=list(var="site", col=paste0("gray", c(100, 80, 60))),
		Points=list(pch=16, col=asRGB("black", .5), cex=1.25),
		MeanLine=list(var=c("int", "site"), col=c("black", "orange"), lwd=c(2,2)),
		Mean=list(col="cyan", cex=1.25, lwd=2), las=1,
		YLabel=list(text="Concentation [mg/dL]", las=0, line=3, font=2, cex=1.25),
		Title=list(main="Multi-Site Evaluation on dataEP05A3_MS_1", cex.main=1.75),
		VCnam=list(font=2, cex=1.5),
		VarLab=list(list(cex=1.5, font=2), list(cex=1.25, font=2)))

## ----Multi_Site_Design_1_Fit, echo=TRUE---------------------------------------
# fit 3 x 5 x 1 x 5 model to data
fit.MS1 <- fitVCA(y~site/day, dataEP05A3_MS_1, method="REML")
fit.MS1

## ----Multi_Site_Design_2, echo=TRUE-------------------------------------------
# simulate fit 3 x 5 x 2 x 3 model to data
set.seed(23)
dat.MS2 <- data.frame( 	y=50 +
						# 3 random effects for sites
						rep(rnorm(3,0,2.5), rep(30, 3)) +
						# 15 random effects for days
						rep(rnorm(15,0, 2), rep(6, 15)) +
						# 30 random effects for runs
						rep(rnorm(30,0, 1), rep(3, 30)) +
						# residual error (repeatability)
						rnorm(90,0,1.5),
						site = gl(3, 30, labels=paste0("Site_", 1:3)),
						day = gl(5, 6, 90),
						run =gl(2, 3, 90)
					)

## ----Multi_Site_Design_2_Plot, echo=TRUE, fig.height=7, fig.width=10----------
varPlot(y~site/day/run, dat.MS2,
	BG=list(var="site", col=paste0("gray", c(100, 80, 60))),
	Points=list(pch=16, col=asRGB("black", .5), cex=1.25),
	MeanLine=list( var=c("int", "site", "day"),
	col=c("black", "orange", "blue"),
	lwd=c(2,2,2)),
	Mean=list(col="cyan", cex=1.25, lwd=2), las=1,
	YLabel=list(text="Concentation [mg/dL]", las=0, line=3, font=2, cex=1.25),
	Title=list(main="3 x 5 x 2 x 3 Multi-Site Evaluation", cex.main=1.75),
	VCnam=list(font=2, cex=1.5),
	# controls for which variable vertical lines are added between levels
	# and how these are plotted
	VLine=list(var="day", col="gray75"),
	VarLab=list(list(cex=1.5), list(cex=1.25), list(cex=1.25)))

## ----Multi_Site_Design_2_Fit, echo=TRUE---------------------------------------

# fit 3 x 5 x 2 x 3 model to data (ANOVA is default)
fit.MS2 <- fitVCA(y~site/day/run, dat.MS2)
print(fit.MS2, digits=4)

## ----Multi_Lot_Multi_Site_Assignment, echo=FALSE------------------------------

mat <- matrix("X", 3, 3)
rownames(mat) <- c("ReagentLot_1", "ReagentLot_2", "ReagentLot_3")
colnames(mat) <- c("Site_1", "Site_2", "Site_3")
print(mat, quote=FALSE)

mat2 <- mat
mat2[1,3] <- mat2[2,1] <- mat2[3,2] <- NA
print(mat2, quote=FALSE, na.print="")

mat3 <- mat
mat3[2,1] <- mat3[3,1:2] <- NA
print(mat2, quote=FALSE, na.print="")

## ----Multi_Lot_Multi_Site_Fit_correct, echo=TRUE------------------------------
fit.MSML <- fitVCA(y~(device+lot)/day/run, datS5)
print(fit.MSML, digits=4)

## ----Multi_Lot_Multi_Site_Fit_wrong, echo=TRUE--------------------------------
fit.MSML2 <- fitVCA(y~device/lot/day/run, datS5)
print(fit.MSML2, digits=4)

## ----ConfidenceIntervals, echo=TRUE-------------------------------------------
inf.MSML 	<- VCAinference(fit.MSML,  VarVC=TRUE)
inf.MSML2	<- VCAinference(fit.MSML2, VarVC=TRUE)

# print CI for CV, other options are "all", "VC", "SD", and "VCA" 
print(inf.MSML,  what="CV", digits=2)
print(inf.MSML2, what="CV", digits=2)

