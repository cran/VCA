# TODO: unit-test functions for function 'remlMM'
# 
# Author: schueta6
###############################################################################


cat("\n\n**************************************************************************")
cat("\nVariance Component Analysis (VCA) - test cases defined in runit.remlMM.R.")
cat("\n**************************************************************************\n\n")

### load all testdata

data(dataEP05A2_1)
data(dataEP05A2_2)
data(dataEP05A2_3)

data(dataEP05A3_MS_1)
data(dataEP05A3_MS_2)
data(dataEP05A3_MS_3)

data(dataRS0003_1)
data(dataRS0003_2)
data(dataRS0003_3)

data(dataRS0005_1)
data(dataRS0005_2)
data(dataRS0005_3)

data(VCAdata1)

datS2 <- VCAdata1[VCAdata1$sample==2, ]


# basic check of equivalence to ANOVA Type-1 implementation for balanced data without zero variance estimates

TF001.anovaVCA_vs_remlMM <- function()
{
	data(dataEP05A2_3)
	fit0 <- anovaVCA(y~day/run, dataEP05A2_3)
	fit1 <- remlMM(y~(day)/run, dataEP05A2_3)
	
	checkEquals(as.numeric(fit0$aov.tab[,"VC"]), fit1$aov.tab[,"VC"], tolerance=1e-7)
	checkEquals(as.numeric(fit0$aov.tab[c("total", "error"), "DF"]), fit1$aov.tab[c("total", "error"), "DF"], tolerance=1e-7)
}


# test covariance parameter estimates, fixed effects, and covariance matrix of VCs

TF002.anovaMM.balanced1 <- function()
{
	
	fit   <- remlMM(y~(lot+device)/(day)/(run), datS2)
	
	checkEquals(as.numeric(round(fit$aov.tab[-1, "VC"], c(4,5,5))), c(0.3147, 0.05382, 0.04409))	# VCs
	
	checkEquals(as.numeric(round(fixef(fit)[,1], 4)), c(23.2393, 0.3655, -0.5570, 0, 1.1730, 0.5925, 0))
	
}

# checks whether both function yield identical results when fitting a random model, 
# once setting negative variance estimates equal to zero, and once allowing negative 
# variance estimates

TF003.anovaMM.NegVC.balanced <- function()
{
	data(dataEP05A2_1)
	
	test.dat <- dataEP05A2_1
	set.seed(123)
	test.dat$y <- test.dat$y + rnorm(40,,2.5)						# add something which yields negative estimates
	
	resRM1 <- remlVCA(y~day/run, test.dat)							# constrain VC to 0
	resMM1 <- remlMM(y~(day)/(run), test.dat)
	
	checkEquals(resRM1$aov.tab, resMM1$aov.tab)
	
	resRM2 <- remlVCA(y~day/run, test.dat)			# allowing negative VCs
	resMM2 <- remlMM(y~(day)/(run), test.dat)
	
	checkEquals(resRM1$aov.tab, resMM1$aov.tab)
}

# test againts SAS PROC MIXED:
#
# proc mixed data=ep5_2 method=type1 asycov cl=wald;
#   class day run;
#   model y = day;
#   random day*run/solution;
# run;

TF004.anovaMM.ProcMixed.balanced1 <- function()
{
	data(dataEP05A2_2)
	res <- remlMM(y~day/(run), dataEP05A2_2)
	
	checkEquals(round(as.numeric(res$aov.tab[-1, "VC" ]),4), c(2.8261, 3.7203))			# VCs
	checkEquals(round(as.numeric(res$FixedEffects[c("day10", "day15", "day20"),]), 4), c(1.8259, 4.4052, 0))		# some fixed effets
		
}



TF005.anova_vs_reml.balanced <- function()
{		
	data(dataEP05A2_2)
	fit0 <- anovaMM(y~day/(run), dataEP05A2_2)
	fit1 <- remlMM(y~day/(run), dataEP05A2_2)
	
	checkEquals(as.numeric(fit0$aov.tab[,"VC"]), fit1$aov.tab[,"VC"], tolerance=1e-7)
	checkEquals(as.numeric(fit0$aov.tab[c("total", "error"), "DF"]), fit1$aov.tab[c("total", "error"), "DF"], tolerance=1e-7)
	
	fit2 <- anovaMM(y~(day)/(run), dataEP05A2_2)
	fit3 <- remlMM(y~(day)/(run), dataEP05A2_2)
	
	checkEquals(as.numeric(fit2$aov.tab[,"VC"]), fit3$aov.tab[,"VC"], tolerance=1e-7)
	checkEquals(as.numeric(fit2$aov.tab[c("total", "error"), "DF"]), fit3$aov.tab[c("total", "error"), "DF"], tolerance=1e-7)
}


TF006.anova_vs_reml.balanced <- function()
{		
	data(dataEP05A2_2)
	fit0 <- anovaMM(y~day/(run), dataEP05A2_2)
	fit1 <- remlMM(y~day/(run), dataEP05A2_2)
	
	checkEquals(as.numeric(fixef(fit0)), as.numeric(fixef(fit1)), tolerance=1e-7)
	checkEquals(as.numeric(vcovVC(fit0)), as.numeric(vcovVC(fit1)), tolerance=1e-7)
}


TF007.anova_vs_reml.regression.balanced <- function()
{		
	data(Orthodont)
	Ortho <- Orthodont
	Ortho$age2 <- Ortho$age - 11
	
	anova.fit <- anovaMM( distance~Sex*age2+(Subject)*age2-1, Ortho)
	reml.fit  <- remlMM(  distance~Sex*age2+(Subject)*age2-1, Ortho, cov=FALSE)
	
	checkEquals(as.numeric(anova.fit$aov.tab[, "VC"]), as.numeric(reml.fit$aov.tab[, "VC"]), tolerance=1e-7)
	checkEquals(as.numeric(anova.fit$aov.tab[c("total", "error"), "DF"]), as.numeric(reml.fit$aov.tab[c("total", "error"), "DF"]), tolerance=1e-7)
	checkEquals(as.numeric(vcovVC(anova.fit)), as.numeric(vcovVC(reml.fit)), tolerance=1e-6)
	checkEquals(as.numeric(fixef(anova.fit)), as.numeric(fixef(reml.fit)), tolerance=1e-7)
}


TF008.remlMM.exception_handling <- function()
{
	checkException(remlMM())                                                               # no input at all 
	checkException(remlMM(Data=1))                                                                                               
	checkException(remlMM(Data=data.frame()))                 
	checkException(remlMM(Data=data.frame(y=1:10)))
	checkException(remlMM(z~day/run, Data=data.frame(y=1:10)))
	checkException(remlMM(y~day/run, Data=data.frame(y=1:10, day=1:10)))
}


TF009.remlMM.Percent_Total_results <- function()
{ 
	res <- remlMM(y~(day)/run, Data=dataEP05A2_1)              
	checkEquals(round(as.numeric(res$aov.tab[,"VC"])*100/sum(as.numeric(res$aov.tab[-1, "VC"])),6), round(as.numeric(res$aov.tab[,"%Total"]), 6))           # exclude total variance in the sum  
	
	res <- remlMM(y~(day)/run, Data=dataEP05A2_2)              
	checkEquals(round(as.numeric(res$aov.tab[,"VC"])*100/sum(as.numeric(res$aov.tab[-1, "VC"])),6), round(as.numeric(res$aov.tab[,"%Total"]), 6)) 
	
	res <- remlMM(y~(day)/run, Data=dataEP05A2_3)              
	checkEquals(round(as.numeric(res$aov.tab[,"VC"])*100/sum(as.numeric(res$aov.tab[-1, "VC"])),6), round(as.numeric(res$aov.tab[,"%Total"]), 6)) 
}


# testcase checks results against SAS PROC MIXED results with REML-estimation
# and centered values of the regression variable, below "sleep" is the original
# sleepstudy dataset imported to SAS
#
# data sleep2;		
#   set sleep;
#   days2 = days - 4.5;
# run;
# 
# proc mixed data=sleep2 method=reml;
#   class subject;
#   model reaction = days2/solution;
#   random subject subject*days2;
# run;

TF010.remlMM.regression.sleepstudy <- function()
{
	data(sleepstudy)
	sleep2 <- sleepstudy
	sleep2$days2 <- sleep2$Days - median(unique(sleep2$Days))
	fit.mm <- remlMM(Reaction~days2*(Subject), sleep2, cov=FALSE)
	
	checkEquals(as.numeric(fit.mm$aov.tab[-1,"VC"]), c(1408.73006360414, 35.0716604983165, 654.941027072368), tolerance=1e-6)
}


TF011.remlMM.zeroVariance <- function()
{
	data(dataEP05A2_3)
	dat1 <- dataEP05A2_3
	dat1$y <- dat1[1,"y"]
	dat1$cov <- rnorm(nrow(dat1),15,3)
	
	fit1 <- remlMM(y~day+cov+day:(run), dat1)
	checkEquals(as.numeric(fit1$aov.tab[,"VC"]), rep(0,3))
	fit2 <- remlMM(y~day/(run), dat1)
	checkEquals(as.numeric(fit2$aov.tab[,"VC"]), rep(0,3))
	fit3 <- remlMM(y~(day)/(run), dat1)
	checkEquals(as.numeric(fit3$aov.tab[,"VC"]), rep(0,4))
}


