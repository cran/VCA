# Test cases for function 'remlVCA' of R-package VCA
# 
# Author: André Schützenmeister
#
# Reference results were obtained by application of SAS 9.2 P PROC MIXED method=reml
# and by checking against R-function 'remlVCA', which is considered to be well tested.
#
#####################################################################################


cat("\n\n***********************************************************************")
cat("\nVariance Component Analysis (VCA) - test cases for function 'remlVCA'.")
cat("\n***********************************************************************\n\n")

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


### check results of function 'remlVCA' ###


# check whether function 'remlVCA' throws exceptions as expected

TF001.remlVCA.exception <- function()
{
	checkException(remlVCA())                                                               # no input at all 
	checkException(remlVCA(Data=1))                                                                                               
	checkException(remlVCA(Data=data.frame()))                 
	checkException(remlVCA(Data=data.frame(y=1:10)))
	checkException(remlVCA(z~day/run, Data=data.frame(y=1:10)))
	checkException(remlVCA(y~day/run, Data=data.frame(y=1:10, day=1:10)))
}

## EP05-A2 20/2/2 Within-Lab Precision Experiments (Single-Site Experiment in EP05-A3) - balanced
## check remlVCA against remlVCA with 4 significant digits
TF002.remlVCA.EP05_A2_intermediate_precision.balanced <- function()
{        
	res11 <- remlVCA(y~day/run, Data=dataEP05A2_1)                                            # call function    
	res01 <- anovaVCA(y~day/run, Data=dataEP05A2_1)
	checkEquals(as.numeric(res11$aov.tab[,"VC"]), as.numeric(res01$aov.tab[,"VC"]), tolerance=1e-7)  
	checkEquals(as.numeric(res11$aov.tab[c("total", "error"),"DF"]), as.numeric(res01$aov.tab[c("total", "error"),"DF"]), tolerance=1e-7)   
	
	res12 <- remlVCA(y~day/run, Data=dataEP05A2_2)                                            # call function    
	res02 <- anovaVCA(y~day/run, Data=dataEP05A2_2)
	checkEquals(as.numeric(res12$aov.tab[,"VC"]), as.numeric(res02$aov.tab[,"VC"]), tolerance=1e-7)     
	checkEquals(as.numeric(res12$aov.tab[c("total", "error"),"DF"]), as.numeric(res02$aov.tab[c("total", "error"),"DF"]), tolerance=1e-7) 
	
	res13 <- remlVCA(y~day/run, Data=dataEP05A2_3)                                            # call function    
	res03 <- anovaVCA(y~day/run, Data=dataEP05A2_3)
	checkEquals(as.numeric(res13$aov.tab[,"VC"]), as.numeric(res03$aov.tab[,"VC"]), tolerance=1e-7)    
	checkEquals(as.numeric(res13$aov.tab[c("total", "error"),"DF"]), as.numeric(res03$aov.tab[c("total", "error"),"DF"]), tolerance=1e-7) 
}

# unbalanced case, check against SAS PROC MIXED REML-estimator
# below, 'data3' corresponds to dataEP05A2_3
#
# data data3ub;
# 	set data3;
# 	Nobs = _N_;
# 	if Nobs in (1,6, 7, 36, 61, 62, 63, 64, 65) then delete;
# run;
# 
# proc mixed data=data3ub;
# 	class day run;
# 	model y =;
# 	random day day*run;
# 	ods output covparms=covparms;
# run;

TF003.remlVCA.EP05_A2_intermediate_precision.unbalanced <- function()
{ 
	res <- remlVCA(y~day/run, Data=dataEP05A2_1[-c(11, 12, 17, 37, 45, 56, 57, 68),])              
	checkEquals(as.numeric(res$aov.tab[-1,"VC"]), c(0.28354702238093, 0.90581093706202, 2.06774293438701), tolerance=3e-4)     
	
	res <- remlVCA(y~day/run, Data=dataEP05A2_2[-c(2, 12, 22, 23, 24, 55, 56, 71),])              
	checkEquals(as.numeric(res$aov.tab[-1,"VC"]), c(1.6266257270871, 3.24912440247053, 3.85555348892268 ), tolerance=3e-4)  
	
	res <- remlVCA(y~day/run, Data=dataEP05A2_3[-c(1,6,7,36,61:65),])              
	checkEquals(as.numeric(res$aov.tab[-1,"VC"]), c(10.8147804686742, 8.1948752793538, 15.1256464181605), tolerance=3e-4) 
}


TF004.remlVCA.EP05_A3_Reproducibility.balanced <- function()
{      
	res <- remlVCA(y~site/day, Data=dataEP05A3_MS_1)
	checkEquals(round(as.numeric(res$aov.tab[,"VC"]), 6), c( 3.635232, 1.017401, 0.345882, 2.271949))
	
	res <- remlVCA(y~site/day, Data=dataEP05A3_MS_2)
	checkEquals(round(as.numeric(res$aov.tab[,"VC"]), 6), c(5.765516, 1.010798, 1.028309, 3.726408))
	
	res <- remlVCA(y~site/day, Data=dataEP05A3_MS_3)
	checkEquals(round(as.numeric(res$aov.tab[,"VC"]), 4), c(38.4389, 17.9908, 5.4237, 15.0245))
}   

## RS0003 21 replications, testing only the residual error component
## WRI = within run imprecision

TF005.remlVCA.WRI <- function()
{ 
	res <- remlVCA(y~1, Data=dataRS0003_1)
	checkEquals( round(as.numeric(res$aov.tab[2,c("SS", "MS")]),5), c(22.44452, 1.12223) )
	
	res <- remlVCA(y~1, Data=dataRS0003_2)
	checkEquals( round(as.numeric(res$aov.tab[2,c("SS", "MS")]),4), c(367.3480, 18.3674) )
	
	res <- remlVCA(y~1, Data=dataRS0003_3)
	checkEquals( round(as.numeric(res$aov.tab[2,c("SS", "MS")]),2), c(2210.59, 110.53) )
} 

## RS0005 - Confirmation of internal data by external labs - balanced
## BDI = between day imprecision

TF006.remlVCA.BDI.external_labs.balanced <- function()
{     
	res <- remlVCA(y~day, Data=dataRS0005_1)
	checkEquals( round(as.numeric(res$aov.tab[,"VC"]), 6), c(5.193341, 1.818128, 3.375212))
	
	res <- remlVCA(y~day, Data=dataRS0005_2)
	checkEquals( round(as.numeric(res$aov.tab[,"VC"]), 6), c(6.611960, 0.055541, 6.556419))
	
	res <- remlVCA(y~day, Data=dataRS0005_3)
	checkEquals( round(as.numeric(res$aov.tab[,"VC"]), 5), c(75.69922, 22.45171, 53.24751))
} 

# check numerical equivalence of a crossed-nested design (balanced)
#
# check vs. SAS PROC MIXED results using this SAS-code:
#
#   proc mixed data=sample1 method=type1 cl;
#       class lot device day run;
#       model y=;
#       random lot device lot*device*day lot*device*day*run;
#   run;

TF007.remlVCA.crossed_nested.balanced <- function()
{
	data(VCAdata1)
	sample1 <- VCAdata1[which(VCAdata1$sample==1),]
	sample1$device <- gl(3,28,252)                                      # add device variable
	set.seed(505)
	sample1$y <- sample1$y + rep(rep(rnorm(3,,.25), c(28,28,28)),3)     # add error component according to 
	
	# write.table(sample1, file="sample1.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")      # export data, import to SAS and apply SAS-code given above
	
	res1 <- remlVCA(y~lot+device+(lot:device:day)/run, sample1)
	
	checkEquals(round(res1$aov.tab[2, "VC"], 5), 0.01552)                        # round to precision of SAS PROC MIXED output
	checkEquals(round(res1$aov.tab[3, "VC"], 5), 0.06214)
	checkEquals(round(res1$aov.tab[4, "VC"], 5), 0.01256)
	checkEquals(round(res1$aov.tab[5, "VC"], 5), 0.05074)
	checkEquals(round(res1$aov.tab[6, "VC"], 6), 0.001152)
	
}


# use subset "sample_8" of VCAdata1

TF008.remlVCA.crossed_nested.balanced <- function()
{
	data(VCAdata1)
	sample8 <- VCAdata1[which(VCAdata1$sample==8),]
	sample8$device <- gl(3,28,252)                                      # add device variable
	set.seed(903211)
	sample8$y <- sample8$y + rep(rep(rnorm(3,,.75), c(28,28,28)),3)     # add error component according to 
	
	# write.table(sample8, file="sample8.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")      # export data, import to SAS and apply SAS-code given above
	
	res8 <- remlVCA(y~lot+device+(lot:device:day)/run, sample8)
	
	checkEquals(round(res8$aov.tab[2, "VC"], 4), 7.3463)                        # round to precision of SAS PROC MIXED output
	checkEquals(round(res8$aov.tab[3, "VC"], 4), 2.2716)
	checkEquals(round(res8$aov.tab[4, "VC"], 4), 3.4588)
	checkEquals(round(res8$aov.tab[5, "VC"], 4), 2.0302)
	checkEquals(round(res8$aov.tab[6, "VC"], 4), 1.2219)
	
}

TF009.remlVCA.zeroVariance <- function()
{
	data(dataEP05A2_3)
	dat1 <- dataEP05A2_3
	dat1$y <- dat1[1,"y"]
	dat1$cov <- rnorm(nrow(dat1),15,3)
	
	fit1 <- remlVCA(y~day, dat1)
	checkEquals(as.numeric(fit1$aov.tab[,"VC"]), rep(0,3))
	fit2 <- remlVCA(y~day/run, dat1)
	checkEquals(as.numeric(fit2$aov.tab[,"VC"]), rep(0,4))
}


TF010.remlVCA.byProcessing <- function()
{
	# reproducible example taken from ticket 12818
	set.seed(42)
	dat <- expand.grid(group=c("group1", "group2"), day=rep(c("day1","day2"), times=3))
	dat$value <- rnorm(nrow(dat))
	
	dat <- subset(dat, !(group=="group2" & day=="day1")) #dataset with valid vca input for group 1 but invalid input for group 2
	
#	VCA::anovaVCA(value~day, Data = subset(dat, group=="group1"), by="group") #valud result (as epcected)
#	VCA::anovaVCA(value~day, Data = subset(dat, group=="group2"), by="group") #invalid result (as expected)
	
	# should generate a warning if quiet=FALSE
	checkTrue(	tryCatch(
					remlVCA(value~day, Data = dat, by = "group", quiet=FALSE),
					warning=function(w) TRUE) )
	# should not generate a warning if quiet=TRUE
	res <- remlVCA(value~day, Data = dat, by = "group", quiet=TRUE)
	checkTrue(	class(res) == "list" && 
					length(res) == 2 &&  
					all.equal(c("VCA", "try-error"), as.character(sapply(res, class))))
}


