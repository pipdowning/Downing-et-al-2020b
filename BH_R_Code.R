
# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #
#																					  					 #
#	Supplementary R code for:														  					 #
#																					  					 #
#				THE BENEFITS OF HELP IN COOPERATIVE BIRDS - NON-EXISTENT OR DIFFICULT TO DETECT?			 #
#																					  					 #
#	Philip A. Downing, Ashleigh S. Griffin and Charlie K. Cornwallis				  					 	 #
#																					  					 #
#	contact: philip.downing@biol.lu.se												  					 #
#																					  					 #
# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #

# packages

library(ape)
library(MCMCglmm)
library(metafor)


# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #

### DATA MANIPULATION ###

## polyandry ##
zrData <- read.csv("enter file path here")

## phylogenetic trees ##
trees <- read.nexus("enter file path here")
is.ultrametric(trees[[1]])

## trim the trees ##
dropTip <- trees[[1]]$tip.label[which(is.na(match(trees[[1]]$tip.label, zrData$animal)))]
zrTrees <- lapply(trees, drop.tip, dropTip, trim.internal=T)
zrTrees <- lapply(zrTrees, makeNodeLabel, method = "number")
zrData$animal[which((zrData$animal %in% zrTrees[[1]]$tip.label) == FALSE)]
zrTrees[[1]]$tip.label[which((zrTrees[[1]]$tip.label %in% zrData$animal) == FALSE)]
# reduced to 19 species


# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #

### WITHIN-STUDY VARIANCE ###

## need this in order to calculate I^2 values
## see Nakagawa & Santos (2012) Methodological issues and advances in biological meta-analysis. Evol. Ecol 26: 1253- 1274
## within-study variance (this is known)
wj <- 1/zrData$Zvar										# inverse measurement error / sampling error variance
k <- length(zrData$Zr)									# number of studies
within <- sum(wj * (k-1)) / (sum(wj)^2 - sum(wj^2))		# 0.046

### OTHER VARIANCE PARAMETERS FROM MCMCglmm ###
# between <- mean(model$VCV[,X]) i.e. the between-study variance (estimated from the data)
# phylo <- mean(model$VCV[,X]) i.e. phylogenetic variance (estimated from the data)


# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #

### PART 1: HELPER BENEFITS WITHOUT PHYLOGENY ###

## A. Mean effect size ##

# metafor #
meanModel1a <- rma(Zr ~ 1, vi=Zvar, data=zrData)
summary(meanModel1a)   				# B = 0.36, lwr = 0.20, upr = 0.51	(Q = 32.8, p = 0.02)
forest(meanModel1a)					# 3/19 studies < 0 (i.e. helpers have a negative effect on r.s.)
trimfill(meanModel1a)				# 3 studies missing from the left-hand side
funnel(trimfill(meanModel1a))		# funnel plot with trim and fill
zrData$invSE <- 1/sqrt(zrData$Zvar)
summary(lm(scale(Zr) ~ invSE, data=zrData))		# Egger's test: intercept = 0.43 (p = 0.36), slope = -0.10 (p = 0.29)
(exp(2*0.3561)-1) / (exp(2*0.3561)+1)			# r = 0.34 (medium effect size)
# between study-variance = 0.0465
# I2.between = 50.06 % (or 0.0465 / (0.0465 + within))


# MCMCglmm #
priorA <- list(R=list(V=1, nu=0.002))
MEV <- zrData$Zvar
meanModel1b <- MCMCglmm(Zr ~ 1, data=zrData, pl=TRUE, slice=TRUE, nitt=1100000, thin=1000, burnin=100000, mev=MEV, prior=priorA, verbose=FALSE)
summary(meanModel1b)
posterior.mode(meanModel1b$Sol)						# Zr = 0.34
HPDinterval(meanModel1b$Sol)							# lwr = 0.20, upr = 0.51
between <- mean(meanModel1b$VCV[,2])				# 0.04
I2.between <- between / (between + within) * 100		# 48.4 %

# run the above three times
meanModel1b.1 <- meanModel1b
meanModel1b.2 <- meanModel1b
meanModel1b.3 <- meanModel1b
save(meanModel1b.1, file="enter file path here")
save(meanModel1b.2, file="enter file path here")
save(meanModel1b.3, file="enter file path here")

# chain convergence
hist(meanModel1b.1$Liab)
plot(meanModel1b.1$VCV) 	   	# units close to 0
plot(meanModel1b.1$Sol)		   	# intercept estimate well mixed
autocorr(meanModel1b.1$VCV)   	# correlation between successive samples < 0.1 for units
autocorr(meanModel1b.1$Sol)   	# correlation between successive samples < 0.1 for intercept
meanModel1bSols <- mcmc.list(list(meanModel1b.1$Sol, meanModel1b.2$Sol, meanModel1b.3$Sol))
plot(meanModel1bSols)
gelman.diag(meanModel1bSols)      # upper CI = 1.00 for intercept suggesting convergence
heidel.diag(meanModel1b.1$VCV)    # units passed halfwidth
heidel.diag(meanModel1b.1$Sol)    # intercept passed halfwidth



## B. The effect of study design ##

# metafor #
sdModel1a <- rma(Zr ~ study.class-1, vi=Zvar, data=zrData)		# note that suppressing the intercept loses R^2 value
summary(sdModel1a)					# Qm = 76.2 (p < 0.001)
# paired comparisons = -0.17 < 0.04 < 0.24 (+/- 0.11 se), n = 9
# multiple regression = 0.33 < 0.45 < 0.57 (+/- 0.06 se), n = 5
# removal experiments = 0.23 < 0.48 < 0.72 (+/- 0.13 se), n = 3
# other controls = 0.19 < 0.49 < 0.80 (+/- 0.16 se), n = 2
# the 5 lowest values including all three negative ones were from paired comparisons (as were the two highest)
# between study-variance = 0.000 (se = 0.013)
# I2.between = 0.00 % (or 0.000 / (0.000 + within))


# MCMCglmm #
# use the prior and MEV term from meanModel1b
sdModel1b <- MCMCglmm(Zr ~ study.class-1, data=zrData, pl=TRUE, slice=TRUE, nitt=1100000, thin=1000, burnin=100000, mev=MEV, prior=priorA, verbose=FALSE)
summary(sdModel1b)
posterior.mode(sdModel1b$Sol)
HPDinterval(sdModel1b$Sol)
# paired comparisons = -0.18 < 0.05 < 0.27
# multiple regression = 0.27 < 0.47 < 0.65
# removal experiments = 0.17 < 0.51 < 0.77
# other controls = 0.13 < 0.46 < 0.88
between <- mean(sdModel1b$VCV[,2])					# 0.019
I2.between <- between / (between + within) * 100		# 28.8 %

# run the above three times
sdModel1b.1 <- sdModel1b
sdModel1b.2 <- sdModel1b
sdModel1b.3 <- sdModel1b
save(sdModel1b.1, file="enter file path here")
save(sdModel1b.2, file="enter file path here")
save(sdModel1b.3, file="enter file path here")

# chain convergence
hist(sdModel1b.1$Liab)
plot(sdModel1b.1$VCV)     		# units close to 0
plot(sdModel1b.1$Sol)		   	# four fixed effect estimates well mixed
autocorr(sdModel1b.1$VCV)   		# correlation between successive samples < 0.1 for units
autocorr(sdModel1b.1$Sol)  	 	# correlation between successive samples < 0.1 for all components
sdModel1bSols <- mcmc.list(list(sdModel1b.1$Sol, sdModel1b.2$Sol, sdModel1b.3$Sol))
plot(sdModel1bSols)
gelman.diag(sdModel1bSols)      # upper CI = 1.00 for x4 fixed effect estimates suggesting convergence
heidel.diag(sdModel1b.1$VCV)    # units passed halfwidth
heidel.diag(sdModel1b.1$Sol)    # fixed effect estimates passed halfwidth


# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #

### PART 2: HELPER BENEFITS WITH PHYLOGENY ###

## A. Mean effect size ##

# metafor #
zrData$ID <- zrData$animal
birdCor <- vcv.phylo(zrTrees[[sample(1:1300, 1)]], cor=TRUE)
meanModel2a <- rma.mv(Zr ~ 1, V=Zvar, random = list(~ 1 | ID, ~ 1 | animal), R = list(animal=birdCor), data=zrData)
summary(meanModel2a)   						# B = 0.36, lwr = 0.20, upr = 0.51	(Q = 32.8, p = 0.02)
(exp(2*0.3561)-1) / (exp(2*0.3561)+1)		# r = 0.34 (medium effect size)
# between study-variance = 0.0465
# phylogenetic variance = 0.0000
I2.between <- 0.0465 / (0.0465 + within + 0.000)		# 50.1 %
I2.phylo <- 0.000 / (0.0465 + within + 0.000)		# 0.0 %


# MCMCglmm #
priorB <- list(R=list(V=1,nu=0.002), G=list(G1=list(V=1,nu=0.002)))
MEV <- zrData$Zvar
INtree <- inverseA(zrTrees[[1]], nodes="TIPS")
interceptModel.start <- MCMCglmm(Zr ~ 1, random = ~ animal, data=zrData, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=0, mev=MEV, prior=priorB, ginverse=list(animal=INtree$Ainv), verbose=FALSE)
interceptModel <- interceptModel.start
for(i in 1:1300){
  INtree <- inverseA(zrTrees[[i]], nodes="TIPS")
  start <- list(Liab=interceptModel$Liab[1,], R=interceptModel$VCV[1,3], G=list(G1=interceptModel$VCV[1,1]))
  interceptModel <- MCMCglmm(Zr ~ 1, random = ~ animal, data=zrData, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=999, mev=MEV, prior=priorB, ginverse=list(animal=INtree$Ainv), start=start, verbose=FALSE)
  if(i>300){
    interceptModel.start$VCV[i-300,] <- interceptModel$VCV[1,]
    interceptModel.start$Sol[i-300,] <- interceptModel$Sol[1,]
    interceptModel.start$Liab[i-300,] <- interceptModel$Liab[1,]
  }
}
meanModel2b.1 <- interceptModel.start
meanModel2b.2 <- interceptModel.start
meanModel2b.3 <- interceptModel.start
save(meanModel2b.1, file="enter file path here")
save(meanModel2b.2, file="enter file path here")
save(meanModel2b.3, file="enter file path here")

# chain convergence
hist(meanModel2b.1$Liab)
plot(meanModel2b.1$VCV)     		# animal and units close to 0
plot(meanModel2b.1$Sol)		   	# intercept estimate well mixed
autocorr(meanModel2b.1$VCV)   	# correlation between successive samples < 0.1 for all components
autocorr(meanModel2b.1$Sol)   	# correlation between successive samples < 0.1 for all components
meanModel2bSols <- mcmc.list(list(meanModel2b.1$Sol, meanModel2b.2$Sol, meanModel2b.3$Sol))
plot(meanModel2bSols)
gelman.diag(meanModel2bSols)      # upper CI = 1.00 for intercept suggesting convergence
heidel.diag(meanModel2b.1$VCV)    # units passed halfwidth
heidel.diag(meanModel2b.1$Sol)    # intercept passed halfwidth

# model parameters
summary(meanModel2b.1)
posterior.mode(meanModel2b.1$Sol)		# intercept = 0.36
HPDinterval(meanModel2b.1$Sol)			# intercept = 0.11 to 0.56
(exp(2*0.3547)-1) / (exp(2*0.3547)+1)	# r = 0.34 (medium effect size)

# I^2 values
between <- mean(meanModel2b.1$VCV[,3])
animal <- mean(meanModel2b.1$VCV[,1])
I2.between <- between / (between + within + animal) * 100	# 28.6 %
I2.animal <- animal / (between + within + animal) * 100		# 31.6 %



## B. The effect of study design ##

# metafor #
sdModel2a <- rma.mv(Zr ~ study.class-1, V=Zvar, random = list(~ 1 | ID, ~ 1 | animal), R = list(animal=birdCor), data=zrData)
summary(sdModel2a)  					# Qm = 76.2 (p < 0.001)
# paired comparisons = -0.17 < 0.04 < 0.24 (+/- 0.11 se), n = 9
# multiple regression = 0.33 < 0.45 < 0.57 (+/- 0.06 se), n = 5
# removal experiments = 0.23 < 0.48 < 0.72 (+/- 0.13 se), n = 3
# other controls = 0.19 < 0.49 < 0.80 (+/- 0.16 se), n = 2
# between study-variance = 0.000
# phylogenetic variance = 0.000
I2.between <- 0.000 / (0.000 + within + 0.000)		# 0.0 %
I2.phylo <- 0.000 / (0.000 + within + 0.000)		# 0.0 %


# MCMCglmm #
# use the prior and MEV term from meanModel2b
INtree <- inverseA(zrTrees[[1]], nodes="TIPS")
fixedModel.start <- MCMCglmm(Zr ~ study.class-1, random = ~ animal, data=zrData, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=0, mev=MEV, prior=priorB, ginverse=list(animal=INtree$Ainv), verbose=FALSE)
fixedModel <- fixedModel.start
for(i in 1:1300){
  INtree <- inverseA(zrTrees[[i]], nodes="TIPS")
  start <- list(Liab=fixedModel$Liab[1,], R=fixedModel$VCV[1,3], G=list(G1=fixedModel$VCV[1,1]))
  fixedModel <- MCMCglmm(Zr ~ study.class-1, random = ~ animal, data=zrData, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=999, mev=MEV, prior=priorB, ginverse=list(animal=INtree$Ainv), start=start, verbose=FALSE)
  if(i>300){
    fixedModel.start$VCV[i-300,] <- fixedModel$VCV[1,]
    fixedModel.start$Sol[i-300,] <- fixedModel$Sol[1,]
    fixedModel.start$Liab[i-300,] <- fixedModel$Liab[1,]
  }
}
sdModel2b.1 <- fixedModel.start
sdModel2b.2 <- fixedModel.start
sdModel2b.3 <- fixedModel.start
save(sdModel2b.1, file="enter file path here")
save(sdModel2b.2, file="enter file path here")
save(sdModel2b.3, file="enter file path here")

# chain convergence
hist(sdModel2b.1$Liab)
plot(sdModel2b.1$VCV)    	 	# animal and units close to 0
plot(sdModel2b.1$Sol)		   	# four fixed effect estimates well mixed
autocorr(sdModel2b.1$VCV)   		# correlation between successive samples < 0.1 for all components
autocorr(sdModel2b.1$Sol)   		# correlation between successive samples < 0.1 for all components
sdModel2bSols <- mcmc.list(list(sdModel2b.1$Sol, sdModel2b.2$Sol, sdModel2b.3$Sol))
plot(sdModel2bSols)
gelman.diag(sdModel2bSols)      # upper CI = 1.00 for x4 fixed effect estimates suggesting convergence
heidel.diag(sdModel2b.1$VCV)    # animal and units failed halfwidth
heidel.diag(sdModel2b.1$Sol)    # paired comparisons failed halfwidth

# model parameters
summary(sdModel2b.1)
posterior.mode(sdModel2b.1$Sol)
HPDinterval(sdModel2b.1$Sol)
# paired comparisons = -0.20 < 0.13 < 0.36
# multiple regression = 0.20 < 0.43 < 0.68
# removal experiments = 0.18 < 0.49 < 0.86
# other controls = 0.06 < 0.59 < 0.87
table(sdModel2b.1$Sol[, 1] > sdModel2b.1$Sol[, 2]) / length(sdModel2b.1$Sol[, 1])	# F = 0.98, T = 0.02
table(sdModel2b.1$Sol[, 1] > sdModel2b.1$Sol[, 3]) / length(sdModel2b.1$Sol[, 1])	# F = 0.96, T = 0.04
table(sdModel2b.1$Sol[, 1] > sdModel2b.1$Sol[, 4]) / length(sdModel2b.1$Sol[, 1])	# F = 0.98, T = 0.02
# paired comparisons significanlty lower than all the others

# I^2 values
between <- mean(sdModel2b.1$VCV[,3])
animal <- mean(sdModel2b.1$VCV[,1])
I2.between <- between / (between + within + animal) * 100	# 22.3 %
I2.animal <- animal / (between + within + animal) * 100		# 22.3 %


# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #
#                                                                                     					 #
#                            		END - thanks for reading! p. 										 #
#                                                                                     					 #  
# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #