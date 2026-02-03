########################################################################
#
#           --- Simulation: low-dimensional setting  ---
#
########################################################################
#
# NOTE: WE DID NOT SET SEEDS, SO SIMULATION RESULTS MAY SLIGHTLY DIFFER
# FROM THOSE IN THE MANUSCRIPT.
#
# Script produces simulation results in the manuscript, including
# plots with marginal symmetric, conditional symmetric and conditional asymmetric
# Shapley values (estimated by SHAP) for genes G, Disease status D, and confounders C1 & C2. 
# Also, global feature importances by SAGE-R2 decomposition for all three
# versions. Finally, it shows the comparative results for the asymmetric,
# refitting based Shapley values (instead of SHAP). 

# Run entire script by entering (will take a few minutes)
#source("SimExample_lowdim.R")

#set working directory
setwd("C:/Synchr/Rscripts/ShapleyDependency/GeneShap/")

#libraries
library(shapr)
library(data.table) # for shapr
library(mgcv) #for fitting gam
library(ggplot2)
library(ggpubr)

# Additional functions for computing asymmetric Shapley values
source("GeneShap_source.R")

#Set up simulation
alpha0 <- alpha3 <- 0
alpha1 <-  alpha2 <- 1
alpha4 <- 2
beta1 <- 2

n <- 1200 #train + test
epsilon <- rnorm(n)

G <- rnorm(n)
C0 <- rnorm(n)
S0 <- rnorm(n)

#for creating correlation between D and C2
U2 <- rnorm(n)

D <- (S0 + beta1*G + U2)/sqrt(6)
C1 <- rnorm(n) #uncorrelated with all other covariates
C2 <- (C0 + 2*U2)/sqrt(5)

#outcome model: regression with squared effect for Disease state D
Y <- alpha0 + alpha1*G  + alpha2*C1 + alpha3*C2 + alpha4*(D^2-mean(D^2)) + epsilon

#training-test split
dat <- data.frame(Y, C1, C2, D, G)
xy_train <- dat[1:1000,]
xy_explain <- dat[-(1:1000),]

#Model fit
print("Fit linear regression model with spline for D (default gam)")
cat("\n")
lmY <- gam(Y ~ G + C1 + C2 + s(D), data = xy_train)
phi0 <- mean(predict(lmY))


#######################################################################
#
#                       Shapley values based on SHAP
#
#######################################################################
print("Computing Shapley values")
Genes <- c("G")
Status <- c("D")
Confounders <- c("C1","C2")

Q1 <- R_D_Matrix(Genes, Status, Confounders,
                 Ordering_between = list("Genes", "Status"),
                 verbose = FALSE)

print("Asymmetric Shapley weight matrix")
print(Q1)
cat("\n")

rowSums_pos <- apply(Q1, 1, function(x) sum(x[x > 0], na.rm = TRUE)) # check if properly normalized (positive weights)
rowSums_neg <- apply(Q1, 1, function(x) sum(x[x < 0], na.rm = TRUE)) # check if properly normalized (negative weights)

print("Check normalization")
print(rowSums_pos) #sum positive weights
print(rowSums_neg) #sum negative weights
cat("\n")

print("Compute Asymmetric Shapley values; ~1 min")
groups <- list("C1", "C2", "D", "G")
shap_names <- c("C1", "C2", "D", "G")
names(groups) <- shap_names

coalition_list <- parse_coalitions(shap_names, names(Q1))

ShapAsymSim <- explain_manual(lmY,
                           x_explain = xy_explain[,-1],
                           x_train = xy_train[,-1],
                           approach = "gaussian",
                           phi0 = phi0,
                           coalition_list = coalition_list,
                           R_D = as.matrix(Q1),
                           group=groups
)

print("Compute Symmetric Shapley values; ~1 min")
ShapSymSim <- explain(lmY,
                           x_explain = xy_explain[,-1],
                           x_train = xy_train[,-1],
                           approach = "gaussian",
                           phi0 = phi0
)

print("Compute Marginal Shapley values (independence); ~1 min")
ShapIndSim <- explain(lmY,
                   x_explain = xy_explain[,-1],
                   x_train = xy_train[,-1],
                   approach = "independence",
                   phi0 = phi0
)

save(ShapIndSim, ShapAsymSim, ShapSymSim, lmY, xy_explain, xy_train, file="ShapLowDimSim.Rdata" )

aShap0 <- copy(ShapAsymSim$shap)
aShap0[,explain_id:=.I]
aShapSim <- ShapSymSim
aShapSim$shapley_values_est <- aShap0

print("Plotting Shapley values...")
windows()
par(mfrow=c(3,1))
mi <- plot(ShapIndSim, plot_type="scatter",scatter_hist=F) + ylim(-2,3) + xlab("Symm-marginal")
mg <- plot(ShapSymSim, plot_type="scatter",scatter_hist=F) + ylim(-2,3) + xlab("Symm-conditional") + ylab("")
mo <- plot(aShapSim, plot_type="scatter",scatter_hist=F) + ylim(-2,3) + xlab("Asymmetric") + ylab("")
ggarrange(mi,mg,mo,nrow=1)

print("Compute SAGE-R2")
#Sage-R2, conditional, Asymmetric
vS <- ShapAsymSim$nu
Y_explain <- xy_explain[,1]
R2 <- function(row) {
  #row <- vS[2,]
  rij <- as.numeric(row)[-1]
  r2 <- 1 - sum((Y_explain - rij)^2)/sum((Y_explain - mean(Y_explain))^2)
  return(r2)
}

R2s_coal <- apply(vS, 1,  R2)
ASAGER2 <- round(as.matrix(Q1) %*% matrix(R2s_coal,ncol=1),3)
ASAGER2b <- ASAGER2[c(1,4,5,3,2),] #reshuffle column order to align with below

#Sage-R2, conditional, symmetric
nm <- ShapSym$internal$objects$feature_specs$labels
W <- ShapSym$internal$objects$W
vS <- ShapSym$internal$output$dt_vS
rownames(W) <- c("Intercept",nm)
R2s_coal <- apply(vS, 1,  R2)
SAGER2 <-  round(as.matrix(W) %*% matrix(R2s_coal,ncol=1),3)

#Sage-R2, marginal (independence)
nm <- ShapInd$internal$objects$feature_specs$labels
W <- ShapInd$internal$objects$W
vS <- ShapInd$internal$output$dt_vS
rownames(W) <- c("Intercept",nm)
R2s_coal <- apply(vS, 1,  R2)
SAGEIndR2 <-  round(as.matrix(W) %*% matrix(R2s_coal,ncol=1),3)

SAGER2all <- data.frame(SAGER2marg=SAGEIndR2, SAGER2symm = SAGER2, SAGER2asym = ASAGER2b) 
Total <- colSums(SAGER2all)
SAGER2all <- rbind(SAGER2all, Total=Total)
print("SAGE-R2s")
print(SAGER2all)

########################################################################
#
#                             Refitting
#
########################################################################

print("Refitting version of Shapley")
#refit prediction model for all allowed coalitions
allpred <- c()
mod0 <- lm(Y ~ 1,data=xy_train)
allpred <- cbind(allpred, pred0=predict(mod0, newdata=xy_explain))
modG <- lm(Y ~ G,data=xy_train)
allpred <- cbind(allpred, predG=predict(modG, newdata=xy_explain))
modC1 <- lm(Y ~ C1,data=xy_train)
allpred <- cbind(allpred, predC1 =predict(modC1, newdata=xy_explain))
modC2 <- lm(Y ~ C2,data=xy_train)
allpred <- cbind(allpred, predC2 = predict(modC2, newdata=xy_explain))
modGD<- gam(Y ~ G + s(D),data=xy_train)
allpred <- cbind(allpred, predGD = predict(modGD, newdata=xy_explain))
modGC1 <- lm(Y ~ C1 + G,data=xy_train)
allpred <- cbind(allpred, predGC1 = predict(modGC1, newdata=xy_explain))
modGC2 <- lm(Y ~ C2 + G,data=xy_train)
allpred <- cbind(allpred, predGC2 = predict(modGC2, newdata=xy_explain))
modC1C2 <- lm(Y ~ C1 + C2,data=xy_train)
allpred <- cbind(allpred, predC1C2 = predict(modC1C2, newdata=xy_explain))
modGDC1<- gam(Y ~ G + C1 + s(D),newdat=xy_train)
allpred <- cbind(allpred, predGDC1 = predict(modGDC1, newdata=xy_explain))
modGDC2<- gam(Y ~ G + C2 + s(D),data=xy_train)
allpred <- cbind(allpred, predGDC2 = predict(modGDC2, newdata=xy_explain))
modGC1C2 <- lm(Y ~ G + C1 + C2,data=xy_train)
allpred <- cbind(allpred, predGC1C2 = predict(modGC1C2, newdata=xy_explain))
modall<- gam(Y ~ G + C1 + C2 + s(D),data=xy_train)
allpred <- cbind(allpred, predall = predict(modall, newdata=xy_explain))

#Compute refitting version of asymmetric Shapley
refitshap <- as.matrix(Q1) %*% t(allpred)

#Plot asymmetric Shapley values based on refitting, showing that this 
#does not work, as it is not able to retrieve the quadratic
#relationships when D is not present in the model
print("Plotting refitting-based asymmetric Shapley values")
windows()
par(mfrow=c(2,2), mar=c(2,2,2,1))
plot(xy_explain$C1,refitshap[4,],main="C1", xlab="",ylab="", ylim=c(-2,3)) 
plot(xy_explain$C2,refitshap[5,],main="C2", xlab="",ylab="", ylim=c(-2,3)) 
plot(xy_explain$D,refitshap[3,],main="D", xlab="",ylab="", ylim=c(-2,3)) 
plot(xy_explain$G,refitshap[2,],main="G", xlab="",ylab="", ylim=c(-2,3)) 


