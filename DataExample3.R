########################################################################
## --- Data Example, part III: Illustrations for explaining blockForest model
## --- 1. Asymmetric global feature importance: SAGE-C and SAGE-AUC
## --- 2. Comparing Asymmetric with Symmetric 
## --- 3. Global inference with local Shapley values
## --- 4. Mediation effect of disease state
########################################################################

# NOTE: AS RESULTS FROM TREE ENSEMBLE LEARNERS (like blockForest) ARE NOT 
# DETERMINISTIC, RESULTS MAY DEVIATE SOMEWHAT FROM THOSE IN THE MANUSCRIPT. 
# QUALITATIVELY, HOWEVER, THEY SHOULD BE VERY SIMILAR.

#'set working directory
setwd("C:/Synchr/Rscripts/ShapleyDependency/GeneShap/")

#'new GeneShap functions
source("GeneShap_source.R")

# load libraries
library(survival)
library(shapr)
library(data.table) #for shapr
library(future) #for shapr
library(progressr) #for shapr
library(blockForest)

# for performance evaluation on survival data (C-index and time-dep AUC)
library(survC1)
library(survivalROC)

# for non-parametric inference
library(coin)

#######

#load data objects (see DataExample1.R)
#show(load("data.Rdata"))
load("data.Rdata")

#load objects containing Shapley values for blockForest model (see DataExample2.R)
#show(load("Shapleys.Rdata"))
load("Shapleys.Rdata")

#asymmetric and symmetric Shapley values for all features and samples
aShapbf <- ShapAsym$shap
Shapbf <- ShapSym$shap

#####################################################
#
#         Global feature importance
#
#####################################################

# Global feature importance by SAGE. SAGE applies an predictive 
# performance evaluation function (like C-index) to all coalitions. 
# For that it uses the same weights as Shapley, and the predictions for each 
# coalition, which are stored in the Shapley object.

print("ASSESSING GLOBAL FEATURE IMPORTANCE BY SAGE")
cat("\n")

# Retrieve marginalized predictions for all coalitions
vS <- as.matrix(ShapAsym$nu[,-1])
rownames(vS) <- colnames(R_D)

# Compute performance for all coalitions
Cinds_coal <- apply(vS, 1, function(row) Est.Cval(cbind(Y_explain[,1], Y_explain[,2], row), 
                                                  tau = 8, nofit=T)$Dhat)
AUCs_coal <- apply(vS, 1, function(row) survivalROC.C(Stime = Y_explain[,1], status = Y_explain[,2], 
                                                      marker = row, predict.time = 5)$AUC)

# Compute SAGE, by weighing coalitions according to Shapley weigths
SAGEAsymC <- as.matrix(R_D) %*% matrix(Cinds_coal,ncol=1)
SAGEAsymAUC <- as.matrix(R_D) %*% matrix(AUCs_coal,ncol=1)

SAGEAsym <- data.frame(Cindex=SAGEAsymC,tAUC = SAGEAsymAUC)
total <- colSums(SAGEAsym)
SAGEAsym <- rbind(SAGEAsym,total=total)

#print(round(SAGEAsym,3))

#Same, SAGE for Symmetric
vS <- as.matrix(ShapSym$nu[,-1])
rownames(vS) <- colnames(R_D_Sym)
Cinds_coal <- apply(vS, 1, function(row) Est.Cval(cbind(Y_explain[,1], Y_explain[,2], row), 
                                                  tau = 8, nofit=T)$Dhat)
AUCs_coal <- apply(vS, 1, function(row) survivalROC.C(Stime = Y_explain[,1], status = Y_explain[,2], 
                                                      marker = row, predict.time = 5)$AUC)

SAGESymC <- as.matrix(R_D_Sym) %*% matrix(Cinds_coal,ncol=1)
SAGESymAUC <- as.matrix(R_D_Sym) %*% matrix(AUCs_coal,ncol=1)

SAGESym <- data.frame(Cindex=SAGESymC,tAUC = SAGESymAUC)
total <- colSums(SAGESym)
SAGESym <- rbind(SAGESym,total=total)

# Compare Asymmetric and Symmetric SAGE values
SAGEs <- data.frame(SAGEAsym,SAGESym)
cns <- colnames(SAGEAsym)
colnames(SAGEs) <- c(paste(cns,"_A",sep=""), paste(cns,"_S",sep=""))
print("Asymmetric (A) and Symmetric (S) SAGE values: decomposition of C-index and tAUC")
print(round(SAGEs,3)) #Asym and Sym SAGE values
cat("\n")

#####################################################
#
#                      Inference 
#
#####################################################

# Global importance: likelihood ratios
# As Shapley is a linear decomposition of the log-cumulative Hazard (for block Forest) 
# a (partial) likelihood can be computed for the full model
# and compared to a likelihood for a model with a feature left out. 


# Conventional likelihood-ratio test for all features (except Intercept)
print("      INFERENCE BASED ON SHAPLEY VALUES   ")
cat("\n")
for(asymmetric in c(TRUE,FALSE)){
if(asymmetric){
  dat <- data.frame(aShapbf)
  print("       Test results for ASYMMETRIC setting     ")
  cat("\n")
} else {
  dat <- data.frame(Shapbf)
  print("       Test results for SYMMETRIC setting      ")  
  cat("\n")
}
nc <- ncol(dat)
cn <- colnames(dat)
llr <- c(NA)
ps <- c(NA)

#intercept is excluded for inference
for(k in 2:nc){
  #print(cn[k])
  cph <- coxph(Y_explain ~ ., data = dat)
  datk <- dat[,-k]
  cph0 <- coxph(Y_explain ~ ., data = datk)
  aov <- anova(cph,cph0)
#  print(aov)
  newp <- aov$`Pr(>|Chi|)`[2]
  ll <- aov$loglik[1] - aov$loglik[2]
  llr <- c(llr,ll)
  ps <- c(ps,newp)
}

names(ps) <- colnames(dat)
print("log-lik ratios and p-values lik ratio test")
print(rbind(llr,pval= ps)) #log-lik and p-values lik-ratio test
cat("\n")

#LRT for all three gene features together
print("genes+CMS+predStage")
cph <- coxph(Y_explain ~ ., data = dat)
datk <- dat[,-(2:4)]
cph0 <- coxph(Y_explain ~ ., data = datk)
aov <- anova(cph,cph0)
#print(aov)
newp <- aov$`Pr(>|Chi|)`[2]
ll <- aov$loglik[1] - aov$loglik[2]
print("log-lik ratio and p-values lik ratio test: all gene features")
print(signif(c(ll,newp),3)) #log-lik ratio + p-values gene features
cat("\n")

#LRT for two gene summaries together
print("CMS+predStage")
cph <- coxph(Y_explain ~ ., data = dat)
datk <- dat[,-(3:4)]
cph0 <- coxph(Y_explain ~ ., data = datk)
aov <- anova(cph,cph0) 
#print(aov)
newp <- aov$`Pr(>|Chi|)`[2]
ll <- aov$loglik[1] - aov$loglik[2]
print("log-lik ratio and p-values lik ratio test: two gene summaries")
print(signif(c(ll,newp),3)) #log-lik ratio + p-values two  gene summary features
cat("\n")

#### Alternative test: conditional independence test
#conditional independence test for all features except intercept
nc <- ncol(dat)
cn <- colnames(dat)
npps <- c(NA)
for(k in 2:nc){
# print(cn[k])
#Matching: creates 50 blocks on the basis of the 
#sum of Shapley values of other variables
  offsets <- rowSums(dat[,-c(1,k)])
  qs <- quantile(offsets, p=seq(0,1,by=0.02))
  whblock <- function(x) which(qs>x)[1]
  blocks <- sapply(offsets, whblock)
  dattest <- data.frame(Y_explain, vartest = dat[,k], blocks=factor(blocks))
#Conditional independence test based on permutation. 
#All samples within one block are permuted 
  newp <- pvalue(independence_test(Y_explain ~ vartest | blocks, data = dattest))
  npps <- c(npps,newp)
}
names(npps) <- colnames(dat)
print("p-values nonparametric cond indep test")
print(signif(npps,4)) #p-values nonparametric cond indep test
cat("\n")

# All three gene groups together
print("genes+CMS+predStage")
offsets <- rowSums(dat[,-c(1,2:4)])
qs <- quantile(offsets, p=seq(0,1,by=0.02))
whblock <- function(x) which(qs>x)[1]
blocks <- sapply(offsets, whblock)
dattest <- data.frame(Y_explain, vartest = rowSums(dat[,2:4]), blocks=factor(blocks))
newp <- pvalue(independence_test(Y_explain ~ vartest | blocks, data = dattest))

print("p-values nonparametric cond indep test: all gene features")
print(signif(newp,3)) #p-value nonparametric cond indep test: all gene features
cat("\n")

# Two gene summaries together
print("CMS+predStage")
offsets <- rowSums(dat[,-c(1,3:4)])
qs <- quantile(offsets, p=seq(0,1,by=0.02))
whblock <- function(x) which(qs>x)[1]
blocks <- sapply(offsets, whblock)
dattest <- data.frame(Y_explain, vartest = rowSums(dat[,3:4]), blocks=factor(blocks))
newp <- pvalue(independence_test(Y_explain ~ vartest | blocks, data = dattest))

print("p-values nonparametric cond indep test: two gene summaries")
print(signif(newp,3)) #p-value nonparametric cond indep test: two gene summaries
cat("\n")
}

############################################################################
#
#                         Mediation effect
#
############################################################################

# Plot Symmetric and Asymmetric Shapley values for Genes by true disease state
print("     MEDIATION EFFECT     ")
print("Plotting Symmetric and Asymmetric Shapley values for Genes by true disease state")
cat("\n")
par(mfrow=c(1,2),mar=c(3,2,2,2))

#Asymmetric
gall <- aShapbf$genes + aShapbf$CMS + aShapbf$predStage
shaps <- data.frame(G_Asym = gall, D = Clin_explain$Stage)
boxplot(G_Asym ~ D, data = shaps, xlab="", main="Asymmetric Shapley", ylim =c(-0.6,1.2),
        names=c("I","II","III","IV"), col=1:4)
Dint <- as.numeric(Clin_explain$Stage) 
points(Dint,gall,col="darkgrey",pch=4,cex=1)
# Test difference between gene Shapley values 
print("Test difference between Asymmetric gene Shapley values:")
print(kruskal_test(gall ~ Clin_explain$Stage)) #Asymmetric

#Symmetric
gall <- Shapbf$genes + Shapbf$CMS + Shapbf$predStage
shaps <- data.frame(G_Sym = gall, D = Clin_explain$Stage)
boxplot(G_Sym ~ D, data = shaps, xlab="", main="Symmetric Shapley", ylim =c(-0.6,1.2),
        names=c("I","II","III","IV"), col=1:4)
Dint <- as.numeric(Clin_explain$Stage)
points(Dint,gall,col="darkgrey",pch=4,cex=1)
# Test difference between gene Shapley values 
print("Test difference between Symmetric gene Shapley values:")
print(kruskal_test(gall ~ Clin_explain$Stage)) #Symmetric






