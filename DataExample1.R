########################################################################
## --- Data Example, part I: data processing and prediction models  ---
########################################################################
#
# NOTE: AS RESULTS FROM TREE ENSEMBLE LEARNERS AND CROSS-VALIDATION (FOR TUNING PARAMETERS) 
# ARE NOT DETERMINISTIC, RESULTS MAY DEVIATE SOMEWHAT FROM THOSE IN THE MANUSCRIPT. 
# QUALITATIVELY, HOWEVER, THEY SHOULD BE VERY SIMILAR. 

# set working directory
setwd("C:/Synchr/Rscripts/ShapleyDependency/GeneShap/")

library(randomForestSRC) #for creating low-dim disease state specific gene summary
library(blockForest)
library(fusedTree)
library(rpart) #used for fusedTree
library(glmnet) #fits ridge regression

# for performance evaluation on survival data (C-index and time-dep AUC)
library(survC1)
library(survivalROC)

print("      INITIAL DATA PROCESSING    ")
cat("\n")

########## Load en Define the Data ##########
#############################################
load("500_CRC.Rdata")
Clinical <- as.data.frame(Clinical)
Omics <- as.data.frame(scale(Omics)) # standardize omics (important for ridge fit and pca)
colnames(Omics) <- paste0("x", seq(1,ncol(Omics)))
Tot <- cbind.data.frame(Omics, Clinical)

### Training and test data ###
set.seed(48) #for reproducibility
ids <- sample(1:nrow(Clinical),size = 0.2*nrow(Clinical),replace = F)
Y_train <- Response[-ids,]; Y_explain <- Response[ids,]
Genes_train <- Omics[-ids,]; Genes_explain <- Omics[ids,]
Clin_train <- Clinical[-ids,]; Clin_explain <- Clinical[ids,]
X_Tot_train <- Tot[-ids,]; X_Tot_explain <- Tot[ids,]

############# Dim Red 0: regressing DS on genes #####
# Create a disease state specific summary of the genes. As this summary is used by the 
# prediction algorithms (see manuscript) it needs to be created before applying these algorithms

#rfsrc requires a data set containing both the response variable and the omics variables
trainlab <- Clin_train$Stage
traindat <- Genes_train
databoth <- data.frame(trainlab,traindat)
#set.seed(455)
rfres <- rfsrc(trainlab ~ ., var.used="all.trees",ntree=1000, data = databoth, importance="none")
# 
# predicted disease state (train)
predstage <- data.frame(rfres$predicted.oob[,-4])
colnames(predstage) <- paste("Pred", colnames(predstage), sep="")
# 
# predicted disease state (test)
predictrf <- data.frame((predict(rfres,newdata=Genes_explain)$predicted)[,-4])
colnames(predictrf) <- colnames(predstage)

#Extend low-dimensional variables by summary
Clin_train <- cbind(Clin_train,predstage)
Clin_explain <- data.frame(Clin_explain,predictrf)

#save processed data, as these are also need for Shapley value computations
save(Clin_train, Clin_explain, Omics, Genes_train, Genes_explain, Y_train, Y_explain, file="data.Rdata")

print("      TRAINING PREDICTION MODELS    ")
cat("\n")
print("NOTE: as results from tree ensemble learners and parameter tuning by cross-validation are not deterministic, results may deviate somewhat from those in the manuscript. Qualitatively, however, they should be very similar.")
cat("\n")

##### --- 1. blockForest --- #####
################################

#Fitting blockForest on training data
X_tot_train_new <- data.frame(Clin_train,Genes_train)

#Two covariate blocks: low-dimensional one and omics. First block is always selected
bloks <- list(clin=1:ncol(Clin_train),gen = ((ncol(Clin_train)+1):ncol(X_tot_train_new)))

#Fitting; takes 10-15 minutes!
print("Training blockForest; may take several minutes")
cat("\n")
bf <- blockfor(X_tot_train_new,Y_train,blocks=bloks,always.select.block = 1,num.trees=200)

# Evaluating predictive performance on test data
X_tot_explain_new <- data.frame(Clin_explain,Genes_explain)

Ypredbf <- predict(bf$forest, data=X_tot_explain_new)
Ypredbf_cumH <- rowSums(Ypredbf$chf)

Cindex <- Est.Cval(cbind(Y_explain[,1], Y_explain[,2], Ypredbf_cumH), tau = 8, nofit=T)$Dhat
AUC <- survivalROC.C(Stime = Y_explain[,1], status = Y_explain[,2], marker = Ypredbf_cumH, predict.time = 5)$AUC

print("Performance blockForest, C-index, tAUC:")
print(Cindex)
print(AUC)
cat("\n")

print("Fit blockForest saved")
cat("\n")
save(bf,AUC,Cindex,file="blockForestFit.Rdata")

##### --- 2. FusedTree --- #####
################################

# Fitting fusedTree on training data. First, decision tree on low dimensional covariates
print("Training fusedTree; may take some minutes")
cat("\n")
dat <- cbind.data.frame(Y_train, Clin_train)
set.seed(10)
rp <- rpart(Y_train ~ .,data = dat, control = rpart.control(xval = 5, minbucket = 30, cp = 0),
            model = T)
minerr <- which.min(rp$cptable[,"xerror"])
bestcp <- rp$cptable[5,"CP"]
Treefit <- prune(rp, cp = bestcp)

# Second: fused ridge regression on omics, using tree fit as offset
# a. tuning fused ridge penalties; takes 3-5 minutes
#set.seed(3683)
folds <- CVfoldsTree(Y = Y_train, Z = Clin_train, Tree = Treefit, model = "cox",
                     nrepeat = 3, kfold = 5)
optPenalties <- PenOpt(Y = Y_train, Z = Clin_train, X = as.matrix(Genes_train),
                       Tree = Treefit, model = "cox", LinVars = TRUE,
                       lambdaInit = 1000, alphaInit = 1000,
                       maxIter = 50)
Lam_tuned <- optPenalties[1]
alf_tuned <- optPenalties[2]

#b. Actual fitting using tuned penalties (and Tree Fit)
fT <- fusedTree(Y = Y_train, Z = Clin_train, X = as.matrix(Genes_train),
                   Tree = Treefit, model = "cox", LinVars = TRUE,
                   lambda = Lam_tuned, alpha = alf_tuned, maxIter = 100)

# Test sample performance
Ypred <- predict(fT, newX = as.matrix(Genes_explain), newZ = Clin_explain, newY = Y_explain)$LinPred[,1]

Cindex <- Est.Cval(cbind(Y_explain[,1], Y_explain[,2], Ypred), tau = 8, nofit=T)$Dhat
AUC <- survivalROC.C(Stime = Y_explain[,1], status = Y_explain[,2], marker = Ypred, predict.time = 5)$AUC


print("Performance fusedTree, C-index, tAUC:")
print(Cindex)
print(AUC)

save(fT, Cindex,AUC, file = "fusedTreeFit.Rdata")

##### --- 3. ridge0 --- #####
################################
# Fitting ridge0 on training data. ridge0 is a ridge regression model that 
# does not penalize the low-dimensional variables

# Fitting: tuning penalty parameter + fit
#ridge regression requires dummy coding for low-dimensional variables
Clin_train_num <- model.matrix(~., Clin_train)[,-1] 
Clin_explain_num <- model.matrix(~., Clin_explain)[,-1] 

p_omics <- ncol(Genes_train)
p_clin <- ncol(Clin_train_num)

print("Training ridge regression; may take some minutes")
Las <- cv.glmnet(as.matrix(cbind(Clin_train_num, Genes_train)),
                 Y_train, alpha = 0,  nfolds = 5,family = "cox",
                 standardize = TRUE, penalty.factor = c(rep(0, p_clin), rep(1, p_omics)))$lambda.min

ridge0 <- glmnet(as.matrix(cbind(Clin_train_num, Genes_train)),
                   Y_train, alpha = 0, lambda = Las, family = "cox",
                   standardize = TRUE, penalty.factor = c(rep(0,p_clin),rep(1,p_omics)))

# Evaluation on test set
LP <- (as.matrix(cbind(Clin_explain_num, Genes_explain)) %*% ridge0$beta)[,1]
Cindex <- Est.Cval(cbind(Y_explain[,1],Y_explain[,2],LP), tau = 8, nofit=T)$Dhat
AUC <- survivalROC.C(Stime = Y_explain[,1], status = Y_explain[,2], marker = LP, predict.time = 5)$AUC


print("Performance ridge0, C-index, tAUC:")
print(Cindex)
print(AUC)

save(ridge0, Cindex, AUC, file="ridge0Fit.Rdata")

##########################