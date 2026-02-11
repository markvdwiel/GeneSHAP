########################################################################
## --- Data Example, part II: (Asymmetric) Shapley values for
## --- best performing model, blockForest
########################################################################

#set working directory
#setwd("C:/Synchr/Rscripts/ShapleyDependency/GeneShap/")

#new GeneShap functions
source("GeneShap_source.R")

# Install specific version of the shapr package
library(remotes)
#remotes::install_github("NorskRegnesentral/shapr", ref = "hacks_for_Mark")
library(shapr)

# load other libraries
library(survival)
library(data.table) #for shapr
library(future) #for shapr
library(progressr) #for shapr
library(blockForest)

#######

#load data objects (see DataExample1.R)
show(load("data.Rdata"))
load("data.Rdata")

#####################################################
#
#      Unsupervised Dimension reduction
#
#####################################################

#PCA for high-dimensional objetct (Omics). These are used for
#dependency modelling when computing conditional Shapley values

pca_dim <- 10
pca_genes <- prcomp(Omics, center = FALSE, scale. = FALSE, rank. = pca_dim)

# Transform the training and explanation data
x_train_genes_pca <- predict(pca_genes, newdata = Genes_train)
x_explain_genes_pca <- predict(pca_genes, newdata = Genes_explain)

# Combine the PCA-reduced genes with the other features
x_train_reduced <- cbind(x_train_genes_pca, Clin_train)
x_explain_reduced <- cbind(x_explain_genes_pca, Clin_explain)

#####################################################
#
#      Shapley values for block forest model
#
#####################################################
# blockForest model was best performing model in terms of prediction
# (as compared to fusedTree and ridge0). Therefore. we focus
# on explaining the blockForest model

print("      COMPUTING SHAPLEY VALUES FOR BLOCKFOREST MODEL    ")
cat("\n")

#loads blockForest model (bf)
load("blockForestFit.Rdata")

#model list need for Shapley value computation
model_list <- list(model = bf,
                   pca_model = pca_genes,
                   x_genes = Omics,
                   pca_dim = pca_dim,
                   get_original_genes = get_original_genes)
class(model_list) <- "custom"

# Define feature group structures
other_list <- as.list(names(x_explain_reduced)[(pca_dim +1):(pca_dim + 3)])
names(other_list) <- names(x_explain_reduced)[(pca_dim +1):(pca_dim + 3)]
group <- c(
  list(genes = c(names(x_explain_reduced)[1:pca_dim]),CMS="CMS", predStage= colnames(Clin_train)[-(1:5)]),
  other_list,
  list(Stage = "Stage")
)

print("Feature groups for Shapley:")
print(group)
cat("\n")

#Feature types
Genes <- c("genes", "CMS","predStage")
Status <- c("Stage")
Confounders <- c("gender", "age","site")

#Compute exact asymmetric Shapley weights, given an ordering Genes -> Status
R_D <- R_D_Matrix(Genes, Status, Confounders,
                  Ordering_between = list("Genes", "Status"),
                  verbose = FALSE)
print("Dimensions Asymmetric Shapley weight matrix")
print(dim(R_D))
cat("\n")

#Asymmetric Shapley weight matrix; features in rows, coalitions in columns
#print(R_D)

#preparing for Shapley value computation
shap_names <- names(group)
coalition_list <- parse_coalitions(shap_names, names(R_D))

#mean prediction for training set
X_tot_train_new <- data.frame(Clin_train,Genes_train)
pred_train_org <- predict_bf(bf, X_tot_train_new)

print("Computing Asymmetric Shapley values for blockForest; takes some time (3-5 min)")
cat("\n")
ShapAsym <- explain_manual(model = model_list,
                            x_explain = x_explain_reduced,
                            x_train = x_train_reduced,
                            predict_model = predict_with_reduced_bf,
                            approach = "ctree",
                            phi0 = mean(pred_train_org),
                            coalition_list = coalition_list,
                            R_D = as.matrix(R_D),
                            group = group
)

### Same but now symmetric Shapley values
#no ordering between feature groups
R_D_Sym <- R_D_Matrix(Genes, Status, Confounders, Ordering_between = NULL,
                  verbose = FALSE)


shap_names <- names(group)
coalition_list <- parse_coalitions(shap_names, names(R_D_Sym))

#Exact symmetric Shapley values; takes some time (10-15 min)
print("Computing Symmetric Shapley values for blockForest; takes some time (3-5 min)")

ShapSym <- explain_manual(model = model_list,
                           x_explain = x_explain_reduced,
                           x_train = x_train_reduced,
                           predict_model = predict_with_reduced_bf,
                           approach = "ctree",
                           phi0 = mean(pred_train_org),
                           coalition_list = coalition_list,
                           R_D = as.matrix(R_D_Sym),
                           group = group
)




# Quick assessment of global importance of features by sum of absolute values
# of (asymmetric) Shapley values (not SAGE!)
print("Sum of absolute values of Asymmetric Shapley values")
print(round(colMeans(abs(ShapAsym$shap)),3)) #Asymmetric
print("Sum of absolute values of Symmetric Shapley values")
print(round(colMeans(abs(ShapSym$shap)),3)) #Symmetric
cat("\n")

#save for further use
save(ShapAsym, ShapSym, R_D, R_D_Sym, file="Shapleys.Rdata")



