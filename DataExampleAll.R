###############################################################
#
#                        DATA EXAMPLE
#
###############################################################

# This script produces the results of the colorectal cancer
# data example in the manuscript. NOTE: as results from tree ensemble learners
# and cross-validation (for tuning parameters) are not deterministic,
# results may deviate somewhat from those in the manuscript.
# Qualitatively, however, they should be very similar.

# For detailed coding info, we refer to the annotation in the
# separate DataExampleX.R scripts, X=1,2,3.

# Set working directory
#setwd("C:/Synchr/Rscripts/ShapleyDependency/GeneShap/")


# libraries
# for model fitting
library(randomForestSRC) #for creating low-dim disease state specific gene summary
library(blockForest)
library(fusedTree)
library(rpart) #used for fusedTree
library(glmnet) #fits ridge regression

# for performance evaluation on survival data (C-index and time-dep AUC)
library(survC1)
library(survivalROC)

# For Shapley values
# Install specific version of the shapr package
library(remotes)
#remotes::install_github("NorskRegnesentral/shapr", ref = "hacks_for_Mark")
library(shapr)

library(data.table)
library(future)
library(progressr)

# for (non-parametric) inference
library(coin)


#####
#Part 1: Initial colorectal cancer data processing and
#model fitting: blockForest, fusedTree, ridge0
system.time(source("DataExample1.R")) #~800sec

#Part 2: Computation of Asymmetric and Symmetric Shapley values for blockForest model
system.time(source("DataExample2.R")) #~1100sec

#Part 3: SAGE, inference and results (Tables, Figure)
system.time(source("DataExample3.R")) #~5sec