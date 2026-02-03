#############################################################################
#
#                             IMPORTANCE SAMPLING     
#
#############################################################################

# Aim: checking whether importance sampling algorithm for asymmetric Shapley 
# weights works accurately. Checks are explained in the manuscript

# Run entire script by entering (will take less than a minute)
#source("ImpSampling_check.R")

# set working directory
setwd("C:/Synchr/Rscripts/ShapleyDependency/GeneShap/")

# new GeneShap functions, including those for Imp Sampling
source("GeneShap_source.R")

# Check 1: verifies correctness of the IS approximation of the asymmetric 
# Shapley weights for a setting in which exact weights are time-consuming:
# one G, one D, 10 confounders C1, ..., C10. 
print("       SETTING 1: many confouders        ")
nconf <- 10 #number of confounders to add
Conf <- paste("C",1:nconf,sep="")

print("Compute approximate weights based on 5000 importance samples")
cat("\n")
system.time(Qapprox <- QIS(Genes = c("G"), Status = c("D"),Confounders = Conf, nsamp = 5000, 
               list("Genes","Status"),Ordering_within=NULL,verbose=FALSE))

#extract approximate weights, and its positive and negative contributions
Qapp <- Qapprox[[1]]
Qpos <- Qapprox[[2]]
Qneg <- Qapprox[[3]]

print("Check: for all features, means of positive and negative contributions should be close to 1 (for intercept only the positive contribution)")
print(data.frame(rbind(pos=round(apply(Qpos,1,mean),2),
                 neg=round(apply(Qneg,1,mean),2))))
cat("\n")
# Check 2: compare Qapp with real Q for large number of samples; 
# are the coalitions sampled with the correct proportions? For a small setting
# with only one G, one D and two confounders C1, C2.
print("     SETTING 2: few confounders, comparison with exact weights  ")
nconf <- 2
Conf <- paste("C",1:nconf,sep="")

print("Compute approximate weights based on 5000 importance samples")
Qapprox <- QIS(Genes = c("G"), Status = c("D"),Confounders = Conf, nsamp = 5000, 
               list("Genes", "Status"), Ordering_within=NULL,verbose=FALSE)

Qapp <- Qapprox[[1]]
Qpos <- Qapprox[[2]]
Qneg <- Qapprox[[3]]

print("Compute exact weights for all subsets")
Qreal <- R_D_Matrix2(Genes = c("G"), Status = c("D"),Confounders = Conf,
                     list("Genes", "Status"),Ordering_within=NULL,verbose=FALSE)
print("Exact weights")
print(Qreal)
cn <- colnames(Qreal)
cnapp <- colnames(Qapp)
nsamp <- ncol(Qapp)

Qappder <- sapply(cn, function(ss){
  #cn = "G"
  wh <- which(cnapp==ss)
  wei <- length(wh)/nsamp
  return(Qapp[,wh[1]]*wei)
})

print("Differences between exact and approximated weigths")
Qreal-Qappder 
