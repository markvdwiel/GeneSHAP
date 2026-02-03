#################################################################
#
#                   Source functions GeneShap 
#
#################################################################

# Source functions for Geneshap. Containing
# 1. Modification of the explain function in shapr 
# 2. Prediction using reduced model 
# 3. Exact Asymmetric Shapley weights 
# 4. Approximation by Importance Sampling 

############################################################################
#
#                     explain_manual (modified from shapr) 
#
############################################################################

# Manual implementation of the explain function in shapr, using a range of internal functions from the shapr package.
# This manual implementation allows us to pass
# coalition_list: A list of the coalitions to evaluate
# R_D: The weight matrix, which we multiply by the v(S) vector/matrix to compute the Shapley values 
# (i.e. R_D[i+1,j] is the weight for coalition j for feature i). May be computed with function R_D_Matrix
explain_manual <- function(
    model,
    x_explain,
    x_train,
    approach,
    phi0,
    group = NULL,
    n_MC_samples = 1e3,
    predict_model = NULL,
    seed = 1,
    coalition_list,
    R_D, # R_D in (7) of Aas et al. (2021)
    coalition_approach_dt = NULL, # EXPERIMENTAL: data.table with columns "coalitions_str" and "approach_new"
    ...
){
  
  internal <- shapr:::setup(
    x_train = x_train,
    x_explain = x_explain,
    approach = approach,
    phi0 = phi0,
    max_n_coalitions = NULL,
    group = group,
    n_MC_samples = n_MC_samples,
    seed = seed,
    feature_specs = shapr:::get_feature_specs(NULL, model),
    verbose = NULL,
    iterative = FALSE,
    ...
  )
  
  predict_model <- shapr:::get_predict_model(predict_model = predict_model, model = model)
  
  # Overwrite internals
  
  internal <- setup_approach(internal, model = model, predict_model = predict_model)
  
  # Manual analogue to shapley_setup
  iter <- 1
  
  dt_coalitions <- data.table(coalitions = coalition_list)
  
  X <- shapr:::create_coalition_table(m = internal$parameters$n_shapley_values,
                                      exact = TRUE,
                                      approach0 = internal$parameters$approach,
                                      coal_feature_list = internal$objects$coal_feature_list,
                                      dt_valid_causal_coalitions = dt_coalitions)
  
  
  # EXPERIMENTAL: Overwriting the set approach using the coalition_approach_dt
  # which specifies the approach to use for each specific coalition
  if(!is.null(coalition_approach_dt)){
    X <- merge(X,coalition_approach_dt,by="coalitions_str",all.x=TRUE,all.y=FALSE)
    
    X[!is.na(approach_new),approach:=approach_new]
    
    X[,approach_new:=NULL]
    
    setorder(X,id_coalition)
    setkey(X,coalition_size)
  }
  
  ## Get feature matrix ---------
  S <- shapr:::coalition_matrix_cpp(
    coalitions = X[["features"]],
    m = internal$parameters$n_features
  )
  
  
  internal$iter_list[[iter]]$X <- X
  internal$iter_list[[iter]]$W <- R_D
  internal$iter_list[[iter]]$S <- S
  internal$iter_list[[iter]]$S_batch <- shapr:::create_S_batch(internal)
  
  # Compute the vS
  vS_list <- compute_vS(internal, model, predict_model)
  
  processed_vS_list <- shapr:::postprocess_vS_list(
    vS_list = vS_list,
    internal = internal
  )
  
  # Compute the Shapley values
  nu <- processed_vS_list$dt_vS
  dt_shapley_est <- shapr:::compute_shapley(internal, nu )
  shap <- dt_shapley_est[]
  #rownames(nu) <- rownames(R_D)
  colnames(shap) <- rownames(R_D)
  #return(shap)
  return(c(shap=list(shap),nu=list(nu)))
}

parse_coalitions <- function(shap_names, coalition_str) {
  result <- vector("list", length(coalition_str))
  
  for (i in seq_along(coalition_str)) {
    if (coalition_str[i] == "Ø") {
      result[[i]] <- numeric()
    } else {
      components <- strsplit(coalition_str[i], "_")[[1]]
      result[[i]] <- match(components, shap_names)
    }
  }
  
  return(result)
}

############################################################################
#
#                     Predict with reduced model 
#
############################################################################
## This module provides functions that map a low-dimensional summary (like PCA), 
## used to represent dependencies when computing SHAP,
## back to the high-dimensional vector, wich can then be used for 
## predicting the outcome when computing conditinional Shapley values by SHAP
##
## Function to map from pca to original genes ##
## nearest neighbour
get_original_genes <- function(pca_values, pca_model, x_genes) {
  # Find the nearest neighbor in the PCA space
  nn_index <- apply(pca_values, 1, function(x) {
    which.min(colSums((t(pca_model$x) - x)^2))
  })
  
  # Return the original gene values corresponding to the nearest neighbors
  return(x_genes[nn_index, ])
}

#for blockForest
predict_with_reduced_bf <- function(model_list, newdata){
  model <- model_list$model
  pca_model <- model_list$pca_model
  x_genes <- model_list$x_genes
  pca_dim <- model_list$pca_dim
  get_original_genes <- model_list$get_original_genes
  
  # Extract the PCA values for the genes
  pca_values <- as.matrix(newdata[, 1:pca_dim]) # newdata should therefore start with PCA scores
  
  # Get the original gene values based on the nearest neighbor in the PCA space
  original_genes <- get_original_genes(pca_values, pca_model, x_genes)
  
  x_combined <- cbind.data.frame(newdata[, (pca_dim + 1):ncol(newdata)], original_genes)
  
  # Predict using the model
  Ypred0 <- predict(model$forest, data = x_combined)
  Ypred <- log(rowSums(Ypred0$chf)) #log cum Haz scale
  return(Ypred) # returns 
}

#for blockForest
predict_bf <- function(model, newdata){
  Ypred0 <- predict(model$forest, data = newdata)
  Ypred <- log(rowSums(Ypred0$chf)) #log cum Haz scale
  # Predict using the model
  return(Ypred) # returns 
}

#for fusedTree
predict_with_reduced_fT <- function(model_list, newdata){
  model <- model_list$model
  Tree <- model$Tree #FusedTree
  Ests <- model$Effects #FusedTree
  pca_model <- model_list$pca_model
  x_genes <- model_list$x_genes
  pca_dim <- model_list$pca_dim
  get_original_genes <- model_list$get_original_genes
  
  # Extract the PCA values for the genes
  pca_values <- as.matrix(newdata[, 1:pca_dim]) # newdata should therefore start with PCA scores
  
  # Get the original gene values based on the nearest neighbor in the PCA space
  original_genes <- get_original_genes(pca_values, pca_model, x_genes)
  
  Clin <- newdata[, (pca_dim + 1):ncol(newdata)]
  
  Dat_test <- fusedTree::Dat_Tree(Tree = Tree, 
                                  X = as.matrix(original_genes), 
                                  Z = Clin, LinVars = TRUE)
  Xtest    <- Dat_test$Omics
  Ztest    <- Dat_test$Clinical
  
  # ---- Align predictors with model coefficients ----
  keep_x   <- colnames(Xtest) %in% names(Ests)
  Xtest    <- Xtest[, keep_x, drop = FALSE]
  
  keep_ids <- names(Ests) %in% c(colnames(Ztest), colnames(Xtest))
  Ests     <- Ests[keep_ids]
  
  # ---- Compute linear predictor ----
  Ypred <- as.numeric(cbind(Ztest, Xtest) %*% Ests)
  
  return(Ypred) # returns linear predictor
}

#for ridge
predict_with_reduced_ridge <- function(model_list, newdata){
  #newdata <- x_train_reduced
  model <- model_list$model
  pca_model <- model_list$pca_model
  x_genes <- model_list$x_genes
  pca_dim <- model_list$pca_dim
  get_original_genes <- model_list$get_original_genes
  
  # Extract the PCA values for the genes
  pca_values <- as.matrix(newdata[, 1:pca_dim]) # newdata should therefore start with PCA scores
  
  # Get the original gene values based on the nearest neighbor in the PCA space
  original_genes <- as.matrix(get_original_genes(pca_values, pca_model, x_genes))
  
  Clin <- as.matrix(newdata[, (pca_dim + 1):ncol(newdata)])
  
  LP <- (cbind(Clin, original_genes) %*% model$beta)[,1]
  return(LP) # returns linear predictor
}

############################################################################
#
#                     Exact Asymmetric Shapley Matrix
#
############################################################################

R_D_Matrix <- function(Genes,
                       Status,
                       Confounders,
                       Ordering_between = NULL,
                       Ordering_within = NULL,
                       verbose = TRUE) {
  
  
  ## Control Statements ##
  valid_groups <- c("Genes", "Status", "Confounders")
  
  if (!is.null(Ordering_between)) {
    if (!is.list(Ordering_between)) stop("`Ordering_between` must be a list or NULL.")
    if (length(Ordering_between) < 2) {
      stop("Ordering should contain at least one level of hierarchy,
           i.e. the list should have at least two elements")
    }
    all_elements <- unlist(Ordering_between)
    if (!all(all_elements %in% valid_groups)) {
      stop("`Ordering_between` can only contain 'Genes', 'Status', and 'Confounders'.")
    }
    if (any(duplicated(all_elements))) {
      dup_elems <- unique(all_elements[duplicated(all_elements)])
      stop(paste0("`ordering` contains duplicated group(s): ",
                  paste(dup_elems, collapse = ", "), "."))
    }
  }
  
  ## define the allowed subsets ##
  all_values <- c(Genes, Status, Confounders)
  group_map <- list(Genes = Genes, Status = Status, Confounders = Confounders)
  
  all_subsets <- function(x) {
    unlist(lapply(0:length(x), function(k) combn(x, k, simplify = FALSE)), recursive = FALSE)
  }
  
  subsets_all <- all_subsets(all_values)
  
  valid_subset <- function(subset) {
    if (is.null(Ordering_between)) return(TRUE)
    for (k in seq_along(Ordering_between)) {
      current_groups <- unlist(group_map[Ordering_between[[k]]])
      if (any(current_groups %in% subset)) {
        required_prior <- unlist(group_map[unlist(Ordering_between[seq_len(k - 1)])])
        if (length(required_prior) > 0 && !all(required_prior %in% subset))
          return(FALSE)
      }
    }
    TRUE
  }
  
  if (is.null(Ordering_between)) {
    Subsets_Constrained <- subsets_all
  } else {
    Ordering_between <- lapply(Ordering_between, as.character)
    Subsets_Constrained <- Filter(valid_subset, subsets_all)
  }
  
  subset_names <- sapply(Subsets_Constrained, function(s) {
    if (length(s) == 0) return("Ø")
    paste(s, collapse = "_")
  })
  
  ## define empty R_D matrix ##
  Q_df <- as.data.frame(matrix(NA, nrow = length(all_values), ncol = length(Subsets_Constrained),
                               dimnames = list(all_values, subset_names)))
  
  convert_map <- c(Genes = "M1", Status = "M2", Confounders = "M3")
  
  
  ## Compute Q-matrix ##
  for (var in all_values) {
    
    if (verbose) {
      message("\n===============================")
      message("Variable of interest: ", var)
      message("===============================")
    }
    
    for (j in seq_along(Subsets_Constrained)) {
      subset <- Subsets_Constrained[[j]]
      subset_name <- if (length(subset) == 0) "Ø" else paste(subset, collapse = "_")
      
      subset_minus_var <- setdiff(subset, var)
      subset_plus_var  <- union(subset, var)
      remaining_minus_var <- setdiff(all_values, union(subset, var))
      
      # check validity of subset_minus_var and subset_plus_var (weight -> 0 if invalid)
      valid_minus <- valid_subset(subset_minus_var)
      valid_plus  <- valid_subset(subset_plus_var)
      
      if (!valid_minus || !valid_plus) {
        if (verbose) {
          message("Subset ", subset_name, " → invalid (subset violates ordering), Q = 0")
        }
        Q_df[var, j] <- 0
        next
      }
      
      # Build M1/M2/M3 sets
      M1_in <- intersect(subset_minus_var, Genes)
      M2_in <- intersect(subset_minus_var, Status)
      M3_in <- intersect(subset_minus_var, Confounders)
      
      M1_out <- intersect(remaining_minus_var, Genes)
      M2_out <- intersect(remaining_minus_var, Status)
      M3_out <- intersect(remaining_minus_var, Confounders)
      
      # Ordering adjustment
      ordering_in <- NULL
      ordering_out <- NULL
      if (!is.null(Ordering_between)) {
        ordering_in <- Filter(function(grp) any(unlist(group_map[grp]) %in% subset_minus_var), Ordering_between)
        ordering_out <- Filter(function(grp) any(unlist(group_map[grp]) %in% remaining_minus_var), Ordering_between)
        ordering_in <- lapply(ordering_in, function(x) convert_map[x])
        ordering_out <- lapply(ordering_out, function(x) convert_map[x])
      }
      
      # Compute Tot_orders
      val_in <- .Tot_orders(M1_in, M2_in, M3_in, ordering = ordering_in)
      val_out <- .Tot_orders(M1_out, M2_out, M3_out, ordering = ordering_out)
      
      # Print debug info
      if (verbose) {
        message("Subset: {", subset_name, "}")
        message("  Inside subset  (subset_minus_var): {", paste(subset_minus_var, collapse = ", "), "}")
        message("  Remaining subset (remaining_minus_var): {", paste(remaining_minus_var, collapse = ", "), "}")
        message("  .Tot_orders in subset  = ", val_in)
        message("  .Tot_orders out subset = ", val_out)
        message("  => Product = ", val_in * val_out)
      }
      
      Q_df[var, j] <- val_in * val_out
      if (!(var %in% subset)) {
        Q_df[var, j] <- -1 * Q_df[var, j]
      }
    }
  }
  ordering_Convert <- lapply(Ordering_between, function(x) convert_map[x])
  TotOrders <- .Tot_orders(M1 = Genes, M2 = Status, M3 = Confounders, ordering = ordering_Convert)
  Q_df <- Q_df/TotOrders # normalize to weights
  Q_df <- rbind(intercept = c(1, rep(0, ncol(Q_df) - 1)), Q_df) # include intercept term to include average across all predictions
  
  return(Q_df)
}


#############################################################################
#
#                             IMPORTANCE SAMPLING     
#
#############################################################################

#Compute approximate Shapley weight matrix by Importance sampling;
QIS <- function(Genes, Status, Confounders, nsamp = 5000, 
                Ordering_between = NULL,Ordering_within=NULL,verbose=FALSE){
  res <- list()
  for(k in 1:nsamp) res <- c(res,list(sampleSet(Genes, Status, Confounders, Ordering_between,Ordering_within)))
  Slist <- c();iswei <- c()
  for(k in 1:nsamp) {
    iswei <- c(iswei,res[[k]][[1]])
    Slist <- c(Slist,list(res[[k]][[2]]))
  }
  Qsamp <- R_D_Matrix2(Genes,Status,Confounders,Ordering_between,Ordering_within, mysubsets = Slist,verbose=FALSE) 
  Qpos <- Qsamp * (Qsamp>=0)
  Qneg <- -Qsamp * (Qsamp<=0)
  ISpos <- t((1/iswei)*t(Qpos))
  ISneg <- t((1/iswei)*t(Qneg))
  Q_approx <-  t((1/iswei)*t(Qpos-Qneg))
  return(list(Q_approx,ISpos,ISneg)) #returns approximation and positive and negative contribution as well
}


#Samples the subsets according to the importance samples
sampleSet <- function(Genes, Status, Confounders, 
                      Ordering_between = NULL, 
                      Ordering_within = NULL){
  #Genes <- c("G1", "G2"); Status <- c("D");Confounders <- c("C1","C2");Ordering_between <- list("Genes", "Status","Confounders") 
  allfeat <- c(Genes,Status,Confounders)
  group_map <- list(Genes = Genes, Status = Status, Confounders = Confounders)
  nG <- length(Genes)
  nS <- length(Status)
  
  nfeat <- length(allfeat) 
  convert_map <- c(Genes = "M1", Status = "M2", Confounders = "M3")
  if(length(Ordering_between)==2 & length(Ordering_between[[2]])==1){  #G-> D;Ordering_between = list("Genes", "Status")
    samall0 <- sample(allfeat) 
    whG <- match(Genes,samall0)
    whS <- match(Status,samall0)
    whGS <- sort(c(whG,whS))
    samall <- samall0
    samall[whGS[1:nG]] <- samall0[whG]
    samall[whGS[(nG+1):(nG+nS)]] <- samall0[whS]
    #samall
  }
  
  if(length(Ordering_between)==2 & length(Ordering_between[[2]])==2){ #G->D, G->C; Ordering_between = list("Genes", c("Status", "Confounders"))
    sam1 <- sample(Genes)
    DC <- c(Status,Confounders)
    sam2 <- sample(DC)
    samall <- c(sam1,sam2)  
  }
  
  if(length(Ordering_between)==3){ #G->D, D->C; Ordering_between = list("Genes", "Status", "Confounders")
    samall <- c(sample(Genes),sample(Status),sample(Confounders))  
  }
  
  if(is.null(Ordering_between)){ #no ordering
    samall <- sample(allfeat)   
  }
  
  orderingall <- lapply(Ordering_between, function(x) convert_map[x])
  valall <- .Tot_orders(Genes,Status,Confounders, ordering = orderingall )
  
  samplecut <- sample(0:nfeat,1)
  samplecut
  if(samplecut==0) actset <- c() else actset <- samall[1:samplecut]
  actset <- sort(actset,decreasing=TRUE)
  inactset <- setdiff(allfeat, actset)
  
  M1_in <- intersect(actset, Genes)
  M2_in <- intersect(actset, Status)
  M3_in <- intersect(actset, Confounders)
  
  M1_out <- intersect(inactset, Genes)
  M2_out <- intersect(inactset, Status)
  M3_out <- intersect(inactset, Confounders)
  
  # Ordering adjustment
  ordering_in <- NULL
  ordering_out <- NULL
  if (!is.null(Ordering_between)) {
    ordering_in <- Filter(function(grp) any(unlist(group_map[grp]) %in% actset), Ordering_between)
    ordering_out <- Filter(function(grp) any(unlist(group_map[grp]) %in% inactset), Ordering_between)
    ordering_in <- lapply(ordering_in, function(x) convert_map[x])
    ordering_out <- lapply(ordering_out, function(x) convert_map[x])
  }
  
  # Compute Tot_orders
  val_in <- .Tot_orders(M1_in, M2_in, M3_in, ordering = ordering_in)
  val_out <- .Tot_orders(M1_out, M2_out, M3_out, ordering = ordering_out)
  wei <- (val_in*val_out/valall)/(nfeat+1)
  return(list(w=wei, S= actset))
}


#Computes exact Shapley weights for sampled subsets
R_D_Matrix2 <- function(Genes, Status, 
                        Confounders, 
                        Ordering_between = NULL, 
                        Ordering_within = NULL,
                        mysubsets = NULL, verbose = TRUE) {
  valid_groups <- c("Genes", "Status", "Confounders")
  
  if (!is.null(Ordering_between)) {
    if (!is.list(Ordering_between)) stop("`Ordering_between` must be a list or NULL.")
    if (length(Ordering_between) < 2) {
      stop("Ordering should contain at least one level of hierarchy,
           i.e. the list should have at least two elements")
    }
    all_elements <- unlist(Ordering_between)
    if (!all(all_elements %in% valid_groups)) {
      stop("`Ordering_between` can only contain 'Genes', 'Status', and 'Confounders'.")
    }
    if (any(duplicated(all_elements))) {
      dup_elems <- unique(all_elements[duplicated(all_elements)])
      stop(paste0("`ordering` contains duplicated group(s): ",
                  paste(dup_elems, collapse = ", "), "."))
    }
  }
  
  ## define the allowed subsets ##
  all_values <- c(Genes, Status, Confounders)
  group_map <- list(Genes = Genes, Status = Status, Confounders = Confounders)
  
  all_subsets <- function(x) {
    unlist(lapply(0:length(x), function(k) combn(x, k, simplify = FALSE)), recursive = FALSE)
  }
  valid_subset <- function(subset) {
    if (is.null(Ordering_between)) return(TRUE)
    for (k in seq_along(Ordering_between)) {
      current_groups <- unlist(group_map[Ordering_between[[k]]])
      if (any(current_groups %in% subset)) {
        required_prior <- unlist(group_map[unlist(Ordering_between[seq_len(k - 1)])])
        if (length(required_prior) > 0 && !all(required_prior %in% subset))
          return(FALSE)
      }
    }
    TRUE
  }
  
  if(is.null(mysubsets)){
    subsets_all <- all_subsets(all_values)
    
    
    
    if (is.null(Ordering_between)) {
      Subsets_Constrained <- subsets_all
    } else {
      Ordering_between <- lapply(Ordering_between, as.character)
      Subsets_Constrained <- Filter(valid_subset, subsets_all)
    }} else { #mysubsets != NULL
      if (!is.null(Ordering_between)) Ordering_between <- lapply(Ordering_between, as.character)
      Subsets_Constrained <- mysubsets
    }
  
  subset_names <- sapply(Subsets_Constrained, function(s) {
    if (length(s) == 0) return("Ø")
    s <- sort(s,decreasing=TRUE) #added to make column names consistent
    paste(s, collapse = "_")
  })
  names(Subsets_Constrained) <- subset_names
  
  ## define empty R_D matrix ##
  Q_df <- as.data.frame(matrix(NA, nrow = length(all_values), ncol = length(Subsets_Constrained),
                               dimnames = list(all_values, subset_names)))
  
  convert_map <- c(Genes = "M1", Status = "M2", Confounders = "M3")
  
  
  jempty <- c();
  for (j in seq_along(Subsets_Constrained)) {
    subset <- Subsets_Constrained[[j]]
    if (length(subset) == 0)  jempty <- c(jempty,j) #to keep track where we should add the intercept 
  }
  
  ## Compute Q-matrix ##
  for (var in all_values) {
    
    if (verbose) {
      message("\n===============================")
      message("Variable of interest: ", var)
      message("===============================")
    }
    for (j in seq_along(Subsets_Constrained)) {
      #j <-1
      subset <- Subsets_Constrained[[j]]
      subset_name <- subset_names[j]
      
      subset_minus_var <- setdiff(subset, var)
      subset_plus_var  <- union(subset, var)
      remaining_minus_var <- setdiff(all_values, union(subset, var))
      
      # check validity of subset_minus_var and subset_plus_var (weight -> 0 if invalid)
      valid_minus <- valid_subset(subset_minus_var)
      valid_plus  <- valid_subset(subset_plus_var)
      
      if (!valid_minus || !valid_plus) {
        if (verbose) {
          message("Subset ", subset_name, " → invalid (subset violates ordering), Q = 0")
        }
        Q_df[var, j] <- 0
        next
      }
      
      # Build M1/M2/M3 sets
      M1_in <- intersect(subset_minus_var, Genes)
      M2_in <- intersect(subset_minus_var, Status)
      M3_in <- intersect(subset_minus_var, Confounders)
      
      M1_out <- intersect(remaining_minus_var, Genes)
      M2_out <- intersect(remaining_minus_var, Status)
      M3_out <- intersect(remaining_minus_var, Confounders)
      
      # Ordering adjustment
      ordering_in <- NULL
      ordering_out <- NULL
      if (!is.null(Ordering_between)) {
        ordering_in <- Filter(function(grp) any(unlist(group_map[grp]) %in% subset_minus_var), Ordering_between)
        ordering_out <- Filter(function(grp) any(unlist(group_map[grp]) %in% remaining_minus_var), Ordering_between)
        ordering_in <- lapply(ordering_in, function(x) convert_map[x])
        ordering_out <- lapply(ordering_out, function(x) convert_map[x])
      }
      
      # Compute Tot_orders
      val_in <- .Tot_orders(M1_in, M2_in, M3_in, ordering = ordering_in)
      val_out <- .Tot_orders(M1_out, M2_out, M3_out, ordering = ordering_out)
      
      # Print debug info
      if (verbose) {
        message("Subset: {", subset_name, "}")
        message("  Inside subset  (subset_minus_var): {", paste(subset_minus_var, collapse = ", "), "}")
        message("  Remaining subset (remaining_minus_var): {", paste(remaining_minus_var, collapse = ", "), "}")
        message("  .Tot_orders in subset  = ", val_in)
        message("  .Tot_orders out subset = ", val_out)
        message("  => Product = ", val_in * val_out)
      }
      
      Q_df[var, j] <- val_in * val_out
      if (!(var %in% subset)) {
        Q_df[var, j] <- -1 * Q_df[var, j]
      }
    }
  }
  ordering_Convert <- lapply(Ordering_between, function(x) convert_map[x])
  TotOrders <- .Tot_orders(M1 = Genes, M2 = Status, M3 = Confounders, ordering = ordering_Convert)
  Q_df <- Q_df/TotOrders # normalize to weights
  interc <- rep(0, ncol(Q_df))
  interc[jempty] <- 1 #1 for empty sets
  Q_df <- rbind(intercept = interc, Q_df) # include intercept term to include average across all predictions
  return(Q_df)
}


######################## -- Helper Functions -- ########################

.Tot_orders <- function(M1, M2, M3, ordering = NULL) {
  
  valid_groups <- c("M1", "M2", "M3")
  group_sizes <- list(M1 = length(M1), M2 = length(M2), M3 = length(M3))
  p_tot <- sum(unlist(group_sizes))
  
  if (p_tot == 0) return(1)  # if no variables, only 1 ordering possible
  
  if (!is.null(ordering)) {
    if (!is.list(ordering)) stop("`ordering` must be a list or NULL.")
    all_elements <- unlist(ordering)
    if (!all(all_elements %in% valid_groups))
      stop("`ordering` can only contain 'M1', 'M2', and 'M3'.")
    if (any(duplicated(all_elements))) {
      dup_elems <- unique(all_elements[duplicated(all_elements)])
      stop(paste0("`ordering` contains duplicated group(s): ",
                  paste(dup_elems, collapse = ", "), "."))
    }
    ordering <- lapply(ordering, as.character)
  }
  
  # Case 0: No ordering restrictions
  if (is.null(ordering)) {
    return(exp(lfactorial(p_tot)))
  }
  
  # Case 1: Simple hierarchy (e.g., list("M1","M2"))
  if (length(ordering) == 2 &&
      all(sapply(ordering, function(x) length(x) == 1))) {
    sizes <- unlist(group_sizes[unlist(ordering)])
    return(exp(lfactorial(sizes[1]) + lfactorial(sizes[2]) +
                 lfactorial(p_tot) - lfactorial(sum(sizes))))
  }
  
  # Case 2: Full hierarchy (e.g., list("M1","M2","M3"))
  if (length(ordering) == 3 &&
      all(sapply(ordering, function(x) length(x) == 1))) {
    sizes <- unlist(group_sizes[unlist(ordering)])
    return(exp(sum(lfactorial(sizes))))
  }
  
  # Case 3: Partial block hierarchy (combined)
  if (length(ordering) == 2 &&
      any(sapply(ordering, length) > 1) &&
      any(sapply(ordering, length) == 1)) {
    block_sizes <- c(
      sum(unlist(group_sizes[ordering[[1]]])),
      sum(unlist(group_sizes[ordering[[2]]]))
    )
    return(exp(sum(lfactorial(block_sizes))))
  }
  
  # Default fallback (unrecognized ordering)
  return(exp(lfactorial(p_tot)))
}






