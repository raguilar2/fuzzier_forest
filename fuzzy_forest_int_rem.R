# computes X values for up to three way interactions and returns as data frame
interactions <- function(X, n, var_order = NULL, three_way = FALSE){
  if(n > ncol(X)) {
    stop("n must be less then or equal to ncol(X)")
  }
  
  if(n < 2) {
    stop("n must be greater than 1")
  }
  
  if(n < 3 & three_way == TRUE) {
    stop("n must be greater than 2 when looking for 3 way interaction")
  }
  
  if(!is.null(var_order)) {
    X <- X[, var_order]
  }
  
  if(three_way == TRUE) {
    X_int3 <- as.data.frame(matrix(0, nrow = nrow(X), ncol=choose(n, 3)))
    i <- 1
    for (j in 1:(n-2)) {
      for (k in (j + 1):(n-1)) {
        for (l in (k + 1):n) {
          X_int3[, i] <- X[, j] * X[, k] * X[, l]
          names(X_int3)[i] <- paste(names(X)[j], names(X)[k],
                                    names(X)[l], sep = "x")
          i <- i + 1 
        }
      }
    }
  }
  
  
  X_int2 <- as.data.frame(matrix(0, nrow = nrow(X), ncol=choose(n, 2)))
  i <- 1
  for (j in 1:(n-1)) {
    for (k in (j + 1):(n)) {
      X_int2[, i] <- X[, j] * X[, k]
      names(X_int2)[i] <- paste(names(X)[j], names(X)[k],
                                sep = "x")
      i <- i + 1
    }
  }
  
  if(three_way == TRUE) {
    cbind(X_int2, X_int3)
  } else {
    X_int2
  }
  
}


# subsets original data frame to only features in 
surv_data <- function(survivor_df, df) {
  df[, names(df) %in% survivor_df[, 1]]
}

# scales VIM by proportion within module
vi_prop <- function(x){
  vi_scaled <- x[,2]/sum(x[,2])
  x <- cbind(x, vi_scaled, 1:nrow(x))
}

fuzzy_forest_int_rem <- function(feature_list, final_rf, module_membership,
                         WGCNA_object=NULL, survivor_list, selection_list,
                         initial_feature_list,  initial_final_rf) {
  out <- list()
  out[[1]] <- feature_list
  out[[2]] <- final_rf
  out[[3]] <- module_membership
  out[[4]] <- module_membership
  out[[5]] <- survivor_list
  out[[6]] <- selection_list
  out[[7]] <- initial_feature_list
  out[[8]] <- initial_final_rf
  names(out) <- c("feature_list", "final_rf", "module_membership",
                  "WGCNA_object", "survivor_list", "selection_list",
                  "initial_feature_list", "initial_final_rf")
  class(out) <- "fuzzy_forest_int_rem"
  return(out)
}

# fuzzy forest with adjustments for within-module interactions
# NOTE: adjustments are highlighted below by bars of '#'
ff_int_rem <- function(X, y, Z=NULL, module_membership,
                          screen_params = screen_control(min_ntree=500),
                          select_params = select_control(min_ntree=500),
                          final_ntree = 5000,
                          num_processors=1, nodesize, test_features=NULL,
                          test_y=NULL, 
                          ############## added args below ############
                          interaction = FALSE, m = NULL, three_way = FALSE, 
                          across = FALSE, r = r, ...) { 
  CLASSIFICATION <- is.factor(y)
  if ( !((mode(y)=="numeric") || is.factor(y)) ) {
    stop("y must be a numeric vector or factor")
  }
  if( (!CLASSIFICATION) && (length(unique(y)) < 5) ) {
    warning("y has 5 or fewer unique values.  In this case, we recommend
            classification instead of regression.  For classification,
            y must be a factor.")
  }
  if(!is.data.frame(X)) {
    stop("X must be a data.frame.")
  }
  if(!is.null(Z)) {
    if (!is.data.frame(Z)) {
      stop("Z must be a data.frame.")
    }
  }
  if(CLASSIFICATION == TRUE) {
    if(missing(nodesize)){
      nodesize <- 1
    }
  }
  if(CLASSIFICATION == FALSE) {
    if(missing(nodesize)){
      nodesize <- 5
    }
  }
  screen_control <- screen_params
  select_control <-  select_params
  module_list <- unique(module_membership)
  if(num_processors > 1) {
    #set up parallel backend
    cl <- parallel::makeCluster(num_processors)
    parallel::clusterCall(cl, library, package = "randomForest", character.only = TRUE)
    doParallel::registerDoParallel(cl)
    #close parallel backend on exit
    on.exit(try(parallel::stopCluster(cl), silent=TRUE))
  }
  survivors <- vector('list', length(module_list))
  drop_fraction <- screen_control$drop_fraction
  mtry_factor <- screen_control$mtry_factor
  ntree_factor <- screen_control$ntree_factor
  min_ntree <- screen_control$min_ntree
  keep_fraction <- screen_control$keep_fraction
  if(ncol(X)*keep_fraction < select_control$number_selected){
    warning(c("ncol(X)*keep_fraction < number_selected", "\n",
              "number_selected will be set to floor(ncol(X)*keep_fraction)"))
    select_control$number_selected <- max(floor(ncol(X)*keep_fraction), 1)
  }
  
  for (i in 1:length(module_list)) {
    module <- X[, which(module_membership == module_list[i]), drop=FALSE]
    num_features <- ncol(module)
    #TUNING PARAMETER mtry_factor
    if(CLASSIFICATION == TRUE) {
      mtry <- min(ceiling(mtry_factor*sqrt(num_features)), num_features)
      if(missing(nodesize)){
        nodesize <- 1
      }
    }
    if(CLASSIFICATION == FALSE) {
      mtry <- min(ceiling(mtry_factor*num_features/3), num_features)
      if(missing(nodesize)){
        nodesize <- 5
      }
    }
    #TUNING PARAMETER ntree_factor
    ntree <- max(num_features*ntree_factor, min_ntree)
    #TUNING PARAMETER keep_fraction
    target <- ceiling(num_features * keep_fraction)
    while (num_features >= target){
      if(num_processors > 1) {
        rf <- foreach(ntree = rep(ntree/num_processors, num_processors),
                      .combine = combine, .packages = 'randomForest') %dopar% {
                        randomForest(module, y, ntree = ntree, mtry = mtry,
                                     importance = TRUE, scale = FALSE, nodesize=nodesize) }
      }
      if(num_processors == 1) {
        rf <- randomForest::randomForest(module, y, ntree = ntree, mtry = mtry,
                                         importance = TRUE, scale = FALSE,
                                         nodesize = nodesize)
      }
      var_importance <- randomForest::importance(rf, type=1, scale=FALSE)[, 1]
      var_importance <- var_importance[order(var_importance,
                                             decreasing=TRUE)]
      reduction <- ceiling(num_features*drop_fraction)
      if(num_features - reduction > target) {
        trimmed_varlist <- var_importance[1:(num_features - reduction)]
        features <- names(trimmed_varlist)
        module <- module[, which(names(module) %in% features)]
        num_features <- length(features)
        if(CLASSIFICATION == TRUE) {
          mtry <- min(ceiling(mtry_factor*sqrt(num_features)), num_features)
        }
        if(CLASSIFICATION == FALSE) {
          mtry <- min(ceiling(mtry_factor*num_features/3), num_features)
        }
        ntree <- max(num_features*ntree_factor, min_ntree)
      }
      else {
        num_features <- target - 1
        mod_varlist <- var_importance[1:target]
        features <- names(var_importance)[1:target]
        survivors[[i]] <- cbind(features, mod_varlist)
        row.names(survivors[[i]]) <- NULL
        survivors[[i]] <- as.data.frame(survivors[[i]])
        survivors[[i]][, 1] <- as.character(survivors[[i]][, 1])
        survivors[[i]][, 2] <- as.numeric(as.character(survivors[[i]][, 2]))
      }
    }
  }
  survivor_list <- survivors
  names(survivor_list) <- module_list
  survivors <- do.call('rbind', survivors)
  survivors <- as.data.frame(survivors, stringsAsFactors = FALSE)
  survivors[, 2] <- as.numeric(survivors[, 2])
  names(survivors) <- c("featureID", "Permutation VIM")
  
  X_surv <- X[, names(X) %in% survivors[, 1]]
  
  ########################### Adjustments ###############################
  
  if(interaction == TRUE & across == FALSE){
    X_surv_mod <- lapply(survivor_list, surv_data, df = X)
    X_surv_int <- do.call(cbind, lapply(X_surv_mod, interactions, 
                                        n = m, three_way = three_way))
    X_surv <- cbind(X_surv, X_surv_int)
  }
  
# testing scaling  
#  if(interaction == TRUE & across == TRUE){
#    surv_list_prop <- lapply(survivor_list, vi_prop)
#    all_surv_list_prop <- do.call(rbind, surv_list_prop)
#    var_order <- all_surv_list_prop[order(all_surv_list_prop[, 4], -all_surv_list_prop[, 3]), 1]
#    X_surv_int <- interactions(X_surv, n = m, var_order = var_order, 
#                               three_way = three_way)
#    X_surv <- cbind(X_surv, X_surv_int)
#  }
 
#  original across mthods   
  if(interaction == TRUE & across == TRUE){
    var_order <- survivors[order(survivors[, 2], decreasing = TRUE), 1]
    X_surv_int <- interactions(X_surv, n = m, var_order = var_order, 
                               three_way = three_way)
    X_surv <- cbind(X_surv, X_surv_int)
  }
  
  if(!is.null(Z)) {
    X_surv <- cbind(X_surv, Z, stringsAsFactors=FALSE)
  }
  #######################################################################
  ############### for loop added to remove top r features ###############
  if (r < 0) {
    stop("r must be greater than or equal to 0")
  } else if (r > 0) {
    l <- 2
  } else {
    l <- 1
  }
  
  final_lists <- list()
  final_rfs <- list()
  
  for (i in 1:l) {
    
    if(i == 2){
      if (r > nrow(final_list)) {
        stop("r > nrow(final_list)")
      }
      if (r >= ncol(X_surv)) {
        "not enough remaining features to remove r of them"
      }
      X_surv <- X_surv[, !(names(X_surv) %in% final_list[(1:r), 1])]
    }
    
    select_args <- list(X_surv, y, num_processors, nodesize)
    select_args <- c(select_args, select_control)
    names(select_args)[1:4] <- c("X", "y", "num_processors", "nodesize")
    select_results <- do.call("select_RF", select_args)
    final_list <- select_results[[1]][, 1, drop=F]
    selection_list <- select_results[[2]]
    row.names(final_list) <- NULL
    colnames(final_list) <- c("feature_name")
    final_list <- as.data.frame(final_list, stringsAsFactors=FALSE)
    #VIMs from last tree in recursive feature elimination should be
    #replaced.
    final_list <- cbind(final_list,
                        matrix(rep(".", 2*dim(final_list)[1]), ncol=2),
                        stringsAsFactors=F)
    ##################### changed X to X_surv #############################
    final_X <- X_surv[, names(X_surv) %in% final_list[, 1], drop=FALSE]
    #Some selected features may be from Z
    if(!is.null(Z)) {
      final_X <- cbind(final_X, Z[, names(Z) %in% final_list[, 1], drop=FALSE],
                       stringsAsFactors=FALSE)
    }
    current_p <- dim(final_X)[2]
    if(CLASSIFICATION == TRUE) {
      final_mtry <- min(ceiling(select_control$mtry_factor*sqrt(current_p)),
                        current_p)
    }
    if(CLASSIFICATION == FALSE) {
      final_mtry <- min(ceiling(select_control$mtry_factor*current_p/3),
                        current_p)
    }
    if(!is.null(test_features)) {
      test_features <- test_features[, which(names(test_features) %in%
                                               names(final_X))]
    }
    
    final_rf <- randomForest::randomForest(x=final_X, y=y, mtry=final_mtry, ntree=final_ntree,
                                           importance=TRUE, nodesize=nodesize,
                                           xtest=test_features, ytest=test_y)
    final_importance <- randomForest::importance(final_rf, type=1, scale = F)
    final_list[, 1] <- row.names(final_importance)
    final_list[, 2] <- final_importance[, 1]
    #Now it's very important to associate the right module to the right
    #feature.  The ordering must be correct.  This is made trickier by
    #by the fact that when Z is not null, there exist elements in the
    #the VIM list that aren't in X.
    
    #select_X is a vector with selected features in order of X.
    select_X <- names(X)[which(names(X) %in% final_list[, 1])]
    #select_mods is a vector with associated module memberships in order of X.
    select_mods <- module_membership[which(names(X) %in% final_list[, 1])]
    #select_order is a vector with selected features given according to
    #the order returned by randomForest.
    select_order <- final_list[, 1][which(final_list[,1] %in% names(X))]
    #select_mods is a vector with module memberships reordered according
    #to the order returned by randomForest
    select_mods <- select_mods[match(select_order, select_X)]
    #Here is where the module membership is entered into the table.
    #Note that for elements of Z, the entry will be "."
    final_list[, 3][final_list[, 1] %in% names(X)] <- select_mods
    names(final_list)[2:3] <- c("variable_importance", "module_membership")
    #Reorder vims so that they are in decreasing order.
    final_list <- final_list[order(final_list[, 2], decreasing=T), ]
    final_lists[[i]] <- final_list
    final_rfs[[i]] <- final_rf
  }
  
  module_membership <- as.data.frame(cbind(names(X), module_membership),
                                     stringsAsFactors=FALSE)
  names(module_membership) <- c("feature_name", "module")
  
  if (r == 0) {
    initial_final_list <- final_list <- final_lists[[1]]
    initial_final_rf <- final_rf <- final_rfs[[1]]
  } else {
    initial_final_list <- final_lists[[1]] 
    final_list <- final_lists[[2]]
    initial_final_rf <- final_rfs[[1]] 
    final_rf <- final_rfs[[2]]
  }
  
  out <- fuzzy_forest_int_rem(final_list, final_rf, module_membership,
                      survivor_list=survivor_list,
                      selection_list=selection_list,
                      initial_feature_list = initial_final_list, 
                      initial_final_rf = initial_final_rf)
  return(out)
}

# wff function but adjust to include arguments for interactions
wff_int_rem <- function(X, y, Z=NULL, WGCNA_params=WGCNA_control(power=6),
                           screen_params=screen_control(min_ntree=500),
                           select_params=select_control(min_ntree=500),
                           final_ntree=5000, num_processors=1, nodesize,
                           test_features=NULL, test_y=NULL, 
                           ############# interaction, m, three_way#################
                           interaction = FALSE, m = NULL, three_way = FALSE, 
                           across = FALSE, r = 0, ...) {
  if (!requireNamespace("WGCNA", quietly = T)) {
    stop("WGCNA must be installed.")
  }
  else{
    if(class(X) != "data.frame"){
      stop("X must be a data.frame")
    }
    if((!is.null(Z)) && (class(Z) != "data.frame")){
      stop("Z must be a data.frame")
    }
    numeric_test <- sapply(X, is.numeric)
    if (sum(numeric_test) != dim(X)[2]) {
      stop("To carry out WGCNA, all columns of X must be numeric.")
    }
    CLASSIFICATION <- is.factor(y)
    if(CLASSIFICATION == TRUE) {
      if(missing(nodesize)){
        nodesize <- 1
      }
    }
    if(CLASSIFICATION == FALSE) {
      if(missing(nodesize)){
        nodesize <- 5
      }
    }
    WGCNA_control <- WGCNA_params
    screen_control <- screen_params
    select_control <-  select_params
    WGCNA_args <- list(X, WGCNA_control$power)
    WGCNA_args <- c(WGCNA_args, WGCNA_control$extra_args)
    names(WGCNA_args) <- c("datExpr", "power", names(WGCNA_control$extra_args))
    bwise <- do.call("blockwiseModules", WGCNA_args)
    module_membership <- bwise$colors
    screen_drop_fraction <- screen_control$drop_fraction
    screen_keep_fraction <- screen_control$keep_fraction
    screen_mtry_factor <- screen_control$mtry_factor
    screen_ntree_factor <- screen_control$ntree_factor
    screen_min_ntree <- screen_control$min_ntree
    out <- ff_int_rem(X, y, Z, module_membership,
                         screen_control, select_control, final_ntree,
                         num_processors, nodesize=nodesize,
                         test_features=test_features, test_y=test_y,
                         ########## added args for everything below ############
                         interaction = interaction, m = m,
                         three_way = three_way,
                         across = across,
                         r = r) 
    out$WGCNA_object <- bwise
    return(out)
  }
}