##########################################
## Loading Packages
##########################################

library(ggplot2)
library(rpart)
library(dplyr)
library(ranger)
library(party)
library(glmnet)
library(randomForest)

##########################################
## Function to subsample various models
##########################################
bag.s <- function(X, y, base.learner = "rpart", ntree, k, mtry = ncol(X),
                  form = as.formula("y~."),
                  alpha = if(base.learner == "lm") 1 else NULL ,
                  lambda = if(base.learner == "lm") 1 else NULL, ranger = F){
  
  y <- unlist(y) # vectorizing y
  #resp_name <- strsplit(as.character(form), split = "")[[1]][1]
  resp_name <- as.character(form)[2]
  D <- data.frame(X, y)
  names(D)[ncol(X) + 1] <- resp_name

  if(!(base.learner %in% c("rpart", "ctree", "rtree", "lm"))){
    stop("Base learner should be one of 'rpart', 'ctree', 'lm' or 'rtree' ")
  }

  if(nrow(D) == 0) stop("data (D) has 0 rows")

  if(base.learner == "rpart"){
    c.sim <- rpart.control(minsplit=10, maxcompete=0, maxsurrogate=0, usesurrogate=0, cp = 0)
    fun = function(){
      ind <- sample(1:nrow(D), size=k, replace=FALSE)
      rpart(form, data = D[ind,], control = c.sim)
    }
  }

  if(base.learner == "ctree"){
    if(ranger){
      fun = function(){ranger(form, num.trees = 1, mtry = mtry, importance = 'none', data = D,
             min.node.size = 1, replace = F, sample.fraction = k/nrow(D), splitrule = "maxstat")}
    }
    else{
    c.sim <- ctree_control(minsplit=1,
                           maxsurrogate=0, mtry = mtry,
                           mincriterion = .3)
    fun = function(){
      ind <- sample(1:nrow(D), size=k, replace=FALSE)
      ctree(form, controls = c.sim, data = D[ind,])
    }
    }
  }
    


  if(base.learner == "rtree") {
    if(ranger){
      fun =  function(){ranger(form, num.trees = 1, mtry = mtry, importance = 'none', data = D,
             min.node.size = 1, replace = F, sample.fraction = k/nrow(D))
      }
    }
    else{
    fun = function(){
      randomForest(x = X, y = y, ntree = 1, mtry = mtry,
                   replace = FALSE,
                   sampsize = floor(k))
    }
    }

  }
  
  if(base.learner == "lm"){
    fun = function(){
      ind <- sample(1:nrow(D), size=k, replace=FALSE)
      features <- model.matrix(form, data = D[ind,])
      glmnet(x = features, y = y[ind], alpha = alpha, lambda = lambda)      
    }
  }

  rfs <- list()
  for(i in 1:ntree){
    rfs[[i]] <- fun()
  }
  if(ranger & base.learner %in% c("rtree", "ctree")){
    class(rfs) <- paste("bagged", base.learner, "ranger", sep = "_")
    
  }
  else{
  class(rfs) <- paste("bagged", base.learner, sep = "_")
  }
  rfs
}


##########################################
## Function to run the test.
##########################################
MSE_Test.default <- function(X, y,  X.test = FALSE, y.test = FALSE,
                     var,  NTest = nrow(X.test), resp = names(y),
                     B = 1000, NTree = 500, p = 1/2, base.learner = "rpart",
                     mtry = ncol(X), importance = T, alpha = if(base.learner == "lm") 1,
                     lambda = if(base.learner == "lm") 1, ranger = F){

  if(mtry > ncol(X)){
    warning("mtry is greater than number of columns in design matrix \n setting mtry = ncol(X)")
    mtry = ncol(X)
  }
  y <- unlist(y)
  N <- nrow(X)
  if(is.logical(X.test)) {
    test <- sample(nrow(X), NTest)
    X.test <- X[test,]
    y.test <- y[test]
    X.train <- X[-test,]
    X.train.pm <- X.train
    y.train <- y[-test]
  }
  else{
    X.train <- X
    X.train.pm <- X
    y.train <- y
  }

  form.resp <- as.formula(paste(resp, "~.", sep = ""))

  for(v in var){
    X.train.pm[,v] <- X.train.pm[sample(nrow(X.train)), v]
  }
  rf_og <- bag.s(X = X.train, y = y.train, base.learner = base.learner,
                 ntree = NTree, k = ceiling(N^p), mtry = mtry, form.resp, 
                 alpha = alpha, lambda = lambda, ranger = ranger)
  rf_pm <- bag.s(X = X.train.pm, y = y.train, base.learner = base.learner,
                 ntree = NTree, k = ceiling(N^p), mtry = mtry, form.resp,
                 alpha = alpha, lambda = lambda, ranger = ranger)
  
  pred_f <- function(model){
    if(class(model) %in% c("bagged_rtree", "bagged_ctree", "bagged_rpart")){
      P <- data.frame(lapply(model, FUN = function(x) predict(x, newdata = X.test)))
    }
    else if(class(model) %in% c("bagged_rtree_ranger","bagged_ctree_ranger")){
      P <- data.frame(lapply(model, FUN = function(x){ P <- predict(x, data = X.test); P["predictions"] }))
    }
    else{
      if(is.null(y.test)){y.dummy <- rep(0, nrow(X.test))}
      else{y.dummy <- y.test}
      form.resp <- as.formula("y~.")
      P <- data.frame(lapply(model, predict, newx = as.matrix(X.test)))
    }
    P
  }
  # if(base.learner == "lm"){
  #   #print(head(data.frame(X.test, "y" = y.test)))
  #   #print(head(model.matrix(form.resp, data = data.frame(X.test, "y" = y.test))))
  #   P <- data.frame(lapply(rf_og, predict, newx = model.matrix(form.resp, data = data.frame(X.test, "y" = y.test))))
  #   PR <- data.frame(lapply(rf_pm, predict, newx = model.matrix(form.resp, data = data.frame(X.test, "y" = y.test))))
  # }
  # else{
  # P <- data.frame(lapply(rf_og, predict, newdata = X.test))
  # PR <- data.frame(lapply(rf_pm, predict, newdata = X.test))
  # }
  P <- pred_f(rf_og)
  PR <- pred_f(rf_pm)
  Pred_0 <- apply(P, FUN = mean, MARGIN = 1)
  Pred_R_0 <- apply(PR, FUN = mean, MARGIN = 1)

  MSE_0 <- mean((Pred_0 - y.test)^2)
  MSE_R_0 <- mean((Pred_R_0 - y.test)^2)
  diff.0 <- MSE_R_0 - MSE_0
  Pool <- cbind(P, PR)
  MSE <- data.frame("Full_MSE" = c(0), "Reduced_MSE" = c(0)) # Intialization
  for(i in 1:B){
    samps <- sample(1:(2*NTree), NTree, replace = F)
    P_t <- Pool[,samps]
    PR_t <- Pool[,-samps]
    Pred_t <- apply(P_t, FUN = mean, MARGIN = 1)
    Pred_R_t <- apply(PR_t, FUN = mean, MARGIN = 1)
    MSE[i,] <-c(mean((Pred_t - y.test)^2), mean((Pred_R_t - y.test)^2))
  }
  MSE.diffs <- MSE$Reduced_MSE - MSE$Full_MSE
  Pvals <- c("Full Model P" = mean(c(diff.0 < MSE.diffs, 1)))

  if(importance){
    sdimp <- (diff.0 - mean(MSE.diffs))/sd(MSE.diffs)
    zimp <- pnorm(sdimp)
    importances <- c("Standard Deviation Importance" = sdimp,
                     "Standard Normal Importance" = zimp)
  }
  else{ importances = NULL}

  result <- list("variables" = var,
       "originalStat" = c("Original MSE" = MSE_0,
                         "Permuted MSE" = MSE_R_0),
         "PermDiffs" = MSE.diffs,
       "Importance" = importances,
        "Pvalue" = Pvals,
       "test_pts" = X.test,
       "weak_learner" = base.learner,
       "model_original" = rf_og,
       "model_permuted" = rf_pm,
       "test_stat" = "MSE",
       "call" = match.call())
  class(result) <- "MSE_Test"
  result
}

MSE_Test <- function(X, y, ...) UseMethod("MSE_Test")



####################################
## Function to compare two pre-trained models
####################################
MSE_compare <- function(m1, m2, X.test, y.test = NULL, B = 1000, return.preds = F, test_stat = if(is.null(y.test)) "KS" else "MSE"){
  #ms <- list(m1, m2)
  if(length(m1) != length(m2)){
    warning("Models should have the same number of base models")
  }
  pred_f <- function(model){
    if(class(model) %in% c("bagged_rtree", "bagged_ctree", "bagged_rpart")){
      P <- data.frame(lapply(model, FUN = function(x) predict(x, newdata = X.test)))
    }
    else if(class(model) %in% c("bagged_rtree_ranger","bagged_ctree_ranger")){
      P <- data.frame(lapply(model, FUN = function(x){ P <- predict(x, data = X.test); P["predictions"] }))
      #P <- lapply(model, FUN = function(x){ P <- predict(x, data = X.test); P["predictions"] })
    }
    else{
      if(is.null(y.test)){y.dummy <- rep(0, nrow(X.test))}
      else{y.dummy <- y.test}
      form.resp <- as.formula("y~.")
      P <- data.frame(lapply(model, predict, newx = as.matrix(X.test)))
    }
    P
  }
  #print(c(length(m1), length(m2)))
  P0 <- pred_f(m1)
  P1 <- pred_f(m2)
  
  
  nt1 <- length(m1)
  nt2 <- length(m2)
  Pool <- cbind(P0, P1)
  
  if(class(test_stat) == "character"){
    if(test_stat == "MSE"){
    f <- function(test, p1, p2){
      #print(c("test_l" = length(test), "p1_l" = length(p1), 
       #       "p2_l" = length(p2)))
      MSE1 <- mean((p1 - test)^2)
      MSE2 <- mean((p2 - test)^2)
      return(c("Original MSE" = MSE1,
               "Permuted MSE" = MSE2))
    }
  }
  
  else if(test_stat == "KS"){
    f <- function(test, p1, p2){
      stat <- max(abs(p1 - p2))
      return(stat)
    }
  }
  }
  else{f <- test_stat}
  
  Pred_0 <- apply(P0, FUN = mean, MARGIN = 1)
  Pred_R_0 <- apply(P1, FUN = mean, MARGIN = 1)
  if(is.null(y.test)){
    y.test <- rep(0, nrow(X.test))
  }
  TS0 <- f(y.test, Pred_R_0, Pred_0)
  TS_temp <- TS0
  if(class(test_stat) == "character"){
    if(test_stat == "MSE") TS0 <- TS0[2] - TS0[1]
  }
  
  perm_stats <- c()
  for(i in 1:B){
        samps <- sample(1:ncol(Pool), nt1, replace = F)
        P_t <- Pool[,samps]
        PR_t <- Pool[,-samps]
        Pred_t <- apply(P_t, FUN = mean, MARGIN = 1)
        Pred_R_t <- apply(PR_t, FUN = mean, MARGIN = 1)
        if(class(test_stat) == "character"){
          if(test_stat == "MSE"){
            temp <- f(y.test, Pred_R_t, Pred_t)
            perm_stats[i] <- temp[2] - temp[1]
            }
          else{
          temp <- NULL
          perm_stats[i] <- f(y.test, Pred_R_t, Pred_t)
          }
        }
        else{
          perm_stats[i] <- f(y.test, Pred_R_t, Pred_t)
        }
  }
  Pvals <- c("Full Model P" = mean(c(TS0 < perm_stats, 1)))

  sdimp <- (TS0 - mean(perm_stats))/sd(perm_stats)
  zimp <- pnorm(sdimp)
  importances <- c("Standard Deviation Importance" = sdimp,
                   "Standard Normal Importance" = zimp)
  
  result <- list("variables" = NULL,
                 "originalStat" = TS_temp,
                 "PermDiffs" = perm_stats,
                 "Importance" = importances,
                 "Pvalue" = Pvals,
                 "test_pts" = X.test,
                 "weak_learner" = c(class(m1), class(m2)),
                 "model_original" = m1,
                 "model_permuted" = m2,
                 "test_stat" = test_stat,
                 "call" = match.call())
  class(result) <- "MSE_Test"
  result
  }

####################################
### Defining a plot method
####################################

plot.MSE_Test <- function(obj, return_plot = F){
  attach(obj)
  if(test_stat == "MSE"){
    diff.0 <- originalStat[2] - originalStat[1]
    if(!is.null(variables)){
      xtitle <- "Reduced MSE - Full MSE"
    }
    else{
    xtitle <- paste("MSE(", weak_learner[2],") - MSE(", weak_learner[1],")", sep = "")
    }
    }
  else{
    diff.0 <- originalStat
    if(test_stat == "KS"){
      xtitle <- expression(paste("max|",Model[1] - Model[2],"|"))
    }
    else{xtitle <- NULL}
  }
  if(!is.null(variables)){
    var.title <- paste(variables, collapse = ", ")
    title.null.MSE <- paste("MSE Differences\n", "Testing for", var.title, sep = " ")
  }
  else{
    title.null.MSE <- paste("Model Comparison between \n", weak_learner[1], "and", weak_learner[2])
  }


  M0 <- min(PermDiffs, diff.0) - 1.5*sd(PermDiffs)
  M1 <- max(PermDiffs, diff.0) + 1.5*sd(PermDiffs)

  g <- ggplot(aes(x = dif), data = data.frame(dif = PermDiffs)) +
    geom_histogram(aes(y = ..density..), binwidth = sd(PermDiffs)/3, col = '#1F49B8', fill = '#6FD6FF') +
    geom_hline(yintercept = 0) + ylab("Density") +
    stat_function(fun = dnorm, args = list(mean = mean(PermDiffs), sd = sd(PermDiffs)),
                  size = 1.075, alpha = .75, col = '#1F49B8',
                  xlim = c(min(PermDiffs, diff.0) - 1.5*sd(PermDiffs),
                           max(PermDiffs) + 1.5*sd(PermDiffs)),
                  n = 750) +
    scale_x_continuous(name = xtitle,
                       limits = c(M0, M1)) +
    ggtitle(title.null.MSE, sub = "") + theme_classic() +
    geom_vline(aes(xintercept = diff.0), col = 'tomato', lwd = 1.45, alpha = .8) +
    theme(plot.title = element_text(size = 21, hjust = .5), axis.title = element_text(size = 15),
          axis.text.y = element_text(size = 18),
          axis.text.x = element_text(size = 18),
          legend.background = element_rect(colour = "black", fill = "white"),
          legend.text = element_text(size = 12), legend.title = element_text(size = 14),
          legend.spacing.y = unit(2, "cm"), plot.subtitle = element_text(hjust = .5),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white"))

  detach(obj)
  plot(g)
  if(return_plot){
    return(g)
  }
}

####################################
### Definining a formula method
####################################
MSE_Test.formula <- function(formula, data=list(), ...)
{
  formstring <- unlist(strsplit(format(formula), split = " "))
  resp <- formstring[1]

  y <- data[resp]
  data[resp] <- NULL

  est <- MSE_Test.default(X = data, y = y, ...)
  est$call <- match.call()
  est$formula <- formula
  est
}


####################################
### Defining a summary method
####################################

summary.MSE_Test <- function(obj){
  cat("Call: \n")
  print(obj$call)
  cat("\n")
  cat("Statistic used:", obj$test_stat, "\n")
  if(!is.null(obj$variables)) cat("Variables Tested:", obj$variables, "\n")
  cat("Number of Test Points:", nrow(obj$test_pts), "\n")
  cat("Base Learner(s):", obj$weak_learner, "\n")
  cat("P-value:", obj$Pvalue, "\n")
  if(!is.null(obj$Importance)) cat("SD Importance:", obj$Importance[1])
}


############################################
### Defining an Importance Measure Function
############################################

permtestImp <- function(X, y, X.test = FALSE, y.test = FALSE, single_forest = T,
                        NTest = nrow(X.test), Nbtree = 30, verbose = F,
                        keep_forest = F,
                        base.learner = "rpart", mtry = ncol(X), p = .5,
                        ...){
  vars <- names(X)
  if(length(Nbtree) == 1) Nbtree = rep(Nbtree, length(vars))
  if(single_forest){
   y <- unlist(y)
   N <- nrow(X)
   if(is.logical(X.test)){
     test <- sample(nrow(X), NTest)
     X.test <- X[test,]
     y.test <- y[test]
     X.train <- X[-test,]
     X.train.pm <- X.train
     y.train <- y[-test]
   }
   else{
     X.train <- X
     X.train.pm <- X
     y.train <- y
   }
   m0 <- bag.s(X = X.train, y = y.train, base.learner = base.learner, ntree = max(Nbtree), mtry = mtry, k = floor(nrow(X)^p), ... )
   fun.ret <- function(var, ntree){
    if(verbose) cat("Testing for ", var, " with ", ntree, " trees ...")
     X.temp <- X.train.pm
      for(v in var){
      X.temp[,v] <- X.temp[sample(nrow(X.train)), v]
      }
     m1 <- bag.s(X = X.temp, y = y.train, base.learner = base.learner, ntree = ntree, mtry = mtry, k = floor(nrow(X)^p), ... )
     obj <- MSE_compare(m0, m1, X.test = X.test, y.test = y.test, B = 1000)
     if(verbose) cat(" done! \n")
     list("SDImp" = obj$Importance[1],
          "MSEImp" = obj$originalStat[2] - obj$originalStat[1],
          "Pval" = obj$Pvalue)
   }
 }
  else{
  fun.ret <- function(var, ntree){
    if(verbose) cat("Testing for ", var, " with ", ntree, " trees ...")
    obj <- MSE_Test(X = X, y = y, X.test = X.test, y.test = y.test,
                    NTest = NTest, NTree = ntree,
                    base.learner = base.learner, mtry = mtry, var = var,
                    ...)
    if(verbose) cat(" done! \n")
    list("SDImp" = obj$Importance[1],
         "MSEImp" = obj$originalStat[2] - obj$originalStat[1],
               "Pval" = obj$Pvalue)
  }
  }
  
  out <- mapply(fun.ret, vars, Nbtree)
  if(keep_forest){
    out <- list("Importance_Table" = out, "TrainedModel" = m0)
  }
  class(out) <- "permtestImp"
  return(out)
}


### Feature Holdout Forest Method

f_holdoutRF <- function(X, y, X.test = FALSE, y.test = FALSE, 
                        NTest = nrow(X.test), mintree = 30, max.trees = 5*ncol(X)*mintree, verbose = F,
                        mtry = ncol(X)/3, p = .5, keep_forest = F, ...){
  if(mtry == ncol(X)) stop("mtry should be less than total number of columns to use holdout procedure")
  # Declaring variable names
  vars <- names(X)
  # Gathering the test points
  if(is.logical(X.test)){
    test <- sample(nrow(X), NTest)
    X.test <- X[test,]
    y.test <- y[test]
    X.train <- X[-test,]
    y.train <- y[-test]
  }
  else{
    X.train <- X
    y.train <- y
  }
  # Vector to keep track of how often each variable has been used
  ntree_with_var <- rep(0, length(vars))
  names(ntree_with_var) <- vars
  # Matrix keeping track of which trees have which variables
  tree_var_mat <- matrix(nrow = 0, ncol = length(vars))
  # Counter for total # of trees
  B <- 1
  # Initializing the ensemble
  vars_temp <- sample(vars, mtry, replace = F)
  X_temp <- X.train[vars_temp]
  mods <- bag.s(X = X_temp, y = y.train, ntree = 1, k = nrow(X.train)^p, mtry = ncol(X_temp), ...)
  c0 <- class(mods)
  ntree_with_var[vars_temp] <- ntree_with_var[vars_temp] + 1
  tree_var_mat <- rbind(tree_var_mat, vars %in% vars_temp)  # Training the ensemble
  while(any(ntree_with_var < mintree) & B < max.trees){
    vars_temp <- sample(vars, mtry, replace = F)
    X_temp <- X.train[vars_temp]
    mods <- c(mods, bag.s(X = X_temp, y = y.train, ntree = 1, k = nrow(X.train)^p, mtry = ncol(X_temp), ...))
    B <- B + 1
    ntree_with_var[vars_temp] <- ntree_with_var[vars_temp] + 1
    tree_var_mat <- rbind(tree_var_mat, vars %in% vars_temp)
    #print(B)
    #print(ntree_with_var)
  }

  tree_var_mat <- data.frame(tree_var_mat)
  names(tree_var_mat) <- vars
  if(any(ntree_with_var < mintree)) warning("Some models may not have enough base learners, consider running again with more models")
  # Comparing all models with variable X against those without variable X
  out_list <- vector(mode  = "list", length = length(vars))
  names(out_list) <- vars
  for(v in vars){
    flags <- tree_var_mat[[v]]
    mod_v <- mods[flags]
    mod_no_v <- mods[!flags]
    #print(c("var In" = length(mod_v), "var Out" = length(mod_no_v)))
    
    # Making sure forests are of equal size
    mod_v <- sample(mod_v, min(length(mod_no_v), length(mod_v)))
    mod_no_v <- sample(mod_no_v, min(length(mod_no_v), length(mod_v)))
    
    class(mod_v) <- c0
    class(mod_no_v) <- c0

    # Comparing the MSE's of each model
    obj <- MSE_compare(m1 = mod_no_v, m2 = mod_v, X.test = X.test, y.test = y.test)
   # results <- cbind(results, c("SDImp" = obj$Importance[1],
                   #      "MSEImp" = obj$originalStat[2] - obj$originalStat[1],
                  #       "Pval" = obj$Pvalue))
    out_list[[v]] <-  c("SDImp" = obj$Importance[1],
                             "MSEImp" = obj$originalStat[2] - obj$originalStat[1],
                               "Pval" = obj$Pvalue)
  }
  out <- do.call(cbind, out_list)
  rownames(out) <- c("SDImp", "MSEImp", "Pval")
  if(keep_forest){
    out <- list("Importance_Table" = out, "TrainedModel" = mods)
  }
  else{
    out <- list("Importance_Table" = out)
  }
  class(out) <- "permtestImp"
  return(out)
}


### Definining an Importance Plot Method
plot.permtestImp <- function(obj, ImpType = "SDImp", col_blind = F){
  obj <- obj$Importance_Table
  title.perm <- "Permutation Test Importance Measures "
  if(ImpType == "SDImp"){
  df <- data.frame("Variable" = colnames(obj),
                   "Importance" = unlist(obj[1,], use.names = F),
                   "Pvalues" = unlist(obj[3,], use.names = F))
  df <- df %>% mutate(Variable = reorder(Variable, Importance)) %>% 
    mutate(signifcode = cut(Pvalues, breaks = c(0,.01, .05, .1,1), labels =
                              c(".01", ".05", ".1", "Not Significant")))
  subtitle <- "Standard Deviation Importance"
  }
  else{
    df <- data.frame("Variable" = colnames(obj),
                     "Importance" = unlist(obj[2,], use.names = F),
                     "Pvalues" = unlist(obj[3,], use.names = F))
    df <- df %>% mutate(Variable = reorder(Variable, Importance)) %>% 
      mutate(signifcode = cut(Pvalues, breaks = c(0,.01, .05, .1,1), labels =
                                c(".01", ".05", ".1", "Not Significant")))
    subtitle <- "MSE Importance"
  }
  # Reducing # of variables shown for ease of plotting
  if(nrow(df) > 15){
    df <- df[1:15,]
  }
  if(!col_blind){
  g <- ggplot(aes(x = Variable), data = df) +
    geom_point(aes(y = Importance, col = signifcode), size = 3.5) + 
    coord_flip() + theme_classic() + ggtitle(title.perm) + ylab(subtitle) +
    scale_colour_manual(values = rev(c("#000000", "#FF8075", "#FF5C55", "#FF3834")), 
                        name = "Significance \nLevel")+
    theme(plot.title = element_text(size = 21, hjust = .5), axis.title = element_text(size = 15),
          axis.text.y = element_text(size = 18),
          axis.text.x = element_text(size = 18),
          legend.background = element_rect(colour = "black", fill = "white"),
          legend.text = element_text(size = 12), legend.title = element_text(size = 14, hjust = .5),
          legend.spacing.y = unit(.3, "cm"), plot.subtitle = element_text(hjust = .5),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white"),
          panel.grid.major.y = element_line(colour = 'gray', linetype = 2))
  }
  else{
    g <- ggplot(aes(x = Variable), data = df) +
      geom_point(aes(y = Importance, shape = signifcode), fill = '#6FD6FF', size = 3.5) + 
      coord_flip() + theme_classic() + ggtitle(title.perm) + ylab(subtitle) +
      scale_shape_discrete(name = "Significance \nLevel")+
      theme(plot.title = element_text(size = 21, hjust = .5), axis.title = element_text(size = 15),
            axis.text.y = element_text(size = 18),
            axis.text.x = element_text(size = 18),
            legend.background = element_rect(colour = "black", fill = "white"),
            legend.text = element_text(size = 12), legend.title = element_text(size = 14, hjust = .5),
            legend.spacing.y = unit(.3, "cm"), plot.subtitle = element_text(hjust = .5),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            plot.background = element_rect(fill = "white"),
            panel.background = element_rect(fill = "white"),
            panel.grid.major.y = element_line(colour = 'gray', linetype = 2))
  }
    plot(g)
}


