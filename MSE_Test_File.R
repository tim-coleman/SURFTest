##########################################
## Loading Packages
##########################################

library(ggplot2)
library(rpart)
library(dplyr)
library(party)
library(glmnet)
library(randomForest)

##########################################
## Function to subsample various trees
##########################################
bag.s <- function(X, y, base.learner = "rpart", ntree, k, mtry = ncol(X),
                  form = if(base.learner %in% c("rpart", "ctree", "lm")) as.formula("y~."),
                  alpha = if(base.learner == "lm") 1,
                  lambda = if(base.learner == "lm") 1){

  D <- data.frame(X, y)

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
    c.sim <- ctree_control(minsplit=1,
                           maxsurrogate=0, mtry = mtry,
                           mincriterion = .3)
    fun = function(){
      ind <- sample(1:nrow(D), size=k, replace=FALSE)
      ctree(form, controls = c.sim, data = D[ind,])
    }
  }


  if(base.learner == "rtree") {
    fun = function(){
      randomForest(x = X, y = y, ntree = 1, mtry = mtry,
                   replace = FALSE,
                   sampsize = floor(k))
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
  rfs
}


##########################################
## Function to run the test.
##########################################
MSE_Test.default <- function(X, y,  X.test = FALSE, y.test = FALSE,
                     var,  NTest = nrow(X.test), resp = names(y),
                     B = 1000, NTree = 500, p = 1/2, base.learner = "rpart",
                     mtry = ncol(X), importance = T){

  if(mtry > ncol(X)){
    warning("mtry is greater than number of columns in design matrix \n setting mtry = ncol(X)")
    mtry = ncol(X)
  }

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
  #print(form.resp)

  for(v in var){
    X.train.pm[,v] <- X.train.pm[sample(nrow(X.train)), v]
  }
  rf_og <- bag.s(X = X.train, y = y.train, base.learner = base.learner,
                 ntree = NTree, k = ceiling(N^p), mtry = mtry, form.resp)
  rf_pm <- bag.s(X = X.train.pm, y = y.train, base.learner = base.learner,
                 ntree = NTree, k = ceiling(N^p), mtry = mtry, form.resp)
  if(base.learner == "lm"){
    #print(head(data.frame(X.test, "y" = y.test)))
    #print(head(model.matrix(form.resp, data = data.frame(X.test, "y" = y.test))))
    P <- data.frame(lapply(rf_og, predict, newx = model.matrix(form.resp, data = data.frame(X.test, "y" = y.test))))
    PR <- data.frame(lapply(rf_pm, predict, newx = model.matrix(form.resp, data = data.frame(X.test, "y" = y.test))))
  }
  else{
  P <- data.frame(lapply(rf_og, predict, newdata = X.test))
  PR <- data.frame(lapply(rf_pm, predict, newdata = X.test))
  }
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
       "originalMSE" = c("Original MSE" = MSE_0,
                         "Permuted MSE" = MSE_R_0),
         "PermDiffs" = MSE.diffs,
       "Importance" = importances,
        "Pvalue" = Pvals,
       "test_pts" = X.test,
       "weak_learner" = base.learner,
       "call" = match.call())
  class(result) <- "MSE_Test"
  result
}

MSE_Test <- function(X, y, ...) UseMethod("MSE_Test")

####################################
### Defining a plot method
####################################

plot.MSE_Test <- function(obj){
  attach(obj)
  diff.0 <- originalMSE[2] - originalMSE[1]
  var.title <- paste(variables, collapse = ", ")
  title.null.MSE <- paste("MSE Differences\n", "Testing for", var.title, sep = " ")


  M0 <- min(PermDiffs, diff.0) - 1.5*sd(PermDiffs)
  M1 <- max(PermDiffs, diff.0) + 1.5*sd(PermDiffs)

  g <- ggplot(aes(x = dif), data = data.frame(dif =PermDiffs)) +
    geom_histogram(aes(y = ..density..), binwidth = sd(PermDiffs)/3, col = '#1F49B8', fill = '#6FD6FF') +
    geom_hline(yintercept = 0) + ylab("Density") +
    stat_function(fun = dnorm, args = list(mean = mean(PermDiffs), sd = sd(PermDiffs)),
                  size = 1.075, alpha = .75, col = '#1F49B8',
                  xlim = c(min(PermDiffs, diff.0) - 1.5*sd(PermDiffs),
                           max(PermDiffs) + 1.5*sd(PermDiffs)),
                  n = 750) +
    scale_x_continuous(name = "Reduced MSE - Full MSE",
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
}

####################################
### Definining a formula method
####################################
MSE_Test.formula <- function(formula, data=list(), ...)
{
  formstring <- unlist(strsplit(format(formula), split = " "))
  resp <- formstring[1]

  y <- data[,resp]
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
  cat("Variables Tested:", obj$variables, "\n")
  cat("Number of Test Points:", nrow(obj$test_pts), "\n")
  cat("Base Learner:", obj$weak_learner, "\n")
  cat("P-value:", obj$Pvalue, "\n")
  if(!is.null(obj$Importance)) cat("SD Importance:", obj$Importance[1])
}


############################################
### Defining an Importance Measure Function
############################################

permtestImp <- function(X, y, X.test = FALSE, y.test = FALSE,
                        NTest = nrow(X.test), Nbtree = 30, verbose = F,
                        base.learner = "rpart", mtry = ncol(X),
                        ...){
  vars <- names(X)
  if(length(Nbtree) == 1) Nbtree = rep(Nbtree, length(vars))

  fun.ret <- function(var, ntree){
    if(verbose) cat("Testing for ", var, " with ", ntree, " trees ...")
    obj <- MSE_Test(X = X, y = y, X.test = X.test, y.test = y.test,
                    NTest = NTest, NTree = ntree,
                    base.learner = base.learner, mtry = mtry, var = var,
                    ...)
    if(verbose) cat(" done! \n")
    list("SDImp" = obj$Importance[1],
               "Pval" = obj$Pvalue)
  }

  out <- mapply(fun.ret, vars, Nbtree)

  class(out) <- "permtestImp"
  out
}


### Definining an Importance Plot Method
plot.permtestImp <- function(obj, col_blind = F){
  title.perm <- "Permutation Test Importance Measures "
  df <- data.frame("Variable" = colnames(obj),
                   "SDImportance" = unlist(obj[1,], use.names = F),
                   "Pvalues" = unlist(obj[2,], use.names = F))
  df <- df %>% mutate(Variable = reorder(Variable, SDImportance)) %>% 
    mutate(signifcode = cut(Pvalues, breaks = c(0,.01, .05, .1,1), labels =
                              c(".01", ".05", ".1", "Not Significant")))
  if(!col_blind){
  g <- ggplot(aes(x = Variable), data = df) +
    geom_point(aes(y = SDImportance, col = signifcode), size = 3.5) + 
    coord_flip() + theme_classic() + ggtitle(title.perm) + ylab("Standard Deviation Importance") +
    scale_colour_manual(values = rev(c("#000000", "#FF8075", "#FF5C55", "#FF3834")), 
                        name = "Significance \nLevel")+
    theme(plot.title = element_text(size = 21, hjust = .5), axis.title = element_text(size = 15),
          axis.text.y = element_text(size = 18),
          axis.text.x = element_text(size = 18),
          legend.background = element_rect(colour = "black", fill = "white"),
          legend.text = element_text(size = 12), legend.title = element_text(size = 14, hjust = .5),
          legend.spacing.y = unit(2, "cm"), plot.subtitle = element_text(hjust = .5),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white"),
          panel.grid.major.y = element_line(colour = 'gray', linetype = 2))
  }
  else{
    g <- ggplot(aes(x = Variable), data = df) +
      geom_point(aes(y = SDImportance, shape = signifcode), fill = '#6FD6FF', size = 3.5) + 
      coord_flip() + theme_classic() + ggtitle(title.perm) + ylab("Standard Deviation Importance") +
      scale_shape_discrete(name = "Significance \nLevel")+
      theme(plot.title = element_text(size = 21, hjust = .5), axis.title = element_text(size = 15),
            axis.text.y = element_text(size = 18),
            axis.text.x = element_text(size = 18),
            legend.background = element_rect(colour = "black", fill = "white"),
            legend.text = element_text(size = 12), legend.title = element_text(size = 14, hjust = .5),
            legend.spacing.y = unit(2, "cm"), plot.subtitle = element_text(hjust = .5),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            plot.background = element_rect(fill = "white"),
            panel.background = element_rect(fill = "white"),
            panel.grid.major.y = element_line(colour = 'gray', linetype = 2))
  }
    plot(g)
}


