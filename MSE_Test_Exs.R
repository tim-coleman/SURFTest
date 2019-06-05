
## Local Loading
#source('C:/Users/drain/Box Sync/Mentch Coleman Shared/RF Permutation Tests/R Stuff/MSE_Test_File.R')
source('~/Documents/RF Perm/R Stuff/MSE_Test_File.R')
library(dplyr)
set.seed(1500)

## Remote Loading
library(devtools)
source_url("https://raw.githubusercontent.com/tim-coleman/SURFTest/master/MSE_Test_File.R")
#^ may not work on government machines, try local loading :) 

## Generating dummy data

N <- 1250
Nvar <- 10
N_test <- 150
name_vec <- paste("X", 1:(2*Nvar), sep = "")

# training data:
X <- data.frame(replicate(Nvar, runif(N)), 
                replicate(Nvar, cut(runif(N), 3, 
                                      labels = as.character(1:3)))) %>%
  mutate(Y = 5*(X3) + .5*X2^2 +  rnorm(N, sd = .5))
names(X) <- c(name_vec, "Y")

# some testing data:
X.t1 <- data.frame(replicate(Nvar, runif(N_test)), 
                   replicate(Nvar, cut(runif(N_test), 3, 
                                       labels = as.character(1:3)))) %>%
  mutate(Y = 5*(X3) + .5*X2^2 +  rnorm(N_test, sd = .5))
names(X.t1) <- c(name_vec, "Y")


## Trying each base learner
b.rpart <- bag.s(X = X %>% select(-Y), y = X %>% select(Y),
                 base.learner = "rpart", ntree = 10, k = N^.85, mtry = 10, form = Y~.)
b.ctree <- bag.s(X = X %>% select(-Y), y = X %>% select(Y),
                 base.learner = "ctree", ntree = 10, k =N^.95, mtry = 2)
b.rf <- bag.s(X = X %>% select(-Y), y = X %>% select(Y),
              base.learner = "rtree", ntree = 10, k = N^.95, mtry = 2, Y~., ranger = F)
b.glmnet <-  bag.s(X = X %>% select(-Y), y = X %>% select(Y),
                   base.learner = "lm", ntree = 10, k = N^.95, mtry = 2)

## Generating an error
bag.s(X = X %>% select(-Y), y = X$Y,
      base.learner = "Not a real model!", ntree = 1, k = 100^.95, mtry = 2, Y~.)


### Trying out the testing method
# Not specifying test points:
M_no_test <- MSE_Test(X = X %>% dplyr::select(-Y), y = X$Y, 
                      base.learner = "rtree", NTest = 100, NTree = 150, B = 1000, var = c( "X3"),
                      p = .85, glm_cv = F)

plot(M_no_test)

# Specifying test points:
M_test <- MSE_Test(X = X %>% select(-Y), y = X$Y, X.test = X.t1 %>% select(-Y), y.test = X.t1$Y,
                      base.learner = "ctree", NTree = 250, B = 1000, var = c( "X2"),
                      p = .85)
plot(M_test)

### Formula implementation
M.formula <- MSE_Test(Y~., data = X, 
                      base.learner = "ctree", NTest = 100, NTree = 250, B = 1000, var = c( "X2", "X3"),
                      p = .85, mtry = 3, ranger = T)
plot(M.formula)

### Formula with prespecified test points
M_test <- MSE_Test(Y~., data = X, X.test = X.t1 %>% dplyr::select(-Y), y.test = X.t1$Y,
                   base.learner = "rtree", NTree = 250, B = 1000, var = c( "X2"),
                   p = .85)

# GLM implementation
M.form.glm <- MSE_Test(Y~., data = X, glm_cv = "external",
                       base.learner = "lm", NTest = 100, NTree = 250, B = 1000, var = c( "X2"),
                       p = .85, mtry = 3)

M.no_form.glm <- MSE_Test(X = X %>% select(-Y), y = X$Y, glm_cv = "external",
                        base.learner = "lm", NTest = 100, NTree = 250, B = 1000, var = c( "X2"),
                        p = .85, mtry = 3)

# High subsample size implementation
M_test_high_k <- MSE_Test(Y~., data = X, X.test = X.t1 %>% dplyr::select(-Y), y.test = X.t1$Y,
                   base.learner = "ctree", NTree = 250, B = 1000, var = c( "X2"),
                   p = log(.63*N, base = N))
plot(M_test_high_k)

# We can also get a summary of a particular test object
# > summary(M.noform.glm)
# Call: 
#   MSE_Test.default(X = X %>% select(-Y), y = X$Y, var = c("X2"), 
#                    NTest = 100, B = 1000, NTree = 250, p = 0.85, base.learner = "lm", 
#                    mtry = 3, glm_cv = "external")
# 
# Statistic used: MSE 
# Variables Tested: X2 
# Number of Test Points: 100 
# Base Learner(s): lm 
# P-value: 0.000999001 
# SD Importance: 19.83285


###########################
# Importance Implementation
###########################

pm1 <- permtestImp(X = X %>% dplyr::select(-Y), y = X$Y, single_forest = F,
                   base.learner = "rtree", mtry = 5, NTest = 30, Nbtree = 10, B = 1000,
                   p = .875, verbose = T)


plot(pm1)

pm.null <- permtestImp(X = X %>% dplyr::select(-Y), y = rnorm(N), single_forest = F,
                       base.learner = "rtree", mtry = 3, NTest = 100, Nbtree = 100, B = 1000,
                       p = .95, verbose = T)

plot(pm.null)


### Trying with the single forest option
pm_single_forest <- permtestImp(X = X %>% dplyr::select(-Y), y = X$Y, single_forest = T, keep_forest = T,
                   base.learner = "rtree", mtry = 5, NTest = 30, Nbtree = 100,  verbose = T, ranger = T)
plot(pm_single_forest)

### Holdout Forest example
# Not specifying test points
hf_no_test <- f_holdoutRF(X = X %>% select(-Y), y = X$Y, mintree = 50, max.trees = 500, verbose = T,
                   NTest = 15, ranger = F, keep_forest = F)
plot(hf_no_test)
# Specifying test points
hf_test <- f_holdoutRF(X = X %>% select(-Y), y = X$Y, 
                       X.test = X.t1 %>% select(-Y), y.test = X.t1$Y,
                       mintree = 50, max.trees = 500, verbose = T,
                       ranger = F, keep_forest = F)
plot(hf_test)


###### Comparing two pre-trained models

m_rpart <- bag.s(X = X %>% select(-Y), y = X %>% select(Y),
                       base.learner = "rpart", ntree = 100, k = 100^.85, mtry = 10, form = Y~., ranger = F)
m_ctree <- bag.s(X = X %>% select(-Y), y = X %>% select(Y),
                 base.learner = "ctree", ntree = 100, k = 100^.85, mtry = 10, form = Y~., ranger = F)
m_ctree_ranger <- bag.s(X = X %>% select(-Y), y = X %>% select(Y),
                    base.learner = "ctree", ntree = 100, k = 100^.85, mtry = 10, form = Y~., ranger = T)
m_rtree <- bag.s(X = X %>% select(-Y), y = X %>% select(Y),
                 base.learner = "rtree", ntree = 100, k = 100^.85, mtry = 10, form = Y~., ranger = F)
m_rtree_ranger <- bag.s(X = X %>% select(-Y), y = X %>% select(Y),
                        base.learner = "rtree", ntree = 100, k = 100^.85, mtry = 10, ranger = T)
m_glm  <- bag.s(X = X %>% select(-Y), y = X %>% select(Y),
                                 base.learner = "lm", ntree = 100, k = 100^.85, mtry = 10, ranger = F)


# Comparing these badboys
MSE_comp_rp_ct <- MSE_compare(m_rpart, m_ctree, X.test = X.t1, y.test = X.t1$Y)
MSE_comp_rp_glm <- MSE_compare(m_rpart, m_glm, X.test = X.t1, y.test = X.t1$Y)
MSE_comp_rtr_glm <- MSE_compare(m_rtree_ranger, m_glm, X.test = X.t1)
MSE_comp_rp_ctr <- MSE_compare(m_rpart, m_ctree_ranger, X.test = X.t1, test_stat = "diff")



# Feature importance usage
X_pm <- X[sample(nrow(X), replace = F),which(names(X) != "Y")]
m_reduced <- bag.s(X = X_pm, y = X %>% select(Y),
                  base.learner = "rpart", ntree = 100, k = 100^.85, mtry = 10)
m_full <- bag.s(X = X %>% select(-Y), y = X %>% select(Y),
                base.learner = "rpart", ntree = 100, k = 100^.85, mtry = 10)
full_vs_red <- MSE_compare(m_reduced, m_full, X.test = X.t1, y.test = X.t1$Y)





