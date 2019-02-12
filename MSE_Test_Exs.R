
## Local Loading
source('C:/Users/drain/Box Sync/Mentch Coleman Shared/RF Permutation Tests/R Stuff/MSE_Test_File.R')
library(dplyr)
set.seed(1500)

## Remote Loading
library(devtools)
source_url("https://raw.githubusercontent.com/tim-coleman/SURFTest/master/MSE_Test_File.R")

## Generating dummy data

N <- 1250
Nvar <- 10
N_test <- 150

# training data:
X <- data.frame(replicate(Nvar, runif(N))) %>%
  mutate(Y = 5*(X3) + .5*X2^2 + ifelse(X6 > 10*X1*X8*X9, 1, 0) +  rnorm(N, sd = .05))

# some testing data:
X.t1 <- data.frame(replicate(Nvar, runif(N_test))) %>%
  mutate(Y = 5*(X3) + .5*X2^2 + ifelse(X6 > 10*X1*X8*X9, 1, 0) +  rnorm(N_test, sd = .05))

## Trying each base learner
b.rpart <- bag.s(X = X %>% select(-Y), y = X %>% select(Y),
                 base.learner = "rpart", ntree = 10, k = 100^.85, mtry = 10, form = Y~.)
b.ctree <- bag.s(X = X %>% select(-Y), y = X %>% select(Y),
                 base.learner = "ctree", ntree = 10, k = 100^.95, mtry = 2, form = Y~.)
b.rf <- bag.s(X = X %>% select(-Y), y = X$Y,
              base.learner = "rtree", ntree = 10, k = 100^.95, mtry = 2, Y~., ranger = F)
b.glmnet <-  bag.s(X = X %>% select(-Y), y = X$Y,
                   base.learner = "lm", ntree = 10, k = 100^.95, mtry = 2, lambda = .5)

## Generating an error
bag.s(X = X %>% select(-Y), y = X$Y,
      base.learner = "error_generator", ntree = 1, k = 100^.95, mtry = 2, Y~.)


### Trying out the testing method
# Not specifying test points:
M_no_test <- MSE_Test(X = X %>% dplyr::select(-Y), y = X$Y, resp = "y",
                      base.learner = "rpart", NTest = 100, NTree = 150, B = 1000, var = c( "X3"),
                      p = .85)

plot(M_no_test)

# Specifying test points:
M_test <- MSE_Test(X = X %>% select(-Y), y = X$Y, resp = "y", X.test = X.t1 %>% select(-Y), y.test = X.t1$Y,
                      base.learner = "rtree", NTree = 250, B = 1000, var = c( "X2"),
                      p = .85)
plot(M_test)

### Formula implementation
M.formula <- MSE_Test(Y~., data = X, 
                      base.learner = "rtree", NTest = 100, NTree = 250, B = 1000, var = c( "X2"),
                      p = .85, mtry = 3, ranger = T)
plot(M.formula)

# GLM implementation
M.form.glm <- MSE_Test(Y~., data = X, #resp = "Y",
                       base.learner = "lm", NTest = 100, NTree = 250, B = 1000, var = c( "X2"),
                       p = .85, mtry = 3)

### 8-Tree Implementation


start.time <- Sys.time()
M.8tree <- MSE_Test(Y~., data = X, resp = "y",
                    base.learner = "rtree", NTest = 1000, NTree = 8, B = 1000, var = c("X2"),
                    p = .95, mtry = 1000)
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)
plot(M.8tree)


### High Dimensional Regression Implementation
N <- 1500
Nvar <- 500
n_non_zero <- 5
N_test <- 150

# training data:
coefs <- c(rep(1, n_non_zero), rep(0, Nvar - n_non_zero))
Xlm <- data.frame(replicate(Nvar, runif(N)))
ylm <- as.matrix(Xlm)%*%coefs + rnorm(N, sd = 1)


# some testing data:
Xlm.t1 <- data.frame(replicate(Nvar, runif(N_test)))
ylm.t1 <- as.matrix(Xlm.t1)%*%coefs + rnorm(N_test, sd = 1)

start.time <- Sys.time()
MSE_high_dim <- MSE_Test(X = Xlm, y = ylm,  X.test = Xlm.t1, y.test = ylm.t1, mtry = Nvar^.75,
                         base.learner = "rtree", NTree = 500, B = 1000, var = "X1", resp = "y",
                         p = .85, lambda = 1.5, ranger = T)
time_taken_ranger <- Sys.time() - start.time
start.time <- Sys.time()
MSE_high_dim <- MSE_Test(X = Xlm, y = ylm,  X.test = Xlm.t1, y.test = ylm.t1, mtry = Nvar^.75,
                         base.learner = "rtree", NTree = 500, B = 1000, var = "X1", resp = "y",
                         p = .85, lambda = 1.5, ranger = F)
time_taken_rf <- Sys.time() - start.time
plot(MSE_high_dim)

###########################
# Importance Implementation
###########################

pm1 <- permtestImp(X = X %>% dplyr::select(-Y), y = X$Y, resp = "y", single_forest = F,
                   base.learner = "rtree", mtry = 5, NTest = 30, Nbtree = 10, B = 1000,
                   p = .875, verbose = T)


plot(pm1)

pm.null <- permtestImp(X = X %>% dplyr::select(-Y), y = rnorm(2500), resp = "y", single_forest = F,
                       base.learner = "rtree", mtry = 3, NTest = 100, Nbtree = 100, B = 1000,
                       p = .95, verbose = T)

plot(pm.null)

pm.highdim <- permtestImp(X = Xlm, y = ylm,  X.test = Xlm.t1, y.test = ylm.t1, mtry = Nvar^.5, single_forest = F,
                          base.learner = "rtree", Nbtree = 100, B = 1000, p = .85, lambda = seq(2.5, 0, length.out = 15),
                          alpha = .5, verbose = T)
plot(pm.highdim)


### Trying with the single forest option
pm1 <- permtestImp(X = X %>% dplyr::select(-Y), y = X$Y,single_forest = T, keep_forest = T,
                   base.learner = "rtree", mtry = 5, NTest = 30, Nbtree = 100,  verbose = T, ranger = T)


### Holdout Forest example
hf1 <- f_holdoutRF(X = X %>% select(-Y), y = X$Y, mintree = 50, max.trees = 500, NTest = 15, ranger = F, keep_forest = T)

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
                        base.learner = "rtree", ntree = 100, k = 100^.85, mtry = 10, form = Y~., ranger = T)
m_glm  <- bag.s(X = X %>% select(-Y), y = X %>% select(Y),
                                 base.learner = "lm", ntree = 100, k = 100^.85, mtry = 10, form = Y~., ranger = F)

p1 <- data.frame(lapply(m_ctree_ranger, FUN = function(x) predict(x, data = X.t1)["predictions"]))
p2 <- data.frame(lapply(m_rpart, FUN = function(x) predict(x, newdata = X.t1)))
p_glm <- data.frame(lapply(m_glm, predict, 
                           newx = model.matrix(y~., data = data.frame(X.t1, "y" = rep(0, N_test)))))
source('C:/Users/drain/Box/Mentch Coleman Shared/RF Permutation Tests/R Stuff/MSE_Test_File.R')


MSE_comp_rp_ct <- MSE_compare(m_rpart, m_ctree, X.test = X.t1, y.test = X.t1$Y)
MSE_comp_rp_glm <- MSE_compare(m_rpart, m_glm, X.test = X.t1)
MSE_comp_rtr_glm <- MSE_compare(m_rtree_ranger, m_glm, X.test = X.t1)
MSE_comp_rp_ctr <- MSE_compare(m_rpart, m_ctree_ranger, X.test = X.t1)



