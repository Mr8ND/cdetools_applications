library(rhdf5)
library(RFCDE)
library(FlexCoDE)
library(plyr)

datadir <- "data/"
n_grid <- 200
y_grid <- seq(0, 1, length.out = n_grid)

###
methods <- list()

name <- "FlexZBoost_basis_search"
train <-  function(x_train, y_train, x_val, y_val, y_grid) {
  n_basis <- 40
  system <- "cosine"
  obj <- fitFlexCoDE(xTrain = x_train, zTrain = y_train,
                     xValidation = x_valid, zValidation = y_valid,
                     regressionFunction = regressionFunction.XGBoost,
                     nIMax = n_basis, system = system, zMin = min(y_grid),
                     zMax = max(y_grid))
  return(obj)
}

pred <- function(obj, x_test, y_grid) {
  n_grid <- length(y_grid)
  return(predict(obj, x_test, B = n_grid)$CDE)
}

###
dataf <- paste0(datadir, "processed.hdf5")
x_train <- h5read(dataf, "/x_train")
y_train <- h5read(dataf, "/y_train")

# Training and Validation Creation
n_training <- as.integer(dim(x_train)[1] * 0.8)
n_testing <- dim(x_train)[1] - n_training
perm <- sample(dim(x_train)[1])
train_ids <- perm[1:n_training]
test_ids <- perm[-(1:n_training)][1:n_testing]

x_train_fl <- x_train[train_ids, ]
y_train_fl <- y_train[train_ids]
x_valid <- x_train[test_ids, ]
y_valid <- y_train[test_ids]

x_test <- h5read(dataf, "/x_test")
y_test <- h5read(dataf, "/y_test")

obj <- train(x_train_fl, y_train_fl, x_valid, y_valid, y_grid)
cde <- pred(obj, x_test, y_grid)

fname <- paste0(datadir, name, ".hdf5")
if (file.exists(fname)) {
  file.remove(fname)
}

h5createFile(fname)
h5write(y_grid, file = fname, name = "/y_grid")
h5write(obj$bestI, file = fname, name = "/best_base")
h5write(y_test, file = fname, name = "/y_true")
h5write(cde, file = fname, name = "/cde")