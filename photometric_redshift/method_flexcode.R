library(rhdf5)
library(RFCDE)
library(FlexCoDE)
library(plyr)

datadir <- "data/"
n_grid <- 200
y_grid <- seq(0, 1, length.out = n_grid)

###
methods <- list()

name <- "FlexZBoost"
train <-  function(x_train, y_train, y_grid) {
  n_basis <- 31
  system <- "cosine"
  obj <- fitFlexCoDE(xTrain = x_train, zTrain = y_train,
                     xValidation = x_train, zValidation = y_train,
                     regressionFunction = regressionFunction.XGBoost,
                     nIMax = n_basis, system = system, zMin = min(y_grid),
                     zMax = max(y_grid), chooseSharpen = TRUE)
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
x_test <- h5read(dataf, "/x_test")
y_test <- h5read(dataf, "/y_test")

obj <- train(x_train, y_train, y_grid)
cde <- pred(obj, x_test, y_grid)

fname <- paste0(datadir, name, ".hdf5")
if (file.exists(fname)) {
  file.remove(fname)
}

h5createFile(fname)
h5write(y_grid, file = fname, name = "/y_grid")
h5write(y_test, file = fname, name = "/y_true")
h5write(cde, file = fname, name = "/cde")