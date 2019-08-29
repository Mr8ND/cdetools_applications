library(rhdf5)
library(RFCDE)
library(FlexCoDE)
library(plyr)

datadir <- "data/"
n_grid <- 200
y_grid <- seq(0, 1, length.out = n_grid)

###
methods <- list()

name <- "RFCDE-Limited"
train <- function(x_train, y_train, y_grid) {
  n_trees <- 500
  n_basis <- 31
  mtry <- 1
  nodesize <- 20
  return(RFCDE::RFCDE(x_train = x_train, z_train = y_train,
                      n_trees = n_trees, mtry = mtry, node_size = nodesize,
                      n_basis = n_basis, min_loss_delta = 0.0, fit_oob = FALSE))
}
pred <- function(obj, x_test, y_grid) {
  return(predict(obj, x_test, y_grid, response='CDE'))
}

###
dataf <- paste0(datadir, "processed.hdf5")
x_train <- h5read(dataf, "/x_train")
y_train <- h5read(dataf, "/y_train")
x_test <- h5read(dataf, "/x_test")
y_test <- h5read(dataf, "/y_test")

# Reduce the number of variables
x_train <- x_train[, 1:3]
x_test <- x_test[, 1:3]

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
