library(rhdf5)
library(NNKCDE)
library(plyr)

datadir <- "data/"
n_grid <- 200
y_grid <- seq(0, 1, length.out = n_grid)

###
methods <- list()

name <- "NNKCDE"
train <- function(x_train, y_train) {
  k <- 5
  h <- 0.01
  return(NNKCDE::NNKCDE$new(x_train, y_train, k = k, h = h))
}
pred <- function(obj, x_test, y_grid) {
  return(obj$predict(x_test, y_grid))
}

###
dataf <- paste0(datadir, "processed.hdf5")
x_train <- h5read(dataf, "/x_train")
y_train <- h5read(dataf, "/y_train")
x_test <- h5read(dataf, "/x_test")
y_test <- h5read(dataf, "/y_test")

obj <- train(x_train, y_train)
cde <- pred(obj, x_test, y_grid)

fname <- paste0(datadir, name, ".hdf5")
if (file.exists(fname)) {
  file.remove(fname)
}

h5createFile(fname)
h5write(y_grid, file = fname, name = "/y_grid")
h5write(y_test, file = fname, name = "/y_true")
h5write(cde, file = fname, name = "/cde")
