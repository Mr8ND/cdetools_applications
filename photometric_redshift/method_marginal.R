library(rhdf5)
library(RFCDE)
library(FlexCoDE)
library(plyr)

datadir <- "data/"
n_grid <- 200
y_grid <- seq(0, 1, length.out = n_grid)

###
name <- "Marginal"
train <- function(x_train, y_train, y_grid) {
  return(y_train)
}
pred <- function(obj, x_test, y_grid) {
  density <- ks::kde(as.numeric(obj), eval.points = y_grid)$estimate
  cde <- matrix(density, nrow(x_test), length(y_grid), byrow = TRUE)
  return(cde)
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