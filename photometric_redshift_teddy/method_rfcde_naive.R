library(rhdf5)
library(quantregForest)
library(FlexCoDE)
library(plyr)

#devtools::install_github("tpospisi/RFCDE/r", force=TRUE)
library(RFCDE)

##### Setup Variables for Predicition
datadir <- "data/"
n_grid <- 200
y_grid <- seq(0, 1, length.out = n_grid)

##### Setup Method
methods <- list()

name <- "RFCDE-naive-hscv"
train <- function(x_train, y_train, y_grid) {
  n_trees <- 500
  mtry <- 2
  nodesize <- 20
  return(quantregForest::quantregForest(x_train, y_train, 
                                        mtry = mtry,
                                        nodesize = nodesize,
                                        n_trees = n_trees))
}
pred <- function(obj, x_test, y_grid) {
  kde_density <- function(x) {
    h_opt <- ks::Hscv(x)
    print(h_opt)
    return(ks::kde(x, eval.points = y_grid, H=h_opt)$estimate)
  }
  
  return(predict(obj, x_test, what = kde_density))
}

#### Setup output file and run method
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
