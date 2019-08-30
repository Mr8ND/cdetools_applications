library(rhdf5)
library(NNKCDE)
library(plyr)
library(pracma)
source('auxiliary_funcs.R')


##### Setup Variables for Predicition
datadir <- "data/"
n_grid <- 200
y_grid <- seq(0, 1, length.out = n_grid)

##### Setup Method
methods <- list()

name <- "NNKCDE"
train <- function(x_train, y_train, x_valid, y_valid) {
  k_grid <- c(1,2,5,10,15,20)
  h_grid <- c(10**seq(-4,0))
  
  obj <- NNKCDE::NNKCDE$new(x_train, y_train)
  obj$tune(x_valid, y_valid, k_grid, h_grid)
  return(obj)
}
pred <- function(obj, x_test, y_grid) {
  return(obj$predict(x_test, y_grid))
}

#### Setup output file and run method
dataf <- paste0(datadir, "processed.hdf5")
x_train <- h5read(dataf, "/x_train")
y_train <- h5read(dataf, "/y_train")

# Training and Validation Creation
n_training <- as.integer(dim(x_train)[1] * 0.8)
n_testing <- dim(x_train)[1] - n_training
perm <- sample(dim(x_train)[1])
train_ids <- perm[1:n_training]
test_ids <- perm[-(1:n_training)][1:n_testing]

x_train_nn <- x_train[train_ids, ]
y_train_nn <- y_train[train_ids]
x_valid <- x_train[test_ids, ]
y_valid <- y_train[test_ids]

x_test <- h5read(dataf, "/x_test")
y_test <- h5read(dataf, "/y_test")

obj <- train(x_train_nn, y_train_nn, x_valid, y_valid)
print(c(obj$h, obj$k))
cde <- pred(obj, x_test, y_grid)

if (pracma::trapz(x = y_grid, y = cde[1,]) > 1.1){
  cde <- t(apply(cde, MARGIN = 1, 
                 FUN = function(x) normalize_density(y_grid[2] - y_grid[1], x)))
}

fname <- paste0(datadir, name, ".hdf5")
if (file.exists(fname)) {
  file.remove(fname)
}

h5createFile(fname)
h5write(y_grid, file = fname, name = "/y_grid")
h5write(y_test, file = fname, name = "/y_true")
h5write(cde, file = fname, name = "/cde")
