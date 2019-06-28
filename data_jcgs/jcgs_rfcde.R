library(microbenchmark)
library(rhdf5)
library(plyr)
library(FlexCoDE)

devtools::install_github("tpospisi/RFCDE/r", force=TRUE)
library(RFCDE)

set.seed(42)

load('data/jcgs_data.RData')

## train/test split
n_train_gals <- 2000
train_ids <- sample(dim(z)[1], size = n_train_gals)
x <- x_mat
z <- z
outdir <- "data/"

## Set methods
methods <- list()

methods[["RFCDE"]] <- function(n_trees, nodesize, n_basis, mtry) {
  return(list(name = "Vector",
              
              "train" = function(x_train, z_train, z_grid) {
                return(RFCDE::RFCDE(x_train = x_train, z_train = z_train,
                                    n_trees = n_trees, node_size = nodesize, mtry = mtry,
                                    n_basis = n_basis, min_loss_delta = 0.0, fit_oob = FALSE))
              },
              
              "predict" = function(obj, x_test, z_grid) {
                return(predict(obj, x_test, z_grid, response='CDE'))
              }))
}

methods[["fRFCDE"]] <- function(n_trees, nodesize, n_basis, mtry) {
  return(list(name = "Functional",
              
              "train" = function(x_train, z_train, z_grid) {
                return(RFCDE::RFCDE(x_train = x_train, z_train = z_train, lens = ncol(x_train),
                                     n_trees = n_trees, node_size = nodesize, mtry = mtry,
                                     n_basis = n_basis, min_loss_delta = 0.0, fit_oob = FALSE, flambda = 50))
              },
              
              "predict" = function(obj, x_test, z_grid) {
                return(predict(obj, x_test, z_grid, response='CDE'))
              }))
}

methods[["Flexcode-Spec"]] <- function(n_basis) {
  return(list(name = "Flexcode-Spec",
              
              "train" = function(x_train, z_train, z_grid, valid_size=400) {
                
                perm <- sample(nrow(x_train))
                train_ids <- perm[1:(nrow(x_train) - valid_size)]
                valid_ids <- perm[-(1:(nrow(x_train) - valid_size))][1:valid_size]
                
                return(FlexCoDE::fitFlexCoDE(xTrain = x_train[train_ids, ], 
                                             zTrain = z_train[train_ids, , drop = FALSE],
                                             xValidation = x_train[valid_ids, ],
                                             zValidation = z_train[valid_ids, , drop = FALSE],
                                            regressionFunction = regressionFunction.Series,
                                            nIMax = n_basis, system = 'cosine', zMin = min(z_grid),
                                            zMax = max(z_grid)))
              },
              
              "predict" = function(obj, x_test, z_grid) {
                return(predict(obj, x_test, B = length(z_grid))$CDE)
              }))
}



## Functions to train and predict
run_method <- function(method, x_train, z_train, z_grid, x_test, best_of = 1) {
  train_times <- microbenchmark(times = best_of, {
    method$train(x_train = x_train, z_train = z_train, z_grid = z_grid)
  })$time
  trained <- method$train(x_train, z_train, z_grid)
  
  predict_times <- microbenchmark(times = best_of, {
    method$predict(trained, x_test, z_grid)
  })$time
  cde <- method$predict(trained, x_test, z_grid)
  
  return(list(cde = cde,
              train_times = train_times,
              predict_times = predict_times))
}

run_simulation <- function(n_train, n_test, outfile, methods, n_grid = 1000) {
  outfile <- paste0(outdir, outfile)
  if (file.exists(outfile)) { unlink(outfile) }
  h5createFile(outfile)
  
  perm <- sample(nrow(x))
  train_ids <- perm[1:n_train]
  test_ids <- perm[-(1:n_train)][1:n_test]
  
  x_train <- x[train_ids, ]
  z_train <- z[train_ids, , drop = FALSE]
  x_test <- x[test_ids, ]
  z_test <- z[test_ids, , drop = FALSE]
  h5write(z_test, outfile, "/z_true")
  
  z_min <- min(z_train)
  z_max <- max(z_train)
  z_grid <- seq(z_min, z_max, length.out = n_grid)
  h5write(z_grid, outfile, "/z_grid")
  
  for (method in methods) {
    results <- run_method(method, x_train, z_train, z_grid, x_test)
    ## Save to hdf5 file
    prefix <- paste0("/", method$name)
    h5createGroup(outfile, prefix)
    h5write(results$cde, outfile, paste0(prefix, "/cde"))
    h5write(results$train_times, outfile, paste0(prefix, "/train_times"))
    h5write(results$predict_times, outfile, paste0(prefix, "/predict_times"))
  }
}

## Parameterization
n_trees <- 1000
nodesize <- 20
n_basis <- 31
n_test_gal <- 800

run_simulation(n_train_gals, n_test_gal, "results_with_flexcode.hdf5", list(
  methods[["Flexcode-Spec"]](n_basis = n_basis),
  methods[["RFCDE"]](n_trees = n_trees, nodesize = nodesize,
                     n_basis = n_basis, mtry = 58),
  methods[["fRFCDE"]](n_trees = n_trees, nodesize = nodesize,
                     n_basis = n_basis, mtry = 8)))
