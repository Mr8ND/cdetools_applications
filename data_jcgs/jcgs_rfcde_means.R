library(microbenchmark)
library(rhdf5)
library(plyr)
library(FlexCoDE)
library(ranger)

#devtools::install_github("tpospisi/RFCDE/r", force=TRUE)
library(RFCDE)

set.seed(42)

load('data/jcgs_data_rafael.RData')
outdir <- "data/"

## Set methods
methods <- list()

methods[["RFCDE"]] <- function(n_trees, nodesize, n_basis, mtry) {
  return(list(name = "Vector",
              
              "train" = function(x_train, z_train, x_valid, z_valid, z_grid) {
                x_train = rbind(x_train, x_valid)
                z_train = rbind(z_train, z_valid)
                return(RFCDE::RFCDE(x_train = x_train, z_train = z_train,
                                    n_trees = n_trees, node_size = nodesize, mtry = mtry,
                                    n_basis = n_basis, min_loss_delta = 0.0, fit_oob = FALSE))
              },
              
              "predict" = function(obj, x_test, z_grid) {
                return(predict(obj, x_test, z_grid, response='CDE', bandwidth='cv'))
              }))
}


methods[["RFCDE-naive"]] <- function(n_trees, nodesize, n_basis, mtry) {
  return(list(name = "rfcde-naive",
              "train" = function(x_train, z_train, x_valid, z_valid, z_grid) {
                x_train = rbind(x_train, x_valid)
                z_train = rbind(z_train, z_valid)
                return(quantregForest::quantregForest(x_train, z_train, 
                                                      mtry = mtry,
                                                      nodesize = nodesize,
                                                      n_trees = n_trees))
              },
              
              "predict" = function(obj, x_test, z_grid) {
                kde_density <- function(x) {
                  return(ks::kde(x, eval.points = z_grid)$estimate)
                }
                
                return(predict(obj, x_test, what = kde_density))
              }))
}

methods[["RFCDE-mean"]] <- function(n_trees, nodesize, n_basis, mtry) {
  return(list(name = "Vector-mean",
              
              "train" = function(x_train, z_train, x_valid, z_valid, z_grid) {
                x_train = rbind(x_train, x_valid)
                z_train = rbind(z_train, z_valid)
                return(RFCDE::RFCDE(x_train = x_train, z_train = z_train,
                                    n_trees = n_trees, node_size = nodesize, mtry = mtry,
                                    n_basis = n_basis, min_loss_delta = 0.0, fit_oob = FALSE))
              },
              
              "predict" = function(obj, x_test, z_grid) {
                return(predict(obj, x_test, z_grid, response='mean'))
              }))
}

methods[["rf-mean"]] <- function(n_trees, nodesize, n_basis, mtry) {
  return(list(name = "ranger-rf-mean",
              
              "train" = function(x_train, z_train, x_valid, z_valid, z_grid) {
                x_train = rbind(x_train, x_valid)
                z_train = rbind(z_train, z_valid)
                
                train_df <- as.data.frame(cbind(x_train, z_train))
                names(train_df) <- make.names(c(as.character(seq(dim(train_df)[2] - 1)),
                                    "target"))
                
                return(ranger(target  ~ ., data = train_df, 
                                        num.trees = n_trees, mtry = mtry,
                                        min.node.size = nodesize))
              },
              
              "predict" = function(obj, x_test, z_grid) {
                pred_df <- as.data.frame(x_test)
                names(pred_df) <- make.names(as.character(seq(dim(pred_df)[2])))
                return(predict(obj, pred_df)$predictions)
              }))
}

methods[["fRFCDE"]] <- function(n_trees, nodesize, n_basis, mtry) {
  return(list(name = "Functional",
              
              "train" = function(x_train, z_train, x_valid, z_valid, z_grid) {
                x_train = rbind(x_train, x_valid)
                z_train = rbind(z_train, z_valid)
                return(RFCDE::RFCDE(x_train = x_train, z_train = z_train, lens = ncol(x_train),
                                    n_trees = n_trees, node_size = nodesize, mtry = mtry,
                                    n_basis = n_basis, min_loss_delta = 0.0, fit_oob = FALSE, flambda = 50))
              },
              
              "predict" = function(obj, x_test, z_grid) {
                return(predict(obj, x_test, z_grid, response='CDE', bandwidth='cv'))
              }))
}

methods[["Flexcode-Spec"]] <- function(n_basis) {
  return(list(name = "Flexcode-Spec",
              
              "train" = function(x_train, z_train, x_valid, z_valid, z_grid) {
                print('Training Flexcode')
                return(FlexCoDE::fitFlexCoDE(xTrain = x_train, 
                                             zTrain = z_train,
                                             xValidation = x_valid,
                                             zValidation = z_valid,
                                             regressionFunction = regressionFunction.Series,
                                             nIMax = n_basis, system = 'cosine', zMin = min(z_grid),
                                             zMax = max(z_grid)))
              },
              
              "predict" = function(obj, x_test, z_grid) {
                return(predict(obj, x_test, B = length(z_grid))$CDE)
              }))
}



## Functions to train and predict
run_method <- function(method, x_train, z_train, x_valid, z_valid, z_grid, x_test, best_of = 1) {
  train_times <- microbenchmark(times = best_of, {
    method$train(x_train = x_train, z_train = z_train, 
                 x_valid = x_valid, z_valid = z_valid, z_grid = z_grid)
  })$time
  trained <- method$train(x_train = x_train, z_train = z_train, 
                          x_valid = x_valid, z_valid = z_valid, z_grid = z_grid)
  
  predict_times <- microbenchmark(times = best_of, {
    method$predict(trained, x_test, z_grid)
  })$time
  cde <- method$predict(trained, x_test, z_grid)
  
  return(list(cde = cde,
              train_times = train_times,
              predict_times = predict_times))
}

run_simulation <- function(x_train, x_valid, x_test,
                           z_train, z_valid, z_test, outfile, 
                           methods, n_grid = 1000) {
  outfile <- paste0(outdir, outfile)
  if (file.exists(outfile)) { unlink(outfile) }
  
  h5createFile(outfile)
  z_min <- min(z_train)
  z_max <- max(z_train)
  z_grid <- seq(z_min, z_max, length.out = n_grid)
  h5write(z_grid, outfile, "/z_grid")
  h5write(z_test, outfile, "/z_true")
  
  for (method in methods) {
    results <- run_method(method, x_train, z_train, x_valid, z_valid, z_grid, x_test)
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

run_simulation(x_train, x_valid, x_test,
               z_train, z_valid, z_test, "results_rfcde_cv_full_withnaive.hdf5", list(
                methods[["Flexcode-Spec"]](n_basis = n_basis),
                methods[["fRFCDE"]](n_trees = n_trees, nodesize = nodesize,
                                    n_basis = n_basis, mtry = 8),
                methods[["RFCDE"]](n_trees = n_trees, nodesize = nodesize,
                                   n_basis = n_basis, mtry = 58),
                methods[["RFCDE-mean"]](n_trees = n_trees, nodesize = nodesize,
                                   n_basis = n_basis, mtry = 58),
                methods[["rf-mean"]](n_trees = n_trees, nodesize = nodesize,
                                     n_basis = n_basis, mtry = 58),
                methods[["RFCDE-naive"]](n_trees = n_trees, nodesize = nodesize,
                                      n_basis = n_basis, mtry = 58)))
