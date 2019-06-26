library(rhdf5)

flux_mat <- NULL
z_vec <- NULL
data_dir = "sdss/data/"
for (filename in list.files(data_dir)){
  load(paste0(data_dir, filename))

  if (is.null(flux_mat)){
    flux_mat <- flux
  } else {
    flux_mat <- rbind(flux_mat, flux)
  }
  print(dim(flux_mat))
  
  if (is.null(z_vec)){
    z_vec <- as.matrix(z)
  } else {
    z_vec <- rbind(z_vec, as.matrix(z))
  }
  print(dim(z_vec))
}
save(flux_mat, z_vec, file='sdss/data/SDSS_galaxies_full_2018.Rdata')


## Creating file for DeepCDE ##################################################
n_train_gals <- 20000
n_test_gal <- 5000
datadir_out <- 'sdss/data/'
fname <- paste0(datadir_out, "sdss_galaxies_full.hdf5")
if (file.exists(fname)) {
  file.remove(fname)
}

perm <- sample(nrow(flux_mat))
train_ids <- perm[1:n_train_gals]
test_ids <- perm[-(1:n_train_gals)][1:n_test_gal]

x_train <- flux_mat[train_ids, ]
z_train <- z_vec[train_ids, , drop = FALSE]
x_test <- flux_mat[test_ids, ]
z_test <- z_vec[test_ids, , drop = FALSE]

h5createFile(fname)
h5write(x_train, file = fname, name = "/x_train")
h5write(z_train, file = fname, name = "/y_train")
h5write(x_test, file = fname, name = "/x_test")
h5write(z_test, file = fname, name = "/y_test")

