library(rhdf5)

datadir <- "Teddy/"

process_data <- function(fname) {
  dat <- read.table(fname, skip = 7, header = FALSE)
  x <- as.matrix(dat[, 8:12])
  y <- as.vector(dat[, 7])
  
  return(list(x = x, y = y))
}

train_data <- process_data(paste0(datadir, "teddy_A"))
x_train <- train_data$x
y_train <- train_data$y

test_data <- process_data(paste0(datadir, "teddy_B"))
x_test <- test_data$x
y_test <- test_data$y

datadir_out <- 'data/'
fname <- paste0(datadir_out, "processed.hdf5")
if (file.exists(fname)) {
  file.remove(fname)
}

h5createFile(fname)
h5write(x_train, file = fname, name = "/x_train")
h5write(y_train, file = fname, name = "/y_train")
h5write(x_test, file = fname, name = "/x_test")
h5write(y_test, file = fname, name = "/y_test")
