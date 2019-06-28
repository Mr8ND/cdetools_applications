library(rhdf5)
library(tidyverse)
library(readr)

x <- readr::read_delim(file = 'data/covariatesSpectra.txt', col_names = FALSE, delim = ' ')
x <- dplyr::mutate_all(x, function(x) as.numeric(x))
x_mat <- as.matrix(x)

z <- readr::read_delim(file = 'data/zSpectra.txt', col_names = FALSE, delim = ' ')
z <- as.matrix(z)

save(x_mat, z, file = 'data/jcgs_data.RData')
