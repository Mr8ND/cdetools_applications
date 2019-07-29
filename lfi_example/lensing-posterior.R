library(fields)
library(ggplot2)
library(plyr)
library(rhdf5)
library(latex2exp)
library(ascii)
library(NNKCDE)
library(RFCDE)
library(cdetools)

data_file <- "data/data.csv"
dat <- read.csv(data_file)
n_draw <- nrow(dat) - 1
n_train <- 1000
n_validation <- 1000
n_grid <- 30

accepts <- c(0.2, 0.5, 1.0)

set.seed(42)

## Read in data
x_obs <- as.vector(dat[1, paste0("X", 1:20)])
z_obs <- dat[1, c("omega_m", "sigma_8")]
dat <- dat[-1, ]

constants <- dat[, c("omega_m", "sigma_8")]
xipm <- as.matrix(dat[, paste0("X", 1:20)])

df <- ldply(accepts, function(accept) {
    ## Take abc samples
    abc_ids <- sample(n_draw, (n_train + n_validation) / accept)
    abc_xipm <- xipm[abc_ids, ]
    abc_constants <- constants[abc_ids, ]

    d <- as.vector(rdist(abc_xipm, x_obs))
    cutoff <- sort(d)[n_train + n_validation]
    close_xipm <- abc_xipm[d <= cutoff, ]
    close_constants <- abc_constants[d <= cutoff, ]

    ## Train/validation split
    train_ids <- sample(nrow(close_xipm), n_train)
    x_train <- as.matrix(close_xipm[train_ids, ])
    z_train <- as.matrix(close_constants[train_ids, ])
    x_validation <- as.matrix(close_xipm[-train_ids, ])
    z_validation <- as.matrix(close_constants[-train_ids, ])

    z_grid <- expand.grid(x = seq(0.1, 0.8, length.out = n_grid),
                          y = seq(0.5, 1.0, length.out = n_grid))

    ## ABC estimate
    H <- ks::Hpi(z_train)
    abc_obj <- NNKCDE$new(x_train, z_train, k = n_train, h = H)
    abc_density <- as.vector(abc_obj$predict(x_obs, z_grid))

    ## Kernel-NN estimate
    k_grid <- seq(100, n_train, 100)
    knn_obj <- NNKCDE$new(x_train, z_train, k = n_train, h = H)
    knn_obj$tune(x_validation, z_validation, k_grid = k_grid)
    knn_density <- as.vector(knn_obj$predict(x_obs, z_grid))

    ## RFCDE estimate
    rfcde_obj <- RFCDE(x_train, z_train, n_basis=15, n_trees=100)
    rfcde_density <- as.vector(predict(rfcde_obj, as.matrix(x_obs), as.matrix(z_grid),
                                       bandwidth = diag(0.0005, 2),
                                       response = 'CDE'))


    get_density <- function(name, density) {
        density <- density / sum(density)
        return(data.frame(z_grid, name = name, density = density, 
                          accept = accept))
    }
    
    return(rbind(get_density("ABC", abc_density),
                 get_density("NNKCDE", knn_density),
                 get_density("RFCDE", rfcde_density)))
})

alpha <- .55
accept <- min(accepts)

d <- as.vector(rdist(xipm, x_obs))
cutoff <- sort(d)[n_train + n_validation]
close_xipm <- xipm[d <= cutoff, ]
close_constants <- constants[d <= cutoff, ]

mod <- lm(log(sigma_8)  ~ log(omega_m), data = as.data.frame(close_constants))
degen <- exp(coef(mod)[1])
alpha <- -coef(mod)[2]
df$accept_plot <- factor(paste0((df$accept)*100, '%'), 
                         levels = c('20%', '50%', '100%'))

fig <- ggplot(df, aes(x = x, y = y, z = density, color = name)) +
    geom_contour(aes(alpha = ..level..), lwd = 1) +
    stat_function(fun = function(x) degen * x ^ -alpha, lty = 2) +
    xlim(0.1, 0.8) + ylim(0.5, 1.0) +
    facet_grid(name ~ accept_plot) +
    xlab(TeX("$\\Omega_{M}$")) + ylab(TeX("$\\sigma_{8}$")) +
    guides(color = FALSE, alpha = FALSE) +
    scale_color_brewer(palette = "Dark2") +
    theme_minimal(base_size = 16)

fig +  theme(plot.title = element_blank(),
             axis.title.x = element_text(size=22),
             axis.title.y = element_text(size=22),
             axis.text.x = element_text(size=16,margin=ggplot2::margin(b=5)),
             axis.text.y = element_text(size=16,margin=ggplot2::margin(l=3)),
             strip.text.x = element_text(size=20),
             strip.text.y = element_text(size=20)) + coord_fixed()

## ggsave("../figures/lensing-posterior.pdf", fig, device = cairo_pdf)


### HPD COVERAGE


df_val <- ldply(accepts, function(accept) {
  ## Take abc samples
  abc_ids <- sample(n_draw, (n_train + n_validation) / accept)
  abc_xipm <- xipm[abc_ids, ]
  abc_constants <- constants[abc_ids, ]
  
  d <- as.vector(rdist(abc_xipm, x_obs))
  cutoff <- sort(d)[n_train + n_validation]
  close_xipm <- abc_xipm[d <= cutoff, ]
  close_constants <- abc_constants[d <= cutoff, ]
  
  ## Train/validation split
  train_ids <- sample(nrow(close_xipm), n_train)
  x_train <- as.matrix(close_xipm[train_ids, ])
  z_train <- as.matrix(close_constants[train_ids, ])
  x_validation <- as.matrix(close_xipm[-train_ids, ])
  z_validation <- as.matrix(close_constants[-train_ids, ])
  
  z_grid <- expand.grid(x = seq(0.1, 0.8, length.out = n_grid),
                        y = seq(0.5, 1.0, length.out = n_grid))
  
  ## ABC estimate
  H <- ks::Hpi(z_train)
  abc_obj <- NNKCDE$new(x_train, z_train, k = n_train, h = H)
  abc_density <- abc_obj$predict(x_validation, z_grid)
  
  ## Kernel-NN estimate
  k_grid <- seq(100, n_train, 100)
  knn_obj <- NNKCDE$new(x_train, z_train, k = n_train, h = H)
  knn_obj$tune(x_validation, z_validation, k_grid = k_grid)
  knn_density <- knn_obj$predict(x_validation, z_grid)
  
  ## RFCDE estimate
  rfcde_obj <- RFCDE(x_train, z_train, n_basis=15, n_trees=100)
  rfcde_density <- predict(rfcde_obj, as.matrix(x_validation), 
                                     as.matrix(z_grid),
                                     bandwidth = diag(0.0005, 2),
                                     response = 'CDE')
  
  
  get_coverage <- function(name, density) {
    density <- density / sum(density)
    p_vals <- cdetools::hpd_coverage(as.matrix(density), z_grid, 
                                     z_validation)
    return(data.frame(name = name, p_vals = as.vector(p_vals), 
                      accept = accept))
  }
  
  return(rbind(get_coverage("ABC", abc_density),
               get_coverage("NNKCDE", knn_density),
               get_coverage("RFCDE", rfcde_density)))
})

df_val$vals <- df_val$p_vals/max(df_val$p_vals)
df_val$accept_plot <- factor(paste0((df_val$accept)*100, '%'), 
                             levels = c('20%', '50%', '100%'))

p1_plot <- ggplot(df_val, aes(x = vals, y = ..density..)) +
  geom_histogram() +
  facet_grid(name ~ accept_plot) +
  xlab("HPD Value") + ylab("Density") +
  #xlim(c(0,1)) +
  theme_minimal()

p1_plot + theme(plot.title = element_text(hjust = 0.5, size=26),
                axis.title.x = element_text(size=22),
                axis.text.x = element_text(size=14),
                axis.text.y = element_text(size=14),
                axis.title.y = element_text(size=22),
                strip.text.x = element_text(size = 20),
                strip.text.y = element_text(size = 20))

library(gtable)
library(gridExtra)
library(grid)
library(ggpubr)
library(pracma)


final_plot <- ggarrange(fig +  theme(plot.title = element_blank(),
                       axis.title.x = element_text(size=22),
                       axis.title.y = element_text(size=22),
                       axis.text.x = element_text(size=16,margin=ggplot2::margin(b=5)),
                       axis.text.y = element_text(size=16,margin=ggplot2::margin(l=3)),
                       strip.text.x = element_text(size=20),
                       strip.text.y = element_blank()),
          p1_plot + theme(plot.title = element_text(hjust = 0.5, size=26),
                          axis.title.x = element_text(size=22),
                          axis.text.x = element_text(size=16,margin=ggplot2::margin(b=5)),
                          axis.text.y = element_text(size=16),
                          axis.title.y = element_text(size=22),
                          strip.text.x = element_text(size = 20),
                          strip.text.y = element_text(size = 20)),
  ncol = 2, nrow = 1
)
ggsave("images/joint_lensing_pval.pdf", final_plot, device = cairo_pdf,
       height=10, width=20)


### CDE LOSS

df_loss <- ldply(accepts, function(accept) {
  ## Take abc samples
  abc_ids <- sample(n_draw, (n_train + n_validation) / accept)
  abc_xipm <- xipm[abc_ids, ]
  abc_constants <- constants[abc_ids, ]
  
  d <- as.vector(rdist(abc_xipm, x_obs))
  cutoff <- sort(d)[n_train + n_validation]
  close_xipm <- abc_xipm[d <= cutoff, ]
  close_constants <- abc_constants[d <= cutoff, ]
  
  ## Train/validation split
  train_ids <- sample(nrow(close_xipm), n_train)
  x_train <- as.matrix(close_xipm[train_ids, ])
  z_train <- as.matrix(close_constants[train_ids, ])
  x_validation <- as.matrix(close_xipm[-train_ids, ])
  z_validation <- as.matrix(close_constants[-train_ids, ])
  
  z_grid <- expand.grid(x = seq(0.1, 0.8, length.out = n_grid),
                        y = seq(0.5, 1.0, length.out = n_grid))
  
  ## ABC estimate
  H <- ks::Hpi(z_train)
  abc_obj <- NNKCDE$new(x_train, z_train, k = n_train, h = H)
  abc_density <- abc_obj$predict(x_validation, z_grid)
  
  ## Kernel-NN estimate
  k_grid <- seq(100, n_train, 100)
  knn_obj <- NNKCDE$new(x_train, z_train, k = n_train, h = H)
  knn_obj$tune(x_validation, z_validation, k_grid = k_grid)
  knn_density <- knn_obj$predict(x_validation, z_grid)
  
  ## RFCDE estimate
  rfcde_obj <- RFCDE(x_train, z_train, n_basis=15, n_trees=100)
  rfcde_density <- predict(rfcde_obj, as.matrix(x_validation), 
                           as.matrix(z_grid),
                           bandwidth = diag(0.0005, 2),
                           response = 'CDE')
  
  
  get_coverage <- function(name, density) {
    density <- density / sum(density)
    loss <- cdetools::cde_loss(as.matrix(density), z_grid, 
                z_validation)
    return(data.frame(name = name, loss = loss, 
                      accept = accept))
  }
  
  return(rbind(get_coverage("ABC", abc_density),
               get_coverage("NNKCDE", knn_density),
               get_coverage("RFCDE", rfcde_density)))
})

df_table <- data.table::copy(df_loss)
mult <- 1e5
df_table['CDE Loss (SE)'] <- paste0(round(df_table$loss.loss*mult,3), " (", 
                                    round(df_table$loss.se*mult,3), ")")

tbl <- df_table[, c("name", "accept", "CDE Loss (SE)")]

print(ascii(tbl, include.rownames = FALSE), type = "org")


