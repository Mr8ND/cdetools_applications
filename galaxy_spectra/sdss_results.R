library(ascii)
library(cdetools)
library(plyr)
library(reshape2)
library(rhdf5)
library(tidyverse)
library(gtable)
library(gridExtra)
library(grid)
library(ggpubr)
library(reshape2)

outdir <- "data/"

grade_simulation <- function(resultsfile) {
  z_test <- h5read(resultsfile, "/z_true")
  z_grid <- h5read(resultsfile, "/z_grid")
  
  methods <- h5ls(resultsfile, recursive = FALSE)$name
  methods <- methods[!methods %in% c("z_grid", "z_true")]
  
  return(ldply(methods, function(method) {
    preds <- h5read(resultsfile, paste0(method, "/cde"))
    train_times <- h5read(resultsfile, paste0(method, "/train_times"))
    predict_times <- h5read(resultsfile, paste0(method, "/predict_times"))
    
    return(data.frame(loss = cdetools::cde_loss(preds, z_grid, z_test),
                      method = method,
                      train_time = min(train_times),
                      predict_time = min(predict_times)))
  }))
}

## Grade individual simulations
fname <- paste0(outdir, "results_full.hdf5")
df <- grade_simulation(fname)

df$Method <- df$method

df$loss.loss = round(df$loss.loss, digits = 3)
df$loss.se = round(df$loss.se, digits = 3)

df[, "Train Time (sec)"] <- df$train_time / 1e9
df[, "Predict Time (sec)"] <- df$predict_time / 1e9
df[, "CDE Loss (SE)"] <- paste0(df$loss.loss, " (", df$loss.se, ")")

tbl <- df[, c("Method", "Train Time (sec)", "CDE Loss (SE)")]

print(ascii(tbl, include.rownames = FALSE), type = "org")


# CDE Examples ############################################################

n_plots <- 10
starting_idx <- 10

get_densities <- function(resultsfile) {
  z_test <- h5read(resultsfile, "/z_true")
  z_grid <- h5read(resultsfile, "/z_grid")
  
  methods <- h5ls(resultsfile, recursive = FALSE)$name
  methods <- methods[!methods %in% c("z_grid", "z_true")]
  
  return(ldply(methods, function(method) {
    cde <- h5read(resultsfile, paste0(method, "/cde"))
    
    return(data.frame(Method = method,
                      density = as.vector(cde[c(starting_idx:(starting_idx+n_plots-1)), ]),
                      y = rep(z_grid, each = n_plots),
                      id = rep(seq_len(n_plots), times = length(z_grid))))
  }))
}


fname <- paste0(outdir, "results_full.hdf5")
densities <- get_densities(fname)
z_test <- h5read(fname, "/z_true")
truth  <-  data.frame(z_test = z_test[c(starting_idx:(starting_idx+n_plots-1))],
                      id = seq_len(n_plots))

densities$Method <- factor(densities$Method, levels = c('Functional', 'Vector'))

fig <- ggplot(NULL, aes(x = y, y = density)) +
  geom_line(data = densities) +
  geom_vline(data = truth, aes(xintercept = z_test), lty = 2) +
  facet_grid(Method ~ id) +
  xlab("z") + ylab("p(z | X)") +
  theme_minimal()
ggsave(filename='images/full/cdes_examples.pdf', plot=fig, dpi=300)
print(fig)



# Redshift Prediction function Examples ############################################################

get_predictions <- function(resultsfile) {
  z_test <- h5read(resultsfile, "/z_true")
  z_grid <- h5read(resultsfile, "/z_grid")
  
  methods <- h5ls(resultsfile, recursive = FALSE)$name
  methods <- methods[!methods %in% c("z_grid", "z_true")]
  
  return(ldply(methods, function(method) {
    cde <- h5read(resultsfile, paste0(method, "/cde"))
    
    return(data.frame(Method = method,
                      prediction = apply(cde, MARGIN = 1, FUN = function(x) z_grid[which.max(x)]),
                      z_true = z_test))
  }))
}

fname <- paste0(outdir, "results_full.hdf5")
df_pred <- get_predictions(fname)
fig_points <- ggplot(data=df_pred,aes(z_true,prediction))+
  geom_point(alpha=1) +  facet_grid(Method ~ .) +
  theme_minimal() + xlim(c(0, 0.75)) + ylim(c(0,0.75)) +
  geom_segment(aes(x=0, xend=0.75, y=0, yend=0.75), color='black', linetype='dashed') +
  labs(colour='Method', y='Predicted Redshift', x='True Redshift',
       title='Predicted vs. True Redshift, 5,000 SDSS Test Galaxies, 20,000 Training Galaxies') +
  theme(strip.text.y = element_text(size=18))
ggsave(filename='images/full/photoz_predictions.pdf', plot=fig_points, dpi=300)
print(fig_points)


# P-values Figure #########################################################


pval_simulation <- function(resultsfile) {
  z_test <- h5read(resultsfile, "/z_true")
  z_grid <- h5read(resultsfile, "/z_grid")
  
  methods <- h5ls(resultsfile, recursive = FALSE)$name
  methods <- methods[!methods %in% c("z_grid", "z_true")]
  
  return(ldply(methods, function(method) {
    cde <- h5read(resultsfile, paste0(method, "/cde"))
    
    cdf_pvals <- cdetools::cdf_coverage(cde, z_grid, z_test)
    hpd_pvals <- cdetools::hpd_coverage(cde, z_grid, z_test)
    
    return(data.frame(Method = method,
                      pval = c(cdf_pvals, hpd_pvals),
                      Test = rep(c("PIT", "HPD"), each = nrow(cde))))
    }))
}

fname <- paste0(outdir, "results_full.hdf5")
df_pval <- pval_simulation(fname)

df_pval$Method <- factor(df_pval$Method, levels = c("Functional", "Vector"))

fig1 <- ggplot(df_pval, aes(x = pval, y = ..density..)) +
  geom_histogram() +
  facet_grid(Test ~ Method) + xlim(c(0,1)) +
  xlab("Values") + ylab("Density") +
  theme_minimal()
ggsave(filename='images/full/values_hpd_pit.pdf', plot=fig1, dpi=300)
print(fig1)


## PP-Plot
test_vec <- c('HPD', 'PIT')
methods <- c("Functional", "Vector")
run_params <- expand.grid(test = test_vec,
                          method = methods)

df_pp_plot <- plyr::mdply(run_params, function(test, method){
  df_temp <- df_pval %>% filter(Test==test, Method==method)
  tmp_pvals <- df_temp$pval
  
  n <- length(tmp_pvals)
  theoretical_pval <- (1 : n) / n - 0.5 / n
  empirical_pval <- sort(punif(tmp_pvals))
  
  return(cbind(theoretical_pval, 
               empirical_pval))
})

fig2 <- ggplot(df_pp_plot, aes(x = theoretical_pval, 
                              y = empirical_pval)) +
  geom_line(size=2) +
  geom_segment(aes(x=0, xend=1, y=0, yend=1), color='red', linetype='dashed') +
  ylim(c(0,1)) +
  facet_grid(test ~ method) +
  xlab("Theoretical Coverage") + ylab("Empirical Coverage") +
  theme_minimal()
ggsave(filename='images/full/ppplots_hpd_pit.pdf', plot=fig2, dpi=300)
print(fig2)


ggarrange(fig1 + xlab("Value \n\n (a)") + theme(plot.title = element_blank(),
                                          axis.title.x = element_text(size=22),
                                          axis.text.x = element_text(size=10,
                                                                     margin=margin(b=10)),
                                          axis.text.y = element_text(size=10,
                                                                     margin=margin(l=10)),
                                          axis.title.y = element_text(size=22),
                                          strip.text.x = element_text(size=22),
                                          strip.text.y = element_blank()),
          fig2 +xlab("Theoretical Coverage \n\n (b)") + theme(plot.title = element_blank(),
                                                             axis.title.x = element_text(size=22),
                                                             axis.text.x = element_text(size=10,
                                                                                        margin=margin(b=10)),
                                                             axis.text.y = element_text(size=10),
                                                             axis.title.y = element_text(size=22),
                                                             strip.text.x = element_text(size=22),
                                                             strip.text.y = element_text(size=22)),
          nrow = 1, ncol =2)

#16.35 x 7.86



