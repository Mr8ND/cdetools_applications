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
library(pracma)

outdir <- "data/"
save_plots <- FALSE
file_run <- "results_rfcde.hdf5"

grade_simulation <- function(resultsfile) {
  z_test <- h5read(resultsfile, "/z_true")
  z_grid <- h5read(resultsfile, "/z_grid")
  
  methods <- h5ls(resultsfile, recursive = FALSE)$name
  methods <- methods[!methods %in% c("z_grid", "z_true")]
  methods <- methods[!str_detect(methods, "mean")]
  
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
fname <- paste0(outdir, file_run)
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

n_plots <- 8
sampled_idx <- sample(x = seq(400), size = n_plots)

get_densities <- function(resultsfile, sampled_idx) {
  z_test <- h5read(resultsfile, "/z_true")
  z_grid <- h5read(resultsfile, "/z_grid")
  
  methods <- h5ls(resultsfile, recursive = FALSE)$name
  methods <- methods[!methods %in% c("z_grid", "z_true")]
  methods <- methods[!str_detect(methods, "mean")]
  
  return(ldply(methods, function(method) {
    cde <- h5read(resultsfile, paste0(method, "/cde"))
    print(dim(cde[sampled_idx, ]))
    
    return(data.frame(Method = method,
                      density = as.vector(cde[sampled_idx, ]),
                      y = rep(z_grid, each = n_plots),
                      id = rep(seq_len(n_plots), times = length(z_grid))))
  }))
}


fname <- paste0(outdir, file_run)
densities <- get_densities(fname, sampled_idx)
z_test <- h5read(fname, "/z_true")
truth  <-  data.frame(z_test = z_test[sampled_idx],
                      id = seq_len(n_plots))
methods <- c("Vector", "Functional", "Flexcode-Spec")
densities$Method <- factor(densities$Method, levels = methods)

fig <- ggplot(NULL, aes(x = y, y = density)) +
  geom_line(data = densities) +
  geom_vline(data = truth, aes(xintercept = z_test), lty = 2) +
  facet_grid(Method ~ id) +
  xlab("Redshift") + ylab("p(z | X)") +
  theme_minimal()
if (save_plots){
  ggsave(filename='images/cdes_examples.pdf', plot=fig, dpi=300)  
}

fig + theme(plot.title = element_blank(),
            axis.title.x = element_text(size=22),
            axis.text.x = element_text(size=9,
                                       margin=margin(b=5)),
            axis.text.y = element_text(size=10),
            axis.title.y = element_text(size=22),
            strip.text.x = element_text(size = 22),
            strip.text.y = element_text(size = 14))


# Stacked CDEs ############################################################

get_stacked_pz <- function(resultsfile) {
  z_test <- h5read(resultsfile, "/z_true")
  z_grid <- h5read(resultsfile, "/z_grid")
  
  methods <- h5ls(resultsfile, recursive = FALSE)$name
  methods <- methods[!methods %in% c("z_grid", "z_true")]
  methods <- methods[!str_detect(methods, "mean")]
  
  return(ldply(methods, function(method) {
    cde <- h5read(resultsfile, paste0(method, "/cde"))
    summed_cde <- apply(cde, MARGIN = 2, FUN = sum)
    summed_cde <- summed_cde/ pracma::trapz(x=z_grid, y=summed_cde)
    return(data.frame(Method = method,
                      stacked_cde = summed_cde,
                      z_grid = z_grid))
  }))
}

fname <- paste0(outdir, file_run)
df_pred <- get_stacked_pz(fname)
z_test <- h5read(fname, "/z_true")
dens_obj <- approxfun(density(z_test, bw='SJ'))

fig_points <- ggplot(data=df_pred) +
  geom_line(data=df_pred, aes(x=z_grid, y=stacked_cde)) + theme_minimal() +
  facet_grid(. ~ Method) + 
  geom_line(aes(x=z_grid, y=dens_obj(z_grid)), color='red') +
  labs(y='Stacked p(z)', x='Redshift',
       title='Stacked CDEs, 800 Test Galaxies, 2,000 Training Galaxies - from JCGS')
if (save_plots){
  ggsave(filename='images/stacked_cdes.pdf', plot=fig_points, dpi=300) 
}
print(fig_points)


# Redshift Prediction function Examples ############################################################

get_predictions <- function(resultsfile) {
  z_test <- h5read(resultsfile, "/z_true")
  z_grid <- h5read(resultsfile, "/z_grid")
  
  methods <- h5ls(resultsfile, recursive = FALSE)$name
  methods <- methods[!methods %in% c("z_grid", "z_true")]
  
  return(ldply(methods, function(method) {
    cde <- h5read(resultsfile, paste0(method, "/cde"))
    if (grepl(pattern='mean', x=method)){
      preds <- cde
    } else {
      preds <- apply(cde, MARGIN = 1, FUN = function(x) z_grid[which.max(x)])
    }
    return(data.frame(Method = method,
                      prediction = preds,
                      z_true = z_test))
  }))
}

fname <- paste0(outdir, file_run)
df_pred <- get_predictions(fname)
fig_points <- ggplot(data=df_pred,aes(z_true,prediction))+
  geom_point(alpha=1) +  facet_wrap(. ~ Method, ncol = 4) +
  theme_minimal() + #xlim(c(0, 0.2)) + ylim(c(0,0.2)) +
  geom_segment(aes(x=0, xend=0.9, y=0, yend=0.9), color='red', linetype='dashed') +
  labs(colour='Method', y='Predicted Redshift', x='True Redshift',
       title='Predicted vs. True Redshift, 800 Test Galaxies, 2,000 Training Galaxies - from JCGS') +
  theme(strip.text.y = element_text(size=18))
if (save_plots){
  ggsave(filename='images/photoz_predictions.pdf', plot=fig_points, dpi=300)
}

fig_points + theme(plot.title = element_blank(),
                   axis.title.x = element_text(size=22),
                   axis.title.y = element_text(size=22),
                   axis.text.x = element_text(size=14),
                   axis.text.y = element_text(size=14),
                   strip.text.x = element_text(size=18),
                   strip.text.y = element_text(size=18))


df_pred_2d <- df_pred %>% filter(!Method %in% 'Marginal') %>%
                     filter(!str_detect(Method, "mean"))
df_pred_2d$Method <- factor(df_pred_2d$Method, 
                               levels=c('Vector', 'Functional', 'Flexcode-Spec', 'rfcde-naive'))
fig_points <- ggplot(data=df_pred_2d, aes(x=z_true,y=prediction))+
  geom_abline(slope = 1, intercept = 0, linetype='dashed') +
  geom_density_2d() +  facet_grid(. ~ Method) +
  theme_minimal() + xlim(c(0.1, 0.75)) + ylim(c(0.1,0.75)) +
  labs(colour='Method', y='Predicted Redshift', x='True Redshift') +
  theme(strip.text.y = element_text(size=18)) 
if (save_plots){
  ggsave(filename='images/photoz_predictions_2d.pdf', plot=fig_points, dpi=300)
}
fig_points + coord_fixed() + theme(plot.title = element_blank(),
                   axis.title.x = element_text(size=22),
                   axis.title.y = element_text(size=22),
                   axis.text.x = element_text(size=14),
                   axis.text.y = element_text(size=14),
                   strip.text.x = element_text(size=18),
                   strip.text.y = element_text(size=18))


# P-values Figure #########################################################


pval_simulation <- function(resultsfile) {
  z_test <- h5read(resultsfile, "/z_true")
  z_grid <- h5read(resultsfile, "/z_grid")
  
  methods <- h5ls(resultsfile, recursive = FALSE)$name
  methods <- methods[!methods %in% c("z_grid", "z_true")]
  methods <- methods[!str_detect(methods, "mean")]
  
  return(ldply(methods, function(method) {
    cde <- h5read(resultsfile, paste0(method, "/cde"))
    
    cdf_pvals <- cdetools::cdf_coverage(cde, z_grid, z_test)
    hpd_pvals <- cdetools::hpd_coverage(cde, z_grid, z_test)
    
    return(data.frame(Method = method,
                      pval = c(cdf_pvals, hpd_pvals),
                      Test = rep(c("PIT", "HPD"), each = nrow(cde))))
    }))
}

fname <- paste0(outdir, file_run)
df_pval <- pval_simulation(fname)

df_pval$Method <- factor(df_pval$Method, levels = c("Vector", "Functional", 'Flexcode-Spec'))

fig1 <- ggplot(df_pval, aes(x = pval, y = ..density..)) +
  geom_histogram() +
  facet_grid(Test ~ Method) + xlim(c(0,1)) +
  xlab("Values") + ylab("Density") +
  theme_minimal()
if (save_plots){
  ggsave(filename='images/values_hpd_pit.pdf', plot=fig1, dpi=300)
}

fig1 + theme(plot.title = element_blank(),
             axis.title.x = element_text(size=22),
             axis.text.x = element_text(size=14),
             axis.text.y = element_text(size=14),
             axis.title.y = element_text(size=22),
             strip.text.x = element_text(size=22),
             strip.text.y = element_blank())

for (method in c("Vector", "Functional", "Flexcode-Spec")){
  for (val in c('HPD', 'PIT')){
    temp_df <- df_pval %>% filter(Method == method, Test == val)
    vals_temp <- temp_df$pval
    pval_temp <- ks.test(vals_temp, "punif")[2]$p.value
    print(c(method, val, pval_temp))
  }
}



## PP-Plot
test_vec <- c('HPD', 'PIT')
methods <- c("Vector", "Functional", "Flexcode-Spec")
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
if (save_plots){
  ggsave(filename='images/ppplots_hpd_pit.pdf', plot=fig2, dpi=300)
}
print(fig2)


ggarrange(fig1 + xlab("Value \n\n (a)") + theme(plot.title = element_blank(),
                                          axis.title.x = element_text(size=22),
                                          axis.text.x = element_text(size=10,
                                                                     margin=ggplot2::margin(b=10)),
                                          axis.text.y = element_text(size=10,
                                                                     margin=ggplot2::margin(l=10)),
                                          axis.title.y = element_text(size=22),
                                          strip.text.x = element_text(size=22),
                                          strip.text.y = element_blank()),
          fig2 +xlab("Theoretical Coverage \n\n (b)") + theme(plot.title = element_blank(),
                                                             axis.title.x = element_text(size=22),
                                                             axis.text.x = element_text(size=10,
                                                                                        margin=ggplot2::margin(b=10)),
                                                             axis.text.y = element_text(size=10),
                                                             axis.title.y = element_text(size=22),
                                                             strip.text.x = element_text(size=22),
                                                             strip.text.y = element_text(size=22)),
          nrow = 1, ncol =2)

#16.35 x 7.86



