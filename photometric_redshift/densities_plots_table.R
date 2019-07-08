library(rhdf5)
library(plyr)
library(ggplot2)
library(cdetools)
library(ascii)
library(tidyverse)
library(latex2exp)
library(gtable)
library(gridExtra)
library(grid)
library(ggpubr)
library(pracma)
source('auxiliary_funcs.R')

datadir <- "data/"
n_plots <- 3
methods <- c("RFCDE", "DeepCDE", "FlexZBoost", "NNKCDE", "RFCDE-Limited", "Marginal")
save_plots <- FALSE

# Method Table ##############################################################
df <- ldply(methods, function(method) {
  fname <- paste0(datadir, method, ".hdf5")
  
  y_grid <- h5read(fname, "/y_grid")
  y_true <- h5read(fname, "/y_true")
  cde <- h5read(fname, "/cde")
  if (method == 'NNKCDE'){
    cde <- t(apply(cde, MARGIN = 1, 
                   FUN = function(x) normalize_density(y_grid[2] - y_grid[1], x)))
  }
  loss <- cde_loss(cde, y_grid, y_true)
  
  return(data.frame(Method = method,
                    loss = loss$loss,
                    se = loss$se))
})

df$Method <- factor(df$Method, levels = methods)
df$loss <- round(df$loss, digits = 3)
df$se <- round(df$se, digits = 3)
df[, "CDE Loss (SE)"] <- paste0(df$loss, " (", df$se, ")")

tbl <- df[, c("Method", "CDE Loss (SE)")]

print(ascii(tbl, include.rownames = FALSE), type = "org")

# Figure
method_names <- c("RFCDE", "DeepCDE", "FlexZBoost", "NNKCDE", "RFCDE-Limited",
                  "Marginal")
method_names <- factor(df$Method, levels = rev(method_names))
values <- df$loss
sd_errs <- df$se

df_plot <- tibble(model = method_names,
                  cde_loss = values,
                  sd = sd_errs)

p1 <- ggplot(df_plot, aes(x=model, y=cde_loss)) + 
  geom_point(stat="identity") +
  geom_errorbar(aes(ymin=cde_loss-2*sd, ymax=cde_loss+2*sd), width=.2,
                position=position_dodge(.9)) +
  labs(title='Method Performance on Teddy data', x='Method', 
       y = latex2exp::TeX("CDE Loss (mean $\\pm$ 2 $\\sigma$)")) +
  theme_bw()

p1 + theme(plot.title = element_blank(),
           axis.title.x = element_text(size=26),
           axis.text.x = element_text(size=22, margin = ggplot2::margin(b=15)),
           axis.text.y = element_text(size=18),
           axis.title.y = element_text(size=26))


# CDE EXAMPLES ########################################################
n_plots <- 8
sampled_idx <- sample(x = seq(500), size = n_plots)

densities <- ldply(methods, function(method) {
  fname <- paste0(datadir, method, ".hdf5")
  
  y_grid <- h5read(fname, "/y_grid")
  y_true <- h5read(fname, "/y_true")
  cde <- h5read(fname, "/cde")
  
  return(data.frame(Method = method,
                    density = as.vector(cde[sampled_idx, ]),
                    y = rep(y_grid, each = n_plots),
                    id = rep(seq_len(n_plots), times = length(y_grid))))
})

method <- methods[1]
fname <- paste0(datadir, method, ".hdf5")
y_true <- h5read(fname, "/y_true")
truth  <-  data.frame(y_true = y_true[sampled_idx],
                      id = seq_len(n_plots))

densities$Method <- factor(densities$Method, levels = rev(methods))

fig <- ggplot(NULL, aes(x = y, y = density)) +
  geom_line(data = densities) +
  geom_vline(data = truth, aes(xintercept = y_true), lty = 2) +
  facet_grid(Method ~ id) +
  xlab("Redshift") + ylab("p(z | X)") +
  theme_minimal()
if (save_plots){
  ggsave(filename='images/cdes_examples.pdf', plot=fig, dpi=300)
}

fig + theme(plot.title = element_blank(),
            axis.title.x = element_text(size=22),
            axis.text.x = element_text(size=9,
                                       margin=ggplot2::margin(b=5)),
            axis.text.y = element_text(size=10),
            axis.title.y = element_text(size=22),
            strip.text.x = element_text(size = 22),
            strip.text.y = element_text(size = 14))

# P-values Figure #########################################################
df_pval <- ldply(methods, function(method) {
  fname <- paste0(datadir, method, ".hdf5")
  
  y_grid <- h5read(fname, "/y_grid")
  y_true <- h5read(fname, "/y_true")
  cde <- h5read(fname, "/cde")
  
  cdf_pvals <- cdetools::cdf_coverage(cde, y_grid, y_true)
  hpd_pvals <- cdetools::hpd_coverage(cde, y_grid, y_true)
  
  return(data.frame(Method = method,
                    pval = c(cdf_pvals, hpd_pvals),
                    Test = rep(c("PIT", "HPD"), each = nrow(cde))))
})

df_pval$Method <- factor(df_pval$Method, levels = methods)

fig <- ggplot(df_pval, aes(x = pval, y = ..density..)) +
  geom_histogram() +
  facet_grid(Test ~ Method) +
  xlab("Values") + ylab("Density") + xlim(c(0,1)) +
  theme_minimal()
if (save_plots){
  ggsave(filename='images/values_hpd_pit.pdf', plot=fig, dpi=300)
}
print(fig)


# Two panels

panel1 <- c("Marginal", "RFCDE-Limited", "NNKCDE")
panel1_df <- df_pval %>% filter(Method %in% panel1)
panel1_df$Method <- factor(panel1_df$Method, levels = panel1)
p1_plot <- ggplot(panel1_df, 
                  aes(x = pval, y = ..density..)) +
  geom_histogram() +
  facet_grid(Test ~ Method) +
  xlab("Value") + ylab("Density") +
  ylim(c(0,3)) +
  theme_minimal()
p1_plot + theme(plot.title = element_text(hjust = 0.5, size=26),
           axis.title.x = element_text(size=22),
           axis.text.x = element_text(size=14),
           axis.text.y = element_text(size=14),
           axis.title.y = element_text(size=20),
           strip.text.x = element_text(size=20),
           strip.text.y = element_blank())


panel2 <- c("FlexZBoost", "DeepCDE", "RFCDE")
panel2_df <- df_pval %>% filter(Method %in% panel2)
panel2_df$Method <- factor(panel2_df$Method, levels = panel2)
p2_plot <- ggplot(panel2_df, 
                  aes(x = pval, y = ..density..)) +
  geom_histogram() +
  facet_grid(Test ~ Method) + ylim(c(0,3)) +
  xlab("Value") + ylab("Density") +
  theme_minimal()
p2_plot + theme(plot.title = element_text(hjust = 0.5, size=26),
                axis.title.x = element_text(size=22),
                axis.text.x = element_text(size=14),
                axis.text.y = element_text(size=14),
                axis.title.y = element_blank(),
                strip.text.x = element_text(size = 20),
                strip.text.y = element_text(size = 20))


ggarrange(p1_plot + xlab('Value\n\n (a)') + theme(plot.title = element_text(hjust = 0.5, size=26),
                          axis.title.x = element_text(size=20),
                          axis.text.x = element_text(size=10,
                                                     margin=ggplot2::margin(b=10)),
                          axis.text.y = element_text(size=10,
                                                     margin=ggplot2::margin(l=10)),
                          axis.title.y = element_text(size=20),
                          strip.text.x = element_text(size = 20),
                          strip.text.y = element_blank()),
          p2_plot + xlab('Value\n\n (b)') + theme(plot.title = element_text(hjust = 0.5, size=26),
                          axis.title.x = element_text(size=20),
                          axis.text.x = element_text(size=10,
                                                     margin=ggplot2::margin(b=10)),
                          axis.text.y = element_text(size=10),
                          axis.title.y = element_blank(),
                          strip.text.x = element_text(size = 20),
                          strip.text.y = element_text(size = 20)),
          nrow = 1, ncol =2)

# PPplot Figure ##############################################################

test_vec <- c('HPD', 'PIT')
methods <- c("RFCDE", "RFCDE-Limited", "FlexZBoost", "DeepCDE", "Marginal", "NNKCDE")
run_params <- expand.grid(test = test_vec,
                          method = methods)

df_pp_plot <- plyr::mdply(run_params, function(test, method){
  df_temp <- df_pval %>% filter(Test==test, Method==method)
  tmp_pvals <- df_temp$pval
  
  n <- length(tmp_pvals)
  theoretical_pval <- (1 : n) / n - 0.5 / n
  empirical_pval <- sort(punif(tmp_pvals))
  
  theoretical_pval <- theoretical_pval[seq(1,length(theoretical_pval), length.out = 1000)]
  empirical_pval <- empirical_pval[seq(1,length(empirical_pval), length.out = 1000)]
  
  return(cbind(theoretical_pval, 
               empirical_pval))
})

fig <- ggplot(df_pp_plot, aes(x = theoretical_pval, 
                              y = empirical_pval)) +
  geom_line(size=2) +
  geom_segment(aes(x=0, xend=1, y=0, yend=1), color='red', linetype='dashed') +
  ylim(c(0,1)) +
  facet_grid(test ~ method) +
  xlab("Theoretical Coverage") + ylab("Empirical Coverage") +
  theme_minimal()
if (save_plots){
  ggsave(filename='images/ppplots_hpd_pit.pdf', plot=fig1, dpi=300)
}
print(fig)


# Two panels

panel1 <- c("Marginal", "RFCDE-Limited", "NNKCDE")
panel1_df <- df_pp_plot %>% filter(method %in% panel1)
panel1_df$method <- factor(panel1_df$method, levels = panel1)
p1_plot <- ggplot(panel1_df, 
                  aes(x = theoretical_pval, 
                      y = empirical_pval)) +
  geom_line(size=1.5) +
  geom_segment(aes(x=0, xend=1, y=0, yend=1), color='red', linetype='dashed', size=1) +
  ylim(c(0,1)) +
  facet_grid(test ~ method) +
  xlab("Theoretical Coverage") + ylab("Empirical Coverage") +
  theme_minimal()

p1_plot + theme(plot.title = element_blank(),
                axis.title.x = element_text(size=22),
                axis.text.x = element_text(size=14),
                axis.text.y = element_text(size=14),
                axis.title.y = element_text(size=22),
                strip.text.x = element_text(size=22),
                strip.text.y = element_blank())


panel2 <- c("FlexZBoost", "DeepCDE", "RFCDE")
panel2_df <- df_pp_plot %>% filter(method %in% panel2)
panel2_df$method <- factor(panel2_df$method, levels = panel2)
p2_plot <- ggplot(panel2_df, 
                  aes(x = theoretical_pval, 
                      y = empirical_pval)) +
  geom_line(size=1.5) +
  geom_segment(aes(x=0, xend=1, y=0, yend=1), color='red', linetype='dashed', size=1) +
  ylim(c(0,1)) +
  facet_grid(test ~ method) +
  xlab("Theoretical Coverage") + ylab("Empirical Coverage") +
  theme_minimal()
p2_plot + theme(plot.title = element_blank(),
                axis.title.x = element_text(size=22),
                axis.text.x = element_text(size=14),
                axis.text.y = element_text(size=14),
                axis.title.y = element_blank(),
                strip.text.x = element_text(size=22),
                strip.text.y = element_text(size=22))


ggarrange(p1_plot + xlab("Theoretical Coverage \n\n (c)") + theme(plot.title = element_blank(),
                          axis.title.x = element_text(size=22),
                          axis.text.x = element_text(size=10,
                                                     margin=ggplot2::margin(b=10)),
                          axis.text.y = element_text(size=10,
                                                     margin=ggplot2::margin(l=10)),
                          axis.title.y = element_text(size=22),
                          strip.text.x = element_text(size=22),
                          strip.text.y = element_blank()),
          p2_plot +xlab("Theoretical Coverage \n\n (d)") + theme(plot.title = element_blank(),
                          axis.title.x = element_text(size=22),
                          axis.text.x = element_text(size=10,
                                                     margin=ggplot2::margin(b=10)),
                          axis.text.y = element_text(size=10),
                          axis.title.y = element_blank(),
                          strip.text.x = element_text(size=22),
                          strip.text.y = element_text(size=22)),
          nrow = 1, ncol =2)

#16.35 x 7.86


# Stacked CDEs ############################################################


df_stacked_pz <- ldply(methods, function(method) {
  fname <- paste0(datadir, method, ".hdf5")
  
  y_grid <- h5read(fname, "/y_grid")
  y_true <- h5read(fname, "/y_true")
  cde <- h5read(fname, "/cde")
  summed_cde <- apply(cde, MARGIN = 2, FUN = sum)
  summed_cde <- summed_cde/ pracma::trapz(x=y_grid, y=summed_cde)
  
  return(data.frame(Method = method,
                    stacked_cde = summed_cde,
                    y_grid = y_grid))
})


method <- methods[1]
fname <- paste0(datadir, method, ".hdf5")
y_true <- h5read(fname, "/y_true")
dens_obj <- approxfun(density(y_true, bw='SJ'))

fig_points <- ggplot(data=df_stacked_pz) +
  geom_line(aes(x=y_grid, y=stacked_cde)) + theme_minimal() +
  facet_grid(. ~ Method) + 
  geom_line(aes(x=y_grid, y=dens_obj(y_grid)), color='red') +
  labs(y='Stacked p(z)', x='Redshift',
       title='Stacked CDEs, 800 Test Galaxies, 2,000 Training Galaxies - from JCGS')
if (save_plots){
  ggsave(filename='images/stacked_cdes.pdf', plot=fig_points, dpi=300) 
}
print(fig_points)


# Redshift Prediction function Examples ############################################################

df_pred <- ldply(methods, function(method) {
  fname <- paste0(datadir, method, ".hdf5")
  
  y_grid <- h5read(fname, "/y_grid")
  y_true <- h5read(fname, "/y_true")
  cde <- h5read(fname, "/cde")
  preds <- apply(cde, MARGIN = 1, FUN = function(x) y_grid[which.max(x)])
  
  return(data.frame(Method = method,
                    prediction = preds,
                    y_true = y_true))
})

methods <- c("RFCDE", "DeepCDE", "FlexZBoost", "NNKCDE", "RFCDE-Limited", "Marginal")
df_pred$Method <- factor(df_pred$Method, levels = rev(methods))

fig_points <- ggplot(data=df_pred,aes(y_true,prediction))+
  geom_point(alpha=1) +  facet_grid(. ~ Method) +
  theme_minimal() + xlim(c(0, 0.75)) + ylim(c(0,0.75)) +
  geom_segment(aes(x=0, xend=0.75, y=0, yend=0.75), color='black', linetype='dashed') +
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
                   strip.text.x = element_text(size=22),
                   strip.text.y = element_text(size=22))

# Stacked 2D densities for all of them apart from the Marginal

df_marginal <- df_pred %>% filter(str_detect(Method, "Marginal"))
fig_marginal <- ggplot(data=df_marginal,aes(y_true,prediction))+
  geom_point(alpha=1) + facet_grid(. ~ Method) +
  theme_minimal() + xlim(c(0, 0.7)) + ylim(c(0,0.7)) +
  geom_abline(slope=1, intercept=0, linetype='dashed') +
  labs(colour='Method', y='Predicted Redshift', x='True Redshift') +
  theme(strip.text.y = element_text(size=18)) + coord_fixed() 
print(fig_marginal)


df_pred_2d <- df_pred %>% filter(!Method %in% 'Marginal')
fig_points <- ggplot(data=df_pred_2d, aes(x=y_true,y=prediction))+
  geom_density_2d() +  facet_grid(. ~ Method) +
  theme_minimal() + #xlim(c(0.1, 0.7)) + ylim(c(0.1,0.7)) +
  geom_abline(slope=1, intercept=0, linetype='dashed') +
  labs(colour='Method', y='Predicted Redshift', x='True Redshift') +
  theme(strip.text.y = element_text(size=18)) + coord_fixed() 
print(fig_points)
fig_points + theme(plot.title = element_blank(),
                   axis.title.x = element_text(size=22),
                   axis.title.y = element_text(size=22),
                   axis.text.x = element_text(size=14),
                   axis.text.y = element_text(size=14),
                   strip.text.x = element_text(size=22),
                   strip.text.y = element_text(size=22))

# Photoz Metrics #################################################################
df_pred_metrics <- ldply(methods, function(method) {
  fname <- paste0(datadir, method, ".hdf5")
  
  y_grid <- h5read(fname, "/y_grid")
  y_true <- h5read(fname, "/y_true")
  cde <- h5read(fname, "/cde")
  preds <- apply(cde, MARGIN = 1, FUN = function(x) y_grid[which.max(x)])
  
  return(data.frame(Method = method,
                    prediction = preds,
                    y_true = y_true))
})

methods <- c("RFCDE", "DeepCDE", "FlexZBoost", "NNKCDE", "RFCDE-Limited", "Marginal")
df_pred_metrics$Method <- factor(df_pred_metrics$Method, levels = rev(methods))


sigma_f_func <- function(z, zspec){
  ratios <- (z - zspec)/(1 + zspec)
  sigmaf_val <- mean(ratios**2)
  return(sigmaf_val)
}

sigma_nmad_func <- function(z, zspec){
  ratios <- (abs(z - zspec))/(1 + zspec)
  sigma_nmad_out <- 1.48*median(ratios)
  return(sigma_nmad_out)
}

olf_func <- function(z, zspec){
  ratios <- (abs(z - zspec))/(1 + zspec)
  olf_val <- 100*sum(ratios > 0.15)/length(z)
  return(olf_val)
}


tbl_results <- df_pred_metrics %>% group_by(Method) %>%
  summarize(sigma_f = sigma_f_func(z = prediction, zspec = y_true),
            sigma_nmad = sigma_nmad_func(z = prediction, zspec = y_true),
            olf = olf_func(z = prediction, zspec = y_true))

print(ascii(tbl_results, include.rownames = FALSE, digits = 6), type = "org")



