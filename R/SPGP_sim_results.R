#### Reproducible code for "A Nonstationary Soft Partitioned Gaussian Process Model via Random Spanning Trees"
#### Code for generating figures and tables in Simulation Studies (Section 4)

library(ggplot2)
library(fields)
library(ggpubr)
library(fossil)
library(colorRamps)
library(RColorBrewer)

### Plotting and util functions ------------------------------------------------

gg_ellipse <- function(rx, ry, xc, yc, color = "red", size = 0.3, ...) {
  x = xc + rx * cos(seq(0, pi, length.out=200))
  ymax = yc + ry * sin(seq(0, pi, length.out=200))
  ymin = yc + ry * sin(seq(0, -pi, length.out=200))
  annotate("ribbon", x = x, ymin = ymin, ymax = ymax, color = color, size = size, fill=NA, ...)
}

# function to plot clusters
plotCluster <- function(coords, cluster, coords_ho = NULL, cluster_ho = NULL, psize = 0.3, title = NULL){
  cluster = as.factor(cluster)
  if(!is.null(cluster_ho)) cluster_ho = as.factor(cluster_ho)
  p = ggplot() + 
    geom_point(data = data.frame(coords), aes(lon, lat, color = cluster), size = psize) +
    xlab(expression(s[h])) + ylab(expression(s[v])) +
    guides(color = guide_legend(title = "Cluster")) +
    ggtitle(title) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5, size = 10),
          legend.text = element_text(size = 9, hjust = 0, margin = margin(l = 3)))
  if(!is.null(cluster_ho)) p = p + geom_point(aes(x = lon, y = lat), shape = 2, size = psize,
                                              data = as.data.frame(coords_ho))
  return(p)
}

plotData <- function(coords, Data, psize = 0.3, col_lim = NULL, legend_name = NULL, title = NULL,
                     colors = rainbow(10)){
  if(missing(col_lim)) {col_lim = range(Data)}
  ggplot() + 
    geom_point(data=data.frame(coords), aes(lon, lat, color = Data), size = psize) +
    scale_color_gradientn(colours = colors, limits = col_lim, name = legend_name,
                          guide = guide_colourbar(barwidth = 0.5)) +
    xlab(expression(s[h])) + ylab(expression(s[v])) +
    ggtitle(title) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          # legend.title=element_blank(), 
          plot.title = element_text(hjust = 0.5, size = 10),
          legend.text = element_text(size = 9, hjust = 0, margin = margin(l = 3)))
}

# function to plot kmean and kvar
plotField <- function(coords, Data, col_lim = NULL, legend_name = NULL, title = NULL, colors = rainbow(10)){
  if(missing(col_lim)) {col_lim = range(Data)}
  ggplot() + 
    geom_raster(data = data.frame(coords), aes(lon, lat, fill = Data)) +
    scale_fill_gradientn(colours = colors, limits = col_lim, name = legend_name, na.value = 'white') +
    xlab(expression(s[h])) + ylab(expression(s[v])) +
    ggtitle(title) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          # legend.title=element_blank(), 
          plot.title = element_text(hjust = 0.5, size = 10),
          legend.text = element_text(size = 9, hjust = 0, margin = margin(l = 3)))
}

## check predictive distribution
# function for density plots
plotDensity <- function(val_mcmc, val_true, val_neighbor = NULL, x_lim = NULL, y_lim = NULL, title = NULL) {
  # get hdi
  require(HDInterval)
  hdinterval = as.data.frame(hdi(density(val_mcmc), allowSplit = T))
  
  legend_pos = 'none'
  
  val_mcmc = data.frame(vid = val_mcmc)
  p = ggplot(val_mcmc, aes(x = vid)) + 
    geom_vline(xintercept = val_true, linetype = 2, colour = 'blue', size = 0.6) +
    geom_density() +
    geom_segment(data = hdinterval, aes(y = 0, yend = 0, x = begin, xend = end), 
                 size = 1, colour = 'red') +
    geom_point(data = hdinterval, aes(y = 0, x = begin), size = 2, shape = 4, colour = 'red') +
    geom_point(data = hdinterval, aes(y = 0, x = end), size = 2, shape = 4, colour = 'red') +
    xlab('Prediction') + ylab('Posterior Density') +
    ggtitle(title) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.title = element_blank(), legend.position = legend_pos,
          plot.title = element_text(hjust = 0.5, size = 10))
  
  if(!is.null(val_neighbor)) {
    val_neighbor = data.frame(x = val_neighbor)
    p = p + geom_vline(data = val_neighbor, aes(xintercept = x, colour = 'neighbor',
                                                linetype = 'neighbor'), size = 0.6)
  }
  if(!is.null(x_lim) & is.null(y_lim)) p = p + coord_cartesian(xlim = x_lim)
  if(!is.null(y_lim) & is.null(x_lim)) p = p + coord_cartesian(ylim = y_lim)
  if(!is.null(y_lim) & !is.null(x_lim)) 
    p = p + coord_cartesian(xlim = x_lim, ylim = y_lim)
  return(p)
}

# get TGP cluster boundary
# function adopted from "TGP" R package
# original code: https://github.com/cran/tgp/blob/689168f5e43941e2808c36bc43603329641028db/R/mapT.R
"tgp.plot.parts.2d.2" <-
  function(parts, dx=c(1,2), what=NULL, trans=matrix(c(1,0,0,1), nrow=2),
           col=NULL, lwd=3)
  {
    indices <- seq(1,dim(parts)[1],4)
    rec = matrix(0, nrow = dim(parts)/4, ncol = 4)
    
    j <- 1
    for(i in indices) {
      a <- parts[i,dx[1]]; b <- parts[i+1,dx[1]];
      c <- parts[i,dx[2]]; d <- parts[i+1,dx[2]];
      x <- c(a, b, b, a, a);
      y <- c(c, c, d, d, c);
      xy <- as.matrix(cbind(x,y)) %*% trans
      rec[j, 1:2] = range(xy[, 1])
      rec[j, 3:4] = range(xy[, 2])
      j <- j+1
    }
    
    colnames(rec) = c('xmin', 'xmax', 'ymin', 'ymax')
    return(data.frame(rec))
  }

# function to get TGP MAP cluster
getTGPCluster <- function(tgp_res, coords) {
  part_tgp = tgp::partition(coords, tgp_res)
  n = nrow(coords)
  cluster_map_tgp = rep(0, n)
  for(j in 1:length(part_tgp)) {
    n_j = nrow(part_tgp[[j]])
    # row indices of coords_ho that are in cluster j
    idx_all_j = match(part_tgp[[j]], coords)
    idx_j = idx_all_j[1:n_j]
    cluster_map_tgp[idx_j] = j
  }
  return(cluster_map_tgp)
}

# function to get scorings
getScores <- function(krig_res, true_val, m = 1, n_post, NN_dist = NULL, cdist_mat = NULL) {
  require(scoringRules)
  
  kmean = krig_res$kmean[, 1:(m*n_post)]
  kvar = krig_res$kvar[, 1:(m*n_post)]
  
  n = nrow(kmean)
  
  weights = matrix(1/(m*n_post), nrow = n, ncol = m*n_post)
  if(!is.null(NN_dist)) {
    NN_weights = 1 / NN_dist ^ 4
    NN_weights = NN_weights / rowSums(NN_weights)
    NN_weights = t(apply(NN_weights, 1, rep, each = n_post))
    weights = NN_weights / n_post
  }
  
  crps = crps_mixnorm(true_val, kmean, sqrt(kvar), weights)
  logs = logs_mixnorm(true_val, kmean, sqrt(kvar), weights)
  
  return(c('CRPS' = mean(crps), 'LogS' = mean(logs)))
}

### Generating figures and tables in Section 4 ---------------------------------

### SPGP model selection =======================================================

load('data/SPGP_sim_DIC.RData')
res_dic = cbind(params_all, pred_res[, c(5, 1, 3, 4)])
colnames(res_dic) = c("L", "DIC", "MSPE", "CRPS", "LogS")

### SPGP model results =========================================================

source('R/SPGP_func_iso.R')
source('R/utils.R')
load("data/sim_input.RData")

Y_iso = Y; coords_iso = coords
col_lim = range(Y_iso)

load('data/SPSP_sim_results.RData')

## plot MAP clusters
log_post_out = mcmc_res$log_post_out
cluster_out = mcmc_res$cluster_out
cluster_obs_out = mcmc_res$cluster_obs_out
phi_out = mcmc_res$phi_out
tau_out = mcmc_res$tau_out
lambda_out = mcmc_res$lambda_out
sigmasq_out = mcmc_res$sigmasq_out
latent_out = mcmc_res$latent_out

# find MAP estimates
map_idx = which.max(log_post_out)
cluster_map = cluster_out[map_idx, ]
cluster_obs_map = cluster_obs_out[map_idx, ]
phi_map = phi_out[[map_idx]][cluster_obs_map]
tau_map = tau_out[[map_idx]][cluster_obs_map]
sigmasq_map = sigmasq_out[[map_idx]][cluster_obs_map]

## out-of-sample Rand index
cluster_ho_map = pred_res$cluster_ho[, map_idx]
rand = adj.rand.index(cluster_ho_map, cluster_ho_true)

## in-sample Rand index
rand_in = adj.rand.index(cluster_obs_map, cluster_true)

## scoring rules
scores_mean = colMeans(scores)

## in-sample estimation for log microergodic parameter eta
eta_true = log(sigmasq_uni[cluster_true]) + 5/2 * log(5) - 5 * log(phi_uni[cluster_true])
eta_map = log(sigmasq_map) + log(tau_map) + 5/2 * log(5) - 5 * log(phi_map)
MSE_eta = mean((eta_true - eta_map) ^ 2)

# plot log microergodic parameters as in Figure 2 (a, b)
col_lim_eta = c(-2.214863, 16.110282)
p_eta = plotData(coords_iso, eta_map, psize = 1, col_lim = col_lim_eta, legend_name = "", 
                 title = 'SPGP Log Mircroergodic\n Parameter Estimate', colors = rainbow(10)) +
  gg_ellipse(0.3, 0.3, 0.5, 0.5, size = 0.9)

p_eta_true = plotData(coords_iso, eta_true, psize = 1, col_lim = col_lim_eta, legend_name = "", 
                      title = 'True Log Mircroergodic Parameter\n', colors = rainbow(10)) +
  gg_ellipse(0.3, 0.3, 0.5, 0.5, size = 0.9) + 
  geom_point(aes(x = lon, y = lat), shape = 2, data = data.frame(coords_ho), 
             size = 1.2, alpha = 1, color = "darkgrey")

## Plot predictive distribution as in Figure S4 (a, d)

# plot predictive distribution
x_lim_bnd = c(-1, 2); y_lim_bnd = c(0, 1.6)
x_lim_int = c(-1.2, 0.6); y_lim_int = c(0, 3.8)
Y_ho_resp = pred_res$Y_ho

i = 96
p_pdist_bnd = plotDensity(Y_ho_resp[i, ], Y_ho[i], x_lim = x_lim_bnd, y_lim = y_lim_bnd,
                          title = 'SPGP Boundary Point')

i = 39
p_pdist_int = plotDensity(Y_ho_resp[i, ], Y_ho[i], x_lim = x_lim_int, y_lim = y_lim_int,
                          title = 'SPGP Interior Point')


## Plot predictive surface as in Figure 3 (b, f)
load('data/SPGP_sim_pred_field.RData')

color_field = RColorBrewer::brewer.pal(11, "Spectral")
col_lim_field = c(-1.6825, 2.58176)
col_lim_sd_field = c(0, 0.9596440)

f_pred = plotField(coords_grid, Y_grid_mean, title = 'SPGP Posterior \nPrediction Mean', 
                   colors = color_field, col_lim = col_lim_field)

f_pred_sd = plotField(coords_grid, Y_grid_sd, title = 'SPGP Posterior \nPrediction SD', 
                      colors = brewer.pal(9, "Reds"), col_lim = col_lim_sd_field)

### TGP model results ==========================================================

library(tgp)
load('data/TGP_sim_results.RData')

## get MAP clusters
cluster_map_tgp = getTGPCluster(tgp_res, coords)

## out-of-sample prediction scoring rules
Ym0r1 = tgp:::mean0.range1(Y)
Y_ho_tgp_all = apply(tgp_res$trace$preds$ZZ, 1, tgp:::undo.mean0.range1, Ym0r1$undo)
n_post = nrow(cluster_out)
kmean_tgp = apply(tgp_res$trace$preds$ZZ.km, 1, tgp:::undo.mean0.range1, Ym0r1$undo)
kvar_tgp = apply(tgp_res$trace$preds$ZZ.ks2, 1, tgp:::undo.mean0.range1, Ym0r1$undo,
                 nomean = T, s2 = T)
krig_res_tgp = list('kmean' = kmean_tgp, 'kvar' = kvar_tgp)

## Plot predictive distribution as in Figure S4 (b, e)
i = 96
p_pdist_bnd_tgp = plotDensity(Y_ho_tgp_all[i, ], Y_ho[i], x_lim = x_lim_bnd, y_lim = y_lim_bnd,
                              title = 'TGP Boundary Point')
i = 39
p_pdist_int_tgp = plotDensity(Y_ho_tgp_all[i, ], Y_ho[i], x_lim = x_lim_int, y_lim = y_lim_int,
                              title = 'TGP Interior Point')

## out-of-sample Rand index
# get MAP partition
cluster_ho_map_tgp = getTGPCluster(tgp_res, coords_ho)
rand_tgp = adj.rand.index(cluster_ho_map_tgp, cluster_ho_true)

## in-sample Rand index
rand_in_tgp = adj.rand.index(cluster_map_tgp, cluster_true)

## in-sample GP parameter estimation
n_post = (tgp_res$BTE[2] - tgp_res$BTE[1]) / tgp_res$BTE[3]
tgp_res_in_sample = predict(tgp_res, coords, BTE = c(0, 1, 1), krige = F, pred.n = F, trace = T)
d_tgp = matrix(0, nrow = n, ncol = n_post)  # d = phi / sqrt(2*nu)
s2_tgp = matrix(0, nrow = n, ncol = n_post) # GP variance
for(i in 1:n) {
  d_tgp[i, ] = tgp_res_in_sample$trace$XX[[i]][1, 7]
  s2_tgp[i, ] = tgp_res_in_sample$trace$XX[[i]][1, 3]
}

eta_tgp_all = log(s2_tgp) - 5 * log(d_tgp)
eta_tgp = rowMeans(eta_tgp_all)
MSE_eta_tgp = mean((eta_tgp - eta_true) ^ 2)

# plot log microergodic parameters as in Figure 2 (c)
p_eta_tgp = plotData(coords_iso, eta_tgp, psize = 1, col_lim = col_lim_eta, legend_name = "", 
                     title = 'TGP Log Mircroergodic\n Parameter Estimate', colors = rainbow(10)) +
  gg_ellipse(0.3, 0.3, 0.5, 0.5, size = 0.9)

## Plot predictive surface
load('data/TGP_sim_pred_field.RData')
n_grid = nrow(coords_grid)
Y_grid_mean_tgp = tgp_res_grid$ZZ.mean
Y_grid_sd_tgp = sqrt(tgp_res_grid$ZZ.s2[1:n_grid])

f_pred_tgp = plotField(coords_grid, Y_grid_mean_tgp, title = 'TGP Posterior \nPrediction Mean', 
                       colors = color_field, col_lim = col_lim_field)

f_pred_sd_tgp = plotField(coords_grid, Y_grid_sd_tgp, title = 'TGP Posterior \nPrediction SD', 
                      colors = brewer.pal(9, "Reds"), col_lim = col_lim_sd_field)

### Generate Table 1 ===========================================================

library(knitr)
table_1 = data.frame('SPGP' = c(rand_in, rand, MSE_eta, MSE_resp, scores_mean),
                     'TGP' = c(rand_in_tgp, rand_tgp, MSE_eta_tgp, mse_ho_tgp, scores_tgp))
rownames(table_1) = c("In-sample ARI", "Hold-out ARI", "MSE", "MSPE", "Mean CRPS", "Mean LogS")

## Generate Figure 2 ===========================================================

figure_2 = ggarrange(p_eta_true, p_eta, p_eta_tgp,
                     ncol = 3, nrow = 1, hjust = 0, widths = c(1, 1, 1),
                     common.legend = T, legend = "right", font.label = list(size=10),
                     labels = c("(a)", "(b)", "(c)"))

## Generate Figure 3 ===========================================================

figure_3_mean = ggarrange(f_pred, f_pred_tgp,
                          ncol = 2, nrow = 1, hjust = 0,
                          common.legend = T, legend = "right", font.label = list(size=10),
                          labels = c("(b)", "(c)"))

figure_3_sd = ggarrange(f_pred_sd, f_pred_sd_tgp,
                        ncol = 2, nrow = 1, hjust = 0,
                        common.legend = T, legend = "right", font.label = list(size=10),
                        labels = c("(f)", "(g)"))

## Generate Table S1 ===========================================================

table_S1 = t(res_dic)

## Generate Figure S4 ==========================================================

figure_S4_bnd = ggarrange(p_pdist_bnd, p_pdist_bnd_tgp,
                          ncol = 2, nrow = 1, hjust = 0,
                          common.legend = T, legend = "right", font.label = list(size=10),
                          labels = c("(a)", "(b)"))

figure_S4_int = ggarrange(p_pdist_int, p_pdist_int_tgp,
                          ncol = 2, nrow = 1, hjust = 0,
                          common.legend = T, legend = "right", font.label = list(size=10),
                          labels = c("(d)", "(e)"))
