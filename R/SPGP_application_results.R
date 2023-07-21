#### Reproducible code for "A Nonstationary Soft Partitioned Gaussian Process Model via Random Spanning Trees"
#### Code for generating figures and tables in precipitation data application (Section 5)

library(igraph)
library(ggplot2)
library(fields)
library(ggpubr)
library(RColorBrewer)
library(colorRamps)
library(maps)
library(sf)

### Plotting and util functions ------------------------------------------------

plotGraphData2 <- function(coords, graph, Data, psize = 0.3, col_lim = NULL, colors = rainbow(10)){
  edgelist = get.edgelist(graph) 
  edgedata = data.frame(coords[edgelist[,1], ], coords[edgelist[, 2], ])
  colnames(edgedata) = c("X1", "Y1", "X2", "Y2")
  if(missing(col_lim)) {col_lim = range(Data)}
  coords = data.frame(coords)
  coords_sf = st_as_sf(coords, coords = c("lon", "lat"), crs = 4269)
  ggplot() + 
    geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2), data = edgedata, size = 0.5, colour = "grey") + 
    geom_sf(data = coords_sf, aes(color = Data), size = psize)+
    scale_color_gradientn(colours = colors, limits = col_lim) +
    xlab('Longitude') + ylab('Latitude') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.title = element_blank(), plot.title = element_text(hjust = 0.5, size = 10),
          legend.text = element_text(size = 9, hjust = 0, margin = margin(l = 3)))
}

plotData <- function(coords, Data, psize = 0.3, col_lim = NULL, colors = rainbow(10)){
  if(missing(col_lim)) {col_lim = range(Data)}
  coords = data.frame(coords)
  coords_sf = st_as_sf(coords, coords = c("lon", "lat"), crs = 4269)
  ggplot() + 
    geom_sf(data = coords_sf, aes(color = Data), size = psize) +
    scale_color_gradientn(colours = colors, limits = col_lim) +
    xlab('Longitude') + ylab('Latitude') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.title=element_blank(), plot.title = element_text(hjust = 0.5, size = 10),
          legend.text = element_text(size = 9, hjust = 0, margin = margin(l = 3)))
}

plotMapData <- function(coords, map, Data, psize = 0.3, lsize = 0.2, col_lim = NULL, colors = rainbow(10)){
  if(missing(col_lim)) {col_lim = range(Data)}
  coords = data.frame(coords)
  coords_sf = st_as_sf(coords, coords = c("lon", "lat"), crs = 4269)
  ggplot() + 
    geom_polygon(data = map, aes(x = long, y = lat, group = group), 
                 color = "black", size = lsize, fill = 'grey90') +
    geom_sf(data = coords_sf, aes(color = Data), size = psize) +
    scale_color_gradientn(colours = colors, limits = col_lim) +
    xlab('Longitude') + ylab('Latitude') +
    coord_sf(xlim = range(coords$lon) + c(-1, 1), ylim = range(coords$lat) + c(-2, 0), expand = F) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.title=element_blank(), plot.title = element_text(hjust = 0.5, size = 10),
          legend.text = element_text(size = 9, hjust = 0, margin = margin(l = 3)))
}

# function to plot spatial field with map
plotMapField <- function(coords_sf, map, Data, lsize = 0.2, col_lim = NULL, legend_name = NULL, title = NULL, 
                         colors = rainbow(10)){
  if(missing(col_lim)) {col_lim = range(Data)}
  coords = data.frame(st_coordinates(coords_sf))
  colnames(coords)[1:2] = c('lon', 'lat')
  p = ggplot() + 
    geom_sf(data = coords_sf, aes(fill = Data, colour = Data)) +
    scale_fill_gradientn(colours = colors, limits = col_lim, name = legend_name, na.value = 'white') +
    scale_colour_gradientn(colours = colors, limits = col_lim, name = legend_name, na.value = 'white', guide = F) +
    geom_polygon(data = map, aes(x = long, y = lat, group = group),
                 color = "grey", size = lsize, fill = NA) +
    xlab('Longitude') + ylab('Latitude') +
    coord_sf(xlim = range(coords$lon) + c(-1, 0), ylim = range(coords$lat) + c(-2, 0), expand = F) +
    ggtitle(title) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          # legend.title=element_blank(), 
          plot.title = element_text(hjust = 0.5, size = 10),
          legend.text = element_text(size = 9, hjust = 0, margin = margin(l = 3)))
  
  return(p)
}

# function to plot clusters
plotCluster <- function(coords, cluster, coords_ho = NULL, cluster_ho = NULL, psize = 0.3, title = NULL, map = NULL){
  cluster = as.factor(cluster)
  coords = as.data.frame(coords)
  coords_sf = st_as_sf(coords, coords = c("lon", "lat"), crs = 4269)
  if(!is.null(cluster_ho)) cluster_ho = as.factor(cluster_ho)
  if(!is.null(map)) {
    p = ggplot() + geom_polygon(data = map, aes(x = long, y = lat, group = group), 
                                color = "black", size = 0.2, fill = 'grey90')
  } else {p = ggplot()}
  p = p + 
    geom_sf(data = coords_sf, aes(color = cluster), size = psize) +
    xlab('Longitude') + ylab('Latitude') +
    guides(color = guide_legend(title = "Cluster")) +
    coord_sf(xlim = range(coords$lon) + c(-1, 1), ylim = range(coords$lat) + c(-2, 0)) +
    ggtitle(title) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5, size = 10),
          legend.text = element_text(size = 9, hjust = 0, margin = margin(l = 3)))
  if(!is.null(cluster_ho)) p = p + geom_point(aes(x = lon, y = lat), shape = 2, size = psize,
                                              data = as.data.frame(coords_ho))
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
    vertices = matrix(0, nrow = dim(parts)/4*5, ncol = 2)
    
    j <- 1
    for(i in indices) {
      a <- parts[i,dx[1]]; b <- parts[i+1,dx[1]];
      c <- parts[i,dx[2]]; d <- parts[i+1,dx[2]];
      x <- c(a, b, b, a, a);
      y <- c(c, c, d, d, c);
      xy <- as.matrix(cbind(x,y)) %*% trans
      vertices[j:(j + 4), ] = xy
      j <- j+5
    }
    
    # transform lcc to long/lat
    for(i in 1:nrow(vertices)) {
      idx = which(round(dat$lcc_x, 3) == round(vertices[i, 1], 3))
      if(length(idx) == 0) {
        warning('lcc value not found')
      } else {
        vertices[i, 1] = dat$longitude[idx[1]]
      }
      
      idx = which(round(dat$lcc_y, 4) == round(vertices[i, 2], 4))
      if(length(idx) == 0) {
        warning('lcc value not found')
      } else {
        vertices[i, 2] = dat$latitude[idx[1]]
      }
    }
    
    # plot boundary
    j = 1
    for(i in indices) {
      xy = vertices[(5*(j-1)+1) : (5*j), ]
      if(is.null(col)) { lines(xy, col=j, lty=j, lwd=lwd); }
      else { lines(xy, col=col, lty=1, lwd=lwd); }
      j <- j+1
    }
  }

# function to get TGP MAP cluster
getTGPCluster <- function(tgp_res, coords) {
  part_tgp = partition(coords, tgp_res)
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

### Generating figures and tables in Section 4 ---------------------------------

source('R/SPGP_func_aniso.R')
source('R/utils.R')

# load precipitation data
dat = read.csv('data/CONUS_precip_west_18.csv')

# load SPGP results
load('data/SPGP_app_results.RData')

coords = as.matrix(dat[, c('lcc_x', 'lcc_y')])
colnames(coords) = c('lon', 'lat')
coords_ho = coords[idx_ho, ]
coords = coords[-idx_ho, ]

coords_orig = as.matrix(dat[, c('longitude', 'latitude')])
colnames(coords_orig) = c('lon', 'lat')
coords_orig_ho = coords_orig[idx_ho, ]
coords_orig = coords_orig[-idx_ho, ]

# load CONUS map
usmap = map_data('state')

### Plot log precipitation as in Figure 4(a) ===================================

# indices of hold-out points near Rocky Mountains
idx_bnd = which(coords_orig_ho[, 1] <= -100 & coords_orig_ho[, 1] >= -115)

# plot raw data with boundary hold-out locations
p_raw_precip = plotGraphData2(coords_orig, graph0, Y, psize = 0.5, colors = hcl.colors(10, palette = "viridis", rev = T)) +
  geom_point(data = data.frame(coords_orig_ho[idx_bnd, ]), aes(lon, lat), size = 0.6, shape = 24, color = 'red') +
  scale_x_continuous(breaks = seq(-125, -90, 10)) +
  scale_y_continuous(breaks = seq(25, 50, 5)) +
  coord_sf(xlim = range(coords_orig[, 1]) + c(-1, 1), ylim = range(coords_orig[, 2]) + c(-2, 0), expand = F) +
  ggtitle('Log Precipitation Rate (log mm/day)')

# plot raw data with non-boundary hold-out locations
p_raw_precip_nonbnd = plotGraphData2(coords_orig, graph0, Y, psize = 0.5, colors = hcl.colors(10, palette = "viridis", rev = T)) +
  geom_point(data = data.frame(coords_orig_ho[-idx_bnd, ]), aes(lon, lat), size = 0.6, shape = 24, color = 'red') +
  scale_x_continuous(breaks = seq(-125, -90, 10)) +
  scale_y_continuous(breaks = seq(25, 50, 5)) +
  coord_sf(xlim = range(coords_orig[, 1]) + c(-1, 1), ylim = range(coords_orig[, 2]) + c(-2, 0), expand = F) +
  ggtitle('Log Precipitation Rate (log mm/day)')

### SPGP model =================================================================

log_post_out = mcmc_res$log_post_out
cluster_out = mcmc_res$cluster_out
phi1_out = mcmc_res$phi1_out
phi2_out = mcmc_res$phi2_out
tau_out = mcmc_res$tau_out
lambda_out = mcmc_res$lambda_out
sigmasq_out = mcmc_res$sigmasq_out
beta_out = mcmc_res$beta_out

# find MAP
map_idx = which.max(log_post_out)
cluster_map = cluster_out[map_idx, ]
phi1_map = phi1_out[[map_idx]][cluster_map]
phi2_map = phi2_out[[map_idx]][cluster_map]
tau_map = tau_out[[map_idx]][cluster_map]
sigmasq_map = sigmasq_out[[map_idx]][cluster_map]
beta_map = beta_out[[map_idx]]

## Plot MAP partition
p_cluster = plotCluster(coords_orig, cluster_map, psize = 0.5, map = usmap, title = 'SPGP MAP Partition') +
  scale_x_continuous(breaks = seq(-125, -90, 10)) +
  scale_y_continuous(breaks = seq(25, 50, 5))

## Get prediction metrics at boundary locations

idx_bnd = which(coords_orig_ho[, 1] <= -100 & coords_orig_ho[, 1] >= -115)

# MSE of hold-out prediction
mse_bnd = mean((rowMeans(Y_ho_resp)[idx_bnd] - Y_ho[idx_bnd])^2)
mse3_bnd = mean((rowMeans(Y_ho_resp3)[idx_bnd] - Y_ho[idx_bnd])^2)
mse5_bnd = mean((rowMeans(Y_ho_resp5)[idx_bnd] - Y_ho[idx_bnd])^2)

# scoring rules
scores_mean_bnd = colMeans(scores[idx_bnd, ])
scores3_mean_bnd = colMeans(scores3[idx_bnd, ])
scores5_mean_bnd = colMeans(scores5[idx_bnd, ])

## Get prediction metrics at non-boundary locations

# MSE of hold-out prediction
mse_nonbnd = mean((rowMeans(Y_ho_resp)[-idx_bnd] - Y_ho[-idx_bnd])^2)
mse3_nonbnd = mean((rowMeans(Y_ho_resp3)[-idx_bnd] - Y_ho[-idx_bnd])^2)
mse5_nonbnd = mean((rowMeans(Y_ho_resp5)[-idx_bnd] - Y_ho[-idx_bnd])^2)

# scoring rules
scores_mean_nonbnd = colMeans(scores[-idx_bnd, ])
scores3_mean_nonbnd = colMeans(scores3[-idx_bnd, ])
scores5_mean_nonbnd = colMeans(scores5[-idx_bnd, ])

## Plot predictive fields

load('data/precip_grid_input.RData')
load('data/SPGP_app_pred_field.RData')

col_lim = c(-4.05, 2.33)
col_lim_sd = c(0, 1.628)

p_pred_field = plotMapField(polyg_grid_orig, usmap, Y_grid_mean, col_lim = col_lim,
                            colors = hcl.colors(10, palette = "viridis", rev = T),
                            title = 'SPGP (L = 3) Posterior \nPrediction Mean') + 
  scale_x_continuous(breaks = seq(-125, -90, 10))

p_pred_sd = plotMapField(polyg_grid_orig, usmap, Y_grid_sd, col_lim = col_lim_sd,
                         colors = brewer.pal(9, "Reds"),
                         title = 'SPGP (L = 3) Posterior \nPrediction SD') + 
  scale_x_continuous(breaks = seq(-125, -90, 10))


### TGP model ==================================================================

library(tgp)
load('data/TGP_app_results.RData')

## plot MAP cluster

# get MAP partition
cluster_map_tgp = getTGPCluster(tgp_res, coords)

p_cluster_tgp = plotCluster(coords_orig, cluster_map_tgp, psize = 0.5, map = usmap, title = 'TGP MAP Partition') +
  scale_x_continuous(breaks = seq(-125, -90, 10)) +
  scale_y_continuous(breaks = seq(25, 50, 5))

## Get prediction metrics at boundary locations

# MSE of hold-out prediction
mse_tgp_bnd = mean((tgp_res$ZZ.mean[idx_bnd] - Y_ho[idx_bnd])^2)

# scoring rules
scores_tgp_mean_bnd = colMeans(scores_tgp[idx_bnd, ])

## Get prediction metrics at non-boundary locations

# MSE of hold-out prediction
mse_tgp_nonbnd = mean((tgp_res$ZZ.mean[-idx_bnd] - Y_ho[-idx_bnd])^2)

# scoring rules
scores_tgp_mean_nonbnd = colMeans(scores_tgp[-idx_bnd, ])

## Plot predictive fields

load('data/precip_grid_input.RData')
load('data/TGP_app_pred_field.RData')
n_grid = nrow(coords_grid)

p_pred_field_tgp = plotMapField(polyg_grid_orig, usmap, tgp_res_grid$ZZ.mean, col_lim = col_lim,
                            colors = hcl.colors(10, palette = "viridis", rev = T),
                            title = 'TGP Posterior \nPrediction Mean') + 
  scale_x_continuous(breaks = seq(-125, -90, 10))

p_pred_sd_tgp = plotMapField(polyg_grid_orig, usmap, sqrt(tgp_res_grid$ZZ.s2[1:n_grid]), col_lim = col_lim_sd,
                         colors = brewer.pal(9, "Reds"),
                         title = 'TGP Posterior \nPrediction SD') + 
  scale_x_continuous(breaks = seq(-125, -90, 10))

### Generate Figure 4 ==========================================================

figure_4 = ggarrange(p_raw_precip, p_cluster, p_cluster_tgp,
                     ncol = 3, nrow = 1, hjust = 0,
                     common.legend = F, legend = "right", font.label = list(size=10),
                     labels = c("(a)", "(b)", "(c)"))

### Generate Table 2 ===========================================================

table_2 = data.frame('SPGP (L=1)' = c(mse_bnd, scores_mean_bnd),
                     'SPGP (L=3)' = c(mse3_bnd, scores3_mean_bnd),
                     'SPGP (L=5)' = c(mse5_bnd, scores5_mean_bnd),
                     'TGP' = c(mse_tgp_bnd, scores_tgp_mean_bnd),
                     check.names = F)
rownames(table_2) = c("MSPE", "Mean CRPS", "Mean LogS")

### Generate Figure S7 =========================================================

figure_S7_mean = ggarrange(p_pred_field, p_pred_field_tgp, 
                           ncol = 2, nrow = 1,
                           common.legend = T, legend = "right", font.label = list(size=10),
                           labels = c("(a)", "(b)"))

figure_S7_sd = ggarrange(p_pred_sd, p_pred_sd_tgp,
                         ncol = 2, nrow = 1,
                         common.legend = T, legend = "right", font.label = list(size=10),
                         labels = c("(d)", "(e)"))

### Generate Table S4 ==========================================================

table_S4 = data.frame('SPGP (L=1)' = c(mse_nonbnd, scores_mean_nonbnd),
                      'SPGP (L=3)' = c(mse3_nonbnd, scores3_mean_nonbnd),
                      'SPGP (L=5)' = c(mse5_nonbnd, scores5_mean_nonbnd),
                      'TGP' = c(mse_tgp_nonbnd, scores_tgp_mean_nonbnd),
                      check.names = F)
rownames(table_S4) = c("MSPE", "Mean CRPS", "Mean LogS")
