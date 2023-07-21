#### Reproducible code for "A Nonstationary Soft Partitioned Gaussian Process Model via Random Spanning Trees"
#### Code for fitting SPGP model in precipitation application (Section 5)
#### NOTE: This code can take days to run.

library(igraph)
library(fields)

source('R/SPGP_func_aniso.R')
source('R/utils.R')

### Read and preprocess precipitation data -------------------------------------

dat = read.csv('data/CONUS_precip_west_18.csv')

Y = log(dat$precip)
n = nrow(dat)
X = matrix(1, nrow = n, ncol = 1)

coords = as.matrix(dat[, c('lcc_x', 'lcc_y')])
colnames(coords) = c('lon', 'lat')

# uniform hold-out data
set.seed(1234)
n_ho = 250
idx_ho = sample.int(n, n_ho, replace = F)
Y_ho = Y[idx_ho]
coords_ho = coords[idx_ho, , drop = F]
X_ho = X[idx_ho, , drop = F]

Y = Y[-idx_ho]
coords = coords[-idx_ho, , drop = F]
n = length(Y)
X = X[-idx_ho, , drop = F]

### Set initial parameters for MCMC --------------------------------------------

## construct spatial graph on observations
graph0 = dentrigraph(coords, threshold = 5)
E(graph0)$eid = c(1:ecount(graph0))  # edge id
V(graph0)$vid = c(1:vcount(graph0))  # vertex id

# divide the domain into 4 initial clusters
k = 9
x1_1 = quantile(coords[, 1], 0.33)
x1_2 = quantile(coords[, 1], 0.67)
x2_1 = quantile(coords[, 2], 0.33)
x2_2 = quantile(coords[, 2], 0.67)
cluster = rep(0, n)
cluster[coords[, 1] < x1_1 & coords[, 2] < x2_1] = 1
cluster[coords[, 1] < x1_2 & coords[, 2] < x2_1 & cluster == 0] = 2
cluster[coords[, 1] >= x1_2 & coords[, 2] < x2_1 & cluster == 0] = 3
cluster[coords[, 1] < x1_1 & coords[, 2] < x2_2 & cluster == 0] = 4
cluster[coords[, 1] < x1_2 & coords[, 2] < x2_2 & cluster == 0] = 5
cluster[coords[, 1] >= x1_2 & coords[, 2] < x2_2 & cluster == 0] = 6
cluster[coords[, 1] < x1_1 & coords[, 2] >= x2_2 & cluster == 0] = 7
cluster[coords[, 1] < x1_2 & coords[, 2] >= x2_2 & cluster == 0] = 8
cluster[coords[, 1] >= x1_2 & coords[, 2] >= x2_2 & cluster == 0] = 9

# get an initial spanning tree that can induce this partition
edge_status = getEdgeStatus(cluster, graph0)
mstgraph0 = proposeMST(graph0, edge_status)
graph0 = delete_edge_attr(graph0, 'weight')
mstgraph0 = delete_edge_attr(mstgraph0, 'weight')

# check if mstgraph0 is connected
if(ecount(mstgraph0) != n - 1) error("MST graph not connected")

# set initial values
init_val = list()
init_val[['tree']] = mstgraph0
init_val[['cluster']] = cluster
init_val[['tau']] = rep(1, k)
init_val[['phi1']] = rep(1, k)
init_val[['phi2']] = rep(1, k)
init_val[['theta']] = rep(pi/2, k)
init_val[['beta']] = rep(0, ncol(X))
init_val[['lambda']] = 1

### Set hyperparameters --------------------------------------------------------

# tuning parameters
tunning = c()
tunning['csize_min'] = 30  # minimum cluster size
tunning['k_max'] = round(2 * sqrt(n * log(n)))  # minimum cluster size

# hyper-parameters
hyperpar = c()
hyperpar['c'] = 0.5
hyperpar['a_phi'] = 1; hyperpar['b_phi'] = 1
hyperpar['a_theta'] = 0; hyperpar['b_theta'] = pi
hyperpar['a_tau'] = 2; hyperpar['b_tau'] = 2
hyperpar['a_sigma'] = 2; hyperpar['b_sigma'] = 2
hyperpar['a_lambda'] = 2; hyperpar['b_lambda'] = 2

MCMC = 50000    # total number of MCMC iterations (including burn-in)
BURNIN = 25000  # number of burn-in iterations
THIN = 5        # thin MCMC by keeping results every 5 iterations 

### Run MCMC for SPGP model ----------------------------------------------------

mcmc_res = STGP_aniso_fit(Y, X, coords, graph0, init_val, hyperpar, tunning,
                          MCMC, BURNIN, THIN, 1234, 
                          EB_control = list('parallel' = T, 'method' = 'L-BFGS-B'),
                          sample_latent = F, standardizeY = F)

# save MAP spanning tree only
log_post_out = mcmc_res$log_post_out
map_idx = which.max(log_post_out)
mcmc_res$ST_map = mcmc_res$ST_out[[map_idx]]
mcmc_res$ST_out = NULL

### Prediction from SPGP model -------------------------------------------------

## Point prediction ============================================================

# prediction with L = 1
set.seed(12345)
Y_ho_resp = STGP_aniso_predict(mcmc_res, X, X_ho, coords, coords_ho, pred_response = T, m = 1, return_values = 'all')$Y_ho
# prediction with L = 3
Y_ho_resp3 = STGP_aniso_predict(mcmc_res, X, X_ho, coords, coords_ho, pred_response = T, m = 3, return_values = 'all')$Y_ho
# prediction with L = 5
set.seed(12345)
Y_ho_resp5 = STGP_aniso_predict(mcmc_res, X, X_ho, coords, coords_ho, pred_response = T, m = 5, return_values = 'all')$Y_ho

## out-of-sample probabilistic prediction with scoring rules ===================

library(scoringRules)

# function to get scorings using Gaussian mixture
getScores <- function(krig_res, true_val, m = 1, n_post, 
                      dist = F, coords = NULL, coords_ho = NULL) {
  kmean = krig_res$kmean[, 1:(m*n_post)]
  kvar = krig_res$kvar[, 1:(m*n_post)]
  
  n = nrow(kmean)
  
  if(!dist) {
    # equal weights of NN
    weights = matrix(1/(m*n_post), nrow = n, ncol = m*n_post)
  } else {
    # distance-weighted
    NNdist = get.knnx(coords, coords_ho, k = m)$nn.dist
    NNweights = 1 / NNdist
    NNweights = NNweights / rowSums(NNweights)

    weights = matrix(0, nrow = n, ncol = m*n_post)
    for(i in 1:m) {
        NNweights_i = matrix(rep(NNweights[, i], n_post), ncol = n_post)
        weights[, ((i-1)*n_post+1):(i*n_post)] = NNweights_i / n_post
    }
  }
  
  crps = crps_mixnorm(true_val, kmean, sqrt(kvar), weights)
  logs = logs_mixnorm(true_val, kmean, sqrt(kvar), weights)
  
  return(data.frame('CRPS' = crps, 'LogS' = logs))
}

n_post = nrow(mcmc_res$cluster_out)
krig_res = STGP_aniso_krigmv(mcmc_res, X, X_ho, coords, coords_ho, pred_response = T, m = 5)
# probabilistic prediction with L=1
scores = getScores(krig_res, Y_ho, m = 1, n_post)
# probabilistic prediction with L=3
scores3 = getScores(krig_res, Y_ho, m = 3, n_post)
# probabilistic prediction with L=5
scores5 = getScores(krig_res, Y_ho, m = 5, n_post)

save.image(file = 'data/SPGP_app_results.RData', compress = 'xz')
