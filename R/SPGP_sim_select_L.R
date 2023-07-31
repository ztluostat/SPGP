#### Reproducible code for "A Nonstationary Soft Partitioned Gaussian Process Model via Random Spanning Trees"
#### Code for selection hyperparameter L in SPGP model in Simulation Studies (Section 4)
#### NOTE: This code can take several hours to run.

library(igraph)
library(fields)
library(foreach)

source('R/SPGP_func_iso.R')
source('R/utils.R')

### Load simulation data -------------------------------------------------------

load("data/sim_input.RData")

### Set initial parameters for MCMC --------------------------------------------

## get initial microergodic parameters

# partition domain into cells
x_cutoff = seq(0, 1, length.out = 7)
y_cutoff = seq(0, 1, length.out = 7)
n_cell = (length(x_cutoff) - 1) * (length(y_cutoff) - 1)
cell_id = getGridIdx(coords, x_cutoff, y_cutoff)

coords_cell = matrix(0, nrow = n_cell, ncol = 2)
for (i in 1:n_cell) {
  idx_cell = which(cell_id == i)
  coords_cell[i, ] = colMeans(coords[idx_cell, ])
}
colnames(coords_cell) = c("lon", "lat")

# get corner coords of cells
bottomleft_cell = expand.grid(y_cutoff[-length(y_cutoff)], x_cutoff[-length(x_cutoff)])
bottomleft_cell[, 1:2] = bottomleft_cell[, c(2, 1)]
colnames(bottomleft_cell) = c("lon", "lat")

topright_cell = expand.grid(y_cutoff[-1], x_cutoff[-1])
topright_cell[, 1:2] = topright_cell[, c(2, 1)]
colnames(topright_cell) = c("lon", "lat")

# estimate microergodic using NNGP within each cell
theta_cell = numeric(n_cell)
for (i in 1:n_cell) {
  idx_cell = which(cell_id == i)
  if (length(idx_cell) == 0) next
  gp_res = GpGp::fit_model(Y[idx_cell], coords[idx_cell, ], matrix(1, ncol = 1, nrow = length(idx_cell)),
                           covfun_name = "matern25_isotropic", m_seq = 15, silent = T)
  theta_cell[i] = log(gp_res$covparms[1]) - 5 * log(gp_res$covparms[2])
}


## generate informed reference knots

n_knot_cells = 600 * theta_cell / sum(theta_cell)
coords_knot = matrix(0, ncol = 2, nrow = 0)
cell_id_knot = numeric()
for (i in 1:n_cell) {
  n_knot_cell = round(sqrt(n_knot_cells[i]))
  xmin = bottomleft_cell[i, 1]; ymin = bottomleft_cell[i, 2]
  xmax = topright_cell[i, 1]; ymax = topright_cell[i, 2]
  coords_knot = rbind(coords_knot, genGrid(n_knot_cell, n_knot_cell, xmin, ymin, xmax, ymax))
  cell_id_knot = c(cell_id_knot, rep(i, n_knot_cell ^ 2))
}
colnames(coords_knot) = c('lon', 'lat')
n_knot = nrow(coords_knot)

## construct spatial graph on knots

graph0 = dentrigraph(coords_knot, threshold = 0.2)
E(graph0)$eid = c(1:ecount(graph0))  # edge id
V(graph0)$vid = c(1:vcount(graph0))  # vertex id


## construct initial partition using fused lasso

# fit fused lasso to obtain initial clusters
k = 5
fl_res = genlasso::fusedlasso2d(matrix(theta_cell, nrow = length(y_cutoff) - 1, ncol = length(x_cutoff) - 1))
theta_cell_hat = fl_res$beta[, which(fl_res$df == k)[1]]
cluster_cell = as.integer(as.factor(theta_cell_hat))
cluster = cluster_cell[cell_id_knot]

# get a spanning tree that can induce this partition
edge_status = getEdgeStatus(cluster, graph0)
mstgraph0 = proposeMST(graph0, edge_status)
if('weight' %in% names(vertex_attr(graph0)))
  graph0 = delete_edge_attr(graph0, 'weight')
if('weight' %in% names(vertex_attr(mstgraph0)))
  mstgraph0 = delete_edge_attr(mstgraph0, 'weight')

# check if mstgraph0 is connected
if(ecount(mstgraph0) != n_knot - 1) stop("MST graph not connected")
# check if cluster is contiguous
if(sum(edge_status[E(mstgraph0)$eid] == 'b') >= k) stop("Non-contiguous partition")

## collect all initial values
init_val0 = list()
init_val0[['tree']] = mstgraph0
init_val0[['cluster']] = cluster
init_val0[['tau']] = rep(1, k)
init_val0[['phi']] = rep(1, k)
init_val0[['sigmasq']] = rep(1, k)
init_val0[['beta']] = 0
init_val0[['lambda']] = 1

### Set hyperparameters --------------------------------------------------------

# tuning parameters
tunning = c()
tunning['csize_min'] = 20    # minimum cluster size of knots
tunning['min_n_nngp'] = 250  # minimum cluster size for using NNGP approximation

# hyper-parameters
hyperpar0 = c()
hyperpar0['c'] = 0.5
hyperpar0['a_phi'] = 1; hyperpar0['b_phi'] = 1
hyperpar0['a_tau'] = 2; hyperpar0['b_tau'] = 2
hyperpar0['a_sigma'] = 2; hyperpar0['b_sigma'] = 0.1
hyperpar0['a_lambda'] = 2; hyperpar0['b_lambda'] = 2
hyperpar0['nngp'] = 15

# candidate parameter L for selection
params_all = expand.grid(
  nn = 1:5
)
n_params = nrow(params_all)

MCMC = 50000     # total number of MCMC iterations (including burn-in)
BURNIN = 25000   # number of burn-in iterations
THIN = 5         # thin MCMC by keeping results every 10 iterations 

### Run MCMC for each L and compute DIC ----------------------------------------

cl = parallel::makeForkCluster(5)
doParallel::registerDoParallel(cl)

res_all = foreach(pid = 1:n_params, .combine = rbind) %dopar% {
  
  ### Finalize hyperparameters =================================================
  
  # set hyperparameter L
  hyperpar = hyperpar0
  hyperpar['nn'] = params_all[pid, 1]

  # compute nearest knots and projection probability
  nn = hyperpar['nn']
  nn_res = FNN::get.knnx(coords_knot, coords, k = nn)
  nearest_knots = nn_res$nn.index
  proj_prob = 1 / nn_res$nn.dist
  for (i in 1:nrow(proj_prob)) {
    if (proj_prob[i, 1] == Inf) {
      proj_prob[i, ] = c(1, rep(0, nn - 1))
    } else {
      proj_prob[i, ] = proj_prob[i, ] / sum(proj_prob[i, ])
    }
  }

  # get random projection
  set.seed(1234)
  init_val = init_val0
  proj = randomProject(nearest_knots, proj_prob, nn)
  init_val[['proj']] = proj

  ### Run MCMC for SPGP model ==================================================
  
  mcmc_res = STGP_fit(Y, X, coords, coords_knot, graph0, init_val, hyperpar, tunning, MCMC, BURNIN, THIN, 
                  1234, NNGP = T, sample_latent = F, standardizeY = F, fix_proj = F)

  ### compute DIC for each L ===================================================
  log_post_out = mcmc_res$log_post_out
  log_like_out = mcmc_res$log_like_out
  map_idx = which.max(log_post_out)
  DIC_res = computeDIC(log_like_out, map_idx)

  ## Out-of-sample prediction for each L =======================================

  set.seed(12345)  # seed for prediction
  pred_res = STGP_predict2(mcmc_res, X, X_ho, coords, coords_ho, 
                           pred_response = T, nn = hyperpar['nn'], return_values = 'all',
                           NNGP = T, nn_nngp = hyperpar['nngp'], min_n_nngp = tunning['min_n_nngp'])
  Y_ho_resp = pred_res$Y_ho
  MSE_resp = mean((rowMeans(Y_ho_resp) - Y_ho)^2)
  MAE_resp = mean(abs(rowMeans(Y_ho_resp) - Y_ho))

  library(scoringRules)
  n_post = (MCMC - BURNIN) / THIN
  krig_res = STGP_krigmv2(mcmc_res, X, X_ho, coords, coords_ho, pred_response = T, 
                          nn = hyperpar['nn'], NNGP = T, nn_nngp = hyperpar['nngp'], min_n_nngp = tunning['min_n_nngp'])
  scores = getScores(krig_res, Y_ho, nn = hyperpar['nn'], n_post)

  # save results
  res = c('MSPE' = MSE_resp, 'MAPE' = MAE_resp)
  res = c(res, colMeans(scores), DIC_res)
  res
}

parallel::stopCluster(cl)

### Save results ---------------------------------------------------------------

pred_res = res_all

save(pred_res, params_all,
     MCMC, BURNIN, THIN, hyperpar0, tunning, 
     file = 'data/SPGP_sim_DIC.RData', compress = 'bzip2')
