#### Reproducible code for "A Nonstationary Soft Partitioned Gaussian Process Model via Random Spanning Trees"
#### Code for generating SPGP predictive field in Simulation Studies (Section 4)

source("R/SPGP_func_iso.R")

### Load inputs and fitted model -----------------------------------------------

load("data/sim_input.RData")
load("data/SPSP_sim_results.RData")

### Generate grid points -------------------------------------------------------

n_grid = 3600
coords_grid = expand.grid(seq(0, 1, length.out = sqrt(n_grid)), seq(0, 1, length.out = sqrt(n_grid)))
coords_grid = as.matrix(coords_grid)
colnames(coords_grid) = c('lon', 'lat')
X_grid = matrix(1, ncol = 1, nrow = n_grid)

### Generate predictive surface ------------------------------------------------

set.seed(12345)  # seed for prediction
pred_res_grid = STGP_predict2(mcmc_res, X, X_grid, coords, coords_grid, 
                              pred_response = T, nn = hyperpar['nn'], return_values = 'all',
                              NNGP = T, nn_nngp = hyperpar['nngp'], min_n_nngp = tunning['min_n_nngp'])

Y_grid_resp = pred_res_grid$Y_ho
Y_grid_mean = rowMeans(Y_grid_resp)
Y_grid_sd = apply(Y_grid_resp, 1, sd)

save(Y_grid_mean, Y_grid_sd, coords_grid, file = 'data/SPGP_sim_pred_field.RData', compress = 'bzip2')
