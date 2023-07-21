#### Reproducible code for "A Nonstationary Soft Partitioned Gaussian Process Model via Random Spanning Trees"
#### Code for generating SPGP predcitive field in precipitation application (Section 5)


### Load inputs and fitted model -----------------------------------------------

# load grid points in CONUS
load('data/precip_grid.RData')
coords_grid = as.matrix(coords_grid)
n_grid = nrow(coords_grid)
X_grid = matrix(1, ncol = 1, nrow = n_grid)

# load fitted model
source('R/SPGP_func_aniso.R')
load('data/SPGP_app_results.RData')

### Generate predictive field --------------------------------------------------

set.seed(12345)
pred_res_grid = STGP_aniso_predict(mcmc_res, X, X_grid, coords, coords_grid, pred_response = T, m = 3, return_values = 'all')
Y_grid_resp = pred_res_grid$Y_ho
Y_grid_mean = rowMeans(Y_grid_resp)
Y_grid_sd = apply(pred_res_grid$Y_ho, 1, sd)

save(Y_grid_mean, Y_grid_sd, coords_grid, file = 'data/SPGP_app_pred_field.RData', compress = 'bzip2')
