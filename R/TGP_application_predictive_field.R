#### Reproducible code for "A Nonstationary Soft Partitioned Gaussian Process Model via Random Spanning Trees"
#### Code for generating TGP predcitive field in precipitation application (Section 5)

### Load inputs and fitted model -----------------------------------------------

# load grid points in CONUS
load('data/precip_grid_input.RData')
coords_grid = as.matrix(coords_grid)
n_grid = nrow(coords_grid)
X_grid = matrix(1, ncol = 1, nrow = n_grid)

# load fitted model
library(tgp)
load('data/TGP_app_results.RData')

### Generate predictive field --------------------------------------------------

set.seed(12345)
n_post = (tgp_res$BTE[2] - tgp_res$BTE[1]) / tgp_res$BTE[3]
tgp_res_grid = predict(tgp_res, coords_grid, BTE = c(0, n_post * 10, 10), MAP = F, krige = F, pred.n = F, trace = F)

save(tgp_res_grid, coords_grid, file = 'data/TGP_app_pred_field.RData', compress = 'bzip2')
