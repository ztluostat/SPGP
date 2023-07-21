#### Reproducible code for "A Nonstationary Soft Partitioned Gaussian Process Model via Random Spanning Trees"
#### Code for generating TGP predictive field in Simulation Studies (Section 4)

library(tgp)

### Load fitted model ----------------------------------------------------------

load("data/TGP_sim_pred_field.RData")

### Generate grid points -------------------------------------------------------

n_grid = 3600
coords_grid = expand.grid(seq(0, 1, length.out = sqrt(n_grid)), seq(0, 1, length.out = sqrt(n_grid)))
colnames(coords_grid) = c('lon', 'lat')
X_grid = matrix(1, ncol = 1, nrow = n_grid)

### Generate predictive surface ------------------------------------------------

set.seed(12345)
n_post = (tgp_res$BTE[2] - tgp_res$BTE[1]) / tgp_res$BTE[3]
tgp_res_grid = predict(tgp_res, coords_grid, BTE = c(0, n_post * 10, 10), MAP = F, krige = F, pred.n = F, trace = F)

save(tgp_res_grid, coords_grid, file = 'data/TGP_sim_pred_field.RData', compress = 'bzip2')
