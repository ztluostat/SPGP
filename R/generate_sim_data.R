#### Reproducible code for "A Nonstationary Soft Partitioned Gaussian Process Model via Random Spanning Trees"
#### Code for generating simulation data in Section 4

set.seed(123)

n = 2000
coords = runif(2*n)
coords = matrix(coords, nrow = n, ncol = 2)
colnames(coords) = c('lon', 'lat')

cluster_true = rep(0, n)
idx = ((coords[, 1] - 0.5)^2 + (coords[, 2] - 0.5)^2 >= 0.3^2)
cluster_true[idx] = 1
cluster_true[cluster_true == 0] = 2
k_true = max(cluster_true)  # true number of clusters

# generate hold-out locations
n_ho = 400
prop_bnd = 0.75  # proportion of hold-out samples near the boundary
coords_ho = matrix(0, nrow = n_ho, ncol = 2)
f_ho = numeric(n_ho)
i_bnd = 0; i_nbnd = 0
while(i_bnd + i_nbnd < n_ho) {
  coords_i = runif(2)
  if((coords_i[1] - 0.5)^2 + (coords_i[2] - 0.5)^2 >= 0.2^2 &
     (coords_i[1] - 0.5)^2 + (coords_i[2] - 0.5)^2 <= 0.4^2) {
    # is near boundary
    if(i_bnd < floor(n_ho * prop_bnd)) {
      coords_ho[i_bnd + i_nbnd + 1, ] = coords_i
      i_bnd = i_bnd + 1
    }
  } else {
    # not near boundary
    if(i_nbnd < n_ho - floor(n_ho * prop_bnd)) {
      coords_ho[i_bnd + i_nbnd + 1, ] = coords_i
      i_nbnd = i_nbnd + 1
    }
  }
}
colnames(coords_ho) = c('lon', 'lat')
cluster_ho_true = ifelse((coords_ho[, 1] - 0.5)^2 + (coords_ho[, 2] - 0.5)^2 >= 0.3^2, 1, 2)

# true covariance parameters
phi_uni = c(0.3, 1)
sigmasq_uni = c(1, 0.5)
beta0_uni = c(1, 1)

# compute covariance matrix
coords_all = rbind(coords, coords_ho)
clust_all = c(cluster_true, cluster_ho_true)
dist_mat = as.matrix(dist(coords_all))
cov_mat_list = list()
for(i in 1:k_true) {
  dist_mat_i = dist_mat[clust_all == i, clust_all == i, drop=F]
  scaled_dist_mat_i = dist_mat_i / phi_uni[i]
  cov_mat_i = (1 + sqrt(5) * scaled_dist_mat_i + 5/3 * scaled_dist_mat_i ^ 2) * exp(-sqrt(5) * scaled_dist_mat_i)
  cov_mat_i = sigmasq_uni[i] * cov_mat_i
  cov_mat_list[[i]] = cov_mat_i
}

# generate true mean function from piecewise GP
f_all = numeric(n + n_ho)
for(i in 1:k_true) {
  chol_i = t(chol(cov_mat_list[[i]]))
  f_all[clust_all == i] = chol_i %*% rnorm(sum(clust_all == i))
}
f_all = f_all + beta0_uni[clust_all]

# generate response Y
Y_all = f_all + rnorm(n + n_ho, 0, 0.1)

idx_ho = (n+1):(n+n_ho)
Y_ho = Y_all[idx_ho]
f_ho = f_all[idx_ho]
X_ho = matrix(1, nrow = n_ho, ncol = 1)

idx_train = 1:n
f = f_all[idx_train]
Y = Y_all[idx_train]
X = matrix(1, nrow = n, ncol = 1)
dist_mat = dist_mat[idx_train, idx_train]

# save input
save(coords, coords_ho, X, X_ho, Y, Y_ho, f, f_ho, cluster_true, cluster_ho_true, n, n_ho,
     phi_uni, sigmasq_uni, beta0_uni,
     file = "data/sim_input.RData", compress = "bzip2")