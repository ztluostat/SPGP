#### Reproducible code for "A Nonstationary Soft Partitioned Gaussian Process Model via Random Spanning Trees"
#### Functions for SPGP model with geometric anisotropic covariance functions

library(igraph)
library(fields)
library(optimParallel)
source('R/SPGP_func_iso.R')

splitCluster2 <- function(mstgraph, k, membership, csize, csize_min) { 
  cnt = 0
  while(TRUE) {
    clust_split = sample.int(k, 1, prob = csize - 1)
    edge_cutted = sample.int(csize[clust_split]-1, 1)
    
    mst_subgraph = induced_subgraph(mstgraph, membership == clust_split) 
    ## the above one can also be potentially improved if we also save/update a list with grouping info.   
    ## replace 'membership==clust.split' with 'group.list[[clust.split]]'
    
    mst_subgraph = delete.edges(mst_subgraph, edge_cutted)
    connect_comp = components(mst_subgraph)
    
    # check cluster sizes
    csize_new = connect_comp$csize
    cnt = cnt + 1
    if(min(csize_new) >= csize_min | cnt >= 500) break
  }
  
  cluster_new = connect_comp$membership
  idx_new = cluster_new == 2  # index for vertices belonging to new cluster
  vid_new = (V(mst_subgraph)$vid)[idx_new]  # vid for vertices belonging to new cluster
  vid_old = (V(mst_subgraph)$vid)[!idx_new]  # vid for vertices left in old cluster
  membership[vid_new] = k + 1
  csize[clust_split] = length(vid_old)
  csize[k+1] = length(vid_new)
  
  return(list(cluster = membership, vid_new = vid_new, vid_old = vid_old,
              clust_old = clust_split, csize = csize))
}

# function to get quadratic forms for each cluster
evalQuad <- function(Y, X, beta, dist_mat, cluster, k, phi, tau) {
  n = length(Y); p = ncol(X)
  idx = split(1:n, cluster)
  nu = 0.5
  e = Y - X %*% beta
  YPY_all = numeric(k)   # (Y-mu)'P^{-1}(Y-mu), where P = cov(Y) 
  XPY_all = matrix(0, nrow = k, ncol = p)  # X'P^{-1}Y
  XPX_all = array(0, dim = c(k, p, p))  # X'P^{-1}X
  chol_all = list()
  for(j in 1:k) {
    idx_j = idx[[j]]
    Y_j = Y[idx_j]; e_j = e[idx_j]; n_j = length(Y_j)
    X_j = X[idx_j, , drop = F]
    phi_j = phi[j]; tau_j = tau[j]
    
    # compute covariance matrix
    scaled_dist_j = dist_mat[idx_j, idx_j, drop=F] / phi_j
    P_j = (1 + sqrt(5) * scaled_dist_j + 5 / 3 * scaled_dist_j ^ 2) * exp(-sqrt(5) * scaled_dist_j) 
    P_j = tau_j * P_j
    diag(P_j) = diag(P_j) + 1  # P = Cov(Y | beta)
    
    # compute quadratic forms
    chol_j = t(chol(P_j))
    sol_e = forwardsolve(chol_j, e_j)
    sol_Y = forwardsolve(chol_j, Y_j)
    sol_X = forwardsolve(chol_j, X_j)
    YPY = sum(sol_e ^ 2)
    XPY = crossprod(sol_X, sol_Y)
    XPX = crossprod(sol_X)
    
    YPY_all[j] = YPY
    chol_all[[j]] = chol_j
    XPY_all[j, ] = XPY
    XPX_all[j, , ] = XPX
  }
  return(list('YPY' = YPY_all, 'XPY' = XPY_all, 'XPX' = XPX_all, 'chol' = chol_all))
}

# function to get log posterior density (up to a constant)
evalLogPost <- function(e, cluster, phi1_all, phi2_all, tau_all, sigmasq_all, beta, lambda, chols, k, hyper) {
  a_phi = hyper['a_phi']; b_phi = hyper['b_phi']
  a_theta = hyper['a_theta']; b_theta = hyper['b_theta']
  a_tau = hyper['a_tau']; b_tau = hyper['b_tau']
  a_sigma = hyper['a_sigma']; b_sigma = hyper['b_sigma']
  a_lambda = hyper['a_lambda']; b_lambda = hyper['b_lambda']
  c = hyper['c']
  
  n = length(e)
  if(k > 1) {
    idx = split(1:n, cluster)
  } else {  # k == 1 case
    idx = list(1:n)
  }
  
  # log likelihood
  log_like = 0
  for(j in 1:k) {
    idx_j = idx[[j]]
    e_j = e[idx_j]; n_j = length(e_j)
    sigmasq_j = sigmasq_all[j]
    chol_j = chols[[j]]
    sol = forwardsolve(chol_j, e_j)
    quad = sum(sol ^ 2)
    
    log_like = log_like - sum(log(diag(chol_j))) - n_j/2*log(sigmasq_j)
    log_like = log_like - quad / (2*sigmasq_j)
  }
  
  log_prior = sum(dtplus(c(phi1_all, phi2_all), a_phi, b_phi, log = T))
  log_prior = log_prior - k * log(b_theta - a_theta)
  log_prior = log_prior + sum(dinvgamma(tau_all, a_tau, b_tau, log = T))
  log_prior = log_prior + sum(dinvgamma(sigmasq_all, a_sigma, b_sigma, log = T))
  log_prior = log_prior + sum(dnorm(beta, 0, sqrt(lambda), log = T))
  log_prior = log_prior + dinvgamma(lambda, a_lambda, b_lambda, log = T)
  log_prior = log_prior - lchoose(n-1, k-1) + k*log(1-c)
  
  return(log_like + log_prior)
}

# function to get log posterior for GP (up to a constant) in a cluster
# Y may be interpreted as residuals here
# cov_param = c(log(phi1), log(phi2), tan(theta - pi / 2), log(tau))
evalLogGPPost <- function(cov_param, Y, coords, hyper) {
  phi1 = exp(cov_param[1]); phi2 = exp(cov_param[2])
  theta = atan(cov_param[3]) + pi / 2
  tau = exp(cov_param[4])
  a_phi = hyper['a_phi']; b_phi = hyper['b_phi']
  a_tau = hyper['a_tau']; b_tau = hyper['b_tau']
  a_theta = hyper['a_theta']; b_theta = hyper['b_theta']
  
  n = length(Y)
  dist_mat = MahalDist(coords, theta, phi1, phi2)$dist_mat
  log_like = evalLogLike(Y, dist_mat, rep(1, n), n, 1, 1, tau, 0, hyper)
  log_prior = sum(dtplus(c(phi1, phi2), a_phi, b_phi, log = T))
  log_prior = log_prior + dinvgamma(tau, a_tau, b_tau, log = T)
  log_prior = log_prior - log(b_theta - a_theta)
  return(log_like + log_prior)
}

# function to get gradient of log posterior for GP (up to a constant) in a cluster
# Y may be interpreted as residuals here
# cov_param = c(log(phi1), log(phi2), tan(theta - pi / 2), log(tau))
evalLogGPGrad <- function(cov_param, Y, coords, hyper) {
  phi1 = exp(cov_param[1]); phi2 = exp(cov_param[2])
  theta = atan(cov_param[3]) + pi / 2
  tau = exp(cov_param[4])
  
  a_sigma = hyper['a_sigma']; b_sigma = hyper['b_sigma']
  a_phi = hyper['a_phi']; b_phi = hyper['b_phi']
  a_tau = hyper['a_tau']; b_tau = hyper['b_tau']
  a_theta = hyper['a_theta']; b_theta = hyper['b_theta']
  n = length(Y)
  
  # get mahalanobis distance matrix
  dist_mat_res = MahalDist(coords, theta, phi1, phi2)
  dist_mat = dist_mat_res$dist_mat
  # get pairwise difference matrix of each coordinate
  D_1 = dist_mat_res$D_1
  D_2 = dist_mat_res$D_2
  
  # compute covariance matrix
  exp_dist = exp(-sqrt(5) * dist_mat)
  cov_mat = (1 + sqrt(5) * dist_mat + 5 / 3 * dist_mat ^ 2) * exp_dist
  cov_mat = tau * cov_mat
  diag(cov_mat) = diag(cov_mat) + 1
  
  # derivative of cov_mat wrt log(phi1)
  dcov_dphi1 = 5 * tau / 6 * (1 + sqrt(5) * dist_mat) * exp_dist / phi1 * (D_1 * cos(theta) + D_2 * sin(theta)) ^ 2
  # derivative of cov_mat wrt log(phi2)
  dcov_dphi2 = 5 * tau / 6 * (1 + sqrt(5) * dist_mat) * exp_dist / phi2 * (D_1 * sin(theta) - D_2 * cos(theta)) ^ 2
  # derivative of cov_mat wrt tan(theta - pi / 2)
  dcov_dtheta = 5 * tau / 6 * (1 + sqrt(5) * dist_mat) * exp_dist
  dcov_dtheta = dcov_dtheta * (1 / phi1 - 1 / phi2) / (1 + tan(theta - pi / 2) ^ 2)
  dcov_dtheta = dcov_dtheta * ((D_1 ^ 2 - D_2 ^ 2) * sin(2 * theta) - 2 * D_1 * D_2 * cos(2 * theta))
  
  # inverse cov_mat
  cov_inv = base::chol2inv( base::chol(cov_mat) )
  cov_inv_y = cov_inv %*% Y  # inv(cov_mat) %*% Y
  quad = crossprod(Y, cov_inv_y) # t(Y) %*% inv(cov_mat) %*% Y
  
  # compute gradient wrt phi1
  grad_phi1 = -0.5 * sum(t(cov_inv) * dcov_dphi1)
  grad_phi1 = grad_phi1 + (n / 2 + a_sigma) / (0.5 * quad + b_sigma) * 0.5 * t(cov_inv_y) %*% dcov_dphi1 %*% cov_inv_y
  grad_phi1 = grad_phi1 - (a_phi + 1) * phi1^2 / (a_phi * b_phi^2 + phi1^2)
  
  # compute gradient wrt phi2
  grad_phi2 = -0.5 * sum(t(cov_inv) * dcov_dphi2)
  grad_phi2 = grad_phi2 + (n / 2 + a_sigma) / (0.5 * quad + b_sigma) * 0.5 * t(cov_inv_y) %*% dcov_dphi2 %*% cov_inv_y
  grad_phi2 = grad_phi2 - (a_phi + 1) * phi2^2 / (a_phi * b_phi^2 + phi2^2)
  
  # compute gradient wrt theta
  grad_theta = -0.5 * sum(t(cov_inv) * dcov_dtheta)
  grad_theta = grad_theta + (n / 2 + a_sigma) / (0.5 * quad + b_sigma) * 0.5 * t(cov_inv_y) %*% dcov_dtheta %*% cov_inv_y
  
  # compute gradient of tau
  grad_tau = -0.5 * (n - sum(diag(cov_inv)))
  grad_tau = grad_tau + (n / 2 + a_sigma) / (0.5 * quad + b_sigma) / 2 * (quad - crossprod(cov_inv_y))
  grad_tau = grad_tau - (a_tau + 1) + b_tau / tau
  
  return(c(grad_phi1, grad_phi2, grad_theta, grad_tau))
}

# wrapper function to perform EB optimization (possibly in parallel)
EBOptim <- function(init, Y, coords, hyper, control, cl = NULL) {
  if (!control$parallel) {
    if (control$method == 'L-BFGS-B') {
      EB_res = try(optim(init, fn = evalLogGPPost, gr = evalLogGPGrad, Y = Y, coords = coords, hyper = hyper,
                         method = 'L-BFGS-B', lower = rep(-10, 4), upper = rep(10, 4),
                         control = list(fnscale = -1, maxit = 500)), silent = F)
    } else {
      EB_res = try(optim(init, fn = evalLogGPPost, gr = evalLogGPGrad, Y = Y, coords = coords, hyper = hyper,
                         method = control$method,
                         control = list(fnscale = -1, maxit = 500)), silent = F)
    }

  } else {
    EB_res = try(optimParallel(init, fn = evalLogGPPost, gr = evalLogGPGrad, Y = Y, coords = coords, hyper = hyper,
                               lower = rep(-10, 4), upper = rep(10, 4),
                               control = list(fnscale = -1, maxit = 500),
                               parallel = list(cl = cl)
                              ), 
                 silent = F)
  }
  
  if (is(EB_res, 'try-error')) {
    return(list('failed' = T))
  } else if (EB_res$convergence != 0) {
    return(list('failed' = T))
  } else {
    EB_res$failed = F
    return(EB_res)
  }
}

# function to get log_likelihood for a cluster using empirical Bayes
# Y may be interpreted as residuals here
evalLogLikeEB <- function(Y, coords, hyper, lphi1_init = 0, lphi2_init = 0, theta_tilde_init = 0, ltau_init = 0, control = list()) {
  a_phi = hyper['a_phi']; b_phi = hyper['b_phi']
  a_tau = hyper['a_tau']; b_tau = hyper['b_tau']
  a_theta = hyper['a_theta']; b_theta = hyper['b_theta']

  # setup default control
  if(is.null(control$parallel)) control$parallel = F
  if(is.null(control$method)) control$method = 'L-BFGS-B'
  
  # try given initial values first
  init = c(lphi1_init, lphi2_init, theta_tilde_init, ltau_init)
  EB_res = EBOptim(init, Y, coords, hyper, control = control)
  
  # if it fails, try default initial values: rep(0, 4)
  if(EB_res$failed & any(init != rep(0, 4))) {
    EB_res = EBOptim(rep(0, 4), Y, coords, hyper, control = control)
  }
  
  # if this also fails, return log_like = -Inf
  if(EB_res$failed)
    return(list('log_like' = -Inf, 'phi1_hat' = NA, 'phi2_hat' = NA, 'theta_hat' = NA, 'tau_hat' = NA))
  
  phi1_hat = exp(EB_res$par[1]); phi2_hat = exp(EB_res$par[2])
  theta_hat = atan(EB_res$par[3]) + pi / 2
  tau_hat = exp(EB_res$par[4])
  log_like = EB_res$value - sum(dtplus(c(phi1_hat, phi2_hat), a_phi, b_phi, log = T))
  log_like = log_like - dinvgamma(tau_hat, a_tau, b_tau, log = T)
  log_like = log_like + log(b_theta - a_theta)
  
  return(list('log_like' = log_like, 'phi1_hat' = phi1_hat, 'phi2_hat' = phi2_hat, 'theta_hat' = theta_hat, 'tau_hat' = tau_hat))
}

# function to update log-likelihood using Cholesky factors
# Y may be interpreted as residuals here
updateLogLike <- function(Y, cluster, csize, k, chols, hyper) {
  n = length(Y)
  if(k > 1) {
    idx = split(1:n, cluster)
  } else {  # k == 1 case
    idx = list(1:n)
  }
  a_sigma = hyper['a_sigma']; b_sigma = hyper['b_sigma']
  nu = 0.5
  
  log_like = numeric(k)
  for(j in 1:k) {
    idx_j = idx[[j]]; n_j = csize[j]
    Y_j = Y[idx_j]
    chol_j = chols[[j]]
    sol = forwardsolve(chol_j, Y_j)
    quad = sum(sol ^ 2)
    
    log_like_j = -sum(log(diag(chol_j))) - n_j/2*log(2*pi) + lgamma(n_j/2 + a_sigma)
    log_like_j = log_like_j - (n_j/2 + a_sigma) * log(b_sigma + quad/2)
    log_like[j] = log_like_j
  }
  log_like = log_like + a_sigma*log(b_sigma) - lgamma(a_sigma)
  return(log_like)
}

# function to get pairwise difference of a vector
# return: n*m matrix D, where n = length(x) and m = length(y),
# D[i, j] = x[i] - y[j]
pairDiff <- function(x, y = x) {
  D = sapply(x, function(z) z - y)
  return(t(D))
}

# function to compute Mahalanobis distance matrix D
# return: D[i, j] = (s_i - s_j)' S^{-1} (s_i - s_j)
MahalDist <- function(coords1, theta, phi1, phi2, coords2 = coords1) {
  # get pairwise difference of each coordinate
  D_1 = pairDiff(coords1[, 1], coords2[, 1])
  D_2 = pairDiff(coords1[, 2], coords2[, 2])
  D_v = cbind(c(D_1), c(D_2))  # vectorized form of D_1, D_2
  
  # compute S^{-1}
  S11 = cos(theta)^2 / phi1 + sin(theta)^2 / phi2
  S12 = sin(theta) * cos(theta) * (1/phi1 - 1/phi2)
  S22 = cos(theta)^2 / phi2 + sin(theta)^2 / phi1
  S = matrix(c(S11, S12, S12, S22), nrow = 2, ncol = 2)
  
  dist_mat = sqrt(rowSums(D_v %*% S * D_v))
  dist_mat = matrix(dist_mat, ncol = nrow(coords2), nrow = nrow(coords1))
  return(list('dist_mat' = dist_mat, 'D_1' = D_1, 'D_2' = D_2))
}


########### MCMC Function ###########

STGP_aniso_fit <- function(Y, X, coords, graph0, init_val, hyper, tunning,
                           MCMC, BURNIN, THIN, seed = 1234, EB_control = list(),
                           sample_latent = F, standardizeY = F) {
  set.seed(seed)
  
  # setup default control and parallel computing
  if(is.null(EB_control$parallel)) EB_control$parallel = F
  if(EB_control$parallel) {
    cl = makeCluster(4)
    setDefaultCluster(cl = cl)
    clusterExport(cl, c("MahalDist", "evalLogLike", "dtplus", "dinvgamma", "pairDiff"))
  }
  
  n = length(Y)
  p = ncol(X)
  
  # standardize Y
  if(standardizeY) {
    std_res = standardize(Y)
    Y = std_res$x
    std_par_Y = std_res$std_par
    rm(std_res)
  } else {
    std_par_Y = c('max' = 1, 'min' = 0, 'mean' = 0)
  }
  
  # initial values
  if('name' %in% names(vertex_attr(graph0))) {
    graph0 = delete_vertex_attr(graph0, 'name')
  }
  stgraph = init_val[['tree']]
  cluster = init_val[['cluster']]
  tau_all = init_val[['tau']]
  phi1_all = init_val[['phi1']]
  phi2_all = init_val[['phi2']]
  theta_all = init_val[['theta']]
  beta = init_val[['beta']]
  lambda = init_val[['lambda']]
  k = max(cluster)  # number of clusters
  csize = as.numeric(table(cluster))  # cluster sizes
  e = Y - X %*% beta  # Y - X %*% beta
  
  # hyper-parameter
  c = hyper['c']
  a_tau = hyper['a_tau']; b_tau = hyper['b_tau']
  a_phi = hyper['a_phi']; b_phi = hyper['b_phi']
  a_theta = hyper['a_theta']; b_theta = hyper['b_theta']
  a_sigma = hyper['a_sigma']; b_sigma = hyper['b_sigma']
  a_lambda = hyper['a_lambda']; b_lambda = hyper['b_lambda']
  
  # tunning paramter
  csize_min = tunning['csize_min']  # minimum cluster size
  k_max = tunning['k_max']  # maximum number of clusters
  #rb = 0.425; rd = 0.425
  #rc = 0.1; rhy = 0.05
  rc = 0.19; rhy = 0.01
  
  # compute distance matrix
  dist_mat = matrix(0, nrow = n, ncol = n)
  for(j in 1:k) {
    idx_j = cluster == j
    coords_j = coords[idx_j, , drop = F]
    dist_mat[idx_j, idx_j] = MahalDist(coords_j, theta_all[j], phi1_all[j], phi2_all[j])$dist_mat
  }
  # log_like: vector of length k
  # loh-likelihood in each cluster
  log_like = evalLogLike(e, dist_mat, cluster, csize, k, rep(1, k), tau_all, rep(0, k), hyper)
  # whether an edge in graph0 is within a cluster or bewteen two clusters
  # n*p matrix
  edge_status = getEdgeStatus(cluster, graph0)
  
  ## MCMC results
  phi1_out = list(); phi2_out = list(); theta_out = list()
  tau_out = list()
  sigmasq_out = list()
  lambda_out = numeric((MCMC-BURNIN)/THIN)
  beta_out = matrix(0, nrow = (MCMC-BURNIN)/THIN, ncol = p)
  latent_out = matrix(0, nrow = (MCMC-BURNIN)/THIN, ncol = n)
  cluster_out = matrix(0, nrow = (MCMC-BURNIN)/THIN, ncol = n)
  ST_out = list()
  log_post_out = numeric((MCMC-BURNIN)/THIN)
  
  cnt_acc = rep(0, 3)
  cnt_att = rep(0, 3)
  
  for(iter in 1:MCMC) {
    if(k == 1) {rb = 0.99; rd = 0; rc = 0
    } else if(k == k_max) {rb = 0; rd = 0.8; rc = 0.19
    } else {rb = 0.4; rd = 0.4; rc = 0.19}
    move = sample(4, 1, prob = c(rb, rd, rc, rhy))
    if(move != 4) cnt_att[move] = cnt_att[move] + 1
    
    if(move == 1) { ## Birth move
      # split an existing cluster
      split_res = splitCluster2(stgraph, k, cluster, csize, csize_min)
      membership_new = split_res$cluster
      vid_new = split_res$vid_new
      vid_old = split_res$vid_old
      clust_old = split_res$clust_old
      
      lphi1 = log(phi1_all[clust_old])
      lphi2 = log(phi2_all[clust_old])
      theta_tilde = tan(theta_all[clust_old] - pi / 2)
      ltau = log(tau_all[clust_old])
      
      # compute log-prior ratio
      log_A = log(1-c)
      
      # compute log-proposal ratio
      if(k == k_max - 1) {
        rd_new = 0.8
      } else {rd_new = 0.4}
      log_P = log(rd_new) - log(rb)
      
      # compute log-likelihood ratio
      log_like_new = log_like
      EB_res_old = evalLogLikeEB(e[vid_old], coords[vid_old, , drop=F], hyper, lphi1, lphi2, theta_tilde, ltau, EB_control)
      log_like_new[clust_old] = EB_res_old$log_like
      EB_res_new = evalLogLikeEB(e[vid_new], coords[vid_new, , drop=F], hyper, lphi1, lphi2, theta_tilde, ltau, EB_control)
      log_like_new[k + 1] = EB_res_new$log_like
      log_L = log_like_new[clust_old] + log_like_new[k+1] - log_like[clust_old]
      
      # acceptance probability
      acc_prob = min(0, log_A + log_P + log_L)
      acc_prob = exp(acc_prob)
      if(runif(1) < acc_prob){
        # accept
        cluster = membership_new
        phi1_all[clust_old] = EB_res_old$phi1_hat
        phi1_all[k + 1] = EB_res_new$phi1_hat
        phi2_all[clust_old] = EB_res_old$phi2_hat
        phi2_all[k + 1] = EB_res_new$phi2_hat
        theta_all[clust_old] = EB_res_old$theta_hat
        theta_all[k + 1] = EB_res_new$theta_hat
        tau_all[clust_old] = EB_res_old$tau_hat
        tau_all[k + 1] = EB_res_new$tau_hat
        k = k + 1
        log_like = log_like_new
        edge_status = getEdgeStatus(membership_new, graph0)
        csize = split_res$csize
        
        dist_mat[vid_old, vid_old] = MahalDist(coords[vid_old, , drop=F], EB_res_old$theta_hat, 
                                               EB_res_old$phi1_hat, EB_res_old$phi2_hat)$dist_mat
        dist_mat[vid_new, vid_new] = MahalDist(coords[vid_new, , drop=F], EB_res_new$theta_hat, 
                                               EB_res_new$phi1_hat, EB_res_new$phi2_hat)$dist_mat
        
        cnt_acc[move] = cnt_acc[move] + 1
      }
    }
    
    if(move == 2) { ## Death move
      # merge two existing clusters (c1, c2) -> c2
      merge_res = mergeCluster(stgraph, edge_status, cluster, csize)
      membership_new = merge_res$cluster
      # c1 is merged into c2
      vid_c1 = merge_res$vid_old  
      vid_c2 = merge_res$vid_new  # vids in c2 before merging
      clust_c1 = merge_res$clust_old # cluster indices are the ones before merging
      clust_c2 = merge_res$clust_new
      vid_new = c(vid_c1, vid_c2)  # vids in merged cluster
      
      if(csize[clust_c1] > csize[clust_c2]) {
        lphi1 = log(phi1_all[clust_c1])
        lphi2 = log(phi2_all[clust_c1])
        theta_tilde = tan(theta_all[clust_c1] - pi / 2)
        ltau = log(tau_all[clust_c1])
      } else {
        lphi1 = log(phi1_all[clust_c2])
        lphi2 = log(phi2_all[clust_c2])
        theta_tilde = tan(theta_all[clust_c2] - pi / 2)
        ltau = log(tau_all[clust_c2])
      }
      
      # compute log-prior ratio
      log_A = -log(1-c)
      
      # compute log-proposal ratio
      if(k == 2) {rb_new = 0.99
      } else {rb_new = 0.4}
      log_P = log(rb_new) - log(rd)
      
      # compute log-likelihood ratio
      vid_new = c(vid_c1, vid_c2)  # vids in merged cluster
      log_like_new = log_like
      # log-likelihood of merged cluster
      EB_res_merged = evalLogLikeEB(e[vid_new], coords[vid_new, , drop=F], hyper, lphi1, lphi2, theta_tilde, ltau, EB_control)
      log_like_new_c2 = EB_res_merged$log_like
      log_like_new[clust_c2] = log_like_new_c2
      log_like_new = log_like_new[-clust_c1]
      log_L = log_like_new_c2 - log_like[clust_c1] - log_like[clust_c2]
      
      #acceptance probability
      acc_prob = min(0, log_A + log_P + log_L)
      acc_prob = exp(acc_prob)
      if(runif(1) < acc_prob){
        # accept
        cluster = membership_new
        k = k - 1
        log_like = log_like_new
        edge_status = getEdgeStatus(membership_new, graph0)
        csize = merge_res$csize
        
        phi1_all[clust_c2] = EB_res_merged$phi1_hat
        phi1_all = phi1_all[-clust_c1]
        phi2_all[clust_c2] = EB_res_merged$phi2_hat
        phi2_all = phi2_all[-clust_c1]
        theta_all[clust_c2] = EB_res_merged$theta_hat
        theta_all = theta_all[-clust_c1]
        tau_all[clust_c2] = EB_res_merged$tau_hat
        tau_all = tau_all[-clust_c1]
        
        dist_mat[vid_new, vid_new] = MahalDist(coords[vid_new, , drop=F], EB_res_merged$theta_hat, 
                                               EB_res_merged$phi1_hat, EB_res_merged$phi2_hat)$dist_mat
        
        cnt_acc[move] = cnt_acc[move] + 1
      }
    }
    
    if(move == 3) { ## change move
      # first perform death move: (c1, c2) -> c2
      merge_res = mergeCluster(stgraph, edge_status, cluster, csize)
      # then perform birth move
      split_res = splitCluster2(stgraph, k-1, merge_res$cluster, merge_res$csize, csize_min)
      
      # update parameters for merge
      vid_c1 = merge_res$vid_old  
      vid_c2 = merge_res$vid_new  # vids in c2 before merging
      clust_c1 = merge_res$clust_old # cluster indices are the ones before merging
      clust_c2 = merge_res$clust_new

      if(csize[clust_c1] > csize[clust_c2]) {
        lphi1 = log(phi1_all[clust_c1])
        lphi2 = log(phi2_all[clust_c1])
        theta_tilde = tan(theta_all[clust_c1] - pi / 2)
        ltau = log(tau_all[clust_c1])
      } else {
        lphi1 = log(phi1_all[clust_c2])
        lphi2 = log(phi2_all[clust_c2])
        theta_tilde = tan(theta_all[clust_c2] - pi / 2)
        ltau = log(tau_all[clust_c2])
      }
      
      log_like_new = log_like
      vid_merged = c(vid_c1, vid_c2)  # vids in merged cluster
      EB_res_merged = evalLogLikeEB(e[vid_merged], coords[vid_merged, , drop=F], hyper, lphi1, lphi2, theta_tilde, ltau, EB_control)
      log_like_new_c2 = EB_res_merged$log_like
      log_like_new[clust_c2] = log_like_new_c2
      log_like_new = log_like_new[-clust_c1]
      
      # reject if EB_res_merged$log_like == -Inf
      if(log_like_new_c2 > -Inf) {
      
        phi1_all_new = phi1_all; phi2_all_new = phi2_all
        theta_all_new = theta_all; tau_all_new = tau_all
        phi1_all_new[clust_c2] = EB_res_merged$phi1_hat
        phi1_all_new = phi1_all_new[-clust_c1]
        phi2_all_new[clust_c2] = EB_res_merged$phi2_hat
        phi2_all_new = phi2_all_new[-clust_c1]
        theta_all_new[clust_c2] = EB_res_merged$theta_hat
        theta_all_new = theta_all_new[-clust_c1]
        tau_all_new[clust_c2] = EB_res_merged$tau_hat
        tau_all_new = tau_all_new[-clust_c1]
        
        # update parameters for split
        membership_new = split_res$cluster
        vid_new = split_res$vid_new
        vid_old = split_res$vid_old
        clust_old = split_res$clust_old
        
        lphi1 = log(phi1_all_new[clust_old])
        lphi2 = log(phi2_all_new[clust_old])
        theta_tilde = tan(theta_all[clust_old] - pi / 2)
        ltau = log(tau_all_new[clust_old])
        
        log_like_old = log_like_new[clust_old]  # log-likelihood of the cluster before splitting
        EB_res_old = evalLogLikeEB(e[vid_old], coords[vid_old, , drop=F], hyper, lphi1, lphi2, theta_tilde, ltau, EB_control)
        log_like_new[clust_old] = EB_res_old$log_like
        EB_res_new = evalLogLikeEB(e[vid_new], coords[vid_new, , drop=F], hyper, lphi1, lphi2, theta_tilde, ltau, EB_control)
        log_like_new[k] = EB_res_new$log_like
        
        phi1_all_new[clust_old] = EB_res_old$phi1_hat
        phi1_all_new[k] = EB_res_new$phi1_hat
        phi2_all_new[clust_old] = EB_res_old$phi2_hat
        phi2_all_new[k] = EB_res_new$phi2_hat
        theta_all_new[clust_old] = EB_res_old$theta_hat
        theta_all_new[k] = EB_res_new$theta_hat
        tau_all_new[clust_old] = EB_res_old$tau_hat
        tau_all_new[k] = EB_res_new$tau_hat
        
        # compute log-prior ratio
        log_A = 0
        
        # compute log-proposal ratio
        log_P = 0
        
        # compute log-likelihood ratio
        log_L = log_like_new_c2 - log_like[clust_c1] - log_like[clust_c2]
        log_L = log_L + log_like_new[clust_old] + log_like_new[k] - log_like_old
        
        #acceptance probability
        acc_prob = min(0, log_A + log_P + log_L)
        acc_prob = exp(acc_prob)
        if(runif(1) < acc_prob){
          # accept
          cluster = membership_new
          phi1_all = phi1_all_new
          phi2_all = phi2_all_new
          theta_all = theta_all_new
          tau_all = tau_all_new
          log_like = log_like_new
          edge_status = getEdgeStatus(membership_new, graph0)
          csize = split_res$csize
          
          dist_mat[vid_merged, vid_merged] = MahalDist(coords[vid_merged, , drop=F], EB_res_merged$theta_hat, 
                                                       EB_res_merged$phi1_hat, EB_res_merged$phi2_hat)$dist_mat
          dist_mat[vid_old, vid_old] = MahalDist(coords[vid_old, , drop=F], EB_res_old$theta_hat, 
                                                 EB_res_old$phi1_hat, EB_res_old$phi2_hat)$dist_mat
          dist_mat[vid_new, vid_new] = MahalDist(coords[vid_new, , drop=F], EB_res_new$theta_hat, 
                                                 EB_res_new$phi1_hat, EB_res_new$phi2_hat)$dist_mat
          
          cnt_acc[move] = cnt_acc[move] + 1
        }
        
      }
    }
    
    if(move == 4) { ## Hyper move
      # update spanning tree
      stgraph = proposeMST(graph0, edge_status)
    }
    
    ## Update [sigmasq | tau, phi, theta, partition, tree, beta, lambda, y]
    # compute quadratic terms
    quad_all = evalQuad(Y, X, beta, dist_mat, cluster, k, rep(1, k), tau_all)
    # sample sigmasq
    shape = csize / 2 + a_sigma
    rate = b_sigma + quad_all$YPY / 2
    sigmasq_all = 1 / rgamma(k, shape, rate)
    
    ## Update [beta, lambda | tau, phi, partition, tree, sigmasq, y]
    # sample beta
    Q_beta = apply(quad_all$XPX / sigmasq_all, c(2, 3), sum)
    diag(Q_beta) = diag(Q_beta) + lambda
    b_beta = colSums(quad_all$XPY / sigmasq_all)
    
    # sample beta from N(Q^{-1}b, Q^{-1}), where Q = R'R
    R_beta = chol(Q_beta); rm(Q_beta)
    z_beta = rnorm(p, 0, 1)
    x_beta = forwardsolve(t(R_beta), b_beta)
    beta = backsolve(R_beta, z_beta + x_beta)
    rm(R_beta, b_beta, z_beta, x_beta)
    
    # sample lambda
    lambda = 1 / rgamma(1, a_lambda + p/2, b_lambda + sum(beta^2) / 2)
    # update loglike
    e = Y - X %*% beta
    log_like = updateLogLike(e, cluster, csize, k, quad_all$chol, hyper)
    
    
    ## save result
    if(iter > BURNIN & (iter - BURNIN) %% THIN == 0) {
      # sample (mean-zero) latent Gaussian field
      if(sample_latent) {
        latent_out[(iter-BURNIN)/THIN, ] = sampleLGP(e, dist_mat, cluster, k, rep(1, k), tau_all, sigmasq_all)
      }
      
      phi1_out[[(iter-BURNIN)/THIN]] = phi1_all
      phi2_out[[(iter-BURNIN)/THIN]] = phi2_all
      theta_out[[(iter-BURNIN)/THIN]] = theta_all
      tau_out[[(iter-BURNIN)/THIN]] = tau_all
      beta_out[(iter-BURNIN)/THIN, ] = beta
      lambda_out[(iter-BURNIN)/THIN] = lambda
      sigmasq_out[[(iter-BURNIN)/THIN]] = sigmasq_all
      ST_out[[(iter-BURNIN)/THIN]] = stgraph
      cluster_out[(iter-BURNIN)/THIN, ] = cluster
      
      log_post_out[(iter-BURNIN)/THIN] = evalLogPost(e, cluster, phi1_all, phi2_all, tau_all, sigmasq_all, beta, lambda, 
                                                     quad_all$chol, k, hyper)
    }
    if(iter %% 1000 == 0) cat('Iteration', iter, 'done\n')
  }
  
  mcmc_out = list('phi1_out' = phi1_out, 'phi2_out' = phi2_out, 'theta_out' = theta_out,
                  'tau_out' = tau_out, 'sigmasq_out' = sigmasq_out, 'beta_out' = beta_out,
                  'lambda_out' = lambda_out,
                  'ST_out' = ST_out, 'cluster_out' = cluster_out, 'log_post_out' = log_post_out,
                  'cnt_att' = cnt_att, 'cnt_acc' = cnt_acc,
                  'std_par_Y' = std_par_Y, 'Y_std' = Y)
  if(sample_latent) mcmc_out[['latent_out']] = latent_out

  if(EB_control$parallel) stopCluster(cl)
  
  return(mcmc_out)
}


########### Prediction Functions ###########

# function to do prediction for SPGP model
STGP_aniso_predict <- function(mcmc_res, X, X_ho, coords, coords_ho, pred_response = T, m = 1,
                         return_values = 'mean', nn_pred = 1) {
  require(FNN)
  
  Y = mcmc_res$Y_std
  std_par_Y = mcmc_res$std_par_Y
  cluster_out = mcmc_res$cluster_out
  phi1_out = mcmc_res$phi1_out
  phi2_out = mcmc_res$phi2_out
  theta_out = mcmc_res$theta_out
  tau_out = mcmc_res$tau_out
  lambda_out = mcmc_res$lambda_out
  sigmasq_out = mcmc_res$sigmasq_out
  beta_out = mcmc_res$beta_out
  if(!pred_response) latent_out = mcmc_res$latent_out
  
  # get nearest neighbor of hold-out dataset
  NNarray = get.knnx(coords, coords_ho, k = m)$nn.index
  
  n_ho = nrow(coords_ho); n = length(Y)
  n_post = nrow(beta_out)
  Y_ho = matrix(0, nrow = n_ho, ncol = n_post)
  cluster_ho_all = matrix(0, nrow = n_ho, ncol = n_post)
  for(i in 1:n_post) {
    cluster = cluster_out[i, ]
    # predict cluster_ho probabilistically
    nn_used = sample(m, n_ho, replace = T, prob = rep(1/m, m))
    nn_used = NNarray[cbind(1:n_ho, nn_used)]
    cluster_ho = cluster[nn_used]
    cluster_ho_all[, i] = cluster_ho
    
    k = max(cluster)
    phi1_uni = phi1_out[[i]]
    phi2_uni = phi2_out[[i]]
    theta_uni = theta_out[[i]]
    tau_uni = tau_out[[i]]
    lambda = lambda_out[i]
    sigmasq_uni = sigmasq_out[[i]]
    beta = beta_out[i, ]
    if(!pred_response) latent = latent_out[i, ]
    idx = split(1:n, cluster)
    idx_ho = split(1:n_ho, cluster_ho)
    
    for(j in names(idx_ho)) {
      phi1_j = phi1_uni[as.integer(j)]
      phi2_j = phi2_uni[as.integer(j)]
      theta_j = theta_uni[as.integer(j)]
      tau_j = tau_uni[as.integer(j)]
      sigmasq_j = sigmasq_uni[as.integer(j)]
      idx_j = idx[[j]]; idx_ho_j = idx_ho[[j]]
      Y_j = Y[idx_j]
      X_j = X[idx_j, , drop = F]; X_ho_j = X_ho[idx_ho_j, , drop = F]
      if(!pred_response) latent_j = latent[idx_j]
      n_j = length(idx_j)
      
      # get distance matrix
      dist_mat_j = MahalDist(coords[idx_j, , drop = F], theta_j, phi1_j, phi2_j)$dist_mat
      dist_mat_ho_j = MahalDist(coords_ho[idx_ho_j, , drop = F], theta_j, phi1_j, phi2_j)$dist_mat
      # get cross-distance matrix
      cdist_mat_j = MahalDist(coords_ho[idx_ho_j, , drop = F], theta_j, phi1_j, phi2_j,
                               coords[idx_j, , drop = F])$dist_mat
      
      # use exact covariance matrix
      if(pred_response) {
        # predict responses
        Y_ho[idx_ho_j, i] = predict_resp(Y_j, X_j, X_ho_j, dist_mat_j, cdist_mat_j, dist_mat_ho_j,
                                         1, tau_j, sigmasq_j, beta = beta)
      } else {
        # predict latent GP
        Y_ho[idx_ho_j, i] = predict_LGP(latent_j, dist_mat_j, cdist_mat_j, dist_mat_ho_j,
                                        1, sigmasq_j * tau_j)
      }
    }
  }
  
  if(return_values == 'mean') {
    Y_ho = rowMeans(Y_ho)
    return(list('Y_ho' = unstandardize(Y_ho, std_par_Y), 'cluster_ho' = cluster_ho_all))
  }
  if(return_values == 'all') {
    Y_ho = apply(Y_ho, 2, unstandardize, std_par_Y)
    return(list('Y_ho' = Y_ho, 'cluster_ho' = cluster_ho_all))
  }
}

# function to get kriging mean and variance
STGP_aniso_krigmv <- function(mcmc_res, X, X_ho, coords, coords_ho, pred_response = T, m = 1) {
  require(FNN)
  
  Y = mcmc_res$Y_std
  cluster_out = mcmc_res$cluster_out
  phi1_out = mcmc_res$phi1_out
  phi2_out = mcmc_res$phi2_out
  theta_out = mcmc_res$theta_out
  tau_out = mcmc_res$tau_out
  lambda_out = mcmc_res$lambda_out
  sigmasq_out = mcmc_res$sigmasq_out
  beta_out = mcmc_res$beta_out
  if(!pred_response) latent_out = mcmc_res$latent_out
  std_par_Y = mcmc_res$std_par_Y
  
  # get nearest neighbor of hold-out dataset
  NNarray = get.knnx(coords, coords_ho, k = m)$nn.index
  
  n_ho = nrow(coords_ho); n = length(Y)
  n_post = nrow(cluster_out)
  kmean = matrix(0, nrow = n_ho, ncol = n_post*m)
  kvar = matrix(0, nrow = n_ho, ncol = n_post*m)
  for(i in 1:n_post) {
    cluster = cluster_out[i, ]
    cluster_ho = matrix(cluster[NNarray], nrow = n_ho, ncol = m)
    
    k = max(cluster)
    phi1_uni = phi1_out[[i]]
    phi2_uni = phi2_out[[i]]
    theta_uni = theta_out[[i]]
    tau_uni = tau_out[[i]]
    sigmasq_uni = sigmasq_out[[i]]
    beta = beta_out[i, ]
    if(!pred_response) latent = latent_out[i, ]
    idx = split(1:n, cluster)
    idx_ho = split(1:(n_ho*m), c(cluster_ho))
    
    kmean_i = matrix(0, nrow = n_ho, ncol = m)
    kvar_i = matrix(0, nrow = n_ho, ncol = m)
    
    for(j in names(idx_ho)) {
      phi1_j = phi1_uni[as.integer(j)]
      phi2_j = phi2_uni[as.integer(j)]
      theta_j = theta_uni[as.integer(j)]
      tau_j = tau_uni[as.integer(j)]
      sigmasq_j = sigmasq_uni[as.integer(j)]
      idx_j = idx[[j]]
      Y_j = Y[idx_j]
      X_j = X[idx_j, , drop = F]
      if(!pred_response) latent_j = latent[idx_j]
      n_j = length(idx_j)
      
      # find row ids of cluster_ho that contains j
      idx_ho_j_all = idx_ho[[j]]
      idx_ho_j = unique(idx_ho_j_all %% n_ho)
      idx_ho_j[idx_ho_j == 0] = n_ho
      
      X_ho_j = X_ho[idx_ho_j, , drop = F]
      
      # get distance matrix
      dist_mat_j = MahalDist(coords[idx_j, , drop = F], theta_j, phi1_j, phi2_j)$dist_mat
      dist_mat_ho_j = MahalDist(coords_ho[idx_ho_j, , drop = F], theta_j, phi1_j, phi2_j)$dist_mat
      # get cross-distance matrix
      cdist_mat_j = MahalDist(coords_ho[idx_ho_j, , drop = F], theta_j, phi1_j, phi2_j,
                               coords[idx_j, , drop = F])$dist_mat
      
      # use exact covariance matrix
      if(pred_response) {
        # predict responses
        krig_res = predict_resp(Y_j, X_j, X_ho_j, dist_mat_j, cdist_mat_j, dist_mat_ho_j,
                                1, tau_j, sigmasq_j, beta = beta, kmv = T)
      } else {
        # predict latent GP
        krig_res = predict_LGP(latent_j, dist_mat_j, cdist_mat_j, dist_mat_ho_j,
                               1, sigmasq_j * tau_j, kmv = T)
      }
      
      for(ii in idx_ho_j_all) {
        jj = which(idx_ho_j == ifelse(ii %% n_ho == 0, n_ho, ii %% n_ho))
        # save krigging mean
        kmean_i[ii] = krig_res$kmean[jj]
        # save krigging variance
        kvar_i[ii] = krig_res$kvar[jj]
      }
    }
    kmean[, seq(i, (m-1)*n_post+i, n_post)] = unstandardize(c(kmean_i), std_par_Y)
    kvar[, seq(i, (m-1)*n_post+i, n_post)] = unstandardize(c(kvar_i), std_par_Y, nomean = T, s2 = T)
  }
  
  return(list('kmean' = kmean, 'kvar' = kvar))
}
