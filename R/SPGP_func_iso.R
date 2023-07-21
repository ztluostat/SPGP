#### Reproducible code for "A Nonstationary Soft Partitioned Gaussian Process Model via Random Spanning Trees"
#### Functions for SPGP model with isotropic covariance functions

library(igraph)
library(fields)
source("./BRISC/estimation.R")
source("./BRISC/prediction.R")

# function to get whether an edge is within a cluster or bewteen two clusters
getEdgeStatus <- function(membership, graph0) {
  #edge_status = character(ecount(graph0))
  inc_mat = get.edgelist(graph0, names = F)
  membership_head = membership[inc_mat[, 1]]
  membership_tail = membership[inc_mat[, 2]]
  edge_status = rep('w', ecount(graph0))
  edge_status[membership_head != membership_tail] = 'b'
  return(edge_status)
}

# function to split an existing cluster given MST
splitCluster <- function(mstgraph, edge_status, k, membership) {
  cluster = membership
  # candidate edges to cut
  ecand = E(mstgraph)[edge_status[E(mstgraph)$eid] == 'w']
  edge_cutted = ecand[sample.int(length(ecand), 1)]
  eid_cutted = edge_cutted$eid
  # endpoints of the cutted edge, note v1$vid > v2$vid
  v1 = head_of(mstgraph, edge_cutted)
  #v2 = tail_of(mstgraph, edge_cutted)
  clust_old = cluster[v1$vid]
  # vertices in the same cluster of v1
  idx_clust_old = (cluster == clust_old)
  vset = V(mstgraph)[idx_clust_old]
  mst_subgragh = induced_subgraph(mstgraph, vset)
  # delete edge to split cluster
  mst_subgragh = delete.edges(mst_subgragh, 
                              E(mst_subgragh)[eid == eid_cutted])
  
  connect_comp = components(mst_subgragh)
  cluster_new = connect_comp$membership
  cluster_v1 = cluster_new[V(mst_subgragh)$vid == v1$vid]
  idx_new = (cluster_new == cluster_v1)  # index for vertices belonging to new cluster
  vid_new = vset$vid[idx_new] # vid for vertices belonging to new cluster
  vid_old = vset$vid[!idx_new] # vid for vertices left in old cluster
  cluster[vid_new] = k + 1
  
  return(list(cluster = cluster, vid_new = vid_new,
              clust_old = clust_old, vid_old = vid_old))
}

splitCluster2 <- function(mstgraph, k, membership, csize, csize_obs_min, proj_size) { 
  cnt = 0
  while(TRUE) {
    clust_split = sample.int(k, 1, prob = csize - 1)
    edge_cutted = sample.int(csize[clust_split]-1, 1)
    
    mst_subgraph = induced_subgraph(mstgraph, membership == clust_split) 
    ## the above one can also be potentially improved if we also save/update a list with grouping info.   
    ## replace 'membership==clust.split' with 'group.list[[clust.split]]'
    
    mst_subgraph = delete.edges(mst_subgraph, edge_cutted)
    connect_comp = components(mst_subgraph)
    cluster_new = connect_comp$membership
    idx_new = cluster_new == 2  # index for vertices belonging to new cluster
    vid_new = (V(mst_subgraph)$vid)[idx_new]  # vid for vertices belonging to new cluster
    vid_old = (V(mst_subgraph)$vid)[!idx_new]  # vid for vertices left in old cluster
    
    # check cluster sizes
    csize_obs_new = sum(proj_size[vid_new])
    csize_obs_old = sum(proj_size[vid_old])
    cnt = cnt + 1
    if(min(csize_obs_old, csize_obs_new) >= csize_obs_min | cnt >= 100) break
  }
  
  # update cluster memberships
  membership[vid_new] = k + 1
  csize[clust_split] = length(vid_old)
  csize[k+1] = length(vid_new)
  
  return(list(cluster = membership, vid_new = vid_new, vid_old = vid_old,
              clust_old = clust_split, csize = csize))
}

# function to merge two existing clusters
mergeCluster <- function(mstgraph, edge_status, membership, csize) {
  # candidate edges for merging
  ecand = E(mstgraph)[edge_status[E(mstgraph)$eid] == 'b']
  #eid_merge = ecand$eid[sample.int(length(ecand$eid), 1)]
  edge_merge = ecand[sample.int(length(ecand), 1)]
  # update cluster information
  # endpoints of edge_merge, note v1$vid > v2$vid
  v1 = head_of(mstgraph, edge_merge); v2 = tail_of(mstgraph, edge_merge)
  # clusters that v1, v2 belonging to
  cluster = membership
  c1 = cluster[v1$vid]; c2 = cluster[v2$vid]
  idx_c1 = (cluster == c1)
  idx_c2 = (cluster == c2)
  
  # vid of vertices in c1
  vid_old = V(mstgraph)$vid[idx_c1]
  # vid of vertices in c2
  vid_new = V(mstgraph)$vid[idx_c2]
  
  # now drop c1
  cluster[idx_c1] = c2
  csize[c2] = csize[c1] + csize[c2]
  cluster[cluster > c1] = cluster[cluster > c1] - 1
  csize = csize[-c1]
  
  # clust_old: idx of dropped cluster (before merging)
  # clust_new: idx of cluster that is mergered into (before merging)
  return(list(cluster = cluster, vid_old = vid_old, vid_new = vid_new,
              clust_old = c1, clust_new = c2, csize = csize))
}


# function to propose a new MST
proposeMST <- function(graph0, edge_status) {
  nedge = length(edge_status)
  nb = sum(edge_status == 'b')
  nw = nedge - nb
  weight = numeric(nedge)
  weight[edge_status == 'w'] = runif(nw, 0, 0.5)
  weight[edge_status == 'b'] = runif(nb, 0.5, 1)
  mstgraph = mst(graph0, weights = weight)
  return(mstgraph)
}

# function to evaluate covariance matrix (including nugget)
evalCov <- function(dist_mat, phi, tau, nugget = T) {
  scaled_dist = dist_mat / phi
  cov_mat = (1 + sqrt(5) * scaled_dist + 5 / 3 * scaled_dist ^ 2) * exp(-sqrt(5) * scaled_dist) 
  cov_mat = tau * cov_mat
  if (nugget)
    diag(cov_mat) = diag(cov_mat) + 1
  return(cov_mat)
}

# function to evaluate lower Cholesky of covariance matrix (including nugget)
evalChol <- function(dist_mat, phi, tau, nugget = T) {
  cov_mat = evalCov(dist_mat, phi, tau, nugget)
  chol = t(chol(cov_mat))
  return(chol)
}

# function to get log likelihood (up to a constant) in each cluster
# where Cov(Y) = tau*Sigma(phi) + I + lambda*J
# Y may be interpreted as residuals here
evalLogLike <- function(Y, dist_mat, cluster, csize, k, phi, tau, lambda, hyper) {
  n = length(Y)
  if(k > 1) {
    idx = split(1:n, cluster)
  } else {  # k == 1 case
    idx = list(1:n)
  }
  a_sigma = hyper['a_sigma']; b_sigma = hyper['b_sigma']
  
  log_like = numeric(k)
  for(j in 1:k) {
    idx_j = idx[[j]]; n_j = csize[j]
    Y_j = Y[idx_j]
    phi_j = phi[j]; tau_j = tau[j]; lambda_j = lambda[j]
    
    scaled_dist_j = dist_mat[idx_j, idx_j, drop=F] / phi_j
    cov_mat_j = (1 + sqrt(5) * scaled_dist_j + 5 / 3 * scaled_dist_j ^ 2) * exp(-sqrt(5) * scaled_dist_j) 
    cov_mat_j = tau_j * cov_mat_j + lambda_j
    diag(cov_mat_j) = diag(cov_mat_j) + 1
    chol_j = t(chol(cov_mat_j))
    sol = forwardsolve(chol_j, Y_j)
    quad = sum(sol ^ 2)
    
    log_like_j = -sum(log(diag(chol_j))) - n_j/2*log(2*pi) + lgamma(n_j/2 + a_sigma)
    log_like_j = log_like_j - (n_j/2 + a_sigma) * log(b_sigma + quad/2)
    
    log_like[j] = log_like_j
  }
  
  log_like = log_like + a_sigma*log(b_sigma) - lgamma(a_sigma)
  return(log_like)
}

# function to get log likelihood from log posterior for GP in a cluster
# Y may be interpreted as residuals here
# For half-t prior, a = df, b = scale
likeFromGPPost <- function(log_post, phi, tau, hyper) {
  a_phi = hyper['a_phi']; b_phi = hyper['b_phi']
  a_tau = hyper['a_tau']; b_tau = hyper['b_tau']
  
  log_prior = dtplus(phi, a_phi, b_phi, log = T)
  log_prior = log_prior + dinvgamma(tau, a_tau, b_tau, log = T)
  return(log_post - log_prior)
}

# function to get log posterior for GP (up to a constant) in a cluster
# Y may be interpreted as residuals here
# cov_param = c(log(phi), log(tau))
# For half-t prior, a = df, b = scale
evalLogGPPost <- function(cov_param, Y, dist_mat, hyper, log_like = NULL) {
  phi = exp(cov_param[1]); tau = exp(cov_param[2])
  a_phi = hyper['a_phi']; b_phi = hyper['b_phi']
  a_tau = hyper['a_tau']; b_tau = hyper['b_tau']
  
  n = length(Y)
  if (is.null(log_like))
    log_like = evalLogLike(Y, dist_mat, rep(1, n), n, 1, phi, tau, 0, hyper)
  log_prior = dtplus(phi, a_phi, b_phi, log = T)
  log_prior = log_prior + dinvgamma(tau, a_tau, b_tau, log = T)
  return(log_like + log_prior)
}

# function to get gradient of log posterior for GP (up to a constant) in a cluster
# Y may be interpreted as residuals here
# cov_param = c(log(phi), log(tau))
evalLogGPGrad <- function(cov_param, Y, dist_mat, hyper) {
  phi = exp(cov_param[1]); tau = exp(cov_param[2])
  a_sigma = hyper['a_sigma']; b_sigma = hyper['b_sigma']
  a_phi = hyper['a_phi']; b_phi = hyper['b_phi']
  a_tau = hyper['a_tau']; b_tau = hyper['b_tau']
  n = length(Y)
  
  # compute covariance matrix and deriatives
  scaled_dist = dist_mat / phi
  exp_dist = exp(-sqrt(5) * scaled_dist)
  cov_mat = (1 + sqrt(5) * scaled_dist + 5 / 3 * scaled_dist ^ 2) * exp_dist
  cov_mat = tau * cov_mat
  diag(cov_mat) = diag(cov_mat) + 1
  # deriavtive of cov_mat wrt phi
  dcov_dphi = tau * 5/3 * scaled_dist ^ 2 * (sqrt(5) * scaled_dist + 1) * exp_dist
  
  # inverse cov_mat
  cov_inv = base::chol2inv( base::chol(cov_mat) )
  cov_inv_y = cov_inv %*% Y  # inv(cov_mat) %*% Y
  quad = crossprod(Y, cov_inv_y) # t(Y) %*% inv(cov_mat) %*% Y
  
  # compute gradient wrt phi
  grad_phi = -0.5 * sum(t(cov_inv) * dcov_dphi)
  grad_phi = grad_phi + (n / 2 + a_sigma) / (0.5 * quad + b_sigma) * 0.5 * t(cov_inv_y) %*% dcov_dphi %*% cov_inv_y
  grad_phi = grad_phi - (a_phi + 1) * phi^2 / (a_phi * b_phi^2 + phi^2)
  
  # compute gradient of tau
  grad_tau = -0.5 * (n - sum(diag(cov_inv)))
  grad_tau = grad_tau + (n / 2 + a_sigma) / (0.5 * quad + b_sigma) * 0.5 * (quad - crossprod(cov_inv_y))
  grad_tau = grad_tau - (a_tau + 1) + b_tau / tau
  
  return(c(grad_phi, grad_tau))
}

# function to get log GP posterior for a cluster using empirical Bayes
# Y may be interpreted as residuals here
evalLogPostEB <- function(Y, dist_mat, hyper, lphi_init = 0, ltau_init = 0,
                          NNGP = F, coords = NULL, sigmasq = NULL, min_n_nngp = 300) {
  a_phi = hyper['a_phi']; b_phi = hyper['b_phi']
  a_tau = hyper['a_tau']; b_tau = hyper['b_tau']
  
  if (!NNGP | length(Y) < min_n_nngp) { # use exact calculation
    # try given initial values first
    init = c(lphi_init, ltau_init)
    EB_res = try(optim(init, fn = evalLogGPPost, gr = evalLogGPGrad, Y = Y, dist_mat = dist_mat, hyper = hyper,
                       method = 'L-BFGS-B', lower = rep(-10, 2), upper = rep(10, 2), 
                       control = list(fnscale = -1, maxit = 500)), silent = F)
    
    if(is(EB_res, 'try-error')) {
      failed = T
    } else if (EB_res$convergence != 0) {
      failed = T
    } else {
      failed = F
    }
    
    # if it fails, try default initial values (0, 0)
    if(failed & (lphi_init != 0 | ltau_init != 0)) {
      EB_res = try(optim(c(0, 0), fn = evalLogGPPost, gr = evalLogGPGrad, Y = Y, dist_mat = dist_mat, hyper = hyper,
                         method = 'L-BFGS-B', lower = rep(-10, 2), upper = rep(10, 2), 
                         control = list(fnscale = -1, maxit = 500)), silent = F)
      
      if(is(EB_res, 'try-error')) {
        failed = T
      } else if (EB_res$convergence != 0) {
        failed = T
      } else {
        failed = F
      }
    }
    
    # if this fails, return log_post = -Inf
    if(failed)
      return(list('log_post' = -Inf, 'phi_hat' = NA, 'tau_hat' = NA))
    
    phi_hat = exp(EB_res$par[1]); tau_hat = exp(EB_res$par[2])
    # log_post = EB_res$value - dtplus(phi_hat, a_phi, b_phi, log = T) - dinvgamma(tau_hat, a_tau, b_tau, log = T)
    log_post = EB_res$value
    
    # cat("EB done\n")
    return(list('log_post' = log_post, 'phi_hat' = phi_hat, 'tau_hat' = tau_hat,
                'method' = 'exact'))
    
  } else { # use NNGP approximation
    
    stopifnot(!is.null(coords) & !is.null(sigmasq))
    prior_param = hyper[c('a_phi', 'b_phi', 'a_tau', 'b_tau', 'a_sigma', 'b_sigma')]
    nngp_res = NNGP_estimation(coords, Y, cov.model = "matern52", n.neighbors = hyper['nngp'],
                               sigma.sq = sigmasq * exp(ltau_init), 
                               tau.sq = sigmasq, 
                               phi = 1 / exp(lphi_init),
                               prior = prior_param, verbose = F)
    
    phi_hat = 1 / nngp_res$Theta[2]
    tau_hat = nngp_res$Theta[1]
    prec_param = list('B' = nngp_res$BRISC_Object$B, 'F' = nngp_res$BRISC_Object$F,
                      'nn_indx' = nngp_res$BRISC_Object$nnIndx, 'nn_indx_LU' = nngp_res$BRISC_Object$nnIndxLU,
                      'ord' = nngp_res$ord)
    # cat("EB done\n")
    return(list('log_post' = nngp_res$log_post, 'phi_hat' = phi_hat, 'tau_hat' = tau_hat,
                'prec_param' = prec_param, 'method' = 'NNGP'))
  }
}

# function for removing j-th column and row from t(R) %*% R
cholDel <- function(R, j) {
  n = nrow(R)
  if(j == n) {
    R_new = R[-n, -n, drop = F]
  } else if(j == 1) {
    R_new = t(ramcmc::chol_update(t(R[-1, -1, drop = F]), R[1, -1]))
  } else {
    R_new = matrix(0, nrow = n-1, ncol = n-1)
    R_new[1:(j-1), ] = R[1:(j-1), -j, drop = F]
    R_new[j:(n-1), j:(n-1)] = t(ramcmc::chol_update(t(R[(j+1):n, (j+1):n, drop = F]), 
                                                    R[j, (j+1):n]))
  }
  return(R_new)
}

# function for adding (A1, A2, A3) as j-th column and row to t(R) %*% R
# note that A2 is a scalar
cholAdd <- function(R, j, A2, A1 = NULL, A3 = NULL) {
  n = nrow(R)
  if(j == n+1) {
    R_new = rbind(R, 0)
    S12 = drop(backsolve(R, A1, transpose = T))
    S22 = sqrt(as.numeric(A2 - sum(S12^2)))
    R_new = cbind(R_new, c(S12, S22))
  } else {
    R_new = matrix(0, nrow = n+1, ncol = n+1)
    if(j > 1) {
      R11 = R[1:(j-1), 1:(j-1), drop = F]
      R_new[1:(j-1), 1:(j-1)] = R11
      S12 = backsolve(R11, A1, transpose = T)
      R_new[1:(j-1), j] = S12
      S13 = R[1:(j-1), j:n, drop = F]
      R_new[1:(j-1), (j+1):(n+1)] = S13
      
      S22 = sqrt(as.numeric(A2 - sum(S12^2)))
      R_new[j, j] = S22
      S23 = (t(A3) - crossprod(S12, S13)) / S22
      R_new[j, (j+1):(n+1)] = S23
      S33 = t(ramcmc::chol_downdate(t(R[j:n, j:n, drop = F]), S23))
      R_new[(j+1):(n+1), (j+1):(n+1)] = S33
    } else {
      S22 = sqrt(as.numeric(A2))
      R_new[1, 1] = S22
      S23 = as.numeric(A3) / S22
      R_new[1, 2:(n+1)] = S23
      S33 = t(ramcmc::chol_downdate(t(R), S23))
      R_new[2:(n+1), 2:(n+1)] = S33
    }
  }
  return(R_new)
}

# function to update log-likelihood using Cholesky/NNGP factors
# Y may be interpreted as residuals here
updateLogLike <- function(Y, cluster, csize, k, cov_factors, hyper) {
  n = length(Y)
  if(k > 1) {
    idx = split(1:n, cluster)
  } else {  # k == 1 case
    idx = list(1:n)
  }
  a_sigma = hyper['a_sigma']; b_sigma = hyper['b_sigma']
  constant = a_sigma*log(b_sigma) - lgamma(a_sigma)
  
  log_like = numeric(k)
  for(j in 1:k) {
    idx_j = idx[[j]]; n_j = csize[j]
    Y_j = Y[idx_j]
    
    if (length(cov_factors[[j]]) == 1) {
      # use Cholesky of exact covariance
      chol_j = cov_factors[[j]]$chol
      sol = forwardsolve(chol_j, Y_j)
      quad = sum(sol ^ 2)
      log_det = sum(log(diag(chol_j)))
      
    } else {
      # use NNGP approximation
      Y_j = Y_j[cov_factors[[j]]$ord]
      quad = evalQuadNNGP(cov_factors[[j]][['B']], cov_factors[[j]][['F']], Y_j, Y_j, 
                          cov_factors[[j]]$nn_indx, cov_factors[[j]]$nn_indx_LU)
      log_det = sum(log(cov_factors[[j]][['F']])) / 2
    }
    
    log_like_j = -log_det - n_j/2*log(2*pi) + lgamma(n_j/2 + a_sigma)
    log_like_j = log_like_j - (n_j/2 + a_sigma) * log(b_sigma + quad/2)
    log_like[j] = log_like_j
  }
  log_like = log_like + constant
  return(log_like)
}

# function to get quadratic forms for each cluster using Cholesky/NNGP factors
evalQuad <- function(Y, X, beta, dist_mat, cluster, k, phi, tau, cov_factors) {
  n = length(Y); p = ncol(X)
  idx = split(1:n, cluster)
  nu = 0.5
  e = Y - X %*% beta
  YPY_all = numeric(k)   # (Y-mu)'P^{-1}(Y-mu), where P = cov(Y) 
  XPY_all = matrix(0, nrow = k, ncol = p)  # X'P^{-1}Y
  XPX_all = array(0, dim = c(k, p, p))  # X'P^{-1}X
  
  for(j in 1:k) {
    # cat("Computing quadratic j \n")
    idx_j = idx[[j]]
    Y_j = Y[idx_j]; e_j = e[idx_j]; n_j = length(Y_j)
    X_j = X[idx_j, , drop = F]
    phi_j = phi[j]; tau_j = tau[j]
    
    use_NNGP = length(cov_factors[[j]]) > 1
    if (!use_NNGP) {
      # use Cholesky of exact covariance
      chol_j = cov_factors[[j]]$chol
      
      # compute quadratic forms
      sol_e = forwardsolve(chol_j, e_j)
      sol_Y = forwardsolve(chol_j, Y_j)
      sol_X = forwardsolve(chol_j, X_j)
      YPY = sum(sol_e ^ 2)
      XPY = crossprod(sol_X, sol_Y)
      XPX = crossprod(sol_X)
    } else {
      # use NNGP approximation
      XPY = numeric(p)
      XPX = matrix(0, nrow = p, ncol = p)
      
      # reordering
      e_j = e_j[cov_factors[[j]]$ord]
      Y_j = Y_j[cov_factors[[j]]$ord]
      X_j = X_j[cov_factors[[j]]$ord, , drop=F]
      
      YPY = evalQuadNNGP(cov_factors[[j]][['B']], cov_factors[[j]][['F']], e_j, e_j, 
                         cov_factors[[j]]$nn_indx, cov_factors[[j]]$nn_indx_LU)
      for (m in 1:p) {
        XPY[m] = evalQuadNNGP(cov_factors[[j]][['B']], cov_factors[[j]][['F']], X_j[, m, drop=F], Y_j, 
                              cov_factors[[j]]$nn_indx, cov_factors[[j]]$nn_indx_LU)
        for (m2 in 1:m)
          XPX[m, m2] = evalQuadNNGP(cov_factors[[j]][['B']], cov_factors[[j]][['F']], X_j[, m, drop=F], X_j[, m2, drop=F], 
                                    cov_factors[[j]]$nn_indx, cov_factors[[j]]$nn_indx_LU)
      }
      XPX = (XPX + t(XPX)) / 2
    }
    
    YPY_all[j] = YPY
    XPY_all[j, ] = XPY
    XPX_all[j, , ] = XPX
  }
  return(list('YPY' = YPY_all, 'XPY' = XPY_all, 'XPX' = XPX_all))
}

# function to sample (mean-zero) latent Gaussian field
# Y may be interpreted as residuals here
sampleLGP <- function(Y, dist_mat, cluster, k, phi, tau, sigmasq) {
  n = length(Y)
  idx = split(1:n, cluster)
  latent = numeric(n)
  
  for(j in 1:k) {
    idx_j = idx[[j]]
    Y_j = Y[idx_j]; n_j = length(Y_j)
    dist_mat_j = dist_mat[idx_j, idx_j, drop=F]
    phi_j = phi[j]; tau_j = tau[j]; sigmasq_j = sigmasq[j]
    
    # compute exact covariance matrix
    scaled_dist_j = dist_mat_j / phi_j
    cov_mat_j = (1 + sqrt(5) * scaled_dist_j + 5 / 3 * scaled_dist_j ^ 2) * exp(-sqrt(5) * scaled_dist_j)
    
    Q = base::chol2inv( base::chol(cov_mat_j) ) / tau_j
    diag(Q) = diag(Q) + 1
    L = t(chol(Q)) / sqrt(sigmasq_j)
    b = Y_j / sigmasq_j
    
    # sample from N(Q^{-1}b, Q^{-1}), where Q = LL'
    R = t(L)
    z = rnorm(n_j, 0, 1)
    x = Matrix::solve(L, b, system = 'L')
    latent[idx_j] = Matrix::solve(R, z + x)
  }
  return(latent)
}

# function to compute deviance information criteria
computeDIC <- function(log_like, map_idx, penalty = 2) {
  # compute deviance at posterior mode
  deviance_hat = -2 * log_like[map_idx]
  
  p_DIC1 = 2 * (log_like[map_idx] - mean(log_like))
  p_DIC2 = 2 * var(log_like)
  DIC1 = deviance_hat + penalty * p_DIC1
  DIC2 = deviance_hat + penalty * p_DIC2
  
  return(c('DIC1' = DIC1, 'DIC2' = DIC2, 'p_DIC1' = p_DIC1, 'p_DIC2' = p_DIC2))
}

# function to get log full posterior density (up to a constant)
evalLogFullPost <- function(e, cluster, phi_all, tau_all, sigmasq_all, beta, lambda, cov_factors, k, hyper) {
  a_phi = hyper['a_phi']; b_phi = hyper['b_phi']
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
    
    use_NNGP = length(cov_factors[[j]]) > 1
    if (!use_NNGP) {
      # use Cholesky of exact covariance
      chol_j = cov_factors[[j]]$chol
      # compute quadratic term
      sol = forwardsolve(chol_j, e_j)
      quad = sum(sol ^ 2)
      # compute log determinant
      log_det = sum(log(diag(chol_j)))
    } else {
      # use NNGP approximation
      e_j = e_j[cov_factors[[j]]$ord]
      quad = evalQuadNNGP(cov_factors[[j]][['B']], cov_factors[[j]][['F']], e_j, e_j, 
                          cov_factors[[j]]$nn_indx, cov_factors[[j]]$nn_indx_LU)
      log_det = sum(log(cov_factors[[j]][['F']])) / 2
    }
    
    log_like = log_like - log_det - n_j/2*log(sigmasq_j)
    log_like = log_like - quad / (2*sigmasq_j)
  }
  
  log_prior = sum(dtplus(phi_all, a_phi, b_phi, log = T))
  log_prior = log_prior + sum(dinvgamma(tau_all, a_tau, b_tau, log = T))
  log_prior = log_prior + sum(dinvgamma(sigmasq_all, a_sigma, b_sigma, log = T))
  log_prior = log_prior + sum(dnorm(beta, 0, sqrt(lambda), log = T))
  log_prior = log_prior + dinvgamma(lambda, a_lambda, b_lambda, log = T)
  log_prior = log_prior - lchoose(n-1, k-1) + k*log(1-c)
  
  return(list('log_post' = log_like + log_prior, 'log_like' = log_like))
}

# function to evaluate density of inver gamma
dinvgamma <- function(x, shape, rate = 1, log = FALSE) {
  log_f = dgamma(1/x, shape, rate, log = TRUE) - 2*log(x)
  if(log) return(log_f)
  return(exp(log_f))
}

# function to evaluate density of half-t
# a = df, b = scale
dtplus <- function(x, df, scale = 1, log = FALSE) {
  f = dt(x / scale, df, log = log)
  if(log) return(f + log(2) - log(scale))
  return(2 * f / scale)
}

# function to sample a random projection
randomProject <- function(nearest_knots, proj_prob, nn = 1) {
  n = nrow(nearest_knots)
  nn_idx = apply(proj_prob, 1, function(prob) sample.int(nn, 1, prob = prob))
  proj = nearest_knots[ cbind(1:n, nn_idx) ]
  return(proj)
}

# function to standardize Y to have mean zero and range one
standardize <- function(x) {
  res = c('max' = max(x), 'min' = min(x))
  x = x / (res['max'] - res['min'])
  res['mean'] = mean(x)
  x = x - res['mean']
  
  return(list(x = x, std_par = res))
}

# function to unstandardize Y
unstandardize <- function(x, std_par, nomean = F, s2 = F) {
  if(!nomean) x = x + std_par['mean']
  if(s2) {
    x = x * (std_par['max'] - std_par['min'])^2
  } else {
    x = x * (std_par['max'] - std_par['min'])
  }
  return(x)
}


########### MCMC Function ###########

STGP_fit <- function(Y, X, coords, coords_knot, graph0, init_val, hyper, tunning,
                     MCMC, BURNIN, THIN, seed = 1234, NNGP = F,
                     sample_latent = F, standardizeY = T, save_tree = F, fix_proj = F) {
  set.seed(seed)
  
  n = length(Y)
  n_knot = nrow(coords_knot)
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
  phi_all = init_val[['phi']]
  sigmasq_all = init_val[['sigmasq']]
  lambda = init_val[['lambda']]
  beta = init_val[['beta']]
  k = max(cluster)  # number of clusters
  csize = as.numeric(table(cluster))  # cluster sizes
  e = Y - X %*% beta
  
  # hyper-parameter
  c = hyper['c']
  a_tau = hyper['a_tau']; b_tau = hyper['b_tau']
  a_phi = hyper['a_phi']; b_phi = hyper['b_phi']
  a_sigma = hyper['a_sigma']; b_sigma = hyper['b_sigma']
  a_lambda = hyper['a_lambda']; b_lambda = hyper['b_lambda']
  nn = hyper['nn']
  
  # tunning paramter
  rc = 0.19; rhy = 0.01
  csize_obs_min = tunning['csize_min']  # minimum cluster size
  min_n_nngp = tunning['min_n_nngp']  # minimum cluster size for using NNGP approximation
  
  # compute distance matrix
  dist_mat = as.matrix(dist(coords))
  
  # compute nearest knots and projection probability
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
  proj = init_val[['proj']]
  cluster_obs = cluster[proj]
  csize_obs = as.integer(table(cluster_obs))
  # how many obs. are projected to the i-th knot
  proj_size = as.integer(table(factor(proj, levels = 1:n_knot)))
  
  ## get log-likelihood in each cluster plus the prior for phi and tau 
  ## and Cholesky/NNGP factors of covariance matrix (including nugget) in each cluster
  
  cov_factors = list()
  for (j in 1:k) {
    idx_j = (cluster_obs == j)
    cov_factors[[j]] = list('chol' = evalChol(dist_mat[idx_j, idx_j, drop=F], phi_all[j], tau_all[j]))
  }
  
  # log_post: vector of length k
  log_post = updateLogLike(e, cluster_obs, csize_obs, k, cov_factors, hyper)
  log_post = log_post + dtplus(phi_all, a_phi, b_phi, log = T)
  log_post = log_post + dinvgamma(tau_all, a_tau, b_tau, log = T)
  
  ## whether an edge in graph0 is within a cluster or bewteen two clusters
  # n*p matrix
  edge_status = getEdgeStatus(cluster, graph0)
  
  ## add warm-up steps
  WARMUP = 0
  MCMC = MCMC + WARMUP
  BURNIN = BURNIN + WARMUP
  
  ## MCMC results
  phi_out = list(); tau_out = list()
  sigmasq_out = list()
  lambda_out = numeric((MCMC-BURNIN)/THIN)
  beta_out = matrix(0, nrow = (MCMC-BURNIN)/THIN, ncol = p)
  latent_out = matrix(0, nrow = (MCMC-BURNIN)/THIN, ncol = n)
  cluster_out = matrix(0, nrow = (MCMC-BURNIN)/THIN, ncol = n_knot)
  cluster_obs_out = matrix(0, nrow = (MCMC-BURNIN)/THIN, ncol = n)
  ST_out = list()
  log_post_out = numeric((MCMC-BURNIN)/THIN)
  log_like_out = numeric((MCMC-BURNIN)/THIN)
  
  cnt_acc = rep(0, 4)
  cnt_att = rep(0, 4)
  
  for(iter in 1:MCMC) {
    ## Update [tau, phi, partition, tree | beta, lambda, y]
    
    if(k == 1) {
      rb = 0.99; rd = 0; rc = 0
    } else if(k == n_knot) {
      rb = 0; rd = 0.8; rc = 0.19
    } else {
      rb = 0.4; rd = 0.4; rc = 0.19
    }
    if (iter > WARMUP) {
      move = sample(4, 1, prob = c(rb, rd, rc, rhy))
    } else {
      move = 4
    }
    cnt_att[move] = cnt_att[move] + 1
    
    if(move == 1) { ## Birth move
      # split an existing cluster
      split_res = splitCluster2(stgraph, k, cluster, csize, csize_obs_min, proj_size)
      cluster_new = split_res$cluster
      clust_old = split_res$clust_old
      
      # get memberships of observations
      cluster_obs_new = cluster_new[proj]
      id_new = which(cluster_obs_new == k + 1)  # observations assigned to new cluster
      id_old = which(cluster_obs_new == clust_old)  # observations left in old cluster
      
      lphi = log(phi_all[clust_old])
      ltau = log(tau_all[clust_old])
      sigmasq = sigmasq_all[clust_old]
      
      # reject if we propose a too small cluster
      if(min(split_res$csize) >= csize_obs_min) {
        
        # compute log-prior ratio
        log_A = log(1-c)
        
        # compute log-proposal ratio
        if(k == n_knot - 1) {
          rd_new = 0.8
        } else {
          rd_new = 0.4
        }
        log_P = log(rd_new) - log(rb)
        
        # compute log-likelihood ratio
        log_post_new = log_post
        EB_res_old = evalLogPostEB(e[id_old], dist_mat[id_old, id_old, drop=F], hyper, lphi, ltau,
                                   NNGP, coords[id_old, ], sigmasq, min_n_nngp)
        log_post_new[clust_old] = EB_res_old$log_post
        EB_res_new = evalLogPostEB(e[id_new], dist_mat[id_new, id_new, drop=F], hyper, lphi, ltau,
                                   NNGP, coords[id_new, ], sigmasq, min_n_nngp)
        log_post_new[k + 1] = EB_res_new$log_post
        log_L = log_post_new[clust_old] + log_post_new[k + 1] - log_post[clust_old]
        
        # acceptance probability
        acc_prob = min(0, log_A + log_P + log_L)
        acc_prob = exp(acc_prob)
        if(runif(1) < acc_prob){
          # accept
          cluster = cluster_new
          cluster_obs = cluster_obs_new
          phi_all[clust_old] = EB_res_old$phi_hat
          phi_all[k + 1] = EB_res_new$phi_hat
          tau_all[clust_old] = EB_res_old$tau_hat
          tau_all[k + 1] = EB_res_new$tau_hat
          k = k + 1
          log_post = log_post_new
          edge_status = getEdgeStatus(cluster_new, graph0)
          csize = split_res$csize
          csize_obs[k] = length(id_new)
          csize_obs[clust_old] = length(id_old)
          
          # update Cholesky/NNGP factors
          if (EB_res_old$method == 'exact') {
            cov_factors[[clust_old]] = list('chol' = evalChol(dist_mat[id_old, id_old, drop=F], 
                                                              EB_res_old$phi_hat, EB_res_old$tau_hat))
          } else {
            cov_factors[[clust_old]] = EB_res_old$prec_param
          }
          if (EB_res_new$method == 'exact') {
            cov_factors[[k]] = list('chol' = evalChol(dist_mat[id_new, id_new, drop=F], 
                                                      EB_res_new$phi_hat, EB_res_new$tau_hat))
          } else {
            cov_factors[[k]] = EB_res_new$prec_param
          }
          
          cnt_acc[move] = cnt_acc[move] + 1
        }
      }
    }
    
    if(move == 2) { ## Death move
      # merge two existing clusters (c1, c2) -> c2
      merge_res = mergeCluster(stgraph, edge_status, cluster, csize)
      cluster_new = merge_res$cluster
      # c1 is merged into c2
      clust_c1 = merge_res$clust_old  # cluster indices are the ones before merging
      clust_c2 = merge_res$clust_new
      
      if(csize_obs[clust_c1] > csize_obs[clust_c2]) {
        lphi = log(phi_all[clust_c1])
        ltau = log(tau_all[clust_c1])
        sigmasq = sigmasq_all[clust_c1]
      } else {
        lphi = log(phi_all[clust_c2])
        ltau = log(tau_all[clust_c2])
        sigmasq = sigmasq_all[clust_c2]
      }
      
      # compute log-prior ratio
      log_A = -log(1-c)
      
      # compute log-proposal ratio
      if(k == 2) {
        rb_new = 0.99
      } else {
        rb_new = 0.4
      }
      log_P = log(rb_new) - log(rd)
      
      # compute log-likelihood ratio
      # cluster_obs_new = cluster_new[proj]
      id_new = which(cluster_obs == clust_c1 | cluster_obs == clust_c2)  # observations in merged cluster
      log_post_new = log_post
      # log-likelihood of merged cluster
      EB_res_merged = evalLogPostEB(e[id_new], dist_mat[id_new, id_new, drop=F], hyper, lphi, ltau,
                                    NNGP, coords[id_new, ], sigmasq, min_n_nngp)
      log_post_new_c2 = EB_res_merged$log_post
      log_post_new[clust_c2] = log_post_new_c2
      log_post_new = log_post_new[-clust_c1]
      log_L = log_post_new_c2 - log_post[clust_c1] - log_post[clust_c2]
      
      #acceptance probability
      acc_prob = min(0, log_A + log_P + log_L)
      acc_prob = exp(acc_prob)
      if(runif(1) < acc_prob) {
        # accept
        cluster = cluster_new
        cluster_obs = cluster_new[proj]
        k = k - 1
        phi_all[clust_c2] = EB_res_merged$phi_hat
        phi_all = phi_all[-clust_c1]
        tau_all[clust_c2] = EB_res_merged$tau_hat
        tau_all = tau_all[-clust_c1]
        log_post = log_post_new
        edge_status = getEdgeStatus(cluster_new, graph0)
        csize = merge_res$csize
        csize_obs[clust_c2] = length(id_new)
        csize_obs = csize_obs[-clust_c1]
        
        # update Cholesky/NNGP factors
        if (EB_res_merged$method == 'exact') {
          cov_factors[[clust_c2]] = list('chol' = evalChol(dist_mat[id_new, id_new, drop=F], 
                                                           EB_res_merged$phi_hat, EB_res_merged$tau_hat))
        } else {
          cov_factors[[clust_c2]] = EB_res_merged$prec_param
        }
        cov_factors = cov_factors[-clust_c1]
        
        cnt_acc[move] = cnt_acc[move] + 1
      }
    }
    
    if(move == 3) { ## change move
      # first perform death move: (c1, c2) -> c2
      merge_res = mergeCluster(stgraph, edge_status, cluster, csize)
      # then perform birth move
      split_res = splitCluster2(stgraph, k-1, merge_res$cluster, merge_res$csize, csize_obs_min, proj_size)
      
      # update parameters for merge
      clust_c1 = merge_res$clust_old # cluster indices are the ones before merging
      clust_c2 = merge_res$clust_new
      
      if(csize_obs[clust_c1] > csize_obs[clust_c2]) {
        lphi = log(phi_all[clust_c1])
        ltau = log(tau_all[clust_c1])
        sigmasq = sigmasq_all[clust_c1]
      } else {
        lphi = log(phi_all[clust_c2])
        ltau = log(tau_all[clust_c2])
        sigmasq = sigmasq_all[clust_c2]
      }
      
      log_post_new = log_post
      id_merged = which(cluster_obs == clust_c1 | cluster_obs == clust_c2)  # obs in merged cluster
      EB_res_merged = evalLogPostEB(e[id_merged], dist_mat[id_merged, id_merged, drop=F], hyper, lphi, ltau,
                                    NNGP, coords[id_merged, , drop=F], sigmasq, min_n_nngp)
      log_post_new_c2 = EB_res_merged$log_post
      log_post_new[clust_c2] = log_post_new_c2
      log_post_new = log_post_new[-clust_c1]
      
      # reject if EB_res_merged$log_post == -Inf or we propose a too small cluster
      if(log_post_new_c2 > -Inf & min(split_res$csize) >= csize_obs_min) {
        
        phi_all_new = phi_all; tau_all_new = tau_all
        phi_all_new[clust_c2] = EB_res_merged$phi_hat
        phi_all_new = phi_all_new[-clust_c1]
        tau_all_new[clust_c2] = EB_res_merged$tau_hat
        tau_all_new = tau_all_new[-clust_c1]
        
        csize_obs_new = csize_obs
        csize_obs_new[clust_c2] = length(id_merged)
        csize_obs_new = csize_obs_new[-clust_c1]
        
        # update parameters for split
        cluster_new = split_res$cluster
        clust_old = split_res$clust_old
        
        # get memberships of observations
        cluster_obs_new = cluster_new[proj]
        id_new = which(cluster_obs_new == k)  # observations assigned to new cluster
        id_old = which(cluster_obs_new == clust_old)  # observations left in old cluster
        
        lphi = log(phi_all_new[clust_old])
        ltau = log(tau_all_new[clust_old])
        sigmasq = sigmasq_all[clust_old]
        
        log_post_old = log_post_new[clust_old]  # log-likelihood of the cluster before splitting
        EB_res_old = evalLogPostEB(e[id_old], dist_mat[id_old, id_old, drop=F], hyper, lphi, ltau,
                                   NNGP, coords[id_old, , drop=F], sigmasq, min_n_nngp)
        log_post_new[clust_old] = EB_res_old$log_post
        EB_res_new = evalLogPostEB(e[id_new], dist_mat[id_new, id_new, drop=F], hyper, lphi, ltau,
                                   NNGP, coords[id_new, , drop=F], sigmasq, min_n_nngp)
        log_post_new[k] = EB_res_new$log_post
        
        phi_all_new[clust_old] = EB_res_old$phi_hat
        phi_all_new[k] = EB_res_new$phi_hat
        tau_all_new[clust_old] = EB_res_old$tau_hat
        tau_all_new[k] = EB_res_new$tau_hat
        
        csize_obs_new[k] = length(id_new)
        csize_obs_new[clust_old] = length(id_old)
        
        # compute log-prior ratio
        log_A = 0
        
        # compute log-proposal ratio
        log_P = 0
        
        # compute log-likelihood ratio
        log_L = log_post_new_c2 - log_post[clust_c1] - log_post[clust_c2]
        log_L = log_L + log_post_new[clust_old] + log_post_new[k] - log_post_old
        
        # acceptance probability
        acc_prob = min(0, log_A + log_P + log_L)
        acc_prob = exp(acc_prob)
        if(runif(1) < acc_prob){
          # accept
          cluster = cluster_new
          cluster_obs = cluster_obs_new
          phi_all = phi_all_new
          tau_all = tau_all_new
          log_post = log_post_new
          edge_status = getEdgeStatus(cluster_new, graph0)
          csize = split_res$csize
          csize_obs = csize_obs_new
          
          # update Cholesky/NNGP factors
          if (EB_res_merged$method == 'exact') {
            cov_factors[[clust_c2]] = list('chol' = evalChol(dist_mat[id_merged, id_merged, drop=F], 
                                                             EB_res_merged$phi_hat, EB_res_merged$tau_hat))
          } else {
            cov_factors[[clust_c2]] = EB_res_merged$prec_param
          }
          cov_factors = cov_factors[-clust_c1]
          
          if (EB_res_old$method == 'exact') {
            cov_factors[[clust_old]] = list('chol' = evalChol(dist_mat[id_old, id_old, drop=F], 
                                                              EB_res_old$phi_hat, EB_res_old$tau_hat))
          } else {
            cov_factors[[clust_old]] = EB_res_old$prec_param
          }
          if (EB_res_new$method == 'exact') {
            cov_factors[[k]] = list('chol' = evalChol(dist_mat[id_new, id_new, drop=F], 
                                                      EB_res_new$phi_hat, EB_res_new$tau_hat))
          } else {
            cov_factors[[k]] = EB_res_new$prec_param
          }
          
          cnt_acc[move] = cnt_acc[move] + 1
        }
        
      }
    }
    
    if(move == 4) { ## Hyper move
      ## update spanning tree
      if (iter > WARMUP)
        stgraph = proposeMST(graph0, edge_status)
      
      ## update random projection
      if (!fix_proj) {
        # transfer NNGP to excat Cholesky factors
        idx = split(1:n, cluster_obs)  # obs. idx in each cluster
        for (j in 1:k) {
          if (length(cov_factors[[j]]) == 1) next
          idx_j = idx[[j]]
          cov_factors[[j]] = list('chol' = evalChol(dist_mat[idx_j, idx_j, drop=F], phi_all[j], tau_all[j]))
        }
        
        for (i in 1:n) {
          # check if this is a boundary point
          nn_cluster = cluster[nearest_knots[i, ]]
          if (length(unique(nn_cluster)) == 1)
            next  # interior point
          
          # get log prior projection prob.
          cluster_prior = tapply(proj_prob[i, ], nn_cluster, sum)
          cluster_i = cluster_obs[i]
          if (csize_obs[cluster_i] <= csize_obs_min)
            cluster_prior[as.character(cluster_i)] = -Inf
          
          # compute posterior projection prob.
          idx = split(1:n, cluster_obs)  # obs. idx in each cluster
          nn_cluster_unique = as.integer(names(cluster_prior))
          phi_nn = phi_all[nn_cluster_unique]; tau_nn = tau_all[nn_cluster_unique]
          log_post_nn = log_post[nn_cluster_unique]  # log posterior of neighboring clusters
          log_like_nn = likeFromGPPost(log_post_nn, phi_nn, tau_nn, hyper)  # log likelihood of neighboring clusters
          names(log_like_nn) = nn_cluster_unique
          log_like_nn_sum = sum(log_like_nn)
          log_like_nn_new = numeric(length(nn_cluster_unique))
          cluster_post = numeric(length(nn_cluster_unique))
          chol_new = list()
          
          # update Cholesky for in cluster_i if obs. i is removed
          idx_i = idx[[cluster_i]]
          pos_i = Rfast::binary_search(idx_i, i, index = T)  # index of obs. i to Choleksy factor
          idx_i_new = idx_i[-pos_i]
          chol_new_i = t( cholDel(t(cov_factors[[cluster_i]]$chol), pos_i) )
          log_like_nn_new_i = updateLogLike(e[idx_i_new], rep(1, csize_obs[cluster_i] - 1),  
                                            csize_obs[cluster_i] - 1, 1, list(list('chol' = chol_new_i)), hyper)
          
          for (j in 1:length(nn_cluster_unique)) {
            cluster_j = nn_cluster_unique[j]
            chol_new[[j]] = list()
            if (cluster_j == cluster_i) {
              # obs. i is already in cluster j
              cluster_post[j] = cluster_prior[j] + log_like_nn_sum
            } else {
              # re-assign obs. i to cluster j
              
              # compute covariance between obs. i and obs. in cluster j
              idx_j = idx[[cluster_j]]
              cov_j = evalCov(dist_mat[i, idx_j, drop=F], phi_all[cluster_j], tau_all[cluster_j], nugget = F)
              
              # update Cholesky for cluster j
              pos_j = Rfast::binary_search(idx_j, i, index = T)  # index to insert i to Choleksy factor
              idx_j_new = append(idx_j, i, pos_j - 1)
              if (pos_j == 1) {
                cov_j_1 = NULL; cov_j_3 = as.numeric(cov_j)
              } else if (pos_j == length(cov_j) + 1) {
                cov_j_1 = as.numeric(cov_j); cov_j_3 = NULL
              } else {
                cov_j_1 = cov_j[, 1:(pos_j-1)]; cov_j_3 = cov_j[, pos_j:length(cov_j)]
              }
              chol_new[[j]] = t( cholAdd(t(cov_factors[[cluster_j]]$chol), pos_j, tau_all[cluster_j] + 1, cov_j_1, cov_j_3) )
              log_like_nn_new[j] = updateLogLike(e[idx_j_new], rep(1, csize_obs[cluster_j] + 1), csize_obs[cluster_j] + 1, 
                                                 1, list(list('chol' = chol_new[[j]])), hyper)
              
              
              cluster_post[j] = cluster_prior[j] + log_like_nn_sum - log_like_nn[j] - log_like_nn[as.character(cluster_i)] +
                log_like_nn_new[j] + log_like_nn_new_i
              
            } # end if
          } # end loop for all neighboring clusters 
          
          # sample a new projection for obs. i
          cluster_post_norm = cluster_post - logSumExp(cluster_post)
          cidx_new = sample.int(length(nn_cluster_unique), 1, prob = exp(cluster_post_norm))
          cluster_i_new = nn_cluster_unique[cidx_new]
          nn_idx = which(nn_cluster == cluster_i_new)
          proj_size[proj[i]] = proj_size[proj[i]] - 1
          proj[i] = nearest_knots[i, nn_idx[sample.int(length(nn_idx), 1, prob = proj_prob[i, nn_idx])] ]
          proj_size[proj[i]] = proj_size[proj[i]] + 1
          
          # update stuff
          if (cluster_i_new != cluster_i) {
            cluster_obs[i] = cluster_i_new
            cov_factors[[cluster_i_new]]$chol = chol_new[[cidx_new]]
            cov_factors[[cluster_i]]$chol = chol_new_i
            
            cidx_changed = c(cluster_i, cluster_i_new)
            log_post[cluster_i_new] = log_like_nn_new[cidx_new]
            log_post[cluster_i] = log_like_nn_new_i
            log_post[cidx_changed] = log_post[cidx_changed] + 
              dtplus(phi_all[cidx_changed], a_phi, b_phi, log = T) +
              dinvgamma(tau_all[cidx_changed], a_tau, b_tau, log = T)
            csize_obs[cidx_changed] = csize_obs[cidx_changed] + c(-1, 1)
          }
        } # end loop for all obs
      } # end fix projection if
    }
    
    ## Update [sigmasq | tau, phi, partition, tree, beta, lambda, y]
    # compute t(y-mu) %*% P^{-1} %*% (y-mu), where P is covariance matrix
    quad_all = evalQuad(Y, X, beta, dist_mat, cluster_obs, k, phi_all, tau_all, cov_factors)
    # sample sigmasq
    shape = csize_obs / 2 + a_sigma
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
    log_post = updateLogLike(e, cluster_obs, csize_obs, k, cov_factors, hyper)
    log_post = log_post + dtplus(phi_all, a_phi, b_phi, log = T)
    log_post = log_post + dinvgamma(tau_all, a_tau, b_tau, log = T)
    
    ## save result
    if(iter > BURNIN & (iter - BURNIN) %% THIN == 0) {
      # sample (mean-zero) latent Gaussian field
      if(sample_latent) {
        latent_out[(iter-BURNIN)/THIN, ] = sampleLGP(e, dist_mat, cluster_obs, k, phi_all, tau_all, sigmasq_all)
      }
      
      phi_out[[(iter-BURNIN)/THIN]] = phi_all
      tau_out[[(iter-BURNIN)/THIN]] = tau_all
      beta_out[(iter-BURNIN)/THIN, ] = beta
      lambda_out[(iter-BURNIN)/THIN] = lambda
      sigmasq_out[[(iter-BURNIN)/THIN]] = sigmasq_all
      if (save_tree)  ST_out[[(iter-BURNIN)/THIN]] = stgraph
      cluster_out[(iter-BURNIN)/THIN, ] = cluster
      cluster_obs_out[(iter-BURNIN)/THIN, ] = cluster_obs
      
      log_full_post_res = evalLogFullPost(e, cluster_obs, phi_all, tau_all, sigmasq_all, beta, lambda, 
                                          cov_factors, k, hyper)
      log_post_out[(iter-BURNIN)/THIN] = log_full_post_res$log_post
      log_like_out[(iter-BURNIN)/THIN] = log_full_post_res$log_like
    }
  }
  
  mcmc_out = list('phi_out' = phi_out, 'tau_out' = tau_out, 'sigmasq_out' = sigmasq_out,
                  'lambda_out' = lambda_out, 'beta_out' = beta_out, 'ST_out' = ST_out,
                  'cluster_out' = cluster_out, 'cluster_obs_out' = cluster_obs_out, 
                  'log_post_out' = log_post_out, 'log_like_out' = log_like_out,
                  'cnt_att' = cnt_att, 'cnt_acc' = cnt_acc,
                  'std_par_Y' = std_par_Y, 'Y_std' = Y)
  if(sample_latent) mcmc_out[['latent_out']] = latent_out
  
  return(mcmc_out)
}


########### Prediction Functions ###########

# function to sample from MVN(mu, Sigma)
rMVN <- function(mu, Sigma) {
  n = length(mu)
  chol_S = t(chol(Sigma))
  z = rnorm(n)
  return(mu + chol_S %*% z)
}

# function to construct a parameter matrix X from a vectors x, y, where X[i, j] = fun(x[i], y[j])
getParamMatrix <- function(x, y = x, fun = c('mean', 'product')) {
  n = length(x)
  if(fun == 'mean')
    return(vapply(y, FUN = function(z) (z + x) / 2, FUN.VALUE = numeric(n)))
  if(fun == 'product')
    return(vapply(y, FUN = function(z) z * x, FUN.VALUE = numeric(n)))
}

# function to rbind Y to by X m times
rbind_m <- function(X, Y, m = 1) {
  res = X
  for(i in 1:m) 
    res = rbind(res, Y)
  return(res)
}

# function to do prediction for responses of a GP model
# if kmv == TRUE: only krigging mean and variance are returned
predict_resp <- function(Y, X, X_new, dist_mat, cdist_mat, dist_mat_new, phi, tau, sigmasq,
                         lambda = NULL, beta = NULL, kmv = F) {
  n = length(Y)
  
  # compute covariance matrices
  scaled_dist = dist_mat / phi
  cov_mat = (1 + sqrt(5) * scaled_dist + 5 / 3 * scaled_dist ^ 2) * exp(-sqrt(5) * scaled_dist)
  cov_mat = tau * cov_mat
  diag(cov_mat) = diag(cov_mat) + 1
  
  scaled_cdist = cdist_mat / phi
  ccov_mat = tau * (1 + sqrt(5) * scaled_cdist + 5 / 3 * scaled_cdist ^ 2) * exp(-sqrt(5) * scaled_cdist)
  
  scaled_dist_new = dist_mat_new / phi
  cov_mat_new = (1 + sqrt(5) * scaled_dist_new + 5 / 3 * scaled_dist_new ^ 2) * exp(-sqrt(5) * scaled_dist_new)
  cov_mat_new = tau * cov_mat_new
  diag(cov_mat_new) = diag(cov_mat_new) + 1
  
  if(!is.null(lambda)) {
    XXt = X %*% t(X)
    cov_mat = cov_mat + lambda * XXt
    ccov_mat = ccov_mat + lambda * XXt
    cov_mat_new = cov_mat_new + lambda * XXt
    
    L_cov = t(chol(cov_mat))
    Linv_Y = forwardsolve(L_cov, Y)
    Linv_ccov = forwardsolve(L_cov, t(ccov_mat))
    
    # krigging mean
    kmean = crossprod(Linv_ccov, Linv_Y)
    # krigging covariance matrix
    kcov = (cov_mat_new - crossprod(Linv_ccov, Linv_ccov)) * sigmasq
    
  } else if(!is.null(beta)){
    L_cov = t(chol(cov_mat))
    Linv_Y = forwardsolve(L_cov, Y - X %*% beta)
    Linv_ccov = forwardsolve(L_cov, t(ccov_mat))
    
    # krigging mean
    kmean = crossprod(Linv_ccov, Linv_Y) + X_new %*% beta
    # krigging covariance matrix
    kcov = (cov_mat_new - crossprod(Linv_ccov, Linv_ccov)) * sigmasq
    
  } else {
    stop('Either lambda or beta has to be provided')
  }
  
  if(!kmv) return(rMVN(kmean, kcov))
  if(kmv) return(data.frame('kmean' = kmean, 'kvar' = diag(kcov)))
}

# function to do prediction for responses of a GP model with NNGP approximation
# if kmv == TRUE: only krigging mean and variance are returned
predict_resp_NNGP <- function(Y, X, X_new, coords, coords_new, phi, tau, sigmasq, beta, 
                              n_neighbors_nngp, kmv = F) {
  n_new = nrow(coords_new)
  
  e = Y - X %*% beta
  pred_nngp = NNGP_prediction(coords, e, NULL, coords_new, NULL, n_neighbors_nngp, 
                              sigmasq * tau, sigmasq, 1 / phi, cov.model = "matern52", verbose = F)
  
  # krigging mean
  kmean = pred_nngp$kmean + X_new %*% beta
  # krigging variance
  kvar = pred_nngp$kvar
  
  if(!kmv) return(rnorm(n_new, kmean, sqrt(kvar)))
  if(kmv) return(data.frame('kmean' = kmean, 'kvar' = kvar))
}

# function to do prediction for latent GP
# if kmv == TRUE: only krigging mean and variance are returned
predict_LGP <- function(latent, dist_mat, cdist_mat, dist_mat_new, phi, sigmasq, kmv = F) {
  n = length(latent)
  
  # compute covariance matrices
  scaled_dist = dist_mat / phi_j
  cov_mat = (1 + sqrt(5) * scaled_dist + 5 / 3 * scaled_dist ^ 2) * exp(-sqrt(5) * scaled_dist)
  
  scaled_cdist = cdist_mat / phi_j
  ccov_mat = (1 + sqrt(5) * scaled_cdist + 5 / 3 * scaled_cdist ^ 2) * exp(-sqrt(5) * scaled_cdist)
  
  scaled_dist_new = dist_mat_new / phi_j
  cov_mat_new = (1 + sqrt(5) * scaled_dist_new + 5 / 3 * scaled_dist_new ^ 2) * exp(-sqrt(5) * scaled_dist_new)
  
  L_cov = t(chol(cov_mat))
  Linv_lat = forwardsolve(L_cov, latent)
  Linv_ccov = forwardsolve(L_cov, t(ccov_mat))
  
  # krigging mean
  kmean = crossprod(Linv_ccov, Linv_lat)
  # krigging covariance matrix
  kcov = (cov_mat_new - crossprod(Linv_ccov, Linv_ccov)) * sigmasq
  
  if(!kmv) return(rMVN(kmean, kcov))
  if(kmv) return(data.frame('kmean' = kmean, 'kvar' = diag(kcov)))
}

# function to do prediction for STGP model (not using knots)
STGP_predict2 <- function(mcmc_res, X, X_ho, coords, coords_ho, pred_response = T, nn = 1,
                          return_values = 'mean', NNGP = F, nn_nngp = 15, min_n_nngp = 300) {
  require(FNN)
  
  Y = mcmc_res$Y_std
  cluster_out = mcmc_res$cluster_out
  cluster_obs_out = mcmc_res$cluster_obs_out
  phi_out = mcmc_res$phi_out
  tau_out = mcmc_res$tau_out
  lambda_out = mcmc_res$lambda_out
  sigmasq_out = mcmc_res$sigmasq_out
  beta_out = mcmc_res$beta_out
  if(!pred_response) latent_out = mcmc_res$latent_out
  std_par_Y = mcmc_res$std_par_Y
  
  # get distance matrix
  dist_mat = as.matrix(fields::rdist(coords))
  dist_mat_ho = as.matrix(fields::rdist(coords_ho))
  # get cross-distance matrix
  cdist_mat = as.matrix(fields::rdist(coords_ho, coords))
  
  # get nearest neighbor of hold-out dataset
  nn_res = get.knnx(coords, coords_ho, k = nn)
  nearest_ho = nn_res$nn.index
  proj_prob_ho = 1 / nn_res$nn.dist
  for (i in 1:nrow(proj_prob_ho)) {
    if (proj_prob_ho[i, 1] == Inf) {
      proj_prob_ho[i, ] = c(1, rep(0, nn - 1))
    } else {
      proj_prob_ho[i, ] = proj_prob_ho[i, ] / sum(proj_prob_ho[i, ])
    }
  }
  
  n_ho = nrow(coords_ho); n = length(Y)
  n_post = length(phi_out)
  Y_ho = matrix(0, nrow = n_ho, ncol = n_post)
  cluster_ho_all = matrix(0, nrow = n_ho, ncol = n_post)
  for(i in 1:n_post) {
    cluster = cluster_out[i, ]
    cluster_obs = cluster_obs_out[i, ]
    # predict cluster_ho probabilistically
    proj_ho = randomProject(nearest_ho, proj_prob_ho, nn)
    cluster_ho = cluster_obs[proj_ho]
    cluster_ho_all[, i] = cluster_ho
    
    k = max(cluster)
    phi_uni = phi_out[[i]]
    tau_uni = tau_out[[i]]
    lambda = lambda_out[i]
    sigmasq_uni = sigmasq_out[[i]]
    beta = beta_out[i, ]
    if(!pred_response) latent = latent_out[i, ]
    idx = split(1:n, cluster_obs)
    idx_ho = split(1:n_ho, cluster_ho)
    
    for(j in names(idx_ho)) {
      phi_j = phi_uni[as.integer(j)]
      tau_j = tau_uni[as.integer(j)]
      sigmasq_j = sigmasq_uni[as.integer(j)]
      idx_j = idx[[j]]; idx_ho_j = idx_ho[[j]]
      Y_j = Y[idx_j]
      X_j = X[idx_j, , drop = F]; X_ho_j = X_ho[idx_ho_j, , drop = F]
      if(!pred_response) latent_j = latent[idx_j]
      n_j = length(idx_j)
      
      if(pred_response) {
        # predict responses
        if (!NNGP | length(Y_j) < min_n_nngp) {
          # use exact covariance matrix
          Y_ho[idx_ho_j, i] = predict_resp(Y_j, X_j, X_ho_j, 
                                           dist_mat[idx_j, idx_j, drop = F],
                                           cdist_mat[idx_ho_j, idx_j, drop = F], 
                                           dist_mat_ho[idx_ho_j, idx_ho_j, drop = F],
                                           phi_j, tau_j, sigmasq_j, beta = beta)
        } else {
          # use NNGP approximation to predict responses
          Y_ho[idx_ho_j, i] = predict_resp_NNGP(Y_j, X_j, X_ho_j, coords[idx_j, , drop = F], 
                                                coords_ho[idx_ho_j, , drop = F], 
                                                phi_j, tau_j, sigmasq_j, beta, nn_nngp)
        }
      } else {
        # predict latent GP
        Y_ho[idx_ho_j, i] = predict_LGP(latent_j, dist_mat[idx_j, idx_j, drop = F],
                                        cdist_mat[idx_ho_j, idx_j, drop = F], 
                                        dist_mat_ho[idx_ho_j, idx_ho_j, drop = F],
                                        phi_j, sigmasq_j * tau_j)
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


# function to get kriging mean and variance (not using knots)
STGP_krigmv2 <- function(mcmc_res, X, X_ho, coords, coords_ho, pred_response = T, nn = 1,
                         NNGP = F, nn_nngp = 15, min_n_nngp = 300) {
  require(FNN)
  
  Y = mcmc_res$Y_std
  cluster_out = mcmc_res$cluster_out
  cluster_obs_out = mcmc_res$cluster_obs_out
  phi_out = mcmc_res$phi_out
  tau_out = mcmc_res$tau_out
  lambda_out = mcmc_res$lambda_out
  sigmasq_out = mcmc_res$sigmasq_out
  beta_out = mcmc_res$beta_out
  if(!pred_response) latent_out = mcmc_res$latent_out
  std_par_Y = mcmc_res$std_par_Y
  
  # get distance matrix
  dist_mat = as.matrix(fields::rdist(coords))
  dist_mat_ho = as.matrix(fields::rdist(coords_ho))
  # get cross-distance matrix
  cdist_mat = as.matrix(fields::rdist(coords_ho, coords))
  
  # get nearest neighbor of hold-out dataset
  nn_res = get.knnx(coords, coords_ho, k = nn)
  nearest_ho = nn_res$nn.index
  proj_prob_ho = 1 / nn_res$nn.dist
  for (i in 1:nrow(proj_prob_ho)) {
    if (proj_prob_ho[i, 1] == Inf) {
      proj_prob_ho[i, ] = c(1, rep(0, nn - 1))
    } else {
      proj_prob_ho[i, ] = proj_prob_ho[i, ] / sum(proj_prob_ho[i, ])
    }
  }
  
  n_ho = nrow(coords_ho); n = length(Y)
  n_post = length(phi_out)
  kmean = matrix(0, nrow = n_ho, ncol = n_post*nn)
  kvar = matrix(0, nrow = n_ho, ncol = n_post*nn)
  for(i in 1:n_post) {
    cluster = cluster_out[i, ]
    cluster_obs = cluster_obs_out[i, ]
    cluster_ho = matrix(cluster_obs[nearest_ho], nrow = n_ho, ncol = nn)
    
    k = max(cluster)
    phi_uni = phi_out[[i]]
    tau_uni = tau_out[[i]]
    sigmasq_uni = sigmasq_out[[i]]
    beta = beta_out[i, ]
    if(!pred_response) latent = latent_out[i, ]
    idx = split(1:n, cluster_obs)
    idx_ho = split(1:(n_ho*nn), c(cluster_ho))
    
    kmean_i = matrix(0, nrow = n_ho, ncol = nn)
    kvar_i = matrix(0, nrow = n_ho, ncol = nn)
    
    for(j in names(idx_ho)) {
      phi_j = phi_uni[as.integer(j)]
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
      
      # use exact covariance matrix
      if(pred_response) {
        # predict responses
        if (!NNGP | length(Y_j) < min_n_nngp) {
          # use exact covariance matrix
          krig_res = predict_resp(Y_j, X_j, X_ho_j, 
                                  dist_mat[idx_j, idx_j, drop = F],
                                  cdist_mat[idx_ho_j, idx_j, drop = F], 
                                  dist_mat_ho[idx_ho_j, idx_ho_j, drop = F],
                                  phi_j, tau_j, sigmasq_j, beta = beta, kmv = T)
        } else {
          # use NNGP approximation to predict responses
          krig_res = predict_resp_NNGP(Y_j, X_j, X_ho_j, coords[idx_j, , drop = F], 
                                       coords_ho[idx_ho_j, , drop = F], 
                                       phi_j, tau_j, sigmasq_j, beta, nn_nngp, kmv = T)
        }
      } else {
        # predict latent GP
        krig_res = predict_LGP(latent_j, dist_mat[idx_j, idx_j, drop = F],
                               cdist_mat[idx_ho_j, idx_j, drop = F], 
                               dist_mat_ho[idx_ho_j, idx_ho_j, drop = F],
                               phi_j, sigmasq_j * tau_j, kmv = T)
      }
      
      for(ii in idx_ho_j_all) {
        jj = which(idx_ho_j == ifelse(ii %% n_ho == 0, n_ho, ii %% n_ho))
        # save krigging mean
        kmean_i[ii] = krig_res$kmean[jj]
        # save krigging variance
        kvar_i[ii] = krig_res$kvar[jj]
      }
    }
    kmean[, seq(i, (nn-1)*n_post+i, n_post)] = unstandardize(c(kmean_i), std_par_Y)
    kvar[, seq(i, (nn-1)*n_post+i, n_post)] = unstandardize(c(kvar_i), std_par_Y, nomean = T, s2 = T)
  }
  
  return(list('kmean' = kmean, 'kvar' = kvar, 'proj_prob' = proj_prob_ho))
}

# function to get scorings
getScores <- function(krig_res, true_val, nn = 1, n_post) {
  require(scoringRules)
  
  kmean = krig_res$kmean[, 1:(nn*n_post)]
  kvar = krig_res$kvar[, 1:(nn*n_post)]
  n = nrow(kmean)
  
  if (is.matrix(krig_res$proj_prob)) {
    weights = t(apply(krig_res$proj_prob[, 1:nn, drop=F], 1, rep, each = n_post))
  } else {
    weights = matrix(1, nrow = n, ncol = ncol(kmean))
  }
  weights = weights / rowSums(weights) 
  
  crps = crps_mixnorm(true_val, kmean, sqrt(kvar), weights)
  logs = logs_mixnorm(true_val, kmean, sqrt(kvar), weights)
  
  return(cbind(crps, logs))
}

