#### Reproducible code for "A Nonstationary Soft Partitioned Gaussian Process Model via Random Spanning Trees"
#### Code for fitting TGP model in Simulation Studies (Section 4)
#### NOTE: This code can take several hours to run.

### Load simulation data -------------------------------------------------------

load("data/sim_input.RData")

### Fit TGP model --------------------------------------------------------------

library(tgp)
set.seed(12345)

MCMC = 100000    # total number of MCMC iterations (including burn-in)
BURNIN = 50000   # number of burn-in iterations
THIN = 10        # thin MCMC by keeping results every 10 iterations 

tgp_res = btgp(coords, Y, XX = coords_ho, meanfn = 'constant', BTE = c(BURNIN, MCMC, THIN), 
               corr = 'matern', nu = 2.5, bprior = 'bmzt', trace = T)
Y_ho_tgp = tgp_res$ZZ.km

# compute MSPE
mse_ho_tgp = mean((Y_ho_tgp - Y_ho)^2)

# function to get probabilistic prediction scores
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
  
  return(c('CRPS' = mean(crps), 'LogS' = mean(logs)))
}

# compute probabilistic prediction scores
n_post = (MCMC - BURNIN) / THIN
Ym0r1 = tgp:::mean0.range1(Y)
kmean_tgp = apply(tgp_res$trace$preds$ZZ.km, 1, tgp:::undo.mean0.range1, Ym0r1$undo)
kvar_tgp = apply(tgp_res$trace$preds$ZZ.ks2, 1, tgp:::undo.mean0.range1, Ym0r1$undo,
                 nomean = T, s2 = T)
krig_res_tgp = list('kmean' = kmean_tgp, 'kvar' = kvar_tgp)
scores_tgp = getScores(krig_res_tgp, Y_ho, nn = 1, n_post)

save(tgp_res, mse_ho_tgp, scores_tgp, Y_ho,
     file = 'data/TGP_sim_results.RData', compress = 'xz')
