#### Reproducible code for "A Nonstationary Soft Partitioned Gaussian Process Model via Random Spanning Trees"
#### Code for fitting TGP model in precipitation application (Section 5)
#### NOTE: This code can take hours to run.

library(tgp)
library(scoringRules)

### Read and preprocess precipitation data -------------------------------------

CONUSprecip = read.csv('data/CONUS_precip_west_18.csv')

# uniform hold-out dataset
set.seed(1234)
n_ho = 250; n = nrow(CONUSprecip) - n_ho
idx_ho = sample.int(nrow(CONUSprecip), n_ho, replace = F)
CONUSprecip_ho = CONUSprecip[idx_ho, ]
CONUSprecip = CONUSprecip[-idx_ho, ]

coords = as.matrix(CONUSprecip[, c('lcc_x', 'lcc_y')])
coords_ho = as.matrix(CONUSprecip_ho[, c('lcc_x', 'lcc_y')])
Y = log(CONUSprecip$precip)
Y_ho = log(CONUSprecip_ho$precip)

### Fit TGP model --------------------------------------------------------------

MCMC = 50000     # total number of MCMC iterations (including burn-in)
BURNIN = 25000   # number of burn-in iterations
THIN = 5         # thin MCMC by keeping results every 10 iterations

set.seed(12345)
tgp_res = btgp(coords, Y, XX = coords_ho, meanfn = 'constant', BTE = c(BURNIN, MCMC, THIN), 
               corr = 'expsep', bprior = 'bmzt', trace = T)

# drop trace to save space
tgp_pred = tgp_res$trace$preds
tgp_res$lpost = tgp_res$trace$post$lpost
tgp_res$trace = NULL

# function to get probabilistic prediction scores
getScores <- function(krig_res, true_val, m = 1, n_post) {
  kmean = krig_res$kmean[, 1:(m*n_post)]
  kvar = krig_res$kvar[, 1:(m*n_post)]
  
  n = nrow(kmean)
  
  weights = matrix(1/(m*n_post), nrow = n, ncol = m*n_post)
  crps = crps_mixnorm(true_val, kmean, sqrt(kvar), weights)
  logs = logs_mixnorm(true_val, kmean, sqrt(kvar), weights)
  
  # return(c('CRPS' = mean(crps), 'LogS' = mean(logs)))
  return(data.frame('CRPS' = crps, 'LogS' = logs))
}

Ym0r1 = tgp:::mean0.range1(Y)
Y_ho_tgp_all = apply(tgp_pred$ZZ, 1, tgp:::undo.mean0.range1, Ym0r1$undo)

# compute probabilistic prediction scores

n_post = (MCMC - BURNIN) / THIN
kmean_tgp = apply(tgp_pred$ZZ.km, 1, tgp:::undo.mean0.range1, Ym0r1$undo)
kvar_tgp = apply(tgp_pred$ZZ.ks2, 1, tgp:::undo.mean0.range1, Ym0r1$undo,
                 nomean = T, s2 = T)
krig_res_tgp = list('kmean' = kmean_tgp, 'kvar' = kvar_tgp)
scores_tgp = getScores(krig_res_tgp, Y_ho, m = 1, n_post)

save(tgp_res, tgp_pred, scores_tgp, Y_ho_tgp_all,
     file = 'data/TGP_app_results.RData', compress = 'xz')
