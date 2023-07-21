#include <string>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <stdio.h>
#include <limits>
#include "lbfgs.h"
#include "util.h"

#ifdef _OPENMP
#include <omp.h>
#endif

extern "C" {

    //Global variables

    double *X_nngp;
    double *y_nngp;

    int n_nngp;
    int p_nngp;
    int m_nngp;

    int covModel_nngp;
    int nThreads_nngp;

    double *D_nngp;
    double *d_nngp;
    int *nnIndx_nngp;
    int *nnIndxLU_nngp;
    int *CIndx_nngp;

    int j_nngp;
    double eps_nngp;
    double fix_nugget_nngp;

    double *prior;


    //covmodel = 0; exponential
    //covmodel = 1; spherical
    //covmodel = 2; matern
    //covmodel = 3; gaussian

    //Defining the likelihood (tausq/sigmasq = alphasq; phi = phi; nu = nu):


    //Update B and F:

    double updateBF(double *B, double *F, double *c, double *C, double *D, double *d, int *nnIndxLU, int *CIndx, int n, double *theta, int covModel, int nThreads, double fix_nugget){
        int i, k, l;
        int info = 0;
        int inc = 1;
        double one = 1.0;
        double zero = 0.0;
        char lower = 'L';
        double logDet = 0;
        double nu = 0;
        //check if the model is 'matern'
        if (covModel == 2) {
            nu = theta[2];
        }

        /* 
        theta = (alpha^2, phi)
        */

        double *bk = (double *) R_alloc(nThreads*(static_cast<int>(1.0+5.0)), sizeof(double));


        //bk must be 1+(int)floor(alpha) * nthread
        int nb = 1+static_cast<int>(floor(5.0));
        int threadID = 0;

#ifdef _OPENMP
#pragma omp parallel for private(k, l, info, threadID)
#endif
        for(i = 0; i < n; i++){
#ifdef _OPENMP
            threadID = omp_get_thread_num();
#endif
            //theta[0] = alphasquareIndex, theta[1] = phiIndex, theta[2] = nuIndex (in case of 'matern')

            // c: cross-covariance between s_i and its neighbors (including nugget)
            // C: inverse of covariance matrix of neighbors (including nugget)
            if(i > 0){
                for(k = 0; k < nnIndxLU[n+i]; k++){
                    c[nnIndxLU[i]+k] = theta[0] * spCor(d[nnIndxLU[i]+k], theta[1], nu, covModel, &bk[threadID*nb]);
                    for(l = 0; l <= k; l++){
                        C[CIndx[i]+l*nnIndxLU[n+i]+k] = theta[0] * spCor(D[CIndx[i]+l*nnIndxLU[n+i]+k], theta[1], nu, covModel, &bk[threadID*nb]);
                        if(l == k){
                            // add nugget part
                            C[CIndx[i]+l*nnIndxLU[n+i]+k] += 1.0*fix_nugget;
                        }
                    }
                }
                // invert C by cholesky
                F77_NAME(dpotrf)(&lower, &nnIndxLU[n+i], &C[CIndx[i]], &nnIndxLU[n+i], &info); if(info != 0){error("c++ error: dpotrf failed\n");}
                F77_NAME(dpotri)(&lower, &nnIndxLU[n+i], &C[CIndx[i]], &nnIndxLU[n+i], &info); if(info != 0){error("c++ error: dpotri failed\n");}
                // compute B = inv(C) %*% c
                F77_NAME(dsymv)(&lower, &nnIndxLU[n+i], &one, &C[CIndx[i]], &nnIndxLU[n+i], &c[nnIndxLU[i]], &inc, &zero, &B[nnIndxLU[i]], &inc);
                // compute F = C(s_i, s_i) - t(B) %*% c = C(s_i, s_i) - t(c) %*% inv(C) %*% c
                F[i] = theta[0] + 1.0*fix_nugget - F77_NAME(ddot)(&nnIndxLU[n+i], &B[nnIndxLU[i]], &inc, &c[nnIndxLU[i]], &inc);
            }else{
                B[i] = 0;
                F[i] = theta[0] + 1.0*fix_nugget;
            }
        }
        for(i = 0; i < n; i++){
            logDet += log(F[i]);
        }

        // det(C)
        return(logDet);
    }

    void solve_B_F(double *B, double *F, double *norm_residual_boot, int n, int *nnIndxLU, int *nnIndx, double *residual_boot){

        residual_boot[0] = norm_residual_boot[0] * sqrt(F[0]);
        double sum;
        for (int i = 1; i < n; i++) {
            sum = norm_residual_boot[i];
            for (int l = 0; l < nnIndxLU[n + i]; l++) {
                sum = sum + B[nnIndxLU[i] + l] * residual_boot[nnIndx[nnIndxLU[i] + l]] / sqrt(F[i]);
            }
            residual_boot[i] = sum * sqrt(F[i]);
        }
    }



    void product_B_F(double *B, double *F, double *residual_nngp, int n, int *nnIndxLU, int *nnIndx, double *norm_residual_nngp){
        norm_residual_nngp[0] = residual_nngp[0]/sqrt(F[0]);
        double sum;
        for (int i = 1; i < n; i++) {
            sum = 0.0;
            for (int l = 0; l < nnIndxLU[n + i]; l++) {
                sum = sum - B[nnIndxLU[i] + l] * residual_nngp[nnIndx[nnIndxLU[i] + l]] / sqrt(F[i]);
            }
            norm_residual_nngp[i] = sum + residual_nngp[i] / sqrt(F[i]);
        }
    }


    void processed_output(double *X, double *y, double *D, double *d, int *nnIndx, int *nnIndxLU, int *CIndx, 
                          int n, int p, int m, double *theta, int covModel, int j, int nThreads, double optimized_likelihod, double *B, double *F, 
                          double *beta, double *Xbeta, double *norm_residual, double *theta_fp, double fix_nugget) {

        char const *ntran = "N";
        int nIndx = static_cast<int>(static_cast<double>(1+m)/2*m+(n-m-1)*m);
        double *c =(double *) R_alloc(nIndx, sizeof(double));
        double *C = (double *) R_alloc(j, sizeof(double)); zeros(C, j);

        double logDet;

        int pp = p*p;
        int info = 0;
        const double negOne = -1.0;
        const double one = 1.0;
        const double zero = 0.0;
        const int inc = 1;
        char const *lower = "L";

        /*
        double *tmp_pp = (double *) R_alloc(pp, sizeof(double));
        double *tmp_p = (double *) R_alloc(p, sizeof(double));
        double *tmp_n = (double *) R_alloc(n, sizeof(double));
        double *residual = (double *) R_alloc(n, sizeof(double));
        */

        /* 
        theta = (alpha^2, phi)
        */

        //create B and F
        logDet = updateBF(B, F, c, C, D, d, nnIndxLU, CIndx, n, theta, covModel, nThreads, fix_nugget);

        /*
        int i;
        for(i = 0; i < p; i++){
            tmp_p[i] = Q(B, F, &X[n*i], y, n, nnIndx, nnIndxLU);
            for(j = 0; j <= i; j++){
                tmp_pp[j*p+i] = Q(B, F, &X[n*j], &X[n*i], n, nnIndx, nnIndxLU);
            }
        }

        // invert tmp_pp (in place) using Cholesky factorization
        F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
        F77_NAME(dpotri)(lower, &p, tmp_pp, &p, &info); if(info != 0){error("c++ error: dpotri failed\n");}

        //create beta = tmp_pp %*% tmp_p
        F77_NAME(dsymv)(lower, &p, &one, tmp_pp, &p, tmp_p, &inc, &zero, beta, &inc);

        //create Xbeta = X %*% beta
        F77_NAME(dgemv)(ntran, &n, &p, &one, X, &n, beta, &inc, &zero, tmp_n, &inc);

        // tmp_n = Xbeta
        dcopy_(&n, tmp_n, &inc, Xbeta, &inc);


        //create normalized residual
        F77_NAME(daxpy)(&n, &negOne, y, &inc, tmp_n, &inc);

        for (int s = 0; s < n; s++) {
            residual[s] = negOne * tmp_n[s];
        }

        product_B_F(B, F, residual, n, nnIndxLU, nnIndx, norm_residual);

        */

        //Create complete theta

        /*
        // 1. Create sigma square
        theta_fp[0] = exp( (2.0 * optimized_likelihod - logDet - n * M_LN_2PI ) / n );  // change

        // 2. Create tau square
        theta_fp[1] = theta[0] * theta_fp[0] * fix_nugget; //change

        // 3. Create phi
        theta_fp[2] = theta[1];

        // 4. Create nu in "matern"
        if (covModel == 2) {
            theta_fp[3] = theta[2];
        }
        */
        // 1. Create alpha square = spatial variance / nugget variance
        theta_fp[0] = theta[0];

        // 2. Create phi
        theta_fp[1] = theta[1];

        // 3. Create nu in "matern"
        if (covModel == 2) {
            theta_fp[2] = theta[2];
        }
    }



    //Defining likelihood in terms of theta.
    double likelihood(double *X, double *y, double *D, double *d, int *nnIndx, int *nnIndxLU, int *CIndx, int n, int p, int m, double *theta, 
                      int covModel, int j, int nThreads, double fix_nugget, double *prior) {

        char const *ntran = "N";
        int nIndx = static_cast<int>(static_cast<double>(1+m)/2*m+(n-m-1)*m);
        double *B = (double *) R_alloc(nIndx, sizeof(double));
        double *F = (double *) R_alloc(n, sizeof(double));
        double *c =(double *) R_alloc(nIndx, sizeof(double));
        double *C = (double *) R_alloc(j, sizeof(double)); zeros(C, j);

        double logDet;

        int pp = p*p;
        int info = 0;
        double neg_log_likelihood;
        const double negOne = -1.0;
        const double one = 1.0;
        const double zero = 0.0;
        const int inc = 1;
        char const *lower = "L";


        double *tmp_pp = (double *) R_alloc(pp, sizeof(double));
        double *tmp_p = (double *) R_alloc(p, sizeof(double));
        double *beta = (double *) R_alloc(p, sizeof(double));
        double *tmp_n = (double *) R_alloc(n, sizeof(double));

        /* 
        theta = (alpha^2, phi)
        */

        // get B, F such that inv(C) = t(I - B) %*% inv(F) %*% (I - B), where F is diagonal
        logDet = updateBF(B, F, c, C, D, d, nnIndxLU, CIndx, n, theta, covModel, nThreads, fix_nugget);

        /*
        // tmp_pp: inverse of t(X) %*% inv(C) %*% X
        // tmp_p: t(X) %*% inv(Sigma) %*% y
        int i;
        for(i = 0; i < p; i++){
            tmp_p[i] = Q(B, F, &X[n*i], y, n, nnIndx, nnIndxLU);
            for(j = 0; j <= i; j++){
                tmp_pp[j*p+i] = Q(B, F, &X[n*j], &X[n*i], n, nnIndx, nnIndxLU);
            }
        }

        F77_NAME(dpotrf)(lower, &p, tmp_pp, &p, &info); if(info != 0){error("c++ error: dpotrf failed\n");}
        F77_NAME(dpotri)(lower, &p, tmp_pp, &p, &info); if(info != 0){error("c++ error: dpotri failed\n");}

        // beta = tmp_pp %*% tmp_p 
        F77_NAME(dsymv)(lower, &p, &one, tmp_pp, &p, tmp_p, &inc, &zero, beta, &inc);
        // tmp_n = X %*% beta
        F77_NAME(dgemv)(ntran, &n, &p, &one, X, &n, beta, &inc, &zero, tmp_n, &inc);
        // tmp_n -= y
        F77_NAME(daxpy)(&n, &negOne, y, &inc, tmp_n, &inc);

        // calculates negative log likelihood
        neg_log_likelihood = n * log(Q(B, F, tmp_n, tmp_n, n, nnIndx, nnIndxLU)/n) + logDet;
        */

        /*
        // calculates negative log likelihood
        neg_log_likelihood = n * log(Q(B, F, y, y, n, nnIndx, nnIndxLU)/n) + logDet;
        // log(2 * pi)
        neg_log_likelihood += n * M_LN_2PI;
        neg_log_likelihood /= 2.0;
        */

        // calculates negative log likelihood
        double a_sigma = prior[4], b_sigma = prior[5];
        // compute t(y) %*% inv(C) %*% y = t(y) %*% t(I - B) %*% inv(F) %*% (I - B) %*% y
        double quad = Q(B, F, y, y, n, nnIndx, nnIndxLU);
        neg_log_likelihood = 0.5*logDet + n/2.0*M_LN_2PI - lgammafn(n/2.0 + a_sigma) + (n/2.0 + a_sigma)*log(b_sigma + quad/2.0);
        neg_log_likelihood += lgammafn(a_sigma) - a_sigma*log(b_sigma);

        // Rprintf("phi: %f\n", theta[1]);
        // Rprintf("alphasq: %f\n", theta[0]);

        // Rprintf("Log like: %f\n", -neg_log_likelihood);
        // Rprintf("Quad: %f\n", quad);
        // Rprintf("Log det: %f\n", logDet);

        // add prior
        double alphasq = theta[0], phi = theta[1];
        double a_phi = prior[0], b_phi = prior[1];
        double a_alpha = prior[2], b_alpha = prior[3];
        neg_log_likelihood -= dinvgamma(alphasq, a_alpha, b_alpha, true);
        neg_log_likelihood -= dtplus(1.0/phi, a_phi, b_phi, true);

        return(neg_log_likelihood);
    }



    //Defining likelihood w.r.t unconstrained optimization with alpha, root_phi, root_nu (in case of matern);

    // a. Non-matern models
    double likelihood_lbfgs_non_matern(double alpha, double root_phi, double *X, double *y, double *D, double *d, int *nnIndx, int *nnIndxLU, int *CIndx, 
                                       int n, int p, int m, int covModel, int j, int nThreads, double fix_nugget, double *prior) {
        double *theta = (double *) R_alloc(2, sizeof(double));
        theta[0] = pow(alpha, 2.0);
        theta[1] = pow(root_phi, 2.0);
        double res = likelihood(X, y, D, d, nnIndx, nnIndxLU, CIndx, n, p, m, theta, covModel, j, nThreads, fix_nugget, prior);//some unnecessary checking are happening here, will remove afterwards
        return(res);
    }


    //b. matern models
    double likelihood_lbfgs_matern(double alpha, double root_phi, double root_nu, double *X, double *y, double *D, double *d, int *nnIndx, int *nnIndxLU, int *CIndx, int n, int p, int m, int covModel, int j, int nThreads, double fix_nugget){
        double *theta = (double *) R_alloc(3, sizeof(double));
        theta[0] = pow(alpha, 2.0);
        theta[1] = pow(root_phi, 2.0);
        theta[2] = pow(root_nu, 2.0);
        double res = likelihood(X, y, D, d, nnIndx, nnIndxLU, CIndx, n, p, m, theta, covModel, j, nThreads, fix_nugget, prior);//some unnecessary checking are happening here, will remove afterwards;
        return(res);
    }


    static lbfgsfloatval_t evaluate(
            void *instance,
            const lbfgsfloatval_t *x,
            lbfgsfloatval_t *g,
            const int n,
            const lbfgsfloatval_t step
    )
    {
        int i;
        lbfgsfloatval_t fx = 0.0;

        if (covModel_nngp != 2) {
            for (i = 0;i < n;i += 2) {
                g[i+1] = (likelihood_lbfgs_non_matern(x[i], (x[i+1]  + eps_nngp), X_nngp, y_nngp, D_nngp, d_nngp, nnIndx_nngp, nnIndxLU_nngp, CIndx_nngp, 
                                                      n_nngp, p_nngp, m_nngp, covModel_nngp, j_nngp, nThreads_nngp, fix_nugget_nngp, prior)
                         - likelihood_lbfgs_non_matern(x[i], (x[i+1] - eps_nngp), X_nngp, y_nngp, D_nngp, d_nngp, nnIndx_nngp, nnIndxLU_nngp, CIndx_nngp, 
                                                           n_nngp, p_nngp, m_nngp, covModel_nngp, j_nngp, nThreads_nngp,fix_nugget_nngp, prior) 
                            ) / (2*eps_nngp);

                g[i] = (likelihood_lbfgs_non_matern((x[i] + eps_nngp), x[i+1], X_nngp, y_nngp, D_nngp, d_nngp, nnIndx_nngp, nnIndxLU_nngp, CIndx_nngp, 
                                                    n_nngp, p_nngp, m_nngp, covModel_nngp, j_nngp, nThreads_nngp, fix_nugget_nngp, prior)
                         - likelihood_lbfgs_non_matern((x[i] - eps_nngp), x[i+1], X_nngp, y_nngp, D_nngp, d_nngp, nnIndx_nngp, nnIndxLU_nngp, CIndx_nngp, 
                                                       n_nngp, p_nngp, m_nngp, covModel_nngp, j_nngp, nThreads_nngp,fix_nugget_nngp, prior)
                            ) / (2*eps_nngp);

                fx += likelihood_lbfgs_non_matern(x[i], x[i+1], X_nngp, y_nngp, D_nngp, d_nngp, nnIndx_nngp, nnIndxLU_nngp, CIndx_nngp, 
                                                  n_nngp, p_nngp, m_nngp, covModel_nngp,j_nngp, nThreads_nngp,fix_nugget_nngp, prior);
            }
        } else {
            for (i = 0;i < n;i += 3) {
                g[i+1] = (likelihood_lbfgs_matern(x[i], (x[i+1]  + eps_nngp), x[i+2], X_nngp, ::y_nngp, D_nngp, d_nngp, nnIndx_nngp, nnIndxLU_nngp, CIndx_nngp, n_nngp, p_nngp, m_nngp, covModel_nngp, j_nngp, nThreads_nngp, fix_nugget_nngp) - likelihood_lbfgs_matern(x[i], (x[i+1] - eps_nngp), x[i+2], X_nngp, y_nngp, D_nngp, d_nngp, nnIndx_nngp, nnIndxLU_nngp, CIndx_nngp, n_nngp, p_nngp, m_nngp, covModel_nngp, j_nngp, nThreads_nngp, fix_nugget_nngp))/(2*eps_nngp);

                g[i] = (likelihood_lbfgs_matern((x[i] + eps_nngp), x[i+1], x[i+2], X_nngp, ::y_nngp, D_nngp, d_nngp, nnIndx_nngp, nnIndxLU_nngp, CIndx_nngp, n_nngp, p_nngp, m_nngp, covModel_nngp, j_nngp, nThreads_nngp, fix_nugget_nngp) - likelihood_lbfgs_matern((x[i] - eps_nngp), x[i+1], x[i+2], X_nngp, y_nngp, D_nngp, d_nngp, nnIndx_nngp, nnIndxLU_nngp, CIndx_nngp, n_nngp, p_nngp, m_nngp, covModel_nngp, j_nngp, nThreads_nngp, fix_nugget_nngp))/(2*eps_nngp);

                g[i+2] = (likelihood_lbfgs_matern(x[i], x[i+1], (x[i+2] + eps_nngp), X_nngp, y_nngp, D_nngp, d_nngp, nnIndx_nngp, nnIndxLU_nngp, CIndx_nngp, n_nngp, p_nngp, m_nngp, covModel_nngp, j_nngp, nThreads_nngp, fix_nugget_nngp) - likelihood_lbfgs_matern(x[i], x[i+1], (x[i+2] - eps_nngp), X_nngp, y_nngp, D_nngp, d_nngp, nnIndx_nngp, nnIndxLU_nngp, CIndx_nngp, n_nngp, p_nngp, m_nngp, covModel_nngp, j_nngp, nThreads_nngp, fix_nugget_nngp))/(2*eps_nngp);

                fx += likelihood_lbfgs_matern(x[i], x[i+1], x[i+2], X_nngp, y_nngp, D_nngp, d_nngp, nnIndx_nngp, nnIndxLU_nngp, CIndx_nngp, n_nngp, p_nngp, m_nngp, covModel_nngp, j_nngp, nThreads_nngp, fix_nugget_nngp);
            }
        }
        return fx;
    }

    /**
     * @brief 
     * 
     * @param y_r 
     * @param X_r 
     * @param p_r 
     * @param n_r 
     * @param m_r 
     * @param coords_r 
     * @param covModel_r 
     * @param alphaStarting_r sqrt(spatial variance / nugget variance)
     * @param phiStarting_r sqrt(phi)
     * @param nuStarting_r 
     * @param sType_r 
     * @param nThreads_r 
     * @param verbose_r 
     * @param eps_r 
     * @param fix_nugget_r 
     * @param prior_r 
     * @return SEXP 
     */
    SEXP BRISC_estimatecpp(SEXP y_r, SEXP X_r, SEXP p_r, SEXP n_r, SEXP m_r, SEXP coords_r, SEXP covModel_r, 
                           SEXP alphaStarting_r, SEXP phiStarting_r, SEXP nuStarting_r,
                           SEXP sType_r, SEXP nThreads_r, SEXP verbose_r, SEXP eps_r, SEXP fix_nugget_r,
                           SEXP prior_r) {

        int i, k, l, nProtect=0;

        //get args
        y_nngp = REAL(y_r);
        X_nngp = REAL(X_r);
        p_nngp = INTEGER(p_r)[0];
        n_nngp = INTEGER(n_r)[0];
        m_nngp = INTEGER(m_r)[0];
        eps_nngp = REAL(eps_r)[0];
        fix_nugget_nngp = REAL(fix_nugget_r)[0];
        double *coords = REAL(coords_r);

        covModel_nngp = INTEGER(covModel_r)[0];
        std::string corName = getCorName(covModel_nngp);

        nThreads_nngp = INTEGER(nThreads_r)[0];
        int verbose = INTEGER(verbose_r)[0];



#ifdef _OPENMP
        omp_set_num_threads(nThreads_nngp);
#else
        if(nThreads_nngp > 1){
            warning("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads_nngp);
            nThreads_nngp = 1;
        }
#endif

        if(verbose){
            Rprintf("----------------------------------------\n");
            Rprintf("\tModel description\n");
            Rprintf("----------------------------------------\n");
            Rprintf("BRISC model fit with %i observations.\n\n", n_nngp);
            Rprintf("Number of covariates %i (including intercept if specified).\n\n", p_nngp);
            Rprintf("Using the %s spatial correlation model.\n\n", corName.c_str());
            Rprintf("Using %i nearest neighbors.\n\n", m_nngp);
#ifdef _OPENMP
            Rprintf("\nSource compiled with OpenMP support and model fit using %i thread(s).\n", nThreads_nngp);
#else
            Rprintf("\n\nSource not compiled with OpenMP support.\n");
#endif
        }

        //parameters
        int nTheta;

        if(corName != "matern"){
            nTheta = 2; //alpha = 0, sqrt(phi) = 1
        }else{
            nTheta = 3; //alpha = 0, sqrt(phi) = 1, nu = 2;
        }
        
        //starting
        double *theta = (double *) R_alloc (nTheta, sizeof(double));
        theta[0] = REAL(alphaStarting_r)[0];
        theta[1] = REAL(phiStarting_r)[0];
        if (corName == "matern") {
            theta[2] = REAL(nuStarting_r)[0];
        }

        // prior hyperparameters
        prior = REAL(prior_r);

        //allocated for the nearest neighbor index vector (note, first location has no neighbors).
        int nIndx = static_cast<int>(static_cast<double>(1+m_nngp)/2*m_nngp+(n_nngp-m_nngp-1)*m_nngp);
        SEXP nnIndx_r; PROTECT(nnIndx_r = allocVector(INTSXP, nIndx)); nProtect++; nnIndx_nngp = INTEGER(nnIndx_r);
        SEXP d_r; PROTECT(d_r = allocVector(REALSXP, nIndx)); nProtect++; d_nngp = REAL(d_r);

        SEXP nnIndxLU_r; PROTECT(nnIndxLU_r = allocVector(INTSXP, 2*n_nngp)); nProtect++; nnIndxLU_nngp = INTEGER(nnIndxLU_r); //first column holds the nnIndx index for the i-th location and the second columns holds the number of neighbors the i-th location has (the second column is a bit of a waste but will simplifying some parallelization).

        //make the neighbor index
        if(verbose){
            Rprintf("----------------------------------------\n");
            Rprintf("\tBuilding neighbor index\n");
#ifdef Win32
            R_FlushConsole();
#endif
        }

        if(INTEGER(sType_r)[0] == 0){
            mkNNIndx(n_nngp, m_nngp, coords, nnIndx_nngp, d_nngp, nnIndxLU_nngp);
        }
        if(INTEGER(sType_r)[0] == 1){
            mkNNIndxTree0(n_nngp, m_nngp, coords, nnIndx_nngp, d_nngp, nnIndxLU_nngp);
        }else{
            mkNNIndxCB(n_nngp, m_nngp, coords, nnIndx_nngp, d_nngp, nnIndxLU_nngp);
        }


        SEXP CIndx_r; PROTECT(CIndx_r = allocVector(INTSXP, 2*n_nngp)); nProtect++; CIndx_nngp = INTEGER(CIndx_r); //index for D and C.
        for(i = 0, j_nngp = 0; i < n_nngp; i++){//zero should never be accessed
            j_nngp += nnIndxLU_nngp[n_nngp+i]*nnIndxLU_nngp[n_nngp+i];
            if(i == 0){
                CIndx_nngp[n_nngp+i] = 0;
                CIndx_nngp[i] = 0;
            }else{
                CIndx_nngp[n_nngp+i] = nnIndxLU_nngp[n_nngp+i]*nnIndxLU_nngp[n_nngp+i];
                CIndx_nngp[i] = CIndx_nngp[n_nngp+i-1] + CIndx_nngp[i-1];
            }
        }

        SEXP j_r; PROTECT(j_r = allocVector(INTSXP, 1)); nProtect++; INTEGER(j_r)[0] = j_nngp;

        SEXP D_r; PROTECT(D_r = allocVector(REALSXP, j_nngp)); nProtect++; D_nngp = REAL(D_r);

        for(i = 0; i < n_nngp; i++){
            for(k = 0; k < nnIndxLU_nngp[n_nngp+i]; k++){
                for(l = 0; l <= k; l++){
                    D_nngp[CIndx_nngp[i]+l*nnIndxLU_nngp[n_nngp+i]+k] = dist2(coords[nnIndx_nngp[nnIndxLU_nngp[i]+k]], coords[n_nngp+nnIndx_nngp[nnIndxLU_nngp[i]+k]], coords[nnIndx_nngp[nnIndxLU_nngp[i]+l]], coords[n_nngp+nnIndx_nngp[nnIndxLU_nngp[i]+l]]);
                }
            }
        }

        if(verbose){
            Rprintf("----------------------------------------\n");
            Rprintf("\tPerforming optimization\n");
#ifdef Win32
            R_FlushConsole();
#endif
        }

        int i_0, ret = 0;
        int k_0 = 0;
        lbfgsfloatval_t fx;
        lbfgsfloatval_t *x = lbfgs_malloc(nTheta);
        lbfgs_parameter_t param;

        /* Initialize the variables. */
        for (i_0 = 0;i_0 < nTheta; i_0++) {
            x[i_0] = theta[i_0];
        }

        /* Initialize the parameters for the L-BFGS optimization. */
        lbfgs_parameter_init(&param);
        param.epsilon = 1e-2;
        param.gtol = 0.9;
        /*param.linesearch = LBFGS_LINESEARCH_BACKTRACKING;*/

        /*
         Start the L-BFGS optimization; this will invoke the callback functions
         evaluate() and progress() when necessary.
         */
        ret = lbfgs(nTheta, x, &fx, evaluate, NULL, NULL, &param);

        // Construct output
        double *theta_nngp = (double *) R_alloc(nTheta, sizeof(double));
        for (k_0 = 0; k_0 < nTheta; k_0++){
            theta_nngp[k_0] = pow(x[k_0], 2.0);
        }

        // Clean up
        lbfgs_free(x);

        if(verbose){
            Rprintf("----------------------------------------\n");
            Rprintf("\tProcessing optimizers\n");
            Rprintf("----------------------------------------\n");
#ifdef Win32
            R_FlushConsole();
#endif
        }

        // int nTheta_full = nTheta + 1;
        int nTheta_full = nTheta;

        SEXP B_r; PROTECT(B_r = allocVector(REALSXP, nIndx)); nProtect++; double *B_nngp = REAL(B_r);

        SEXP F_r; PROTECT(F_r = allocVector(REALSXP, n_nngp)); nProtect++; double *F_nngp = REAL(F_r);

        SEXP beta_r; PROTECT(beta_r = allocVector(REALSXP, p_nngp)); nProtect++; double *beta_nngp = REAL(beta_r);


        SEXP Xbeta_r; PROTECT(Xbeta_r = allocVector(REALSXP, n_nngp)); nProtect++; double *Xbeta_nngp = REAL(Xbeta_r);

        SEXP norm_residual_r; PROTECT(norm_residual_r = allocVector(REALSXP, n_nngp)); nProtect++; double *norm_residual_nngp = REAL(norm_residual_r);

        SEXP theta_fp_r; PROTECT(theta_fp_r = allocVector(REALSXP, nTheta_full)); nProtect++; double *theta_fp_nngp = REAL(theta_fp_r);

        SEXP max_log_like_r; PROTECT(max_log_like_r = allocVector(REALSXP, 1)); nProtect++; REAL(max_log_like_r)[0] = -fx;

        processed_output(X_nngp, y_nngp, D_nngp, d_nngp, nnIndx_nngp, nnIndxLU_nngp, CIndx_nngp, n_nngp, p_nngp, m_nngp, theta_nngp, covModel_nngp, j_nngp, nThreads_nngp, fx, B_nngp, F_nngp, beta_nngp, Xbeta_nngp, norm_residual_nngp, theta_fp_nngp, fix_nugget_nngp);

        //return stuff
        SEXP result_r, resultName_r;
        int nResultListObjs = 12;  // change

        PROTECT(result_r = allocVector(VECSXP, nResultListObjs)); nProtect++;
        PROTECT(resultName_r = allocVector(VECSXP, nResultListObjs)); nProtect++;

        SET_VECTOR_ELT(result_r, 0, B_r);
        SET_VECTOR_ELT(resultName_r, 0, mkChar("B"));

        SET_VECTOR_ELT(result_r, 1, F_r);
        SET_VECTOR_ELT(resultName_r, 1, mkChar("F"));

        SET_VECTOR_ELT(result_r, 2, max_log_like_r);
        SET_VECTOR_ELT(resultName_r, 2, mkChar("log_post"));

        SET_VECTOR_ELT(result_r, 3, norm_residual_r);
        SET_VECTOR_ELT(resultName_r, 3, mkChar("norm.residual"));

        SET_VECTOR_ELT(result_r, 4, theta_fp_r);
        SET_VECTOR_ELT(resultName_r, 4, mkChar("theta"));


        SET_VECTOR_ELT(result_r, 5, Xbeta_r);
        SET_VECTOR_ELT(resultName_r, 5, mkChar("Xbeta"));

        SET_VECTOR_ELT(result_r, 6, nnIndxLU_r);
        SET_VECTOR_ELT(resultName_r, 6, mkChar("nnIndxLU"));

        SET_VECTOR_ELT(result_r, 7, CIndx_r);
        SET_VECTOR_ELT(resultName_r, 7, mkChar("CIndx"));

        SET_VECTOR_ELT(result_r, 8, D_r);
        SET_VECTOR_ELT(resultName_r, 8, mkChar("D"));

        SET_VECTOR_ELT(result_r, 9, d_r);
        SET_VECTOR_ELT(resultName_r, 9, mkChar("d"));

        SET_VECTOR_ELT(result_r, 10, nnIndx_r);
        SET_VECTOR_ELT(resultName_r, 10, mkChar("nnIndx"));

        SET_VECTOR_ELT(result_r, 11, j_r);
        SET_VECTOR_ELT(resultName_r, 11, mkChar("Length.D"));


        namesgets(result_r, resultName_r);

        //unprotect
        UNPROTECT(nProtect);


        return(result_r);
    }

    // R interface to computes the quadratic term: t(u) %*% t(I - B) %*% inv(F) %*% (I - B) %*% v
    SEXP evalQuadCpp(SEXP B_r, SEXP F_r, SEXP u_r, SEXP v_r, SEXP nnIndx_r, SEXP nnIndxLU_r) {
        int n = length(u_r);
        double quad = Q(REAL(B_r), REAL(F_r), REAL(u_r), REAL(v_r), n, INTEGER(nnIndx_r), INTEGER(nnIndxLU_r));

        SEXP quad_r;  PROTECT(quad_r = allocVector(REALSXP, 1));
        REAL(quad_r)[0] = quad;
        UNPROTECT(1);
        return quad_r;
    }
}
