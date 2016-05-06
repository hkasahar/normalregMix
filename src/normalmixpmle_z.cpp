#include <R.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Arith.h>

/* Compute the EM step and the matrix of "posterior" probabilities
   for the penalized MLE of a finite mixture of regressions
   Based on normpost.c in mixtools package.  */
extern "C" {
    void normalmixpmle_z(
        int *set, /* vector of n, m, p, ninits, maxit, jpvt */
        double *y, /* n by 1 vector of dependent varaible */
        double *z, /* n by p matrix of explanatory variables wrt gamma */
        double *alphaset, /* current vector of mixing parameters */
        double *muset, /* current vector of means */
        double *sigmaset, /* current vector of component stdevs */
        double *gammaset, /* current vector of regression coefficients */
        double *sigma0, /* stdevs used in the penalty term */
        double *mu0,    /* estimate of mu from m-1 component model */
        double *aan, /* constant in the penalty term */
        int *hh, /* the component that is split into two */
        double *lub, /* lower and upper bound on mu */
        double *work, /* 3*m-vector of workspace, which will be broken into 3 parts */
        double *post, /* n by m matrix of posterior probabilities */
        double *loglikset, /* scalar loglikelihood value */
        double *penloglikset, /* scalar penalized loglikelihood value */
        int *notcg,  /* =1 if not convergent */
        double *ttol, /* tolerance level for exiting iteration */
        double *wxy /* n by (2+p) matrix of zgamma, yhat, and zCopy */
        ) {
            int *nn = set;          /* sample size */
            int *mm = set+1;        /* number of components */
            int *pp = set+2;        /* number of variables in z */
            int *nninits = set+3;   /* number of initial values */
            int *mmaxit = set+4;    /* maximum number of iterations */
            int *jpvt = set+5;      /* pivots used in dgelsy */

            int n=*nn, m=*mm, p=*pp, h=*hh, i, j, ii, minj=0, info, rk;
            int ninits=*nninits, maxit = *mmaxit, emit, sing;
            const int np = n*p;
            double tol = *ttol, oldpenloglik;
            double wmu=0.0, r, rowsum, min=0.0, an=*aan, worksize, ssr_j, s0j;
            double *AlpSigRatio = work+m; /* Second 1/3 of workspace, for frequently used constants */
            double *logAlpSigRatio = work+2*m; /* Third 1/3 of workspace, for frequently used constants */
            double *lb = lub;
            double *ub = lub + m;
            double *zgamma = wxy;
            double *yhat = wxy + n;
            double *zCopy = wxy + n*2;

            /* First, make the query, store the result in 'nwork', and allocate the workspace. */
            double rcond = sqrt((2.2204460492503131e-16)/2);    //  2.22... is Machine epsilon
            const int iOne = 1;
            int query = -1;  // Code for "query mode"
            F77_CALL(dgelsy)(nn, pp, &iOne, z, nn, y, nn, jpvt, &rcond, &rk, &worksize, &query, &info);
            int nwork = static_cast<int>(worksize);
            double *work2;
            work2 = (double *) R_alloc(nwork, sizeof(double));
            if (work2==NULL) error("Cannot allocate workspace\n");

            if (h!=0) {  // If h!=0, compute upper and lower bounds
                mu0[0] = R_NegInf;
                mu0[m] = R_PosInf;
                double *lb0 = mu0;
                double *ub0 = mu0+1;

                for (j=0; j<h; j++) {
                    lb[j] = (lb0[j]+lb0[j+1])/2.0;
                    ub[j] = (ub0[j]+ub0[j+1])/2.0;
                }
                for (j=h; j<m; j++) {
                    lb[j] = (lb0[j-1]+lb0[j])/2.0;
                    ub[j] = (ub0[j-1]+ub0[j])/2.0;
                }
            }
            int jn;
            double penloglik = 0.0;
            double loglik = 0.0;

            for (jn=0; jn<ninits; jn++) {
                double *alpha = alphaset + jn*m;
                double *mu    = muset + jn*m;
                double *sigma = sigmaset + jn*m;
                double *gamma = gammaset + jn*p;

                /* initialize EM iteration */

                oldpenloglik = R_NegInf;
                emit = 0;
                double diff = 1.0;
                sing = 0;

                /* EM loop begins */
                while(1)
                {
                loglik = -(double)n * 0.91893853320467274178; /* n/2 times log(2pi) */
                for (j=0; j<m; j++) { /* store some often-used values to save time later */
                    AlpSigRatio[j] = alpha[j] / sigma[j];
                    logAlpSigRatio[j] = log(AlpSigRatio[j]);
                }

                for (i=0; i<n; i++) {
                    zgamma[i] = 0.0;    // Initialize zgamma
                    for (ii=0; ii<p; ii++) {  // Compute zgamma
                        zgamma[i] += z[i + n*ii]*gamma[ii];
                    }
                    for (j=0; j<m; j++) {
                        r = y[i] - zgamma[i] - mu[j];
                        r = r*r;
                        work[j] = r = r / (2.0 * sigma[j] * sigma[j]);
                        /* Keep track of the smallest standardized squared residual. By dividing
                         everything by the component density with the smallest such residual,
                         the denominator of the posterior is guaranteed to be at least one and
                         cannot be infinite unless the values of alpha or sigma are very large or small.
                         This helps prevent numerical problems when calculating the posteriors.*/
                        if (j==0 || r < min) {
                            minj = j;
                            min = r;
                        }
                    }
                    /* At this stage, work contains the squared st'dized resids over 2 */
                    rowsum = 1.0;
                    for (j=0; j<m; j++) {
                        if (j==minj)
                            work[j] = 1.0;
                        else {
                            work[j] = (AlpSigRatio[j] / AlpSigRatio[minj]) * exp(min - work[j]);
                            rowsum += work[j];
                        }
                    }
                    /* At this stage, work contains the normal density at the ith observation divided by the
                     normal density with the largest st'dized resid Thus, dividing by rowsum gives the  posteriors: */
                    for (j=0; j<m; j++) {
                        post[i + n*j] = work[j] / rowsum;
                    }
                    /* Finally, adjust the loglikelihood correctly */
                    loglik += log(rowsum) - min + logAlpSigRatio[minj];
                } /* end for (i=0; i<n; i++) loop */

                /* Compute the penalized loglik. Note that penalized loglik uses old (not updated) sigma */
                penloglik = loglik;
                for (j=0; j<m; j++) {
                    s0j = sigma0[j]/sigma[j];
                    penloglik += -an*(s0j*s0j - 2.0*log(s0j) -1.0);
                }

                diff = penloglik - oldpenloglik;
                oldpenloglik = penloglik;

                /* Normal exit */
                    if (diff < tol || emit>=maxit) {
                        break;
                    }

                /* If not exit, update alpha, mu, and sigma. */

                /* Update gamma */
                for (i=0; i<n; i++) {
                    wmu = 0.0; // Initialize wmu
                    for (j=0; j<m; j++) {
                        wmu += post[i + n*j]*mu[j];
                    }
                    yhat[i] = y[i] - wmu;
                }

                // Copy z to zCopy because z will be overwritten.
                F77_CALL(dcopy)(&np, z, &iOne, zCopy, &iOne);
                    
                /* initialize jpvt */
                for (ii=0; ii<p; ii++) {
                    jpvt[ii] = 0;
                }

                F77_CALL(dgelsy)(nn, pp, &iOne, zCopy, nn, yhat, nn, jpvt, &rcond, &rk, work2, &nwork, &info);
                if (info!=0) error("Error: info=%d\n",info);

                for (ii=0; ii<p; ii++) {
                    gamma[ii] = yhat[ii];    // Store the value of updated gamma
                }

                for (j=0; j<m; j++) {
                    alpha[j] = 0.0;     // initialize alpha
                    mu[j] = 0.0;     // initialize mu
                    /* Update alpha and compute weighted observations */
                    for (i=0; i<n; i++) {
                        alpha[j] += post[i + n*j] / n;
                        zgamma[i] = 0.0; // Initialize zgamma
                        for (ii=0; ii<p; ii++) {  // Compute zgamma with updated gamma
                            zgamma[i] += z[i + n*ii]*gamma[ii];
                        }
                        mu[j] += post[i + n*j]*(y[i] - zgamma[i]);
                    }
                    mu[j] = mu[j] / (alpha[j] * n);

                    /* If h!=0, impose upper and lower bound */
                    if (h!=0) {
                        mu[j] = fmax(mu[j],lb[j]);
                        mu[j] = fmin(mu[j],ub[j]);
                    }

                    /*  Compute the residuals and newsigma. */
                    ssr_j = 0.0;
                    for (i=0; i<n; i++) {
                        r = y[i] - mu[j] - zgamma[i];
                        ssr_j += r*r*post[i + n*j];
                    }
                    /* Update sigma */
                    sigma[j] = sqrt((ssr_j + 2.0*an*sigma0[j]*sigma0[j])  / (alpha[j] * n + 2.0*an));
                    sigma[j] = fmax(sigma[j],0.01*sigma0[j]);
                }   /* end for j=0; j<m; j++ loop updating alpha, mu and sigma */

                /* If h!=0, update the h and h+1 th element of alpha */
                if (h!=0) {
                    double alphah = (alpha[h-1]+alpha[h])/2;
                    alpha[h-1] = alphah;
                    alpha[h] = alphah;
                }

                /* Check singularity */
                for (j=0; j<m; j++) {
                    if (alpha[j] < 1e-8 || isnan(alpha[j]) || sigma[j] < 1e-8) {
                        sing = 1;
                    }
                }

                /* Exit from the loop if singular */
                if (sing) {
                    notcg[jn] = 1;
                    break;
                }

                emit++;
            } // end while loop

            penloglikset[jn] = penloglik;
            loglikset[jn] = loglik;
            for (j=0; j<m; j++) {
                alphaset[jn*m+j] = alpha[j];
                muset[jn*m+j]    = mu[j];
                sigmaset[jn*m+j] = sigma[j];
            }

        } /* end for (jn=0; jn<ninits; jn++) loop */

        return;
        } /* end normalmixpmle_z()  */

} /* end extern "C" */
