#include <R.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Arith.h>

/* Compute the EM step and the matrix of "posterior" probabilities
   for the penalized MLE of a finite mixture of regressions
   Based on normpost.c in mixtools package.  */
extern "C" {
    void regmixpmle(
        int *set,   /* vector of n, m, k, ninits, maxit, jpvt */
        double *y, /* n by 1 vector of dependent varaible */
        double *x, /* n by k matrix of explanatory variables */
        double *alphaset, /* current vector of mixing parameters */
        double *mubetaset, /* current vector of regression coefficients */
        double *sigmaset, /* current vector of component stdevs */
        double *sigma0, /* stdevs used in the penalty term */
        double *mu0,    /* estimate of mu from m-1 component model */
        double *aan, /* constant in the penalty term */
		double *ttau,   /* the parameter that splits the h-th component */
        int *hh, /* the component that is split into two. If =0, unconstrained PMLE is computed. */
		int *ttaupenoption,  /* the steps in EM update. If =0, no penalty on tau. If >=1, penalty on tau. If =1, constraints are imposed on mu. */
        double *lub, /* lower and upper bound on mu */
        double *work, /* 3*m-vector of workspace, which will be broken into 3 parts */
        double *post, /* n by m matrix of posterior probabilities */
        double *loglikset, /* scalar loglikelihood value */
        double *penloglikset, /* scalar penalized loglikelihood value */
        int *notcg,  /* =1 if not convergent */
        double *ttol, /* tolerance level for exiting iteration */
        
		double *wxy /* n by (1+k) vector of weighted dependent varaible, wy and wx */
        ) {
            int *nn = set;          /* sample size */
            int *mm = set+1;        /* number of components */
            int *kk = set+2;        /* number of explanatory variables */
            int *nninits = set+3;   /* number of initial values */
            int *mmaxit = set+4;    /* maximum number of iterations */
            int *jpvt = set+5;      /* pivots used in dgelsy */

            int n=*nn, m=*mm, k=*kk, h = *hh, taupenoption = *ttaupenoption, i, j, ii, minj=0, info, rk;
            int ninits=*nninits, maxit = *mmaxit, emit, sing;
            double tol = *ttol, oldpenloglik;
            double xtheta=0.0, r, rowsum, min=0.0, an=*aan, alphah, tau=*ttau, tauhat, worksize, s0j, ssr_j;
            double *AlpSigRatio = work+m; /* Second 1/3 of workspace, for frequently used constants */
            double *logAlpSigRatio = work+2*m; /* Third 1/3 of workspace, for frequently used constants */
            double *lb = lub;
            double *ub = lub + m;
            double *wy = wxy;
            double *wx = wxy + n;

            /* First, make the query, store the result in 'nwork', and allocate the workspace. */
            double rcond = sqrt((2.2204460492503131e-16)/2);    //  2.22... is Machine epsilon
            const int iOne = 1;
            int query = -1;  // Code for "query mode"
            F77_CALL(dgelsy)(nn, kk, &iOne, x, nn, y, nn, jpvt, &rcond, &rk, &worksize, &query, &info);
            int nwork = static_cast<int>(worksize);
            double *work2;
            work2 = (double *) R_alloc(nwork, sizeof(double));
            if (work2==NULL) error("Cannot allocate workspace\n");

            if (taupenoption==1) {  // If taupenoption==1, compute upper and lower bounds
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
                double *mubeta = mubetaset + jn*m*k;
                double *sigma = sigmaset + jn*m;

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
                    for (j=0; j<m; j++) {
                        xtheta = 0.0; // Initialize xtheta
                        for (ii=0; ii<k; ii++) {  // Compute xtheta
                            xtheta += x[i + n*ii]*mubeta[ii + k*j];
                        }
                        r = y[i]-xtheta;
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
                penloglik = loglik + log(2.0) + fmin(log(tau),log(1-tau));
                for (j=0; j<m; j++) {
                    s0j = sigma0[j]/sigma[j];
                    penloglik += -an*(s0j*s0j - 2.0*log(s0j) -1.0);
                }

                diff = penloglik - oldpenloglik;
                oldpenloglik = penloglik;

                /* Normal exit */
                    if (diff < tol || emit>=maxit){
                        break;
                    }

                /* If not exit, update alpha, mubeta, and sigma. */
                for (j=0; j<m; j++) {
                    alpha[j] = 0.0;     // initialize alpha
                    /* Update alpha and compute weighted observations */
                    for (i=0; i<n; i++) {
                        alpha[j] += post[i + n*j] / n;
                        wy[i] = sqrt(post[i + n*j])*y[i];
                        for (ii=0; ii<k; ii++) {
                            wx[i + ii*n] = sqrt(post[i + n*j])*x[i + n*ii];
                        }
                    }
                    /* initialize jpvt */
                    for (ii=0; ii<k; ii++) {
                        jpvt[ii] = 0;
                    }

                    /* Update mubeta */
                    F77_CALL(dgelsy)(nn, kk, &iOne, wx, nn, wy, nn, jpvt, &rcond, &rk, work2, &nwork, &info);
                    if (info!=0) error("Error: info=%d\n",info);

                    for (ii=0; ii<k; ii++) {
                        mubeta[ii + k*j] = wy[ii];    // Store the value of updated mubeta
                    }
                    /* If taupenoption==1, impose upper and lower bound */
                    if (taupenoption==1) {
                        mubeta[k*j] = fmax(mubeta[k*j],lb[j]);
                        mubeta[k*j] = fmin(mubeta[k*j],ub[j]);
                    }

                    /*  Compute the residuals and newsigma. Note that 'theta_hat' values are now in 'wy'. */
                    ssr_j = 0.0;
                    for (i=0; i<n; i++) {
                        xtheta = 0.0; // Initialize and reuse xtheta
                        for (ii=0; ii<k; ii++){
                            xtheta += x[i + n*ii]*wy[ii];
                        }
                        r = y[i] - xtheta;
                        ssr_j += r*r*post[i + n*j];
                    }
                    sigma[j] = sqrt((ssr_j + 2.0*an*sigma0[j]*sigma0[j])  / (alpha[j] * n + 2.0*an));
                    sigma[j] = fmax(sigma[j],0.01*sigma0[j]);
                }   /* end for j=0; j<m; j++ loop updating alpha, mubeta and sigma */

				/* if taupenoption!=0, update alpha and/or tau */
				if (taupenoption!=0){
					alphah = (alpha[h-1]+alpha[h]);
					/* If taupenoption!=1, update tau. If taupenoption==1, no update of tau. */
					if (taupenoption!=1) {
						tauhat = alpha[h-1]/(alpha[h-1]+alpha[h]);
						if (tauhat <= 0.5) {
							tau = fmin((alpha[h-1]*n + 1.0)/(alpha[h-1]*n + alpha[h]*n + 1.0), 0.5);
						} else {
							tau = fmax(alpha[h-1]*n /(alpha[h-1]*n + alpha[h]*n + 1.0), 0.5);
						}
					}
					/* Using tau, revise the h and h+1 th element of alpha */
					alpha[h-1] = alphah*tau;
					alpha[h] = alphah*(1-tau);
				}


                /* Check singularity */
                for (j=0; j<m; j++) {
                    if (alpha[j] < 1e-8 || isnan(alpha[j]) || sigma[j] < 1e-8){
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
                sigmaset[jn*m+j] = sigma[j];
            }
            for (j=0; j<k*m; j++) {
                mubetaset[jn*m*k+j] = mubeta[j];
            }

        } /* end for (jn=0; jn<ninits; jn++) loop */

        return;
        } /* end regmixpmle()  */

} /* end extern "C" */
