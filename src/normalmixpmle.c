#include <R.h>
#include <Rmath.h>

/* Compute the EM step and the matrix of "posterior" probabilities
   for the penalized MLE of a finite mixture of regressions
   Based on normpost.c in mixtools package.  */
void normalmixpmle(
    int *set,   /* vector of n, m, ninits, maxit */
    double *y, /* n by 1 vector of data */
    double *alphaset, /* current vector of mixing parameters */
    double *muset, /* current vector of regression coefficients */
    double *sigmaset, /* current vector of component stdevs */
    double *sigma0, /* stdevs used in the penalty term */
    double *mu0,    /* estimate of mu from m-1 component model */
    double *aan, /* constant in the penalty term */
    double *ttau,   /* the parameter that splits the h-th component */
    int *hh, /* the component that is split into two. If ==0, no penalty on tau is imposed */
    int *kk, /* the steps in EM update. If =0, no penalty on tau. If >=1, penalty on tau. If =1, constraints are imposed on mu. */
    double *lub, /* lower and upper bound on mu */
    double *work, /* 3*m-vector of workspace, which will be broken into 3 parts */
    double *post, /* n by m matrix of posterior probabilities */
    double *loglikset, /* scalar loglikelihood value */
    double *penloglikset, /* scalar penalized loglikelihood value */
    int *notcg,  /* =1 if not convergent */
    double *ttol /* tolerance level for exiting iteration */
    ) {
        int *nn = set;          /* sample size */
        int *mm = set+1;        /* number of components */
        int *nninits = set+2;   /* number of initial values */
        int *mmaxit = set+3;    /* maximum number of iterations */

        int n=*nn, m=*mm, h = *hh, k = *kk, i, j, minj=0;
        int ninits=*nninits, maxit = *mmaxit, emit, sing;
        double tol = *ttol, oldpenloglik;
        double r, rowsum, min=0.0, an=*aan, alphah, tau=*ttau, tauhat, s0j, ssr_j;
        double *AlpSigRatio = work+m; /* Second 1/3 of workspace, for frequently used constants */
        double *logAlpSigRatio = work+2*m; /* Third 1/3 of workspace, for frequently used constants */
        double *lb = lub;
        double *ub = lub + m;

        if (k==1) {  // If k==1, compute upper and lower bounds
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
            double *mu = muset + jn*m;
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
                    r = y[i]-mu[j];
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

            /* If not exit, update alpha, mu, and sigma. */
            for (j=0; j<m; j++) {
                alpha[j] = 0.0;     // initialize alpha
                mu[j] = 0.0;     // initialize mu
                /* Update alpha and mu */
                for (i=0; i<n; i++) {
                    alpha[j] += post[i + n*j] / n;
                    mu[j] += post[i + n*j] * y[i];
                }
                mu[j] = mu[j] / (alpha[j] * n);

                /* If k==1, impose upper and lower bound */
                if (k==1) {
                    mu[j] = fmax(mu[j],lb[j]);
                    mu[j] = fmin(mu[j],ub[j]);
                }

                /*  Compute the residuals and newsigma. */
                ssr_j = 0.0;
                for (i=0; i<n; i++) {
                    r = y[i] - mu[j];
                    ssr_j += r*r*post[i + n*j];
                }
                sigma[j] = sqrt((ssr_j + 2.0*an*sigma0[j]*sigma0[j])  / (alpha[j] * n + 2.0*an));
                sigma[j] = fmax(sigma[j],0.01*sigma0[j]);
            }   /* end for j=0; j<m; j++ loop updating alpha, mubeta and sigma */
                
            /* if k!=0, update alpha and/or tau */
            if (k!=0){
                alphah = (alpha[h-1]+alpha[h]);
                /* If k!=1, update tau. If k==1, no update of tau. */
                if (k!=1) {
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
            muset[jn*m+j]    = mu[j];
            sigmaset[jn*m+j] = sigma[j];
        }

    } /* end for (jn=0; jn<ninits; jn++) loop */

    return;
} /* end normalmixpmle()  */

