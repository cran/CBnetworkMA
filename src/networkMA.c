/****************************************************************************************
*
*
* My attempt to produce an MCMC algorithm for the three models we propose
# in the BNP Network meta-analysis paper using the .Call function
*
*
****************************************************************************************/
#include "matrix.h"
#include "Rutil.h"

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>

#include <R_ext/Lapack.h>
#include <R_ext/Utils.h>

#include <math.h>
#include <stdio.h>
#include <time.h>


// inputs
// r - number of successes
// n - number of trials
// tid - treatment id
// sid - study id
// N - number of studies
// K - total number of treatments across all studies
// model - integer indicating which model to fit
//         1 - gaussian
//         2 - DP with gaussian base measure
//         3 - DP with spike and slab base measure
// H - upper bound on the number of components on the random probability measure
// modelPriors - priors values that I will fill out shortly
// mh - vector containing candidate density sd parameters for metropolis updates
// verbose - logical indicating if information should be print to screen
// niter - total number of MCMC iterates that will be sampled
// nburn - number of MCMC iterates to be discarded as burn-in
// nthin - amount of thinning applied to the MCMC chain
//
// outputs
// mu - matrix containing MCMC draws for study-specific intercepts
// delta - matrix containing MCMC draws for component variances
// tau - matrix containing MCMC draws for component weights
// ci - matrix containing MCMC draws for component labels
// lpml - matrix containing MCMC draws for component labels


static void network_ma(
                  double *r_mat, double *n_mat, int *tid_mat, int *tid_logical,
                  int *nobs, int *N, int *K, int *H,
                  int *model,
                  double *modelPriors, int *tau_prior, double *mh,
                  int *verbose,
                  int *niter, int *nburn, int *nthin,
                  double *mu, double *delta, double *tau2, double *d1,
                  int *ci, double *omega, double *sh, double *lpml){
//
//                  double *ispred, int *isordpred, double *ppred, int *predclass, int *ordppred,
//                  double *rbpred, int *rbordpred, double *predclass_prob){



  // i - MCMC index
  // j - study index
  // k - delta index
  // kk - second delta index
  // d -  d1 index
  // h - d1star index (up to *H)
  // ii - save MCMC iterates index


  int i, j, k, kk, h, d, ii;
  int nout = (*niter - *nburn)/(*nthin);
  int npwcomp = factorial(*K)/(factorial(2)*factorial(*K-2));

//  Rprintf("npwcomp = %d\n", npwcomp);
//  RprintVecAsMat("r_mat", r_mat, *N, npwcomp+1);
//  RprintVecAsMat("n_mat", n_mat, *N, npwcomp+1);
//  RprintIVecAsMat("tid_mat", tid_mat, *N, npwcomp+1);
//  RprintIVecAsMat("tid_logical", tid_logical, *N, npwcomp+1);
//  Rprintf("nobs = %d\n", *nobs);
//  Rprintf("N = %d\n", *N);
//  Rprintf("K = %d\n", *K);
//  Rprintf("tau_prior = %d\n", *tau_prior);
//  Rprintf("model = %d\n", *model);
  // =====================================================================================
  //
  // Memory vectors to hold a single MCMC iterate
  //
  // =====================================================================================
  double *_mu = R_VectorInit(*N,0.0);
  double *_delta = R_VectorInit((*N)*(npwcomp+1),-99.0);
  double *_d1 = R_VectorInit(*K, 0.0);
  double *_d1star = R_VectorInit(*H, 0.0);
  double *_V = R_VectorInit(*H, 0.5); _V[*H] = 1.0;
  double *_pi_vec = R_VectorInit(*H, 0);
  double *_sh = R_VectorInit(*H, 1); _sh[0] = 0.1;
  double _tau2 = 0.5*modelPriors[2]*modelPriors[2];
  double _omega = 0.5;

  int _ci[*K]; _ci[0] = 1;

  // Initialize variables
  for(j = 0; j < *N; j++){
//    _mu[j] = rnorm(0,0.5);
    _mu[j] = 0.0;
    _delta[j*(npwcomp+1) + 0] = 0;
    for(k = 1; k < (npwcomp+1); k++){
      if(tid_logical[j*(npwcomp+1)+k] == 1){
        _delta[j*(npwcomp+1) + k] = 0;
      }
    }
  }

  double tmp=1.0;
  _pi_vec[0] = _V[0];
  tmp = 1.0;
  for(h = 1; h < *H; h++){
    tmp = tmp*(1-_V[h-1]);
    _pi_vec[h] = _V[h]*tmp;
    _d1star[h] = rnorm(modelPriors[5], modelPriors[6]);
  }
//  RprintVecAsMat("pi_vec", _pi_vec, 1, *H);

  for(d = 1; d < *K; d++){
    _ci[d] = rbinom(2, 0.25) + 2;
    _d1[d] = rnorm(modelPriors[5],modelPriors[6]);
  }


//  RprintVecAsMat("delta", _delta, *N, (npwcomp+1));
//  RprintIVecAsMat("ci", _ci, 1, *K);

  // ===================================================================================
  //
  // scratch vectors of memory needed to update parameters
  //
  // ===================================================================================

  // stuff that I need to update Si (cluster labels);
  int pwj;
  double llo,lln, llr, uu;

  // stuff that I need to update mu
  double muold, munew, piold, pinew;

  // stuff that I need to update delta
  double dold, dnew, ld;
  double *del_veco = R_VectorInit(*K, 0.0);
  double *del_vecn = R_VectorInit(*K, 0.0);
  double *d1mn_vec = R_VectorInit(*K, 0.0);
  double *S = R_VectorInit((*K)*(*K), 0.0);
  double *scr1 = R_VectorInit((*K)*(*K), 0.0);
  double *scr2 = R_VectorInit((*K)*(*K), 0.0);

  // stuff that I need to update tau2
  double tauo, taun, ldo, ldn;
  double *del_vec = R_VectorInit(*K, 0.0);
  double *So = R_VectorInit((*K)*(*K), 0.0);
  double *Sn = R_VectorInit((*K)*(*K), 0.0);

  // stuff I need to update d1 and d1star (for models 2 and 3)
  int study_used_d1k;
  double d1so, d1sn, d1k, d1_baseline;
  double *d1mn_veco = R_VectorInit(*K, 0.0);
  double *d1mn_vecn = R_VectorInit(*K, 0.0);

  // stuff I need to update ci (for models 2 and 3)
  double lup, maxp, denp, cprob;
  double *pp = R_VectorInit(*H, 0.0);
  double *prob = R_VectorInit(*H, 0.0);

  // stuff I need to update V
  double ccount, cumccount, astar, bstar;

  // stuff I need to update omega
  double sum1, sumv0;

  // Stuff to compute lpml, likelihood, WAIC, and Rao-Blackwellized density values
  double _lpml;
//    double elppdWAIC;
  double *CPOinv = R_VectorInit(*N, 0.0);
//  double *mnlike = R_VectorInit(*N, 0.0);
//  double *mnllike = R_VectorInit(*N, 0.0);


  // ===================================================================================
  //
  // Prior distribution parameter values
  //
  // ===================================================================================

  // priors for mu_{i, b_i} (baseline mean)
  double mb = modelPriors[0]; double sb = modelPriors[1];

  // prior for tau if tau ~ UN(0, tau_max)
  double tau_max = modelPriors[2];

  // prior for tau if tau ~ LN(tau_lm, tau_lsd)
  double tau_lm = modelPriors[3], tau_lsd = modelPriors[4];

  // prior values d1 or d1star for dp_gaussian model;
  double md=modelPriors[5], sd=modelPriors[6];

  // DP scale parameter
  double alpha=modelPriors[7];

  // Prior from the nonlocal density if model = 3;
  double kap0=modelPriors[8], scale=modelPriors[9], nu=modelPriors[10];

  // Prior for the spike component of the spike and slab
  double v0=modelPriors[11];

  // prior for omega (probability in spike and slab model)
  double aw = modelPriors[12], bw=modelPriors[13];

//  Rprintf("kap0 = %f\n", kap0);
//  Rprintf("scale = %f\n", scale);
//  Rprintf("aw = %f\n", aw);
//  Rprintf("bw = %f\n", bw);

  // M-H tuning parameters
  double csigMU=mh[0], csigDEL=mh[1], csigTAU=mh[2], csigD1=mh[3];

  if(*verbose){
    Rprintf("N = %d\n", *N);
    if(*tau_prior==1)Rprintf("Prior for tau ~ UN(0, tau_max) \n");
    if(*tau_prior==2)Rprintf("Prior for tau ~ logN(tau_lsm, tau_lsd) \n");
    Rprintf("Prior values: mb = %.2f, sb = %.2f md = %.2f, sd = %.2f\n", mb, sb, md, sd);
  }




  ii = 0;

  // ===================================================================================
  //
  // Beginning of the MCMC loop
  //
  // ===================================================================================
  double calc_time = 0.0;
  clock_t  begin = clock();


  for(i=0; i<*niter; i++){

//    Rprintf("i ====================================== %d\n", i);

    if(*verbose){
      clock_t ith_iterate = clock();
      calc_time = (ith_iterate - begin)/CLOCKS_PER_SEC;

      Rprintf("Progress:%.1f%%, Time:%.1f seconds\r", ((double) (i+1) / (double) (*niter))*100.0, calc_time);
    }

    //////////////////////////////////////////////////////////////////////////////////
    //
    //      update mu_{i, b_i}
    //
    //////////////////////////////////////////////////////////////////////////////////

    for(j = 0; j < *N; j++){
//      Rprintf("j = %d\n", j);

      muold = _mu[j];
      munew = rnorm(muold, csigMU);

//      Rprintf("muold = %f\n", muold);
//      Rprintf("munew = %f\n", munew);

      llo = 0.0; lln = 0.0;

      for(k = 0; k < (npwcomp+1); k++){ // loop through each pairwise comparison
        if(tid_logical[j*(npwcomp+1)+k] == 1){  // if kth pairwise comparison observed in jth study
//          Rprintf("k = %d\n", k);
//          Rprintf("tid_logical[j*(npwcomp+1)+k] = %d\n", tid_logical[j*(npwcomp+1)+k]);
//          Rprintf("_delta[j*(npwcomp+1)+k] = %f\n", _delta[j*(npwcomp+1)+k]);

          piold = exp(muold + _delta[j*(npwcomp+1)+k])/(1 + exp(muold + _delta[j*(npwcomp+1)+k]));
          pinew = exp(munew + _delta[j*(npwcomp+1)+k])/(1 + exp(munew + _delta[j*(npwcomp+1)+k]));

//          Rprintf("piold = %f\n", piold);
//          Rprintf("pinew = %f\n", pinew);

//         Rprintf("r_mat[j*(npwcomp+1)+k] = %f\n", r_mat[j*(npwcomp+1)+k]);
//         Rprintf("n_mat[j*(npwcomp+1)+k] = %f\n", n_mat[j*(npwcomp+1)+k]);
//         Rprintf("dbinom(r, n, piold, 1) = %f\n", dbinom(r_mat[j*(npwcomp+1)+k], n_mat[j*(npwcomp+1)+k], piold, 1));

          llo = llo + dbinom(r_mat[j*(npwcomp+1)+k], n_mat[j*(npwcomp+1)+k], piold, 1);
          lln = lln + dbinom(r_mat[j*(npwcomp+1)+k], n_mat[j*(npwcomp+1)+k], pinew, 1);

//          Rprintf("llo = %f\n", llo);
//          Rprintf("lln = %f\n", lln);
        }
      }

      llo = llo + dnorm(muold, mb, sb, 1);
      lln = lln + dnorm(munew, mb, sb, 1);

//      Rprintf("llo = %f\n", llo);
//      Rprintf("lln = %f\n", lln);

      llr = lln - llo;
//      Rprintf("llr = %f\n", llr);

      uu = runif(0,1);
//      Rprintf("uu = %f\n", uu);

      if(llr > log(uu)) _mu[j] = munew;
//      Rprintf("_mu[j] = %f\n", _mu[j]);

    }

//    RprintVecAsMat("_mu", _mu, 1, *N);

    //////////////////////////////////////////////////////////////////////////////////
    //
    // update the deltas. Note if the study contains more than one arm
    // (i.e., more than two treatments), then these are modeled with
    // a multivariate normal.
    //
    //////////////////////////////////////////////////////////////////////////////////
    for(j = 0; j < *N; j++){
//      Rprintf("j = %d\n", j);
      // Need to update delta vector for the jth study
      kk = 0;
      // Notice that when using tid_mat[j*(npwcomp+1)+0] to index a vector we must subtract
      // by 1 to accommodate C's storage.
      for(k = 1; k < (npwcomp+1); k++){ // loop through each pairwise comparison skip delta11
        if(tid_logical[j*(npwcomp+1)+k] == 1){
          del_veco[kk] = _delta[j*(npwcomp+1)+k];
          del_vecn[kk] = _delta[j*(npwcomp+1)+k];
          d1mn_vec[kk] = _d1[tid_mat[j*(npwcomp+1)+k]-1] - _d1[tid_mat[j*(npwcomp+1)+0]-1]; // d1k - d1_baseline
          kk = kk+1;
        }
      }
//      Rprintf("kk = %d\n", kk);
//      RprintVecAsMat("del_veco", del_veco, 1, kk);
//      RprintVecAsMat("del_vecn", del_vecn, 1, kk);

      // create the S matrix
      for(k = 0; k < (kk*kk); k++){
        S[k] = _tau2*0.5;
        if(k % (kk + 1) == 0) S[k] = _tau2;
      }

//      RprintVecAsMat("S", S, kk, kk);
      cholesky(S, kk, &ld);
	  inverse_from_cholesky(S, scr1, scr2, kk);
//      RprintVecAsMat("S", S, kk, kk);


      pwj = kk;
      kk = 0;
      for(k = 1; k < (npwcomp+1); k++){ // loop through each pairwise comparison.  Need to skip delta11 so it stays zero;
        if(tid_logical[j*(npwcomp+1)+k] == 1){  // if kth pairwise comparison observed in jth study
//          Rprintf("k = %d\n", k);
          dold = _delta[j*(npwcomp+1)+k];
          dnew = rnorm(dold, csigDEL);

 //         Rprintf("dold = %f\n", dold);
//          Rprintf("dnew = %f\n", dnew);
//          Rprintf("mu[j] = %f\n", _mu[j]);
          piold = exp(_mu[j] + dold)/(1 + exp(_mu[j] + dold));
          pinew = exp(_mu[j] + dnew)/(1 + exp(_mu[j] + dnew));

//          Rprintf("piold = %f\n", piold);
//          Rprintf("pinew = %f\n", pinew);

//          Rprintf("r_mat[j*(npwcomp+1)+k] = %f\n", r_mat[j*(npwcomp+1)+k]);
//          Rprintf("n_mat[j*(npwcomp+1)+k] = %f\n", n_mat[j*(npwcomp+1)+k]);

          llo = dbinom(r_mat[j*(npwcomp+1)+k], n_mat[j*(npwcomp+1)+k], piold, 1);
          lln = dbinom(r_mat[j*(npwcomp+1)+k], n_mat[j*(npwcomp+1)+k], pinew, 1);

//          Rprintf("llo = %f\n", llo);
//          Rprintf("lln = %f\n", lln);

          // Now I need to evaluate the MVN density for delta.

          del_vecn[kk] = dnew;

//          RprintVecAsMat("del_veco", del_veco, 1, pwj);
//          RprintVecAsMat("del_vecn", del_vecn, 1, pwj);
//          RprintVecAsMat("d1mn_vec", d1mn_vec, 1, pwj);


          llo = llo + dmvnorm(del_veco, d1mn_vec, S, pwj, ld, scr1, 1);
          lln = lln + dmvnorm(del_vecn, d1mn_vec, S, pwj, ld, scr1, 1);
//          Rprintf("dmvnorm(del_veco, d1mn_vec, S, pwj, ld, scr1, 1); = %f\n", dmvnorm(del_veco, d1mn_vec, S, pwj, ld, scr1, 1));
//          Rprintf("dmvnorm(del_vecn, d1mn_vec, S, pwj, ld, scr1, 1); = %f\n", dmvnorm(del_vecn, d1mn_vec, S, pwj, ld, scr1, 1));
//          Rprintf("llo = %f\n", llo);
//          Rprintf("lln = %f\n", lln);

          llr = lln - llo;
//          Rprintf("llr = %f\n", llr);

          uu = runif(0,1);
//          Rprintf("uu = %f\n", uu);

          del_vecn[kk] = dold; // Resort to old delta vector unless proposal is accepted.
          if(llr > log(uu)){
             _delta[j*(npwcomp+1)+k] = dnew;
             del_veco[kk] = dnew; // I need to update the delta vector as well.
             del_vecn[kk] = dnew;
          }
//          RprintVecAsMat("del_veco", del_veco, 1, pwj);
//          RprintVecAsMat("del_vecn", del_vecn, 1, pwj);
          kk = kk + 1;
        }
      }
    }
//    RprintVecAsMat("delta", _delta, *N, (npwcomp+1));



    //////////////////////////////////////////////////////////////////////////////////
    //
	//	# update tau2
	//
    //////////////////////////////////////////////////////////////////////////////////
    tauo = sqrt(_tau2);
    taun = rnorm(tauo, csigTAU);
//    Rprintf("tauo = %f\n", tauo);
//    Rprintf("taun = %f\n", taun);
    if(taun > 0){

      llo = 0.0; lln=0.0;
      for(j = 0; j < *N; j++){

//        Rprintf("j = %d\n", j);
        kk = 0;
        for(k = 1; k < (npwcomp+1); k++){ // loop through each pairwise comparison skip delta11
          if(tid_logical[j*(npwcomp+1)+k] == 1){
            del_vec[kk] = _delta[j*(npwcomp+1)+k];

//            Rprintf("k = %d\n", k);
//            Rprintf("tid_mat[j*(npwcomp+1)+0] = %d\n", tid_mat[j*(npwcomp+1)+0]);
//            Rprintf("tid_mat[j*(npwcomp+1)+k] = %d\n", tid_mat[j*(npwcomp+1)+k]);
//            Rprintf("_delta[j*(npwcomp+1)+k] = %f\n", _delta[j*(npwcomp+1)+k]);
//            Rprintf("_d1[tid_mat[j*(npwcomp+1)+0]] = %f\n", _d1[tid_mat[j*(npwcomp+1)+0]-1]);
//            Rprintf("_d1[tid_mat[j*(npwcomp+1)+k]] = %f\n", _d1[tid_mat[j*(npwcomp+1)+k]-1]);

            d1mn_vec[kk] = _d1[tid_mat[j*(npwcomp+1)+k]-1] - _d1[tid_mat[j*(npwcomp+1)+0]-1]; // d1k - d1_baseline
            kk = kk+1;

//            Rprintf("kk = %d\n", kk);
          }
        }

//        RprintVecAsMat("del_vec", del_vec, 1, kk);
//        RprintVecAsMat("d1mn_vec", d1mn_vec, 1, kk);

        for(k = 0; k < (kk*kk); k++){
          So[k] = tauo*tauo*0.5;
          Sn[k] = taun*taun*0.5;
          if(k % (kk + 1) == 0){
            So[k] = tauo*tauo;
            Sn[k] = taun*taun;
          }
        }

//        RprintVecAsMat("So", So, kk, kk);
        cholesky(So, kk, &ldo);
	    inverse_from_cholesky(So, scr1, scr2, kk);

//        RprintVecAsMat("Sn", Sn, kk, kk);
        cholesky(Sn, kk, &ldn);
	    inverse_from_cholesky(Sn, scr1, scr2, kk);

//	      Rprintf("dmvnorm(del_vec, d1mn_vec, So, kk, ldo, scr1, 1)")
//	      Rprintf("dmvnorm(del_vec, d1mn_vec, Sn, kk, ldn, scr1, 1)")

        llo = llo + dmvnorm(del_vec, d1mn_vec, So, kk, ldo, scr1, 1);
        lln = lln + dmvnorm(del_vec, d1mn_vec, Sn, kk, ldn, scr1, 1);

//        Rprintf("llo = %f\n", llo);
//        Rprintf("lln = %f\n", lln);
      }

      if(*tau_prior == 1){
        llo = llo + dunif(tauo, 0.0, tau_max, 1);
        lln = lln + dunif(tauo, 0.0, tau_max, 1);
      }
      if(*tau_prior == 2){
//        Rprintf("dlnorm(tauo, tau_lm, tau_lsd, 1) = %f\n", dlnorm(tauo, tau_lm, tau_lsd, 1));
//        Rprintf("dlnorm(taun, tau_lm, tau_lsd, 1) = %f\n", dlnorm(taun, tau_lm, tau_lsd, 1));
        llo = llo + dlnorm(tauo*tauo, tau_lm, tau_lsd, 1);
        lln = lln + dlnorm(taun*taun, tau_lm, tau_lsd, 1);
      }

//      Rprintf("llo = %f\n", llo);
//      Rprintf("lln = %f\n", lln);
      llr = lln - llo;
//      Rprintf("llr = %f\n", llr);

      uu = runif(0,1);
//      Rprintf("uu = %f\n", uu);

//      Rprintf("tauo = %f\n", tauo);
//      Rprintf("taun = %f\n", taun);

      if(llr > log(uu)) _tau2 = taun*taun;
//      Rprintf("_tau2 = %f\n", _tau2);
    }

    //////////////////////////////////////////////////////////////////////////////////////
    //
    //  # update the the kth entry of d1k.
    //  Note that I use a Metropolis algorithm here
    //  Note that d11 = 0 and is not updated.
    //
    //  Note that the updating of d1 depends on the model
    //  model 1 - Gaussian prior for d1
    //  model 2 - DP prior with a Gaussian base measure
    //  model 3 - DP prior with a spike & slab base measure
    //
    //////////////////////////////////////////////////////////////////////////////////////
    if(*model==1){
      for(d = 1; d < *K; d++){

        d1so = _d1[d];
        d1sn = rnorm(d1so, csigD1);

        llo=0.0, lln=0.0;

        for(j=0; j<*N; j++){

          study_used_d1k=0;

          // First need to check of d1h is the baseline treatment
          // Notice that I compare to h+1 because d11 = 0 always and so it is not updated.
          // Need to determine if treatment h is baseline treatment or not

          kk = 0;
          for(k = 1; k < (npwcomp+1); k++){ // loop through each pairwise comparison skip delta11
            if(tid_logical[j*(npwcomp+1)+k] == 1){
              del_vec[kk] = _delta[j*(npwcomp+1)+k];
              d1mn_veco[kk] = _d1[tid_mat[j*(npwcomp+1)+k]-1] - _d1[tid_mat[j*(npwcomp+1)+0]-1]; // d1k - d1_baseline
              d1mn_vecn[kk] = _d1[tid_mat[j*(npwcomp+1)+k]-1] - _d1[tid_mat[j*(npwcomp+1)+0]-1]; // d1k - d1_baseline
              if(tid_mat[j*(npwcomp+1)+k] == d+1){
                d1mn_vecn[kk] = d1sn - _d1[tid_mat[j*(npwcomp+1)+0]-1];
                study_used_d1k=1;
              }
              if(tid_mat[j*(npwcomp+1)+0] == d+1){
                d1mn_vecn[kk] = _d1[tid_mat[j*(npwcomp+1)+k]-1] - d1sn;
                study_used_d1k=1;
              }
              kk = kk+1;
            }
          }

          if(study_used_d1k){
            for(k = 0; k < (kk*kk); k++){
              S[k] = (_tau2)*0.5;
              if(k % (kk + 1) == 0) S[k] = _tau2;
            }

            cholesky(S, kk, &ld);
	        inverse_from_cholesky(S, scr1, scr2, kk);

            llo = llo + dmvnorm(del_vec, d1mn_veco, S, kk, ld, scr1, 1);
            lln = lln + dmvnorm(del_vec, d1mn_vecn, S, kk, ld, scr1, 1);
          }
        }


        llo = llo + dnorm(d1so, md, sd, 1);
        lln = lln + dnorm(d1sn, md, sd, 1);

        llr = lln - llo;
//        Rprintf("llr = %f\n", llr);

        uu = runif(0,1);
//        Rprintf("uu = %f\n", uu);

        if(llr > log(uu)) _d1[d] = d1sn;
      }

    }


    /////////////////////////////////////////////////////////////////////////////////////
    //
    // Update the d1stars.  This model employs a DP with Gaussian base measure model 2.
    //
    /////////////////////////////////////////////////////////////////////////////////////
    if((*model==2) | (*model==3)){
      for(h = 1; h < *H; h++){
        d1so = _d1star[h];
        d1sn = rnorm(d1so, csigD1);

        llo=0.0, lln=0.0;
        for(j=0; j<*N; j++){

          // Need to check if there is a treatment in the jth study such that
          // treatment k is allocated to component d, i.e., ck = d with k in Ti
          // recall that tid_mat is matrix that identifies which treatments were used
          // in the jth study.
          //
          // This includes the baseline treatment and for that reason
          // we begin the for-loop at 0;


          // If study does have c_k = d, then need to include it when updating _d1star_d;
          study_used_d1k=0;
          kk = 0; // kk becomes the dimension of N(delta_i d1star_i S);
          for(k = 1; k < (npwcomp+1); k++){ // loop through each pairwise comparison skip delta11
            if(tid_logical[j*(npwcomp+1)+k] == 1){
//              Rprintf("h = %d\n", h);
//              Rprintf("j = %d\n", j);
//              Rprintf("k = %d\n", k);
//              RprintIVecAsMat("ci", _ci, 1, *K);
//              RprintVecAsMat("_d1star", _d1star, 1, *H);
//              Rprintf("d1sn = %f\n", d1sn);
              del_vec[kk] = _delta[j*(npwcomp+1)+k];

              d1k = _d1star[_ci[tid_mat[j*(npwcomp+1)+k]-1]-1];
              d1_baseline = _d1star[_ci[tid_mat[j*(npwcomp+1)+0]-1]-1];

//              Rprintf("tid_mat[j*(npwcomp+1)+k] = %d\n",tid_mat[j*(npwcomp+1)+k]);
//              Rprintf("tid_mat[j*(npwcomp+1)+0] = %d\n",tid_mat[j*(npwcomp+1)+0]);
//              Rprintf("_ci[tid_mat[j*(npwcomp+1)+k]-1] = %d\n", _ci[tid_mat[j*(npwcomp+1)+k]-1]);
//              Rprintf("_ci[tid_mat[j*(npwcomp+1)+0]-1]  = %d\n",_ci[tid_mat[j*(npwcomp+1)+0]-1] );
//              Rprintf("d1k = %f\n", d1k);
//              Rprintf("d1baseline = %f\n", d1_baseline);

              d1mn_veco[kk] =  d1k - d1_baseline; // d1k - d1_baseline
              d1mn_vecn[kk] =  d1k - d1_baseline; // d1k - d1_baseline
              if(_ci[tid_mat[j*(npwcomp+1)+k]-1] == h+1){// ci[0]=1 always and so cluster labeling starts at two and hence put h+1
                d1mn_vecn[kk] = d1sn - d1_baseline;
                study_used_d1k=1;
              }
              if(_ci[tid_mat[j*(npwcomp+1)+0]-1] == h+1){// ci[0]=1 always and so cluster labeling starts at two and hence put h+1
                d1mn_vecn[kk] = d1k - d1sn;
                study_used_d1k=1;
              }
              kk = kk+1;
            }
          }
//          Rprintf("study_used_d1k = %d\n", study_used_d1k);
//          RprintVecAsMat("d1mn_veco", d1mn_veco, 1, kk);
//          RprintVecAsMat("d1mn_vecn", d1mn_vecn, 1, kk);
//          Rprintf("llo = %f\n", llo);
//          Rprintf("lln = %f\n", lln);
          if(study_used_d1k){
            for(k = 0; k < (kk*kk); k++){
              S[k] = (_tau2)*0.5;
              if(k % (kk + 1) == 0) S[k] = _tau2;
            }

            cholesky(S, kk, &ld);
	          inverse_from_cholesky(S, scr1, scr2, kk);

            llo = llo + dmvnorm(del_vec, d1mn_veco, S, kk, ld, scr1, 1);
            lln = lln + dmvnorm(del_vec, d1mn_vecn, S, kk, ld, scr1, 1);

//            Rprintf("kk = %d\n", kk);
//            RprintVecAsMat("del_vec", del_vec, 1, kk);
          }
        }

//        Rprintf("llo = %f\n", llo);
//        Rprintf("lln = %f\n", lln);
        if(*model==2){
          llo = llo + dnorm(d1so, md, sd, 1);
          lln = lln + dnorm(d1sn, md, sd, 1);
        }
        if(*model==3){// spike and slabe base measure using a nlp slab and gaussian spike;
//          Rprintf("d1so = %f\n", d1so);
//          Rprintf("kap0 = %f\n", kap0);
//          Rprintf("scale = %f\n", scale);
//          Rprintf("nu = %f\n", nu);
//          Rprintf("nlp density = %f\n", dnlp(d1so, 0, kap0, scale, nu, 1));
          if(_sh[h] == 1){
            llo = llo + dnlp(d1so, 0, kap0, scale, nu, 1);
            lln = lln + dnlp(d1sn, 0, kap0, scale, nu, 1);
          } else {
            llo = llo + dnorm(d1so, 0, v0/3, 1);
            lln = lln + dnorm(d1sn, 0, v0/3, 1);
          }
        }

        llr = lln - llo;
//        Rprintf("llr = %f\n", llr);

        uu = runif(0,1);
//        Rprintf("uu = %f\n", uu);

        if(llr > log(uu)) _d1star[h] = d1sn;
      }

      // Recreate the d1 vector (recall the first entry is always a zero as d11=0)
      for(d = 1; d < *K; d++){
        _d1[d] = _d1star[_ci[d]-1];
      }
//      RprintVecAsMat("d1", _d1, 1, *K);
//      _d1[0] = 0; _d1[1] = 0.306; _d1[2] = 0; _d1[3] = 0; _d1[4] = 0.306; _d1[5] = 0.306;

      ////////////////////////////////////////////////////////////////////////////////////
      //
      // Now I need to update ci for effect 2 through *H
      //
      ////////////////////////////////////////////////////////////////////////////////////
      for(d = 1; d < *K; d++){

        for(h=1; h<*H; h++){ // ci can be any of 2,...,*H clusters (ci[0] = 1 always)
          lup = 0.0;
          for(j=0; j<*N; j++){

            study_used_d1k=0;
            kk = 0; // kk becomes the dimension of N(delta_i d1star_i S);
            for(k = 1; k < (npwcomp+1); k++){ // loop through each pairwise comparison skip delta11
              if(tid_logical[j*(npwcomp+1)+k] == 1){
//              Rprintf("d = %d\n", d);
//              Rprintf("h = %d\n", h);
//              Rprintf("j = %d\n", j);
//              Rprintf("k = %d\n", k);
//              RprintIVecAsMat("ci", _ci, 1, *K);
//              RprintVecAsMat("_d1star", _d1star, 1, *H);

                del_vec[kk] = _delta[j*(npwcomp+1)+k];

                d1k = _d1star[_ci[tid_mat[j*(npwcomp+1)+k]-1]-1];
                d1_baseline = _d1star[_ci[tid_mat[j*(npwcomp+1)+0]-1]-1];

//              Rprintf("tid_mat[j*(npwcomp+1)+k] = %d\n",tid_mat[j*(npwcomp+1)+k]);
//              Rprintf("tid_mat[j*(npwcomp+1)+0] = %d\n",tid_mat[j*(npwcomp+1)+0]);
//              Rprintf("_ci[tid_mat[j*(npwcomp+1)+k]-1] = %d\n", _ci[tid_mat[j*(npwcomp+1)+k]-1]);
//              Rprintf("_ci[tid_mat[j*(npwcomp+1)+0]-1]  = %d\n",_ci[tid_mat[j*(npwcomp+1)+0]-1] );
//              Rprintf("d1k = %f\n", d1k);
//              Rprintf("d1baseline = %f\n", d1_baseline);
                d1mn_vec[kk] =  d1k - d1_baseline; // d1k - d1_baseline

                if(tid_mat[j*(npwcomp+1)+k] == d+1){// determine of jth study includes dth treatment
                  d1mn_vec[kk] = _d1star[h] - d1_baseline; // d1k - d1_baseline
                  study_used_d1k=1;
                }
                if(tid_mat[j*(npwcomp+1)+0] == d+1){ // determine of jth study includes dth treatment
                  d1mn_vec[kk] = d1k - _d1star[h]; // d1k - d1_baseline
                  study_used_d1k=1;
                }
                kk = kk+1;
              }
            }
//            Rprintf("study_used_d1k = %d\n", study_used_d1k);
//            RprintVecAsMat("d1mn_vec", d1mn_vec, 1, kk);
            if(study_used_d1k){
              for(k = 0; k < (kk*kk); k++){
                S[k] = (_tau2)*0.5;
                if(k % (kk + 1) == 0) S[k] = _tau2;
              }

              cholesky(S, kk, &ld);
	          inverse_from_cholesky(S, scr1, scr2, kk);

              lup = lup + dmvnorm(del_vec, d1mn_vec, S, kk, ld, scr1, 1);

            }
          }
          pp[h-1] = lup + log(_pi_vec[h-1]);
        }

//        RprintVecAsMat("pp", pp, 1, *H-1);
        // All these indices are strange because I am skipping the first ci as it is always
        // 1.  So we start with the second which is just a labeling trick.
        maxp = pp[0];
        for(h = 2; h < *H; h++){
          if(pp[h-1] > maxp) maxp=pp[h-1];
        }

        denp = 0.0;
        for(h = 1; h < *H; h++){
          pp[h-1] = exp(pp[h-1] - maxp);
          denp = denp + pp[h-1];
        }

        for(h = 1; h < *H; h++){
          prob[h-1] = pp[h-1]/denp;
        }
//        RprintVecAsMat("prob", prob, 1, *H-1);
        uu = runif(0.0,1.0);
//        Rprintf("uu = %f\n", uu);

        cprob = 0.0;
        for(h = 1; h < *H; h++){
          cprob = cprob + prob[h-1];
          if (uu < cprob){
            _ci[d] = h+1;  // Note that ci begins with 2 as ci[0] = 1 always since d11=0 always
            break;
          }
        }


      }
//      RprintIVecAsMat("ci", _ci, 1, *K);
      // Recreate the d1 vector (recall the first entry is always a zero as d11=0)
      for(d = 1; d < *K; d++){
        _d1[d] = _d1star[_ci[d]-1];
      }
//      RprintVecAsMat("d1", _d1, 1, *K);
//      _d1[0] = 0; _d1[1] = 0.306; _d1[2] = 0; _d1[3] = 0; _d1[4] = 0.306; _d1[5] = 0.306;


      ////////////////////////////////////////////////////////////////////////////////////
      //
      // Update the stickbreaking weight
      //
      // Notice that for-loop starts on 1 for both h and d.  This is because
      // d1star[0] = 0 always and ci[0] = 1 always as d11 = 0 always.  Thus cluster
      // labels start at 2 and we only have H-1 free atoms not H.
      //
      ////////////////////////////////////////////////////////////////////////////////////
      tmp=1.0;
      for(h=1; h<*H; h++){
	    ccount = 0.0;
	    cumccount = 0.0;

		for(d=1; d<*K; d++){
		  if(_ci[d] == h+1){
		    ccount = ccount + 1.0;
		  }
		  if(_ci[d] > h+1){
		    cumccount = cumccount + 1.0;
		  }

		}
        astar = 1.0 + ccount;
        bstar = alpha + cumccount;
//		Rprintf("astar = %f\n", astar);
//		Rprintf("bstar = %f\n", bstar);

        // Note that I store V values starting in slot zero.
        // So V goes from slot 0 to H*-2
        _V[h-1] = rbeta(astar,bstar);

		if(h==(*H-1)) _V[h-1] = 1.0;
		if(h==1) _pi_vec[h-1] = _V[h-1];
		if(h>1){
		  tmp = tmp*(1-_V[h-2]);
		  _pi_vec[h-1] = _V[h-1]*tmp;
		}
	  }
//      RprintVecAsMat("V", _V, 1, *H-1);
//	    RprintVecAsMat("pi_vec", _pi_vec, 1, *H-1);



      // Update sh and omega
      if(*model==3){

        // update sh
        sum1 = 0, sumv0 = 0;
        for(h=1; h < *H; h++){
          pp[0] = dnorm(_d1star[h], 0, v0/3, 1) + log(_omega);
          pp[1] = dnlp(_d1star[h], 0, kap0, scale, nu, 1) + log(1- _omega);
//          RprintVecAsMat("pp", pp, 1, 2);
          maxp = pp[0]; if(pp[1] > pp[0]) maxp = pp[1];

          pp[0] = exp(pp[0] - maxp); pp[1] = exp(pp[1] - maxp);
          prob[0] = pp[0]/(pp[0] + pp[1]); prob[1] = pp[1]/(pp[0] + pp[1]);
//          Rprintf("prob[0] = %f\n", prob[0]);
          if(rbinom(1,prob[1]) == 0){
            _sh[h] = v0;
            sumv0 = sumv0 + 1;
          } else {
            _sh[h] = 1;
            sum1 = sum1 + 1;
          }

        }
//        Rprintf("sum1 = %f\n", sum1);
//        Rprintf("sumv0 = %f\n", sumv0);
//        RprintVecAsMat("_sh", _sh, 1, *H);

        // update omega
        _omega = rbeta(aw+sum1, bw+sumv0);
      }

    }




/*
    ////////////////////////////////////////////////////////////////////////////////////////////
    //
    // in sample prediction to assess model fit
    //
    ////////////////////////////////////////////////////////////////////////////////////////////
    if((i > (*nburn-1)) & (i % (*nthin) == 0)){

      for(j = 0; j < *N; j++){

        mn = _muh[_Si[j]-1];
        if(*model == 2){
          xb = 0.0;
          for(b = 0; b < ncov; b++){
      	    xb = xb + fullXmat[j*ncov+b]*_beta[b];
      	  }

          mn = _muh[_Si[j]-1] + xb;
        }

        _ispred[j] = rnorm(mn, sqrt(_sig2h[_Si[j]-1]));


        /////////////////////////////////////////////
        //
        // Compute the CPO and lpml using the mixture
        //
        /////////////////////////////////////////////

        _like[j] = dnorm(y[j], mn, sqrt(_sig2h[_Si[j]-1]), 0);

        // These are needed for WAIC
        mnlike[j] = mnlike[j] + (_like[j])/(double) nout;
        mnllike[j] = mnllike[j] + log(_like[j])/(double) nout;

        CPOinv[j] = CPOinv[j] + (1/(double) nout)*(1/_like[j]);

      }

    }
*/
    //////////////////////////////////////////////////////////////////////////////////////
    //
    // Store MCMC iterates
    //
    //////////////////////////////////////////////////////////////////////////////////////
    if((i > (*nburn-1)) & ((i+1) % *nthin ==0)){

      tau2[ii] = _tau2;
      omega[ii] = _omega;
      for(j=0; j<*N; j++){
        mu[ii + nout*j] = _mu[j];
      }

      for(j=0; j<(*N)*(npwcomp+1); j++){
        delta[ii + nout*j] =  _delta[j];
      }

//      RprintVecAsMat("sh", _sh, 1, *H);
//      RprintIVecAsMat("ci", _ci, 1, *K);
      for(d=0; d<*K; d++){
        d1[ii + nout*d] = _d1[d];
        ci[ii + nout*d] = _ci[d];
        sh[ii + nout*d] = _sh[_ci[d]-1];
      }

      ii = ii + 1;
    }

  }


  //////////////////////////////////////////////////////////////////////////////////
  // calculate LPML
  //////////////////////////////////////////////////////////////////////////////////
  _lpml=0.0;
  for(j = 0; j < *N; j++){
    _lpml = _lpml + log(1/CPOinv[j]);
  }
  lpml[0] = _lpml;

  ////////////////////////////////////////////////////////////////////////////////////////////
  // Computing WAIC  (see Gelman article in lit review folder)
  ////////////////////////////////////////////////////////////////////////////////////////////
//  elppdWAIC = 0.0;
//  for(j = 0; j < *N; j++){
//    elppdWAIC = elppdWAIC + (2*mnllike[j] - log(mnlike[j]));
//  }
//  waic[0] = -2*elppdWAIC;


}



SEXP NETWORK_MA(SEXP r_mat, SEXP n_mat, SEXP tid_mat, SEXP tid_logical,
                SEXP nobs, SEXP N, SEXP K, SEXP H,
                SEXP model, SEXP modelPriors, SEXP tau_prior, SEXP mh,
                SEXP verbose,
                SEXP niter, SEXP nburn, SEXP nthin) {
  int nprot = 0;

  int _nobs = asInteger(nobs);
  int _N = asInteger(N);
  int _K = asInteger(K);
  int _H = asInteger(H);
  int _model = asInteger(model);
  int _tau_prior = asInteger(tau_prior);
  int _verbose = asInteger(verbose);

  int _niter = asInteger(niter);
  int _nburn = asInteger(nburn);
  int _nthin = asInteger(nthin);

  double nout = (_niter-_nburn)/_nthin;
  int npwcomp = factorial(_K)/(factorial(2)*factorial(_K-2));

  r_mat = PROTECT(coerceVector(r_mat, REALSXP)); nprot++;
  n_mat =  PROTECT(coerceVector(n_mat, REALSXP)); nprot++;
  tid_mat =  PROTECT(coerceVector(tid_mat, INTSXP)); nprot++;
  tid_logical =  PROTECT(coerceVector(tid_logical, INTSXP)); nprot++;
  modelPriors = PROTECT(coerceVector(modelPriors, REALSXP)); nprot++;
  mh = PROTECT(coerceVector(mh, REALSXP)); nprot++;

  SEXP MU = PROTECT(allocMatrix(REALSXP, nout, _N)); nprot++;
  SEXP DELTA = PROTECT(allocMatrix(REALSXP, nout, (_N)*(npwcomp+1))); nprot++;
  SEXP TAU2 = PROTECT(allocMatrix(REALSXP, nout, 1)); nprot++;
  SEXP D1 = PROTECT(allocMatrix(REALSXP, nout, _K)); nprot++;
  SEXP CI = PROTECT(allocMatrix(INTSXP, nout, _K)); nprot++;
  SEXP OMEGA = PROTECT(allocMatrix(REALSXP, nout, 1)); nprot++;
  SEXP SH = PROTECT(allocMatrix(REALSXP, nout, _K)); nprot++;

  SEXP LPML = PROTECT(Rf_allocVector(REALSXP, 1)); nprot++;

  double *MUout, *DELTAout, *TAU2out, *D1out, *OMEGAout, *SHout, *LPMLout;
  int *CIout;

  MUout = REAL(MU);
  DELTAout = REAL(DELTA);
  TAU2out = REAL(TAU2);
  D1out = REAL(D1);
  CIout = INTEGER(CI);
  SHout = REAL(SH);
  OMEGAout = REAL(OMEGA);

  LPMLout = REAL(LPML);


  GetRNGstate();

  network_ma(REAL(r_mat), REAL(n_mat), INTEGER(tid_mat), INTEGER(tid_logical),
                &_nobs, &_N, &_K, &_H,
                &_model, REAL(modelPriors), &_tau_prior, REAL(mh), &_verbose,
                &_niter, &_nburn, &_nthin,
                MUout, DELTAout, TAU2out, D1out, CIout, OMEGAout, SHout, LPMLout);

  PutRNGstate();


  SEXP ans = PROTECT(allocVector(VECSXP, 8)); nprot++;
  SET_VECTOR_ELT(ans, 0, MU);
  SET_VECTOR_ELT(ans, 1, DELTA);
  SET_VECTOR_ELT(ans, 2, TAU2);
  SET_VECTOR_ELT(ans, 3, D1);
  SET_VECTOR_ELT(ans, 4, CI);
  SET_VECTOR_ELT(ans, 5, OMEGA);
  SET_VECTOR_ELT(ans, 6, SH);
  SET_VECTOR_ELT(ans, 7, LPML);


  SEXP nm = allocVector(STRSXP, 8);
  setAttrib(ans, R_NamesSymbol, nm);
  SET_STRING_ELT(nm, 0, mkChar("mu"));
  SET_STRING_ELT(nm, 1, mkChar("delta"));
  SET_STRING_ELT(nm, 2, mkChar("tau2"));
  SET_STRING_ELT(nm, 3, mkChar("d1"));
  SET_STRING_ELT(nm, 4, mkChar("ci"));
  SET_STRING_ELT(nm, 5, mkChar("omega"));
  SET_STRING_ELT(nm, 6, mkChar("sh"));
  SET_STRING_ELT(nm, 7, mkChar("lpml"));

  UNPROTECT(nprot);
  return(ans);
}



