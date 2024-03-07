# Wrapper that executes C code to sample from the posterior distribution
# from network meta analysis model
# mb=0; sb=1; md=0; sd=1; tau_prior = "lognormal"; tau_max=5; tau_lm = -2.34; tau_lsd = 1.62; alpha=1; aw=1; bw=1; v0=0.1; kap=1; scale=1; nu=1; mh=c(0.5, 0.5, 0.1, 0.5); H=20; verbose=FALSE
networkMA <- function(data, model="gaussian",
                        niter=1100, nburn=100, nthin=1,
                        mb=0, sb=1, md=0, sd=1,
                        tau_prior = "uniform", tau_max=5, tau_lm = -2.34, tau_lsd = 1.62,
                        alpha=1, aw=1, bw=1, v0=0.1, scale=1, nu=1,
                        mh=c(0.5, 0.5, 0.1, 0.5), H=20, verbose=FALSE){


  # data must have the following columns .
  # sid - study id (must an integer beginning with 1)
  # tid - treatment id (must be an integer beginning with 1)
  # r - number of "success" associated with treatment x study
  # n - number of "trials" associated with treatment x study

  N <- max(data$sid) # number of studies
  K <- length(unique(data$tid)) # total number of treatments across studies
  nobs <- length(data$sid) # total number of observations


  # Create matrices that indicate comparisons in each study and the number of
  # trials and successes.
  tid_mat <- r_mat <- n_mat <- tid_logical <- matrix(FALSE, nrow=N, ncol=choose(K,2)+1)
  colnames(tid_mat) <- colnames(r_mat) <- colnames(n_mat) <- colnames(tid_logical) <-
           c("baseline", apply(combn(1:K,2), 2, function(x) paste(x, collapse="_")))
  tmp1 <- tapply(data$tid, data$sid, function(x) c(paste(x[1],x[1],sep="_"), paste(x[1], x[-1],sep="_")))
  tmp2 <- tapply(data$r, data$sid, function(x) x)
  tmp3 <- tapply(data$n, data$sid, function(x) x)
  tmp4 <- tapply(data$tid, data$sid, function(x) x)
  for(i in 1:N){
    tid_logical[i,tmp1[[i]][-1]] <- TRUE
    r_mat[i, c("baseline", tmp1[[i]][-1])] <- tmp2[[i]]
    n_mat[i, c("baseline", tmp1[[i]][-1])] <- tmp3[[i]]
    tid_mat[i, c("baseline", tmp1[[i]][-1])] <- tmp4[[i]]
  }
  tid_logical[,1] <- 1



  ## non-local inverse moment prior
  dnlp <- function(x,x0, kap,scale,nu,logd=FALSE){
    ld <- log(kap) + (nu/2)*log(scale) - lgamma(nu/(2*kap)) +
          -0.5*(nu+1)*log((x - x0)^2) - ((x-x0)^2/scale)^(-kap)
    if(x == x0){ld <- -Inf}
    if(logd){out <- ld}
    if(!logd){out <- exp(ld)}
    out
  }


  # Calibrate the spike and slab nlp prior
  xx <- seq(-9*v0, 0, length=100001)
  xx0 <- xx[which.min(abs(dnorm(xx, 0, v0/3) - 0.01))]

  kap0 <- seq(0.05,5, length=10001)
  tmp <- which.min(abs(dnlp(xx0, 0, kap=kap0, scale=1,nu=1) - dnorm(xx0, 0, v0/3)))
  kap0 <- kap0[tmp]

  # priors
  # mb - mean of mu_{i, b_i}
  # sb - standard deviation of mu_{i, b_i}
  # md -  mean of d*_j
  # sd - standard deviation of d*_j
  # tau_max - upper-bound of tau
  # alpha - precision/scale parameter of DP
  # aw - shape 1 parameter for omega
  # bw - shape 2 parameter for omega
  # v0 - constant to which sd is multiplied in the spike.  Determines the similarity of comparisons necessary
  #      to conlcude that two treatments are equal
  # H - truncation of the infinite mixture prior
  modelPriors <- c(mb, sb, tau_max, tau_lm, tau_lsd, md, sd, alpha, kap0, scale, nu, v0, aw, bw)
  if(tau_prior == "uniform")tauprior <- 1
  if(tau_prior == "lognormal") tauprior <- 2

  if(model=="gaussian"){
    run <- .Call("NETWORK_MA",
                  as.double(t(r_mat)), as.double(t(n_mat)), as.integer(t(tid_mat)), as.integer(t(tid_logical)),
                  as.integer(nobs), as.integer(N), as.integer(K), as.integer(H),
                  as.integer(1), as.double(modelPriors), as.integer(tauprior), as.double(mh),
                  as.integer(verbose),
                  as.integer(niter), as.integer(nburn), as.integer(nthin))
    if(tau_prior == "uniform"){
      priors <- round(c(mb, sb, tau_max, md, sd),2)
      names(priors) <- c("mb", "sb","tau_max","md","sd")
    }
    if(tau_prior == "lognormal"){
      priors <- round(c(mb, sb, tau_lm, tau_lsd, md, sd),2)
      names(priors) <- c("mb", "sb","tau_lm","tau_lsd","md","sd")
    }
    run <- run[c(1,2,3,4,8)]

  }

  if(model=="dp_gaussian"){
    run <- .Call("NETWORK_MA",
                  as.double(t(r_mat)), as.double(t(n_mat)), as.integer(t(tid_mat)), as.integer(t(tid_logical)),
                  as.integer(nobs), as.integer(N), as.integer(K), as.integer(H),
                  as.integer(2), as.double(modelPriors), as.integer(tauprior), as.double(mh),
                  as.integer(verbose),
                  as.integer(niter), as.integer(nburn), as.integer(nthin))
    if(tau_prior == "uniform"){
      priors <- round(c(mb, sb, tau_max, md, sd, alpha),2)
      names(priors) <- c("mb", "sb","tau_max","md","sd", "alpha")
    }
    if(tau_prior == "lognormal"){
      priors <- round(c(mb, sb, tau_lm, tau_lsd, md, sd, alpha),2)
      names(priors) <- c("mb", "sb","tau_lm","tau_lsd","md","sd", "alpha")
    }
    run <- run[c(1,2,3,4,5,8)]
  }

  if(model=="dp_spike_slab"){
    # cat("kap0 = ", kap0, "\n")
    run <- .Call("NETWORK_MA",
                  as.double(t(r_mat)), as.double(t(n_mat)), as.integer(t(tid_mat)), as.integer(t(tid_logical)),
                  as.integer(nobs), as.integer(N), as.integer(K), as.integer(H),
                  as.integer(3), as.double(modelPriors), as.integer(tauprior), as.double(mh),
                  as.integer(verbose),
                  as.integer(niter), as.integer(nburn), as.integer(nthin))

    if(tau_prior == "uniform"){
      priors <- round(c(mb, sb, tau_max, alpha, v0, aw, bw, kap0, scale, nu),2)
      names(priors) <- c("mb", "sb","tau_max","alpha", "v0","aw","bw","kap","scale","nu")
    }
    if(tau_prior == "lognormal"){
      priors <- round(c(mb, sb, tau_lm, tau_lsd, alpha, v0, aw, bw, kap0, scale, nu),2)
      names(priors) <- c("mb", "sb","tau_lm","tau_lsd", "alpha", "v0","aw","bw","kap","scale","nu")
    }

  }


  # Creat the order pairwise comparison matrix
  nout <- (niter - nburn)/nthin
  ordmat <- list()
  for(t in 1:nout){
    # Here I create all the pairwise differences
  	if(model=="gaussian"){
  	  if(K > 2) dtilde  <- c(run$d1[t,-1], apply(combn(run$d1[t,-1],2), 2, diff))
  	}
    if(model=="dp_gaussian"){
  	  # Here I create all the pairwise differences
  	  if(K > 2) dtilde  <- c(run$d1[t,-1], apply(combn(run$d1[t,-1],2), 2, diff))
    }
    if(model=="dp_spike_slab"){
      # Here I create all the pairwise differences and adjust for those that are allocatd
  	  # to slab and those that are allocated to spike
  	  if(K > 2){
  	    dtilde  <- c(run$d1[t,-1]*(run$sh[t,-1]==1),
  	                          apply(combn(run$d1[t,-1],2)*(combn(run$sh[t,-1],2)==1), 2, diff))
  	  }
    }
  	# Here I create a pairwise comparison matrix with upper triangle containing zero if treatments
  	# are equal
    if(K==2) dtilde <- run$d1[t,2]
  	tmp <- matrix(0, nrow=K, ncol=K)
  	tmp[lower.tri(tmp)] <- dtilde
  	tmp[tmp < 0] <- -1
  	tmp[tmp > 0] <- 1
  	ordmat[[t]] <- t(tmp)

  }

  # list of possible treatments

  run$ordmat <- ordmat
  run$prior_values <- priors
  run[-(length(run)-2)]

}


if(FALSE){
# This code will be used to run function to debug etc.
  library(TeachingDemos)
  library(NetworkMA)

  coverage <- matrix(NA, nrow=100, ncol=4)
  colnames(coverage) <- c("mu","delta","d1","tau2")

  source("~/Research/BYU/NetworkMetaAnalysis/analysis/dataGeneration.R")
  load("~/Research/BYU/NetworkMetaAnalysis/analysis/ReducedCipriani_dataAnalysis.RData")
  mu_true <- apply(m3$mu,2,mean)
  tau_true <- 0.1
  for(ii in 1:100){
    cat("dataset number = ", ii, "\n")
    set.seed(ii)
    synth.data <- dat.gen.cipriani(ci=c(1,2,2),
                                 mu_true=mu_true,
                                 tmt_keep=c(2,4,5,9,11,12),
                                 d1_true=(c(0,0.306,0,0,0.306,0.306)*mult),
                                 tau2_true=tau_true^2,
                                 tau_max=0.1, mb=0, sb=1, md=0, sd=1)


    out <- networkMA(synth.data$data, model="gaussian", niter=30000, nburn=20000, nthin=10, mb=0, sb=10, md=0, sd=5,
                        tau_prior = "lognormal", tau_max=5, tau_lm = -2.34, tau_lsd = 1.62,
                        alpha=1, aw=1, bw=1, v0=0.1, kap=1, scale=1, nu=1)

    # check mus
    coverage[ii,1] <- mean(apply(t(apply(out$mu, 2, emp.hpd)) - synth.data$mu, 1, prod) < 0)

    # check d1
    coverage[ii,3] <- mean(apply(t(apply(out$d1, 2, emp.hpd))[-1,] - synth.data$d1[-1], 1, prod) < 0)

    # check tau2
    coverage[ii,4] <- prod(apply(out$tau2,2,emp.hpd) - tau_true^2) < 0

    # check deltas
    tmp <- apply(out$delta,2,emp.hpd)
    tmpkeep <- apply(t(tmp),1,diff) != 0
    coverage[ii,2] <- mean(apply(t(tmp)[tmpkeep,] - c(t(synth.data$delta))[tmpkeep],1,prod) < 0)
  }
}
