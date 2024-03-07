\name{networkMA}
\alias{networkMA}
\title{Contrast-Based Bayesian Network Meta-Analysis Model}
\description{
Fits a contrast-based Bayesian network meta-anlaysis model to data that contains a binary reponse.
}
\usage{
networkMA(data,
          model="gaussian",
          niter=1100, nburn=100, nthin=1,
          mb=0, sb=1,
          md=0, sd=1,
          tau_prior = "uniform", tau_max=5,
          tau_lm = -2.34, tau_lsd = 1.62,
          alpha=1,
          aw=1, bw=1,
          v0=0.1, scale=1, nu=1,
          mh=c(0.5, 0.5, 0.1, 0.5),
          H=20,
          verbose=FALSE)

}

\arguments{
\item{data}{Data frame that must consist of the following four columns
  \itemize{
    \item  sid - study id (must be contiguous integers beginning with 1).
    \item tid - treatment id (must be contiguous integers beginning with 1).
    \item r - number of "successes" recorded for each treatment by study combination.
    \item n - number of "trials" recorded for each treatment by study combination.
  }
Note that the reference treatment must be labeled using 1 and treatment labels must be contiguous integers.  The colnames of the data fram must also be ``sid'',``tid'',``r'', and ``n''.

}
\item{model}{Specifies the model that will be fit.  The options are
  \itemize{
\item gaussian - treatment effects are modeled with a Gaussian.
    \item dp_gaussian - treatment effects are modeled with a DPM with a Gaussian base measure.
    \item dp_spike_slab - treatment effects are modeled with a DPM with a spike & slab base measure.
  }
}
\item{niter}{number of MCMC iterates to be collected. default is 1100.}

\item{nburn}{number of MCMC iterates discared as burn-in. default is 100.}

\item{nthin}{number by which the MCMC chain is thinned. default is 1.  Thin must be selected so that it is a multilple of (draws - thin).}

\item{mb}{prior mean for baseline treatment.}

\item{sb}{prior standard deviation for baseline treatment.}

\item{md}{prior mean for d1k.  Only used for gaussian and dp_gaussian models.}

\item{sd}{prior standard deviation for d1k.  Only used for gaussian and dp_gaussian models.}

\item{tau_prior}{prior distribution for \eqn{\tau}.  Options are
  \itemize{
    \item "uniform" which implies that \eqn{\tau \sim UN(0, tau\_max)}.
    \item "lognormal" which implies that \eqn{\tau \sim LogNormal(tau\_lm, tau\_lsd)}.
 }
}
\item{tau_max}{Upper bound on \eqn{\tau}.  Only used when tau_prior = ``uniform''.  Default is 5.}

\item{tau_lm}{mean of \eqn{log(\tau)}.  Only used if tau_prior = ``lognormal''.  Default is -2.34.}

\item{tau_lsd}{standard deviation of \eqn{log(\tau)}.  Only used if tau_prior = lognormal.  Default is 1.62.}

\item{alpha}{Precision parameter of the DP. Only used if model is dp_gaussian or dp_spike_slab. Default value is 1.}

\item{aw}{first shape parameter of omega's beta prior.  only used if model is dp_spike_slab.  Default value is 1.}

\item{bw}{second shape parameter of omega's beta prior.  only used if model is dp_spike_slab.  Default value is 1.}

\item{v0}{Parameter that. Default is 0.1.}

\item{scale}{Parameter that. Default is 1.}

\item{nu}{Parameter that. Default is 2.}

\item{mh}{4-dimensional vector containing values for the standard deviation of the Gaussian proposal associated with the MH update for \eqn{\mu}, \eqn{\delta}, \eqn{\tau}, \eqn{d_{1k}}.}

\item{H}{Truncated upper bound of the stick-breaking representation of the DP.  Only used for the dp_gaussian or dp_spike_slab models.}

\item{verbose}{Logical indicating if information regarding data and MCMC iterate should be printed to screen.}

}

\details{
This function permits the user to fit three different types of binomial models for contrast-based Bayesian network meta-analysis.
}

\value{
The function returns a list containing arrays filled with MCMC iterates that correspond to model parameters. In order to provide more detail, in what follows let

"T" - be the number of MCMC iterates collected,

"N" - be the number of studies,

"K" - be the total number of treatments.

The output list contains the following
\itemize{
\item mu - a matrix  of dimension (T, N) containing MCMC iterates associated with each study's baseline treatment (\eqn{mu_{i, b_i}}).

\item delta - a matrix of dimension (T, N*K(K-1)/2+1) containing MCMC iterates associated \eqn{\delta_{i, b_{i}k}}.  For ease of storage, each study is alloted enough columns to contain all pairwise comparisons (i.e., deltas).  If a comparison is not included in a study the column associated with it is set to -99.  By default the first column for each study is a vector of zeros (corresponding to \eqn{d_{11}})

\item tau2 -  a matrix of dimension (T,1) containing MCMC iterates associated with the variance of delta (common across studies)

\item d1 - a matrix of dimension (T, K) containing MCMC iterates assocated with \eqn{d_{1k}}.  Note that \eqn{d_{11} = 0} so the first column contains zeros.

\item ci - a matrix of dimension (T, K) containing MCMC iterates associated with cluster labels of the K treatments.  The first column corresponds to \eqn{d_{11}} which is always allocated to ``cluster 1''.  The last K-1 columns are cluster labels for the remaining \eqn{d_{1k}} treatments.  This object is provided only if model is "dp_gaussian" or "dp_spike_slab".

\item omega - a matrix of dimension (T,1) containing MCMC iterates for omega0 the weight associated with spike-and-slabe mixture. This object is provided only if model is "dp_spike_slab".

\item sh - a matrix of dimension (T, K) containing MCMC iterates for binary indicator of being allocated to spike (labeled using ``0.1'') or slab (labeled using ``1''). The first column corresponds to \eqn{d_{11}} which is always allocated to spike.  The last K-1 columns correspond to the remaining \eqn{d_{1k}} treatments.  This is object is provided only if model is "dp_spike_slab".

\item ordmat - a list of size T with each entry being a KxK pairwise treatment comparison matrix.

\item prior_values - a vector returning the priors values used.

}
}
\examples{



# This number of MCMC samples is for illustrative purposes only, it may
# be necessary to increase the total
ni <- 10000
nb <- 5000
nt <- 5

dat <- smoking # Use the smoking cessation dataset.

# total number of treatments
K <- length(unique(dat$tid))

# Fit model 1
set.seed(101)
# Fit the Guassian Effects model.
m1 <- networkMA(dat, model="gaussian", niter=ni, nburn=nb, nthin=nt,
                mb=0, sb=10, md=0, sd=1,
                tau_prior = "lognormal", tau_lm = -2.34, tau_lsd = 2,
                mh=c(0.5, 0.5, 0.05, 0.5))

mean(m1$d1[,2])
quantile(m1$d1[,2], c(0.025, 0.975))

# Fit the DP Gaussian base measure model.
m2 <- networkMA(dat, model="dp_gaussian", niter=ni, nburn=nb, nthin=nt,
                mb=0, sb=10, md=0, sd=1,
                tau_prior = "lognormal", tau_lm = -2.34, tau_lsd = 2,
                alpha=1,
                mh=c(0.5, 0.5, 0.05, 0.5))

mean(m2$d1[,2])
quantile(m2$d1[,2], c(0.025, 0.975))


# Fit the DP spike and slab base measure model.
m3 <- networkMA(dat, model="dp_spike_slab", niter=ni, nburn=nb, nthin=nt,
                mb=0, sb=10, md=0, sd=1,
                tau_prior = "lognormal", tau_lm = -2.34, tau_lsd = 2,
                alpha=1, aw=1, bw=1, v0=0.1, scale=1, nu=1,
                mh=c(0.5, 0.5, 0.05, 0.5))

mean(m3$d1[,2])
quantile(m3$d1[,2], c(0.025, 0.975))

# Function that finds the graph corresponding to the posterior samples, and
# graphs for a sequence of threshold probabilities (denoted as gamma in
# the article)

gamma_vec <- c(0.5, 0.75, 0.9, 0.95, 0.99)
networks <- network_graphs(m3[["ordmat"]], gamma=gamma_vec)


# One way of plotting the directed graph based on the output of the function
# above is the following.  The "igraph" package can be used to facilitate
# producing pair-wise graphical model display


oldpar <- par(no.readonly = TRUE)


# Plot network that corresponds to posterior mode
Network = networks[[1]]
out = cbind(from=1:ncol(Network),to=1:ncol(Network),color=0)
for(i in 1:(ncol(Network)-1)){
  for(j in (i+1):ncol(Network)){
    if(Network[i,j]==1) out = rbind(out,c(i,j,2))
    if(Network[i,j]==-1)out = rbind(out,c(j,i,2))
    if(Network[i,j]==0) out = rbind(out,c(i,j,1),c(j,i,1))
  }
}


mynet <- igraph::graph_from_data_frame(out,directed = TRUE)
igraph::V(mynet)$label.cex <- 0.5
igraph::V(mynet)$label.cex <- 1
names <- igraph::V(mynet)$name

igraph::E(mynet)$lty <- igraph::E(mynet)$color
igraph::E(mynet)$arrow.mode <- ifelse(out[,"color"]==0,0,2)
igraph::E(mynet)$color <- 'black'

plot(mynet,margin=c(0.05, 0.05, 0.05, 0.05),
     vertex.shape = "circle",
     vertex.color = "white",
     vertex.label.family = "Helvetica",
     edge.arrow.size = 0.3,
     vertex.label.color = "black",
     vertex.frame.color = "black",
     layout = igraph::layout_with_kk,
     asp = 0, ylim=c(-0.9,0.9), xlim=c(-0.9,0.9),
     main = paste("P[mode Graph|Data] =",networks[[2]]),
     sub = paste("Number of edges = ",nrow(out)-ncol(Network)))


# Or alternatively
coords <- igraph::layout_as_star(mynet, center = igraph::V(mynet)[1])

plot(mynet,margin=c(0.05, 0.05, 0.05, 0.05),
     vertex.shape = "circle",
     vertex.color = "white",
     layout = coords,
     vertex.label.family = "Helvetica",
     edge.arrow.size = 0.3,
     vertex.label.color = "black",
     vertex.frame.color = "black",
     layout = igraph::layout_with_kk,
     asp = 0, ylim=c(-0.9,0.9), xlim=c(-0.9,0.9),
     main = paste("P[mode Graph|Data] =",networks[[2]]),
     sub = paste("Number of edges = ",nrow(out)-ncol(Network)))




# Plot the sequence of graphs based on gamma
network_seq <- networks[[3]]

for(i in 1:length(network_seq)){
  Probpair <- gamma_vec[i]
  Network <- network_seq[[i]]
  # Plot network
  out = cbind(from=1:ncol(Network),to=1:ncol(Network),color=0)
  for(i in 1:(ncol(Network)-1)){
    for(j in (i+1):ncol(Network)){
      if(Network[i,j]==1) out = rbind(out,c(i,j,2))
      if(Network[i,j]==-1)out = rbind(out,c(j,i,2))
      if(Network[i,j]==0) out = rbind(out,c(i,j,1),c(j,i,1))
    }
  }
  # out

  # Compute joint probability
  PointEst = (Network + 10)*(upper.tri(Network) & Network!=-1111)
  prob = mean(sapply(m3[["ordmat"]],
         function(aux){
           sum(abs((aux+10)*(upper.tri(Network)&Network!=-1111)-PointEst))}
           ==0))


  mynet <- igraph::graph_from_data_frame(out,directed = TRUE)
  igraph::V(mynet)$label.cex <- 0.5
  igraph::V(mynet)$label.cex <- 1
  names <- igraph::V(mynet)$name

  igraph::E(mynet)$lty <- igraph::E(mynet)$color
  igraph::E(mynet)$arrow.mode <- ifelse(out[,"color"]==0,0,2)
  igraph::E(mynet)$color <- 'black'


  plot(mynet,margin=c(0.05, 0.05, 0.05, 0.05),
         vertex.shape = "circle",
         vertex.color = "white",
         layout = coords,
         vertex.label.family = "Helvetica",
         edge.arrow.size = 0.3,
         vertex.label.color = "black",
         vertex.frame.color = "black",
         layout = igraph::layout_with_kk,
         asp = 0, ylim=c(-0.9,0.9), xlim=c(-0.9,0.9),
         main = paste("max P[di - dj|Data] >=",Probpair,"and P[Graph|Data] =",prob),
         sub = paste("Number of edges = ",nrow(out)-ncol(Network)))
}




# Extract cliques

ordmat <- m3[["ordmat"]]

clique_extract(ordmat,
               type="Highest_Post_Prob")


clique_extract(ordmat,
               type="Highest_Pairwise_Post_Prob",
               gamma=0.9)

clique_extract(ordmat,
               type="Highest_Pairwise_Post_Prob",
               clique_size=5,
               gamma=0.95)


par(oldpar)







}
