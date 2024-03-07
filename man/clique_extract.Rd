\name{clique_extract}
\alias{clique_extract}
\title{Extract cliques from network}
\description{
Extracts cliques of a certain size from treatment network estimated based on a model fit using the `networkMA' function.
}
\usage{
clique_extract(ordmat,
              type = "Highest_Post_Prob",
              clique_size = NULL,
              gamma = 0.95,
              plot_graph = FALSE)

}

\arguments{
\item{ordmat}{The 'ordmat' object obtained from model fit using 'networMA' function.}
\item{type}{The type of treatment network summary from which cliques are to be extracted.  Options include
  \itemize{
    \item Highest_Post_Prob - Cliques are extracted from treatment network with highest posterior probability.
    \item Highest_Pairwise_Post_Prob - Cliques are extracted from treatment network based on lower bound of pair-wise posterior probability of treatment effects.
  }
}
\item{clique_size}{size of cliques that are to be returned.}
\item{gamma}{lower-bound for pairwise treatment comparison probabilities.  Only used if type = "Highest_Pairwise_Post_Prob"}
\item{plot_graph}{logical - if true a plot of the network is supplied}
}

\details{
This function takes as its input the MCMC iterates of the comparison matrix that is provided as an output to the 'networkMA' function and returns cliques from the the treatment network based on a criteria associated with the posterior distribution of treatment graphs.
}

\value{
A list containing a cliques of a particular size.
}
\examples{
# See the example supplied in the help file for 'networkMA'
}



