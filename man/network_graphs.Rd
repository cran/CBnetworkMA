\name{network_graphs}
\alias{network_graphs}
\title{Plot treatment network}
\description{
Extracts graphs based on a model fit using the `networkMA' function.  The output of this function facilitates producing the network graphical display of treatment pair-wise comparisons.
}
\usage{
network_graphs(ordmat,
               gamma=c(0.5,0.75,0.9,0.95,0.99))

}

\arguments{
\item{ordmat}{The 'ordmat' object obtained from model fit using 'networMA' function.}
\item{gamma}{Sequence of probabilities that represent lower-bound of pairwise comparison probabilities that are included to produce a network graph.}
}

\details{
This function takes as its input the MCMC iterates of the comparison matrix that is provided as an output to the 'networkMA' function and returns the sequence of pairwise comparison graphs associated with the lower bound pobabilities supplied in gamma.
}

\value{
A list with the following entries
\itemize{
\item pairwise comparison matrix that corresponds to the posterior mode.  Note that pairwise matrices have the following
  \itemize{
    \item  1 - implies effect of column treatment is greater than the row treatment.
    \item -1 - implies  effect of row treatment is greater than the column treatment.
    \item 0 - implies column treatment and row treatment have same effect.
    \item -1111 - implies that there is not enough information to make conclusions regarding two treatments.
  }

\item posterior probability associated with the posterior mode pairwise comparison matrix

\item a list, the length of which corresponds to the length of the gamma vector.  Each entry corresponds to a pairwise comparison matrix associated with a particular gamma value.

}
}
\examples{
# See the example supplied in the help file for 'networkMA'
}



