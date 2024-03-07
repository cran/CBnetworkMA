# ordmat: MCMC output
network_graphs = function(ordmat, gamma=c(0.5, 0.75, 0.9, 0.95, 0.99)){

  # gamma is a vector of probability thresholds so that all pairwise edges
  # in a directed graph have posterior probability of at least gamma
  #################################
  # Network with highest post prob
  ##################################


  # Find network with highest prob
  ordmatcollapse = sapply(ordmat, function(aux)paste(aux, collapse = "|"))
  Postordmat = sort(table(ordmatcollapse),decreasing = T)
  prob = Postordmat[1]/sum(Postordmat)
  Network = ordmat[[which(names(prob)==ordmatcollapse)[1]]]

  mode_graph <- Network
  post_prob_mode_graph <- as.vector(prob)

  ##################################
  # Network with pairwise post prob above threshold (Probpair)
  ##################################
  # Find network with highest pairwise post prob
  Network = ordmat[[1]]-ordmat[[1]]
  pairProbs = ordmat[[1]]-ordmat[[1]]
  for(i in 1:(ncol(Network)-1))
  {
    for(j in (i+1):ncol(Network))
    {
      prob = (table(c(sapply(ordmat,function(x) x[i,j]),c(-1,0,1)))-1)/length(ordmat)
      Network[i,j] = as.numeric(names(which.max(prob)))
      pairProbs[i,j] = max(prob)
      # if(max(prob)<=Probpair)
      #   Network[i,j] = -1111
    }
  }

  # Find Network0 which is the coherent graph in the mcmc closest to Network
  weightl1 = sapply(1:length(ordmat), function(j1) sum(abs(ordmat[[j1]]-Network)*pairProbs))
  Network0 = ordmat[[sample(which(weightl1==min(weightl1)),1)]]



  out <- list()
  cnt <- 1
  for(Probpair in gamma){
    # Remove from Network the edges with pairwise probability less than Probpair
    Network = ifelse(pairProbs <= Probpair, -1111, Network0)
    Network = ifelse(lower.tri(Network, diag = T), 0, Network)
    out[[cnt]] <- Network
    cnt <- cnt+1
  }
  names(out) <- gamma
  list(mode_graph= mode_graph, post_prob_mode_graph=post_prob_mode_graph,out)

}

