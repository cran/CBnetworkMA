clique_extract = function(ordmat,
                              type = "Highest_Post_Prob",
                              clique_size = NULL,
                              gamma = 0.95,
                              plot_graph = FALSE){
  K <- nrow(ordmat[[1]])
  if(type == "Highest_Post_Prob"){
    # Find network with highest prob
    ordmatcollapse = sapply(ordmat, function(aux)paste(aux, collapse = "|"))
    Postordmat = sort(table(ordmatcollapse),decreasing = TRUE)
    prob = Postordmat[1]/sum(Postordmat)

      # Plot network
    Network = ordmat[[which(names(prob)==ordmatcollapse)[1]]]
    out = cbind(from=1:ncol(Network),to=1:ncol(Network),color=0)
    for(i in 1:(ncol(Network)-1))
    {
      for(j in (i+1):ncol(Network)){
        if(Network[i,j]==1)
          out = rbind(out,c(i,j,2))
        if(Network[i,j]==-1)
          out = rbind(out,c(j,i,2))
        if(Network[i,j]==0)
        {
          out = rbind(out,c(i,j,1),c(j,i,1))
        }
      }
    }

    # out <- out[,-3]
    # out <- out[!(out[,1] == out[,2]),]

    # Blue, one direction, orange equal (both direction)
    mynet <- igraph::graph_from_data_frame(out,directed = TRUE)
    graph <- mynet
    igraph::V(graph)$label <- 1:K

    # Plot network
    if(plot_graph){
      plot(graph, layout=igraph::layout.circle)
    }
    # Find the NMA path for each large clique
    if(is.null(clique_size)){
      v <- suppressWarnings(igraph::largest.cliques(graph))
    } else {
      # Find the NMA path for smaller clique
      v <- suppressWarnings(igraph::cliques(graph, min = 2, max = clique_size))

    }
  }
  if(type == "Highest_Pairwise_Post_Prob"){
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

    Probpair <- gamma ########## THRESHOLD
    # Remove from Network the edges with pairwise probability less than Probpair
    Network = ifelse(pairProbs <= Probpair, -1111, Network0)
    Network = ifelse(lower.tri(Network, diag = TRUE), 0, Network)

    # Plot network
    out = cbind(from=1:ncol(Network),to=1:ncol(Network),color=0)
    for(i in 1:(ncol(Network)-1))
    {
      for(j in (i+1):ncol(Network)){
        if(Network[i,j]==1)
          out = rbind(out,c(i,j,2))
        if(Network[i,j]==-1)
          out = rbind(out,c(j,i,2))
        if(Network[i,j]==0)
        {
          out = rbind(out,c(i,j,1),c(j,i,1))
        }

      }
    }
    out



    # out <- out[,-3]
    # out <- out[!(out[,1] == out[,2]),]


    mynet <- igraph::graph_from_data_frame(out, directed = TRUE)
    graph <- mynet
    igraph::V(graph)$label <- 1:K
    if(plot_graph){
      plot(graph, layout=igraph::layout.circle)
    }


    # Find the NMA path for each large clique
    if(is.null(clique_size)){
      v <- suppressWarnings(igraph::largest.cliques(graph))
    } else {
      # Find the NMA path for smaller clique
      v <- suppressWarnings(igraph::cliques(graph, min = 2, max = clique_size))
    }

  }
  cl_list <- rep(0, length(v))
  count <- 1
  for(i in 1:length(v)){
    for(j in 1:length(v)){
      if (j != i){
        if(all(v[[i]] %in% v[[j]])){
          cl_list[i] <- 1
          break
        }
      }
    }
  }
  suppressWarnings(lapply(v[cl_list==0], path_NMA, graph = graph))
}
