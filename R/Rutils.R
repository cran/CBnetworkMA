
# Function to find the longest simple path
find_longest_simple_path <- function(graph) {
  longest_path <- NULL
  max_length <- 0

  for (node in igraph::V(graph)) {
    paths <- igraph::all_simple_paths(graph, from=node)
    max_path_length <- max(lengths(paths))

    if (max_path_length > max_length) {
      longest_path <- paths[[which.max(lengths(paths))]]
      max_length <- max_path_length
    }
  }

  return(longest_path)
}

# Function to find the NMA path within a clique
path_NMA <- function(v, graph) {
  # Extract a subnetwork with the given nodes
  subgraph <- igraph::induced_subgraph(graph, v)
  igraph::V(subgraph)$label <- igraph::V(graph)$label[v]

  # Find the longest simple path
  m <- igraph::get.adjacency(subgraph, sparse = FALSE)
  if (is.null(colnames(m))) colnames(m) <- rownames(m) <- igraph::V(subgraph)$label

  largest_path <- names(sort(rowSums(m), decreasing = TRUE))
  size_path <- length(largest_path)

  i <- 1
  aux <- sum(sapply(igraph::all_simple_paths(subgraph, from=which(colnames(m) %in% largest_path[i]),
                                             to=which(colnames(m) %in% largest_path[i+1])),
                    length) == 2) ==
    sum(sapply(igraph::all_simple_paths(subgraph, from=which(colnames(m) %in% largest_path[i+1]),
                                        to=which(colnames(m) %in% largest_path[i])),
               length) == 2)

  if (aux) final_path <- paste(colnames(m)[which(colnames(m) %in% largest_path[i])], "=",
                               colnames(m)[which(colnames(m) %in% largest_path[i+1])])
  if (!aux) final_path <-  paste(colnames(m)[which(colnames(m) %in% largest_path[i])], "<",
                                 colnames(m)[which(colnames(m) %in% largest_path[i+1])])

  if (size_path > 2) {
    for (i in 2:(size_path-1)) {
      aux <- sum(sapply(igraph::all_simple_paths(subgraph, from=which(colnames(m) %in% largest_path[i]),
                                                 to=which(colnames(m) %in% largest_path[i+1])),
                        length) == 2) ==
        sum(sapply(igraph::all_simple_paths(subgraph, from=which(colnames(m) %in% largest_path[i+1]),
                                            to=which(colnames(m) %in% largest_path[i])),
                   length) == 2)

      if (aux) final_path <- paste(final_path, "=", colnames(m)[which(colnames(m) %in% largest_path[i+1])])
      if (!aux) final_path <- paste(final_path, "<",  colnames(m)[which(colnames(m) %in% largest_path[i+1])])
    }
  }
  final_path
}
