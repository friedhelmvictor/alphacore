#' alphaCore ranking of nodes in a complex, directed network
#'
#' \code{alphaCore} returns a node ranking of a graph
#'
#' Iteratively computes a node ranking based on a feature set derived from
#' edge attributes and optionally static node features using the
#' mahalanobis data depth function at the origin.
#'
#' @param input_graph An igraph object of a directed graph
#' @param features A list of selected node features; default is ["all"]
#' @param stepSize defines the stepsize of each iteration as percentage of node count
#'   Higher values (e.g., > 0.1) speed up execution but yield lower ranking resolution, while lower values (e.g., < 0.1)
#'   provide finer ranking but increase runtime.
#' @param startEpsilon the epsilon to start with. Removes all nodes with depth>epsilon at start
#'   Higher values (e.g., > 1.0) remove more nodes early, emphasizing denser cores,
#'   whereas lower values (e.g., < 1.0) allow for a more gradual refinement of the ranking.
#' @param exponentialDecay dynamically reduces the step size, to have high cores with few nodes
#'   which results in higher resolution for the highest cores. Set to TRUE if node ranking precision is more important than speed.
#' @return A dataframe of node name and alpha value indicating the ranking
alphaCore <- function(input_graph,
                      features = "all",
                      stepSize = 0.1,
                      startEpsilon = 1,
                      exponentialDecay = TRUE) {
  
  #1 Extract numerical node features from the graph
  node_features <- extractFeatures(input_graph, features)
  #2 Compute covariance matrix
  cov_mat <- cov(node_features[, -c("node")])
  #Try computing the inverse; fallback to pseudo-inverse for stability
  tryCatch({
    cov_mat_inv <- solve(cov_mat)
  }, error = function(e) {
    message("Covariance matrix is not invertible; using pseudo-inverse for stability.")
    cov_mat_inv <- MASS::ginv(cov_mat)  
  })
  #3
  node_features$depth <- mhdOrigin(node_features[, !c("node")], cov_mat_inv)
  #4
  epsilon <- startEpsilon
  #5
  result <- node_features[, c("node")]
  result$alpha <- 0
  #6
  result$batch <- 0
  #7
  alpha <- 1 - epsilon # modified
  #8
  alpha_prev <- alpha # modified
  #9
  batch_ID <- 0
  #10
  while (vcount(input_graph) > 0) {
    #11
    repeat{
      #12: for each depth geq epsilon
      #13
      nodes <- node_features[depth >= epsilon]$node
      result[node %in% nodes, alpha := alpha_prev]
      #14
      result[node %in% nodes, batch := batch_ID]
      #15
      input_graph <- delete_vertices(input_graph, nodes)
      #16
      batch_ID <- batch_ID + 1
      #17
      node_features <- extractFeatures(input_graph, features)
      #18
      node_features$depth <- mhdOrigin(node_features[, !c("node")], cov_mat_inv)
      
      if(length(nodes) == 0){ # 19 : if there doesn't exist a vertex with depth geq epsilon
        break
      }
    }
    #20
    alpha_prev <- alpha
    #22 # reduce epsilon with strategy
    if(exponentialDecay) {
      localStepSize <- ceiling(vcount(input_graph) * stepSize)
      epsilon <- min(head(node_features[order(depth, decreasing = T)], localStepSize)$depth)
    } else {
      epsilon <- epsilon - stepSize
    }
    #21
    alpha <- 1 - epsilon
    
  }
  return(result)
  #
}


############################## AUX FUNCTIONS ################################
# Extract numerical node features from the graph
extractFeatures <- function(graph, features) {
  all_features <- names(vertex.attributes(graph))
  numeric_features <- all_features[sapply(all_features, function(f) is.numeric(V(graph)[[f]]))]
  
  if (features == "all") {
    if (length(numeric_features) == 0) {
      message("No numerical node features found. Reverting to default AlphaCore node features.")
      return(computeNodeFeatures(graph))
    }
    selected_features <- numeric_features
  } else {
    selected_features <- intersect(features, numeric_features)
    if (length(selected_features) == 0) {
      message("None of the selected features are numerical. Reverting to default AlphaCore node features.")
      return(computeNodeFeatures(graph))
    }
  }
  
  data <- as.data.table(setNames(lapply(selected_features, function(f) V(graph)[[f]]), selected_features))
  data[, node := V(graph)$name]
  return(data)
}


computeNodeFeaturFun <- function(graph) {
  # get node names
  nodes <- V(graph)$name
  # compute indegree
  indegree <- degree(graph, mode = "in")
  # compute sum of incoming weights (a.k.a. "strength")
  strength <- strength(graph, mode="in")
  # combine results
  nodeFeatures <- data.table(node=nodes, indegree=indegree, strength=strength)
  return(nodeFeatures)
}

computeNodeFeatures <- function(graph, features) {
  nodeFeatures <- data.table(node = V(graph)$name)
  if("degree" %in% features) {
    nodeFeatures[, indegree := degree(graph, mode = "all")]
  }
  if("indegree" %in% features) {
    nodeFeatures[, indegree := degree(graph, mode = "in")]
  }
  if("outdegree" %in% features) {
    nodeFeatures[, outdegree := degree(graph, mode = "out")]
  }
  if("strength" %in% features) {
    nodeFeatures[, strength := strength(graph, mode = "all")]
  }
  if("instrength" %in% features) {
    nodeFeatures[, instrength := strength(graph, mode = "in")]
  }
  if("outstrength" %in% features) {
    nodeFeatures[, outstrength := strength(graph, mode = "out")]
  }
  if("triangles" %in% features) {
    nodeFeatures[, triangles := count_triangles(graph)]
  }
  if("neighborhoodsize" %in% features) {
    nodeFeatures[, neighborhoodsize := neighborhood.size(graph, mode = "all", mindist = 1)]
  }
  if("inneighborhoodsize" %in% features) {
    nodeFeatures[, inneighborhoodsize := neighborhood.size(graph, mode = "in", mindist = 1)]
  }
  if("outneighborhoodsize" %in% features) {
    nodeFeatures[, outneighborhoodsize := neighborhood.size(graph, mode = "out", mindist = 1)]
  }
  return(nodeFeatures[])
}

customNodeFeatures <- function(features) {
  return(function(x) {
    return(computeNodeFeatures(x, features = features))
  })
}


mhdOrigin <- function(data, sigma_inv) {
  origin <- rep(0,ncol(data)) # c(0,0,...)
  # We reuse the Mahalanobis distance implementation of the stats package,
  # which returns the squared Mahalanobis distance: (x - μ)' Σ^-1 (x - μ) = D^2
  # To arrive at the Mahalanobis Depth to the origin, we only need to add 1 and
  # take the reciprocal.
  return((1 + pmax(stats::mahalanobis(data, center = origin, sigma_inv, inverted = TRUE), 0))^-1)
}
