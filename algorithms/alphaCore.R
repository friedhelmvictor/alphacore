#' alphaCore ranking of nodes in a complex, directed network
#'
#' \code{alphaCore} returns a node ranking of a graph
#'
#' Iteratively computes a node ranking based on a feature set derived from
#' edge attributes and optionally static node features using the
#' mahalanobis data depth function at the origin.
#'
#' @param input_graph An igraph object of a directed graph
#' @param featureComputeFun A function that converts a node's edges (with
#'   attributes) into node features. Default computes in-degree and strength
#' @param stepSize defines the stepsize of each iteration as percentage of node count
#' @param startEpsilon the epsilon to start with. Removes all nodes with depth>epsilon at start
#' @param exponentialDecay dynamically reduces the step size, to have high cores with few nodes
#' @return A dataframe of node name and alpha value indicating the ranking
alphaCore <- function(input_graph,
                      featureComputeFun = computeNodeFeaturFun,
                      stepSize = 0.1,
                      startEpsilon = 1,
                      exponentialDecay = TRUE) {
  
  #1
  node_features = featureComputeFun(input_graph)
  #2
  cov_mat_inv <- solve(cov(node_features[, -c("node")]), tol = NULL)
  # temporary ----
  cov_mat <- cov(node_features[, -c("node")])
  # temporary ----
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
      node_features = featureComputeFun(input_graph)
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
  return((1 + stats::mahalanobis(data, center = origin, sigma_inv, inverted = TRUE))^-1)
}
