#' alphaCore ranking of nodes in a complex, directed network
#'
#' \code{alphaCore} returns a node ranking of a graph
#'
#' Iteratively computes a node ranking based on a featureset derived from
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
                      stepSize = 1,
                      startEpsilon = 1,
                      exponentialDecay = FALSE) {
  
  #edge_attr(input_graph, "weight") <- normalizeData(edge_attr(input_graph, "weight"))
  
  epsilon <- startEpsilon
  previousEpsilon <- startEpsilon
  result <- list()
  localStepSize <- ceiling(vcount(input_graph) * stepSize)
  
  # internal function that reduces the epsilon (either dynamically or fixed stepSize)
  reduceEpsilon <- function() {
    previousEpsilon <<- epsilon
    if(!exponentialDecay) {
      return(epsilon - stepSize)
    }
    
    if(stepSize < 1) {
      if(exponentialDecay) {
        localStepSize <- ceiling(vcount(input_graph) * stepSize)
      }
      newEpsilon <- min(head(node_features[order(depth, decreasing = T)], localStepSize)$depth)
      if(newEpsilon < epsilon) {
        epsilon <- newEpsilon
      }
    } else {
      stop("Step size should be < 1")
    }
    return(epsilon)
  }
  
  batch = 0
  print("Computing initial features...")
  initialFeatures <- featureComputeFun(input_graph)
  initial_covariance_mat <- cov(initialFeatures[, -c("node")])
  while(epsilon >= 0) {
    print(paste("nodes remaining:", vcount(input_graph), "current epsilon:", epsilon))
    
    if(vcount(input_graph)<=1){ # EXIT condition
      message("Graph has no nodes left.")
      return(do.call(rbind, result))
    }
    
    # compute features
    if(exists("initialFeatures")) {
      node_features <- initialFeatures
      remove(initialFeatures)
    } else {
      print("Computing features...")
      node_features <- featureComputeFun(input_graph) 
    }
    
    # compute depth of each node based on features
    only_features <- node_features[, !c("node")]
    node_features$depth <- round(1/(1 + mahalanobis.origin(only_features, initial_covariance_mat)), 10)
    
    # select nodes with high depth value for deletion
    nodesToBeDeleted <- node_features[depth >= epsilon]
    
    # if there are any nodes that are to be deleted
    if(nrow(nodesToBeDeleted) > 0) {
      batch <- batch + 1
      
      # delete vertices from input_graph
      input_graph <- delete_vertices(input_graph, nodesToBeDeleted$node)
      # find and add isolated nodes
      isolated <- names(which(degree(input_graph) == 0))
      input_graph <- delete_vertices(input_graph, isolated)
      nodesToBeDeleted <- rbind(nodesToBeDeleted, node_features[node %in% isolated])
      
      # check if after deletion only one node would remain
      # if so, delete the remaining one as well, as we can't compute data depth for one node
      if(nrow(node_features) - nrow(nodesToBeDeleted) == 1) {
        nodesToBeDeleted <- node_features
      }
      #record alpha level for deleted nodes prior to removing them
      nodesToBeDeleted[, alpha := 1-previousEpsilon]
      print(paste("previous eps", previousEpsilon))
      nodesToBeDeleted[, batch := batch]
      result[[batch]] <- nodesToBeDeleted
    }
    else { # No nodes to delete
      print(paste("No nodes to be deleted"))
      epsilon <- reduceEpsilon()
    }
  }
  return(do.call(rbind, result))
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


mahalanobis.origin <- function(data, sigma) {
  # calculate the mahalanobis depth with respect to the origin
  sigma.inv <- try(solve(sigma), silent = TRUE)
  if (!is.matrix(sigma.inv)) {
    sv <- svd(sigma)
    distance <- tryCatch(
      {
        sigma.inv <- sv$v %*% diag(1/sv$d) %*% t(sv$u)
        warning("Inverse of sigma computed by SVD")
        distance <- apply(data, 1, FUN=function(x){sqrt(sum((x-0) ^ 2))})
        return(distance)
      },
      error=function(cond) {
        return(rep(0, nrow(data)))
      }
    )
    return(distance)
  }
  apply(data, 1, function(x) x %*% sigma.inv %*% matrix(x))
}