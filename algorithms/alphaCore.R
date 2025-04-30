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
  cov_mat  <- cov(node_features[, -1, with = FALSE])
  cov_mat_inv <- tryCatch(
    solve(cov_mat),
    error = function(e1) {
      message("Covariance singular; falling back to pseudo-inverse.")
      mat <- MASS::ginv(cov_mat)
      if (nrow(mat) == 0 || ncol(mat) == 0) {
        stop("Could not invert covariance: no valid features remain.")
      }
      mat
    }
  )
  #3
  node_features$depth <- mhdOrigin(node_features[, -1, with = FALSE], cov_mat_inv)

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
      input_graph <- delete_vertices(input_graph, nodes)
      #16 
      batch_ID <- batch_ID + 1
      #17
      node_features <- extractFeatures(input_graph, features)
      #18
      node_features$depth <- mhdOrigin(node_features[, -1, with = FALSE], cov_mat_inv)
      
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

#' Compute default node features for alphaCore
#'
#' \code{computeNodeFeaturFun} returns a data.table with indegree and strength
#'
#' @param graph An igraph object
#' @return A data.table with node names, indegree and strength
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


#' Compute specified node features
#'
#' \code{computeNodeFeatures} calculates requested network features for each node
#'
#' @param graph An igraph object
#' @param features A character vector of features to compute
#' @return A data.table with node names and the requested features
computeNodeFeatures <- function(graph, features = c("indegree","strength")) {
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


#' Extract numerical node features from a graph
#'
#' \code{extractFeatures} returns a data.table of node features
#'
#' Extracts numerical vertex attributes from the graph. If specified features
#' don't exist as vertex attributes, they will be computed on demand.
#'
#' @param graph An igraph object
#' @param features Features to extract: "all" for all numeric attributes, or a character vector of specific features
#' @return A data.table with node names and requested features
extractFeatures <- function(graph, features = "all") {
  #If the user explicitly requested a feature list, compute those
  if (!identical(features, "all")) {
    df <- computeNodeFeatures(graph, features = features)
    missing <- setdiff(features, names(df))
    if (length(missing) > 0L) {
      stop(
        "Vertex attribute(s) not found or not numeric: ",
        paste(missing, collapse = ", "),
        call. = FALSE
      )
    }
    return(df)
  }

  #Otherwise, pull *all* numeric vertex attributes off the graph
  attrs <- vertex.attributes(graph)
  numeric_names <- names(attrs)[sapply(attrs, is.numeric)]

  #If none are present, fall back to the original indegree+strength
  if (length(numeric_names) == 0L) {
    message(
      "No numeric vertex attributes found; reverting to indegree+strength."
    )
    return(computeNodeFeaturFun(graph))
  }

  #Build a data.table: one column 'node', then each numeric attribute
  dt <- data.table(node = V(graph)$name)
  for (nm in numeric_names) dt[[nm]] <- attrs[[nm]]
  return(dt[])
}



#' Calculate the Mahalanobis depth to the origin
#'
#' \code{mhdOrigin} calculates the Mahalanobis depth of data points to the origin
#'
#' @param data A matrix or data.frame of multivariate data
#' @param sigma_inv The inverse covariance matrix
#' @return A vector of depth values
mhdOrigin <- function(data, sigma_inv) {
  origin <- rep(0,ncol(data)) # c(0,0,...)
  # We reuse the Mahalanobis distance implementation of the stats package,
  # which returns the squared Mahalanobis distance: (x - μ)' Σ^-1 (x - μ) = D^2
  # To arrive at the Mahalanobis Depth to the origin, we only need to add 1 and
  # take the reciprocal.
  return((1 + pmax(stats::mahalanobis(data, center = origin, sigma_inv, inverted = TRUE), 0))^-1)
}