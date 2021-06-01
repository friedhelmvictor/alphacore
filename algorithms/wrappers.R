library(data.table)
library(igraph)
library(tnet)
library(Rcpp)
sourceCpp("algorithms/weightedKCore.cpp")
source("algorithms/inOutSizeCore.R")
source("algorithms/alphaCore.R")

###########################################################################
###########################################################################
###                                                                     ###
###                         UTILITY FUNCTIONS:                          ###
###        Wrap some measures to always return one named vector         ###
###                                                                     ###
###########################################################################
###########################################################################

page_rank_wrapped <- function(graph) {
  print("Computing page rank")
  return(page_rank(graph)$vector)
}

degree_centrality_wrapped <- function(graph, mode="all") {
  print("Computing degree centrality")
  res <- centr_degree(graph, mode = mode)$res
  names(res) <- V(graph)$name
  return(res)
}

w_degree_centrality_wrapped <- function(graph, mode="all") {
  mapping <- data.table(node=as.numeric(V(graph)), name=V(graph)$name)
  graph_df <- as.data.table(as_edgelist(graph, names=F))
  colnames(graph_df) <- c("i", "j")
  graph_df$w <- edge_attr(graph, "weight")
  graph_df <- graph_df[, list(w = sum(w)), by=list(i, j)]
  if(mode=="all") {
    graph_df <- as.data.frame(symmetrise_w(as.matrix(graph_df)))
  }
  res <- degree_w(net=graph_df, measure=c("alpha"), alpha=0.5, type=mode)
  res <- merge(as.data.frame(res), mapping)
  result <- res$alpha
  names(result) <- res$name
  return(result)
}

betweenness_centrality_wrapped <- function(graph) {
  print("Computing betweenness centrality")
  res <- centr_betw(graph)$res
  names(res) <- V(graph)$name
  return(res)
}

w_betweenness_centrality_wrapped <- function(graph) {
  mapping <- data.table(node=as.numeric(V(graph)), name=V(graph)$name)
  graph_df <- as.data.table(as_edgelist(graph, names=F))
  colnames(graph_df) <- c("i", "j")
  graph_df$w <- edge_attr(graph, "weight")
  graph_df <- graph_df[, list(w = sum(w)), by=list(i, j)]
  res <- betweenness_w(graph_df, alpha=0.5)
  res <- merge(as.data.frame(res), mapping)
  result <- res$betweenness
  names(result) <- res$name
  return(result)
}

closeness_centrality_wrapped <- function(graph) {
  print("Computing closeness centrality")
  res <- centr_clo(graph)$res
  names(res) <- V(graph)$name
  return(res)
}

w_closeness_centrality_wrapped <- function(graph) {
  mapping <- data.table(node=as.numeric(V(graph)), name=V(graph)$name)
  graph_df <- as.data.table(as_edgelist(graph, names=F))
  colnames(graph_df) <- c("i", "j")
  graph_df$w <- edge_attr(graph, "weight")
  graph_df <- graph_df[, list(w = sum(w)), by=list(i, j)]
  res <- closeness_w(graph_df, alpha=0.5)
  res <- merge(as.data.frame(res), mapping)
  result <- res$closeness
  names(result) <- res$name
  return(result)
}

only_data_depth_wrapped <- function(graph, features) {
  edge_attr(graph, "weight") <- normalizeData(edge_attr(graph, "weight"))
  featureDF <- customNodeFeatures(features)(graph)
  initial_covariance_mat <- cov(featureDF[, -c("node")])
  featureDF$depth <- round(1/(1 + mahalanobis.origin(featureDF[, !c("node")], initial_covariance_mat)), 10)
  res <- 1 - featureDF$depth
  names(res) <- featureDF$node
  return(res)
}


alphaCore_wrapped <- function(graph, features=c("inneighborhoodsize", "outneighborhoodsize", "strength"), startEps = 0.1, stepSize = 0.1, expDecay = T) {
  res <- alphaCore(graph,
                   customNodeFeatures(features),
                   startEpsilon = startEps,
                   stepSize = stepSize,
                   exponentialDecay = expDecay)
  namedResult <- res$batch # TODO order by batch AND alpha!
  names(namedResult) <- res$node
  return(namedResult)
}
