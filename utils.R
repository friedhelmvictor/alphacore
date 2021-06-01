library(data.table)
library(RSQLite)
source("algorithms/wrappers.R")

getTokenAddresses <- function(dbLocation) {
  sqlitecon = dbConnect(RSQLite::SQLite(), dbname=dbLocation)
  tokenNetworks <- dbGetQuery(sqlitecon, "SELECT DISTINCT token_address AS address FROM token_transfers")$address
  dbDisconnect(sqlitecon)
  return(tokenNetworks)
}

# function to retrieve a token network as an igraph object
getGraph <- function(networkAddress, dbLocation) {
  con = dbConnect(RSQLite::SQLite(), dbname=dbLocation)
  print(paste("Fetching transfers for token network", networkAddress))
  transfers <- as.data.table(dbGetQuery(con, paste0("SELECT from_address, to_address, value, block_number FROM token_transfers WHERE token_address='",networkAddress,"';"))
  )[,list(source = from_address, target = to_address, amount = as.numeric(value), blockNumber = as.numeric(block_number))]
  transfers <- transfers[amount > 0]
  transfers <- transfers[amount >= max(transfers$amount)/(10^10), list(source, target, weight=amount, blockNumber)]
  graph <- graph_from_data_frame(transfers)
  dbDisconnect(con)
  return(graph)
}

#' Evaluate an algorithm's result
#' 
#' @param nodeRanking A named vector of numerics, where the number corresponds to the rank, and name to the node.
#' @param nodesOfInterest A vector consisting of nodes of interest (should all be in node column of nodeRanking)
#' @param k A vector of integers, corresponding to the scores at the top k results
#' @param dataName A name for this dataset / graph
#' @param algoName A name for describing the algorithm that produced the ranking
#' @return A DataDable   \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' add(10, 1)
scoreResult <- function(nodeRanking, nodesOfInterest, k = c(10,20,50,100), dataName = "newData", algoName = "newAlgo") {
  sortedRankingDF <- data.table(node=names(nodeRanking), k=nodeRanking)[order(-k)]
  vCount <- nrow(sortedRankingDF)
  nodesOfInterestInGraph = length(intersect(nodesOfInterest, sortedRankingDF$node))
  
  res <- do.call(rbind, lapply(k, function(x) {
    data.table(dataName = dataName,
               vCount = vCount,
               nodesOfInterestInGraph = nodesOfInterestInGraph,
               algorithmName = algoName,
               found = length(intersect(head(sortedRankingDF, x)$node, nodesOfInterest)),
               k=x)
  }))
  
  return(res)
}

evaluateAlphaCoreCombinations <- function(graph, dataName, features, startEpsilon, stepSize) {
  getCombinations <- function(features) { # https://stackoverflow.com/a/18718066
    n <- length(features)
    masks <- 2^(1:n-1)
    lapply( 2:2^n-1, function(u) features[ bitwAnd(u, masks) != 0 ] )
  }
  combs <- getCombinations(features)
  
  res <- do.call(rbind, lapply(combs, function(x) {scoreResult(alphaCore_wrapped(graph, x), top100subreddits, dataName = dataName, algoName = paste0("alphaCore(",paste0(x, collapse = ","),")"))}))
  return(res)
}

evaluateAlphaCore <- function(graph, nodesOfInterest, dataName = "newData") {
  res <- rbind(
    scoreResult(alphaCore_wrapped(graph, c("inneighborhoodsize", "outneighborhoodsize")), nodesOfInterest, dataName = dataName, algoName = "alphaCore(inneighborhoodsize,outneighborhoodsize)"),
    scoreResult(alphaCore_wrapped(graph, c("outneighborhoodsize")), nodesOfInterest, dataName = dataName, algoName = "alphaCore(outneighborhoodsize)"),
    scoreResult(alphaCore_wrapped(graph, c("outneighborhoodsize", "instrength", "outstrength")), nodesOfInterest, dataName = dataName, algoName = "alphaCore(outneighborhoodsize,instrength,outstrength)"),
    scoreResult(alphaCore_wrapped(graph, c("inneighborhoodsize", "outneighborhoodsize", "instrength", "outstrength")), nodesOfInterest, dataName = dataName, algoName = "alphaCore(inneighborhoodsize,outneighborhoodsize,instrength,outstrength)"),
    scoreResult(alphaCore_wrapped(graph, c("inneighborhoodsize")), nodesOfInterest, dataName = dataName, algoName = "alphaCore(inneighborhoodsize)"),
    scoreResult(alphaCore_wrapped(graph, c("inneighborhoodsize", "instrength")), nodesOfInterest, dataName = dataName, algoName = "alphaCore(inneighborhoodsize,instrength)")
  )
}

evaluateKCore <- function(graph, nodesOfInterest, dataName = "newData") {
  res <- rbind(
    scoreResult(coreness(simplify(graph), "all"), nodesOfInterest, dataName = dataName, algoName = "k-core(simplified,all)"),
    scoreResult(coreness(simplify(graph), "all"), nodesOfInterest, dataName = dataName, algoName = "k-core(simplified,all)"),
    scoreResult(coreness(simplify(graph), "out"), nodesOfInterest, dataName = dataName, algoName = "k-core(simplified,out)"),
    scoreResult(coreness(simplify(graph), "in"), nodesOfInterest, dataName = dataName, algoName = "k-core(simplified,in)")
  )
}

evaluateWKCore <- function(graph, nodesOfInterest, dataName = "newData") {
  res <- rbind(
    scoreResult(wkCore(graph, "all"), nodesOfInterest, dataName = dataName, algoName = "wk-core(neighborhoodsize,strength)"),
    scoreResult(wkCore(graph, "in"), nodesOfInterest, dataName = dataName, algoName = "wk-core(inneighborhoodsize,instrength)"),
    scoreResult(wkCore(graph, "out"), nodesOfInterest, dataName = dataName, algoName = "wk-core(outneighborhoodsize,outstrength)"),
    scoreResult(inOutSizeCore(graph), nodesOfInterest, dataName = dataName, algoName = "wk-core(outneighborhoodsize,inneighborhoodsize)")
  )
}

evaluatePageRank <- function(graph, nodesOfInterest, dataName = "newData") {
  res <- scoreResult(page_rank_wrapped(graph), nodesOfInterest, dataName = dataName, algoName = "Pagerank")
  return(res)
}

evaluateDegCentr <- function(graph, nodesOfInterest, dataName = "newData") {
  res <- rbind(
    scoreResult(degree_centrality_wrapped(graph, "in"), nodesOfInterest, dataName = dataName, algoName = "InDegree Centrality"),
    scoreResult(degree_centrality_wrapped(graph, "out"), nodesOfInterest, dataName = dataName, algoName = "OutDegree Centrality")
  )
}

evaluateWDegCentr <- function(graph, nodesOfInterest, dataName = "newData") {
  res <- rbind(
    scoreResult(w_degree_centrality_wrapped(graph, "in"), nodesOfInterest, dataName = dataName, algoName = "Weighted InDegree Centrality"),
    scoreResult(w_degree_centrality_wrapped(graph, "out"), nodesOfInterest, dataName = dataName, algoName = "Weighted OutDegree Centrality")
  )
}

evaluateBetwCentr <- function(graph, nodesOfInterest, dataName = "newData") {
  scoreResult(betweenness_centrality_wrapped(graph), nodesOfInterest, dataName = dataName, algoName = "Betweenness Centrality")
}

evaluateWBetwCentr <- function(graph, nodesOfInterest, dataName = "newData") {
  scoreResult(w_betweenness_centrality_wrapped(graph), nodesOfInterest, dataName = dataName, algoName = "Weighted Betweenness Centrality")
}

evaluateCloCentr <- function(graph, nodesOfInterest, dataName = "newData") {
  scoreResult(closeness_centrality_wrapped(graph), nodesOfInterest, dataName = dataName, algoName = "Closeness Centrality")
}

evaluateWCloCentr <- function(graph, nodesOfInterest, dataName = "newData") {
  scoreResult(w_closeness_centrality_wrapped(graph), nodesOfInterest, dataName = dataName, algoName = "Weighted Closeness Centrality")
}

evaluateDataDepth <- function(graph, nodesOfInterest, dataName = "newData") {
  res <- rbind(
    scoreResult(only_data_depth_wrapped(graph, c("inneighborhoodsize", "outneighborhoodsize")), nodesOfInterest, dataName = dataName, algoName = "dataDepth(inneighborhoodsize,outneighborhoodsize)"),
    scoreResult(only_data_depth_wrapped(graph, c("outneighborhoodsize")), nodesOfInterest, dataName = dataName, algoName = "dataDepth(outneighborhoodsize)"),
    scoreResult(only_data_depth_wrapped(graph, c("outneighborhoodsize", "instrength", "outstrength")), nodesOfInterest, dataName = dataName, algoName = "dataDepth(outneighborhoodsize,instrength,outstrength)"),
    scoreResult(only_data_depth_wrapped(graph, c("inneighborhoodsize", "outneighborhoodsize", "instrength", "outstrength")), nodesOfInterest, dataName = dataName, algoName = "dataDepth(inneighborhoodsize,outneighborhoodsize,instrength,outstrength)"),
    scoreResult(only_data_depth_wrapped(graph, c("inneighborhoodsize")), nodesOfInterest, dataName = dataName, algoName = "dataDepth(inneighborhoodsize)"),
    scoreResult(only_data_depth_wrapped(graph, c("inneighborhoodsize", "instrength")), nodesOfInterest, dataName = dataName, algoName = "dataDepth(inneighborhoodsize,instrength)")
  )
}

evaluateAll <- function(graph, nodesOfInterest, dataName = "newData") {
  res <- rbind(
    evaluateAlphaCore(graph, nodesOfInterest, dataName),
    evaluateKCore(graph, nodesOfInterest, dataName),
    evaluateWKCore(graph, nodesOfInterest, dataName),
    evaluatePageRank(graph, nodesOfInterest, dataName),
    evaluateDegCentr(graph, nodesOfInterest, dataName),
    evaluateWDegCentr(graph, nodesOfInterest, dataName),
    evaluateBetwCentr(graph, nodesOfInterest, dataName),
    evaluateWBetwCentr(graph, nodesOfInterest, dataName),
    evaluateCloCentr(graph, nodesOfInterest, dataName),
    evaluateWCloCentr(graph, nodesOfInterest, dataName),
    evaluateDataDepth(graph, nodesOfInterest, dataName)
  )
}

evaluateFastPart <- function(graph, nodesOfInterest, dataName = "newData") {
  res <- rbind(
    evaluateAlphaCore(graph, nodesOfInterest, dataName),
    evaluateKCore(graph, nodesOfInterest, dataName),
    evaluatePageRank(graph, nodesOfInterest, dataName),
    evaluateDegCentr(graph, nodesOfInterest, dataName),
    evaluateWDegCentr(graph, nodesOfInterest, dataName),
    evaluateDataDepth(graph, nodesOfInterest, dataName)
  )
}

evaluateSlowPart <- function(graph, nodesOfInterest, dataName = "newData") {
  res <- rbind(
    evaluateWKCore(graph, nodesOfInterest, dataName),
    evaluateBetwCentr(graph, nodesOfInterest, dataName),
    evaluateWBetwCentr(graph, nodesOfInterest, dataName),
    evaluateCloCentr(graph, nodesOfInterest, dataName),
    evaluateWCloCentr(graph, nodesOfInterest, dataName)
  )
}
