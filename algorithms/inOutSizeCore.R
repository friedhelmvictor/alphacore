inOutSizeCore <- function(inputGr) {
  # simplify the graph, so that the degree corresponds to the neighborhood size
  inputGr <- simplify(inputGr, remove.multiple = T, remove.loops = T)
  k <- 1
  batch <- 1
  alpha <- 0.5
  beta <- 0.5
  nodes_to_delete <- c()
  cores <- c()
  
  result <- list()
  
  # while there are still vertices left in the graph
  while (vcount(inputGr) > 0){ 
    
    outsize_factor <- neighborhood.size(inputGr, mode = "out", mindist = 1)**alpha
    insize_factor <- neighborhood.size(inputGr, mode = "in", mindist = 1)**beta
    scores <- round((outsize_factor * insize_factor) ** (1/(alpha+beta)), digits = 10)
    scores <- data.table(node=V(inputGr)$name, score=scores)
    nodes_to_delete <- scores[score < k]$node
    batch <- batch+1
    
    if (length(nodes_to_delete) == 0) { # if no nodes can be deleted, increase k by 1
      print(paste("Cant delete anything at",k,"(nothing less than that) - next higher is", min(scores$score)))
      k <- floor(min(scores$score) + 1)
      #k <- ceiling(min(scores$score) + 0.000001)
      print(paste("k is",k, "vcount is", vcount(inputGr)))
    } else {
      message("Number of nodes to be deleted: ", length(nodes_to_delete))
      finished <- data.table(node=nodes_to_delete, k=k-1)
      result[[batch]] <- finished
      inputGr <- delete_vertices(inputGr, nodes_to_delete)
      nodes_to_delete <- c()
    }
  }
  resultDF <- do.call(rbind, result)
  res <- resultDF$k
  names(res) <- resultDF$node
  return(res)
}
