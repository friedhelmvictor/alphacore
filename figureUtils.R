library(ggplot2)
library(stringr)

# plotting utils

theme_Publication <- function(base_size=14, base_family="Times New Roman") {
  library(grid)
  library(ggthemes)
  library(extrafont)
  #extrafont::loadfonts()
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(#face = "bold",
      size = rel(1.2), hjust = 0.5),
      text = element_text(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      axis.title = element_text(size = rel(1)),
      axis.title.y = element_text(angle=90,vjust =2),
      axis.title.x = element_text(vjust = -0.2),
      axis.text = element_text(), 
      axis.line = element_line(colour="black"),
      axis.ticks = element_line(),
      panel.grid.major = element_line(colour="#f0f0f0"),
      panel.grid.minor = element_blank(),
      legend.key = element_rect(colour = NA),
      legend.title = element_text(face="italic"),
      plot.margin=unit(c(10,5,5,5),"mm"),
      strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
      strip.text = element_text(size = 14)
    ))
}

createPrecRecallTable <- function(resultDF) {
  resultStats <- resultDF[, list("AP@k" = mean(found/k),
                                  "AR@k" = mean(found/nodesOfInterestInGraph), count = .N),
                           by=list(algorithmName, k)]
  
  at10Res <- resultStats[k == 10]
  at10Res[, k := NULL]
  setnames(at10Res, "AP@k", "AP@10")
  setnames(at10Res, "AR@k", "AR@10")
  at20Res <- resultStats[k == 20]
  at20Res[, k := NULL]
  setnames(at20Res, "AP@k", "AP@20")
  setnames(at20Res, "AR@k", "AR@20")
  at50Res <- resultStats[k == 50]
  at50Res[, k := NULL]
  setnames(at50Res, "AP@k", "AP@50")
  setnames(at50Res, "AR@k", "AR@50")
  at100Res <- resultStats[k == 100]
  at100Res[, k := NULL]
  setnames(at100Res, "AP@k", "AP@100")
  setnames(at100Res, "AR@k", "AR@100")
  
  table <- merge(merge(merge(at10Res, at20Res), at50Res), at100Res)
  table[, count := NULL]
  table$"Input features" <- unname(str_replace(sapply(table$algorithmName, function(x) {tail(strsplit(x, "\\(")[[1]], 1)}), "\\)", ""))
  table$algorithmName <- unname(sapply(table$algorithmName, function(x) {head(strsplit(x, "\\(")[[1]], 1)}))
  setnames(table, "algorithmName", "Algorithm")
  setcolorder(table, c("Algorithm", "Input features", "AP@10", "AP@20", "AP@50", "AP@100", "AR@10", "AR@20", "AR@50", "AR@100"))
  
  table <- table[order(Algorithm, -`AP@10`)]
  return(table)
}
