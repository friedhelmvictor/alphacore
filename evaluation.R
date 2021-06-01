library(parallel)
source("utils.R")
source("figureUtils.R")

###########################################################################
###########################################################################
###                                                                     ###
###                              SECTION 1:                             ###
###                    DATA INPUT AND INITIALIZATION                    ###
###                                                                     ###
###########################################################################
###########################################################################

# Token networks SQLite database
# Please make sure you have downloaded and extracted the database file
# https://mega.nz/file/6WIiEQ5b#s4HmjOkO9WbB_-egCpmjqWoRup9iVUG2RaR6IBtS3UA
tokenDBFile <- "data/tokens/transfers.db"
accountLabelsFile <- "data/tokens/etherscanLabels20200529.csv"

# Reddit crosslinks
# Post_crosslinks_info.tsv is part of http://snap.stanford.edu/conflict/conflict_data.zip
# In the zip file, it can be found in /prediction/detailed_data/
redditCrossLinksFile <- "data/reddit/post_crosslinks_info.tsv"
topSubredditsFile <- "data/reddit/top100_subreddits.csv"

# Get this file from http://opsahl.co.uk/tnet/datasets/openflights.txt
flightsFile <- "data/flights/openflights.txt"
airportRankingFile <- "data/flights/airport_ranking_2011_04_14.csv"

# SoftChainCoin token network
sccNetworkFile <- "data/sccTransfers.csv"

cpuCoreCount <- 8 # Adjust to improve evaluation runtime

out_dir = 'output/'
dir.create(out_dir, showWarnings = TRUE)
### outputs ###

##---------------------------------------------------------------
##     Outputs (only load if you want to skip computation)     --
##---------------------------------------------------------------
# you can read these files if you want to take a look at the outputs without long processing times
tokenResults <- do.call(rbind, lapply(paste0(out_dir, list.files(path = out_dir, pattern="^0x*")), fread))
redditResults <- fread(paste0(out_dir,"reddit.csv"))
flightResults <- fread(paste0(out_dir,"flights.csv"))


##----------------------------------------------------------------
##                Flight routes between airports                --
##----------------------------------------------------------------
flights_g <- graph_from_data_frame(fread(flightsFile)[, list(V1, V2, weight=V3)])
busiest_airports <- fread(airportRankingFile)$RefId

##---------------------------------------------------------------
##                      Reddit crosslinks                      --
##---------------------------------------------------------------
reddit_crosslinks <- fread(redditCrossLinksFile, col.names = c("Source", "Target", "a", "b", "username"))
top100subreddits <- fread(topSubredditsFile)$subreddit
reddit_g <- graph_from_data_frame(reddit_crosslinks[, list(weight = .N), by=list(Source, Target)], directed = T)

##---------------------------------------------------------------
##                      28 Token Networks                      --
##---------------------------------------------------------------
# The database contains the 28 token networks with at least 10 exchange nodes.
# We retrieve all 28 network identifiers from the database, to load them one by one.
tokenNetworks <- getTokenAddresses(tokenDBFile)
exchangeLabels <- unique(fread(accountLabelsFile)[type %in% c("exchange", "dex")]$address)

###########################################################################
###########################################################################
###                                                                     ###
###                              SECTION 2:                             ###
###               COMPARING ALPHACORE WITH OTHER MEASURES               ###
###                                                                     ###
###########################################################################
###########################################################################

### compute ###
redditResults <- evaluateAll(reddit_g, top100subreddits, dataName = "reddit")
fwrite(redditResults, paste0(out_dir,"reddit.csv"), row.names = F)

flightResults <- evaluateAll(flights_g, busiest_airports, dataName = "flights")
fwrite(flightResults, paste0(out_dir,"flights.csv"), row.names = F)

evalFastAndSave <- function(x) {
  res <- evaluateFastPart(getGraph(x, dbLocation = tokenDBFile), exchangeLabels, x)
  fwrite(res, paste0(out_dir,x,"-fast.csv"), row.names = F)
  return(res)
}

evalSlowAndSave <- function(x) {
  res <- evaluateSlowPart(getGraph(x, dbLocation = tokenDBFile), exchangeLabels, x)
  fwrite(res, paste0(out_dir,x,"-slow.csv"), row.names = F)
  return(res)
}

tokenResultsFastPart <- do.call(rbind, mclapply(tokenNetworks, evalFastAndSave, mc.preschedule = F, mc.cores = cpuCoreCount))
# *** WARNING *** the next line will take about 2 days to run through
tokenResultsSlowPart <- do.call(rbind, mclapply(tokenNetworks, evalSlowAndSave, mc.preschedule = F, mc.cores = cpuCoreCount))
# Combine
tokenResults <- rbind(tokenResultsFastPart, tokenResultsSlowPart)

##---------------------------------------------------------------
##                   Aggregate the results                     --
##---------------------------------------------------------------
tokenStats <- createPrecRecallTable(tokenResults)
redditStats <- createPrecRecallTable(redditResults)
colnames(redditStats) <- sapply(colnames(redditStats), function(x) {str_replace(x, "AP@", "rP@")})
colnames(redditStats) <- sapply(colnames(redditStats), function(x) {str_replace(x, "AR@", "rR@")})
flightStats <- createPrecRecallTable(flightResults)
colnames(flightStats) <- sapply(colnames(flightStats), function(x) {str_replace(x, "AP@", "fP@")})
colnames(flightStats) <- sapply(colnames(flightStats), function(x) {str_replace(x, "AR@", "fR@")})

comparisonTable <- merge(merge(tokenStats, redditStats, all = T), flightStats, all = T)[order(`AP@10`, decreasing = T)]
comparisonTable[ , (c(6, 10, 14, 18, 22, 26)) := NULL]

xtable::xtable(comparisonTable[order(Algorithm)])
xtable::xtable(comparisonTableAppendix)

###########################################################################
###########################################################################
###                                                                     ###
###                              SECTION 3:                             ###
###                   ALPHACORE PARAMETER EXPLORATION                   ###
###                                                                     ###
###########################################################################
###########################################################################

evaluateParameters <- function(graph, features, nodesOfInterest, name) {
  scores <- do.call(rbind, apply(expand.grid(startEps=c(0.0001, 0.001, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 1),
                                             stepSize=seq(0.1, 0.5, 0.1)), 1,
                                 function(x) {
                                   res <- scoreResult(alphaCore_wrapped(graph, features, startEps = x[1], stepSize = x[2], expDecay = T),
                                                      nodesOfInterest, dataName = name)
                                   res$algorithmName = paste0("alphaCore(",x[1],"/",x[2],")(",paste0(features, collapse = "+"),")")
                                   res$startEps <- x[1]
                                   res$stepSize <- x[2]
                                   return(res)
                                 }))
  return(scores)
}

pExploreAirports <- evaluateParameters(flights_g, c("outneighborhoodsize"), busiest_airports, "airports") # was inneighborhoodsize
pExploreReddit <- evaluateParameters(reddit_g, c("inneighborhoodsize"), top100subreddits, "reddit")
# *** WARNING *** This will take several minutes
pExploreTokenList <- do.call(rbind, mclapply(tokenNetworks, function(x) {
  evaluateParameters(getGraph(x, dbLocation = tokenDBFile), c("inneighborhoodsize", "outneighborhoodsize"), exchangeLabels, x)
}, mc.preschedule = F, mc.cores = cpuCoreCount)
)
pExploreTokens <- pExploreTokenList[k == 10, list(dataName = "28 token nets", found = mean(found)), by=list(k, startEps, stepSize)]

parameterExplore <- rbind(rbind(pExploreAirports, pExploreReddit), pExploreTokens, fill=T)
parameterExplore$dataName_f <- factor(parameterExplore$dataName, levels = c("28 token nets", "airports", "reddit"))

pePlot <- ggplot(parameterExplore[k==10, list(meanFound = mean(found/k)), by=list(dataName_f, k, startEps, stepSize)]) +
  geom_point(aes(x=as.factor(startEps), y=meanFound, color=as.factor(stepSize), shape=as.factor(stepSize))) +
  facet_grid(dataName_f ~ ., scales = "free_y") +
  labs(x="Start epsilon", y="AP@10", color="Step size", shape="Step size") +
  theme_Publication() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pePlot

ggsave(filename = paste0(out_dir, "pePlot.pdf"), width = 5, height = 5)


###########################################################################
###########################################################################
###                                                                     ###
###                              SECTION 4:                             ###
###                    GRAPH ORDER VS PRECISION PLOT                    ###
###                                                                     ###
###########################################################################
###########################################################################

graphOrderData <- tokenResults[algorithmName == "alphaCore(inneighborhoodsize,outneighborhoodsize)" & k==10]
graphOrderData[, `P@10` := found/k]
ggplot(graphOrderData) +
  geom_point(aes(x=vCount, y=`P@10`)) +
  scale_x_log10() +
  scale_y_continuous(breaks = seq(0,1,0.2)) +
  labs(x="Vertex count in graph", y="P@10") +
  theme_Publication()

ggsave(filename = paste0(out_dir, "graphOrderPlot.pdf"), width = 5, height = 3)


###########################################################################
###########################################################################
###                                                                     ###
###                              SECTION 5:                             ###
###                        SOFTCHAINCOIN EXAMPLE                        ###
###                                                                     ###
###########################################################################
###########################################################################

sccTransfers <- fread(sccNetworkFile)
sccGraph <- graph_from_data_frame(sccTransfers[, list(source, target, weight = normalizeData(amount))])

alphaRes <- alphaCore(sccGraph, featureComputeFun = customNodeFeatures(c("inneighborhoodsize","instrength", "outstrength")),
                      exponentialDecay = T, stepSize = 0.25)
fwrite(alphaRes[, list(id = node, expDecay = paste0("c",alpha))], paste0(out_dir,"sccExpDecay.csv"), row.names = F)
alphaRes <- alphaCore(sccGraph, featureComputeFun = customNodeFeatures(c("inneighborhoodsize","instrength", "outstrength")),
                      exponentialDecay = F, stepSize = 0.25)
fwrite(alphaRes[, list(id = node, linDecay = paste0("c",alpha))], paste0(out_dir,"sccLinDecay.csv"), row.names = F)
