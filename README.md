# AlphaCore

These are the supplementary files for the [AlphaCore KDD paper](https://dl.acm.org/doi/10.1145/3447548.3467322).

## Usage
Until this is packaged, source the file algorithms/alphaCore.R

### Get a graph with named vertices and edge weights:
```
g <- erdos.renyi.game(200, 2/200, directed = T)
E(g)$weight <- 1:ecount(g)
V(g)$name <- paste("v", 1:vcount(g), sep="")
```

### Run alphaCore with default parameters
```
> alphaCore(g)
     node     alpha batch
  1:   v1 0.3281309    19
  2:   v2 0.3281309    30
  3:   v3 0.3724843    34
  4:   v4 0.3281309    23
  5:   v5 0.3281309    19
 ---                     
196: v196 0.3281309    28
197: v197 0.0000000     1
198: v198 0.3724843    34
199: v199 0.3281309    23
200: v200 0.3281309    23
```

### Run alphaCore with builtin edge property functions
```
> alphaCore(g, featureComputeFun = customNodeFeatures(c("indegree", "triangles")))
     node     alpha batch
  1:   v1 0.3060040    14
  2:   v2 0.0000000     6
  3:   v3 0.0000000     5
  4:   v4 0.0000000     5
  5:   v5 0.0000000     6
 ---                     
196: v196 0.6381691    19
197: v197 0.0000000     5
198: v198 0.0000000     5
199: v199 0.6381691    19
200: v200 0.0000000     8
```

### Python implementation (contributed by Jason Zhu!) example:
```
> import networkx as nx
> G = nx.erdos_renyi_graph(n=200, seed=1, p=2/200, directed=True)
> for idx, (u,v,w) in enumerate(G.edges(data=True)):
      w['value'] = idx

> alphaCore(G)
   nodeID 	alpha   batchID
0 	    18 	  0.0         0
1 	    75 	  0.0 	     0
2 	    78 	  0.0 	     0
3 	    25 	  0.3 	     5
4 	    91 	  0.3 	     5
... 	... 	  ... 	   ...
195 	  8 	  0.7 	    27
196 	131 	  0.7 	    27
197 	185 	  0.7 	    27
198 	192 	  0.7 	    27
199 	158 	  0.7 	    28
```

## Evaluation

In order to run, you first need to download the three datasets:

### Token Networks (11.4GB)
The files are hosted at: https://zenodo.org/record/4898412
Store the transfers.db in data/tokens/transfers.db
The matching exchangeLabels.csv file should be placed at data/tokens/exchangeLabels.csv

### Reddit Crosslinks
The reddit crosslinks are part of http://snap.stanford.edu/conflict/conflict_data.zip
In the zip file, they can be found in /prediction/detailed_data/
Place the file at the location data/reddit/post_crosslinks_info.tsv

### Flight routes
Get this file from http://opsahl.co.uk/tnet/datasets/openflights.txt
and store it in data/flights/openflights.txt.

### Running the evaluation
Open up the file evaluation.R from the main directory.

This is the main file to run the evaluation.
Input your number of CPU cores and the path to the extracted database file.

To run the entire evaluation, ideally do so on a server with Rscript, as some algorithms, i.e. some centralities depending on all pairs shortest path computations, and our own implementation of weighted k-core have very long running times. The full executing takes between 1-3 days, depending on hardware.

AlphaCore with an exponentially decaying step size however is quite fast, even though it is only an R implementation.

# Citing
Please use the following BibTeX entry:
```
@inproceedings{10.1145/3447548.3467322,
author = {Victor, Friedhelm and Akcora, Cuneyt G. and Gel, Yulia R. and Kantarcioglu, Murat},
title = {Alphacore: Data Depth Based Core Decomposition},
year = {2021},
isbn = {9781450383325},
publisher = {Association for Computing Machinery},
address = {New York, NY, USA},
url = {https://doi.org/10.1145/3447548.3467322},
doi = {10.1145/3447548.3467322},
booktitle = {Proceedings of the 27th ACM SIGKDD Conference on Knowledge Discovery &amp; Data Mining},
pages = {1625–1633},
numpages = {9},
keywords = {core decomposition, networks, data depth},
location = {Virtual Event, Singapore},
series = {KDD '21}
}
```
