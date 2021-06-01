#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <limits>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
IntegerVector wkCore(SEXP graph, std::string mode="all"){
  
  Environment igraph("package:igraph");
  Function as_adjacency_matrix = igraph["as_adjacency_matrix"];
  Function as_undirected = igraph["as.undirected"];
  Function vertex_attr = igraph["vertex_attr"];
  
  sp_mat weighted_adjacency;
  // get the weighted adjacency matrix. This already sums the weights for each edge
  if(mode == "all") {
    weighted_adjacency = as<sp_mat>(as_adjacency_matrix(Named("graph", as_undirected(graph)), Named("sparse", true), Named("attr", "weight")));
  } else if(mode == "in") {
    weighted_adjacency = as<sp_mat>(as_adjacency_matrix(Named("graph", graph), Named("sparse", true), Named("attr", "weight")));
  } else if(mode == "out") {
    weighted_adjacency = as<sp_mat>(as_adjacency_matrix(Named("graph", graph), Named("sparse", true), Named("attr", "weight"))).t();
  } else {
    throw std::range_error("mode must be one of \"in\", \"out\" or \"all\"");
  }
  
  // Apply normalization by mean weight and divide by minimum weight, so that the minimum weight will be 1
  double meanWeight = mean(nonzeros(weighted_adjacency));
  weighted_adjacency = weighted_adjacency/meanWeight;
  double minWeight = min(nonzeros(weighted_adjacency));
  weighted_adjacency = weighted_adjacency/minWeight;
  
  int vertex_count = weighted_adjacency.n_rows;
  LogicalVector vertex_removed (vertex_count);
  IntegerVector vertex_core (vertex_count);
  
  int degree;
  double strength, score;
  bool no_vertex_removed;
  
  long long k_shell = 1;
  long long old_k_shell = 1;
  long double min_score = std::numeric_limits<long double>::max();
  int vertices_remaining = vertex_count;
  
  Progress p(vertex_count, true);
  
  do {
    no_vertex_removed = true;
    min_score = std::numeric_limits<double>::max();
    
    for(int i=0; i<vertex_count; i++) {
      if(vertex_removed[i])
        continue;
      
      strength = sum(weighted_adjacency.col(i));
      degree = weighted_adjacency.col(i).n_nonzero;
      //score = sqrt(degree * strength); FUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUCK
      score = Rf_fround(sqrt(degree) * sqrt(strength), 10);
      //Rcout << "Vertex " << i << " has score " << score << " d:" << degree << " s:" << strength << " k is " << k_shell << std::endl;
      
      if((score >= k_shell) && (score < min_score)) {
        min_score = score;
      }
      if(score >= k_shell) {
        //vertex_core[i] = k_shell - 1;
      } else {
        vertex_core[i] = k_shell - 1;
        vertices_remaining--;
        vertex_removed[i] = true;
        //Rcout << "removed " << i << std::endl;
        no_vertex_removed = false;
        weighted_adjacency.col(i).zeros();
        weighted_adjacency.row(i).zeros();
      }
    }
    //Rcout << "Remaining: " << vertices_remaining << " vertices." << std::endl;
    if(no_vertex_removed) {
      old_k_shell = k_shell;
      k_shell = std::floor(min_score);
      if(k_shell == old_k_shell) {
        k_shell++;
      }
    }
    p.update(vertex_count - vertices_remaining);
  } while(vertices_remaining > 0);
  
  vertex_core.attr("names") = vertex_attr(graph, "name");
  
  return vertex_core;
}