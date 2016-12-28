#' Spectral clustering
#'
#' Unnormalized spectral clustering function. Uses Partitioning Around Medoids clustering instead of K-means.
#
#' @param W NxN similarity matrix
#' @param k Number of clusters
#'
#' @return Cluster labels
#'
#' @import cluster
#'
#' @examples
#' #install igraph if necessary
#' #install.packages('igraph')
#' #install.packages('cluster')
#'
#' library(igraph)
#' library(cluster)
#' library(anocva)
#'
#' set.seed(2000)
#'
#' #Create a tree graph
#' treeGraph = make_tree(80, children = 4, mode = "undirected")
#'
#' #Visualize the tree graph
#' plot(treeGraph, vertex.size=10, vertex.label=NA)
#'
#' #Get the adjacency matrix of the tree graph
#' adj = as.matrix(get.adjacency(treeGraph))
#'
#' #Cluster the tree graph in to four clusters
#' cluster = spectralClustering(adj, 4)
#'
#' #See the result clustering
#' plot(treeGraph, vertex.size=10, vertex.color = cluster, vertex.label=NA)
#'
#' @export
spectralClustering <- function(W, k){
  n <- ncol(W)
  S = rowSums(W)
  D = diag(S)
  L <- D - W
  U <- (eigen(L)$vectors)[ , ((n-k+1):n)]
  C <- pam(x = U, k = k)
  return(C$clustering)
}
