#' Estimates the optimal number of clusters using the Slope criterion when silhouette statistics already exists.
#' The optimal number of clusters will be verified in the range {2,..., maxClust}.
#'
#' @param silStats Silhouette statistics.
#' @param p Slope adjust parameter.
#' @param maxClust The max number of clusters to be tried.
#'
#' @return The optimal number of clusters.
#'
#' @references
#' Fujita A, Takahashi DY, Patriota AG (2014b) A non-parametric method to
#' estimate the number of clusters. Computational Statistics & Data Analysis 73:27–39
#'
#' @keywords internal
#'
optimalSlope = function(silStats, p = 1, maxClust){
  slope = array(NA, maxClust - 2)
  for (j in seq(maxClust - 2)){
    slope[j] = - (silStats[j+1] - silStats[j]) * (silStats[j]^p)
  }
  slopeIndice = which(slope == max(slope))[1]
  return(slopeIndice + 1)
}

#' Estimates the optimal number of clusters using the Silhouette criterion when silhouette statistics already exists.
#'
#' @param silStats Silhouette statistics for each possible tested number of clusters.
#'
#' @return The optimal number of clusters.
#'
#' @references
#' Rousseeuw PJ (1987) Sihouettes: a graphical aid to the interpretation and validation of cluster
#' analysis. Journal of Computational and Applied Mathematics 20:53–65
#'
#' @keywords internal
#'
optimalSilhouette = function(silStats){
  return(which.max(silStats) + 1)
}

#' Optimal Number of Clusters Estimation
#'
#' Estimates the optimal number of clusters using either Slope or Silhouette criterion. The optimal number of clusters will be
#' verified in the range {2,..., maxClust}.
#'
#' @param meanDist An NxN matrix that represents the distances between the N items of the sample.
#' @param p Slope adjust parameter.
#' @param maxClust The maximum number of clusters to be tried. The default value is 20.
#' @param clusteringFunction The clustering function to be used.
#' @param criterion The criterion that will be used for estimating the number of clusters. The options are "slope" or "silhouette". If not defined, "slope" will be used.
#'
#' @return The optimal number of clusters.
#'
#' @references
#' Fujita A, Takahashi DY, Patriota AG (2014b) A non-parametric method to
#' estimate the number of clusters. Computational Statistics & Data Analysis 73:27–39
#'
#' Rousseeuw PJ (1987) Sihouettes: a graphical aid to the interpretation and validation of cluster
#' analysis. Journal of Computational and Applied Mathematics 20:53–65
#'
#' @examples
#' # Install packages if necessary
#' # install.packages('MASS')
#' # install.packages('cluster')
#'
#' library(MASS)
#' library(cluster)
#' library(anocva)
#'
#' set.seed(2000)
#'
#' # Defines a k-means function that returns cluster labels directly
#' myKmeans = function(dist, k){
#'   return(kmeans(dist, k, iter.max = 50, nstart = 5)$cluster)
#' }
#'
#' # Generate simulated data
#' nitem = 70
#' sigma = matrix(c(0.04, 0, 0, 0.04), 2)
#' simuData = rbind(mvrnorm(nitem, mu = c(0, 0), Sigma = sigma ),
#'              mvrnorm(nitem, mu = c(3,0), Sigma = sigma),
#'              mvrnorm(nitem, mu = c(2.5,2), Sigma = sigma))
#'
#' plot(simuData, asp = 1, xlab = '', ylab = '', main = 'Data for clustering')
#'
#' # Calculate distances and perform {0,1} normalization
#' distMatrix = as.matrix(dist(simuData))
#' distMatrix = checkRange01(distMatrix)
#'
#' # Estimate the optimal number of clusters
#' r = nClust(meanDist = distMatrix, p = 1, maxClust = 10,
#'            clusteringFunction = myKmeans, criterion = "silhouette")
#' sprintf("The optimal number of clusters found was %d.", r)
#'
#' # K-means Clustering
#' labels = myKmeans(distMatrix, r)
#'
#' plot(simuData, col = labels, asp = 1, xlab = '', ylab = '', main = 'K-means clustered data')
#'
#' @export
#'
nClust = function(meanDist, p = 1, maxClust = 20, clusteringFunction, criterion = c("slope", "silhouette")){

  # Slope is the default criterion for estimating the number of cluster
  criterion = match.arg(criterion)

  N = dim(meanDist)[2]

  clusterList = array(NA, c(maxClust - 1, N))
  silStats = array(NA, maxClust - 1)

  for(i in seq(2, maxClust)){
    clusters = clusteringFunction(1 - meanDist, i)
    clusterList[i - 1, ] = clusters
    S = silhouette(clusters, dmatrix = meanDist)[ , 3]
    silStats[i - 1] = mean(S)
  }

  if(criterion == "slope"){
    return(optimalSlope(silStats, p, maxClust))
  }else{
    return(optimalSilhouette(silStats))
  }
}

#' Multiple Samples Optimal Number of Clusters Estimation
#'
#' Estimates the optimal number of clusters for multiple samples using either Slope or Silhouette criterion. The optimal number of clusters will be
#' verified in the range {2,..., maxClust}. Takes the mean of all samples in order to perform the estimation.
#'
#' @param dataDist An  matrix with n subjects. Each subject has the size of NxN and represents the distances between the elements of the sample.
#' @param p Slope adjust parameter.
#' @param maxClust The maximum number of clusters to be tried.
#' @param clusteringFunction The clustering function to be used.
#' @param criterion The criterion that will be used for estimating the number of clusters. The options are "slope" or "silhouette". If not defined, "slope" will be used.
#'
#' @return The optimal number of clusters.
#'
#' @references
#' Fujita A, Takahashi DY, Patriota AG (2014b) A non-parametric method to
#' estimate the number of clusters. Computational Statistics & Data Analysis 73:27–39
#'
#' Rousseeuw PJ (1987) Sihouettes: a graphical aid to the interpretation and validation of cluster
#' analysis. Journal of Computational and Applied Mathematics 20:53–65
#'
#' @examples
#' # Install packages if necessary
#' # install.packages('MASS')
#' # install.packages('cluster')
#'
#' library(anocva)
#' library(MASS)
#' library(cluster)
#'
#' set.seed(5000)
#'
#' # A k-means function that returns cluster labels directly.
#' myKmeans = function(dist, k){
#'   return(kmeans(dist, k, iter.max = 50, nstart = 5)$cluster)
#' }
#'
#' # Number of subjects in each population
#' nsub = 25
#' # Number of items in each subject
#' nitem = 60
#'
#' # Generate simulated data
#' data = array(NA, c(nsub, nitem*2, 2))
#' data.dist = array(NA, c(nsub, nitem*2, nitem*2))
#' meanx = 2
#' delta = 0.5
#' # Covariance matrix
#' sigma = matrix(c(0.03, 0, 0, 0.03), 2)
#' for (i in seq(nsub)){
#'   sub = rbind(mvrnorm(nitem, mu = c(0, 0), Sigma = sigma ),
#'               mvrnorm(nitem, mu = c(meanx,0), Sigma = sigma))
#'   data[i,,] = sub
#'   data.dist[i,,] = as.matrix(dist(data[i,,]))
#' }
#'
#' # Estimate the optimal number of clusters
#' r = nClustMulti(dataDist = data.dist, p = 1, maxClust = 20,
#'                 clusteringFunction = myKmeans, criterion = "slope")
#' sprintf("The optimal number of clusters found was %d.", r)
#'
#' @export
#'
nClustMulti = function(dataDist, p = 1, maxClust = 20, clusteringFunction, criterion = c("slope", "silhouette")){

  #Calculates the mean distance matrix
  meanDist = colMeans(dataDist)

  return(nClust(meanDist, p, maxClust, clusteringFunction, criterion))
}
