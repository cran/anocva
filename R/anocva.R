#' If the number of clusters has not been set, estimates it by using Slope criterion in the range {2, 20}.
#'
#' @param dataDist A matrix with n subjects. Each subject has the size of NxN and represents the distances between the elements of the sample.
#' @param r The number of clusters. NULL if it's unknown.
#' @param p Slope adjust parameter.
#' @param maxClust The maximum number of clusters to be tried if estimating optimal number of clusters. The default value is 20.
#' @param clusteringFunction The clustering function that Slope should use.
#' @param criterion The criterion that will be used for estimating the number of clusters. The options are "slope" or "silhouette". If not defined, "slope" will be used.
#'
#' @keywords internal
#'
checkNClust = function(dataDist, r = NULL, p = 1, maxClust = 20, clusteringFunction, criterion = c("slope", "silhouette")){
  if(is.null(r)){
    optimalN = nClustMulti(dataDist = dataDist, p = p, maxClust = maxClust, clusteringFunction = clusteringFunction, criterion = criterion)
    return(optimalN)
  }
  return(r)
}

#' Check for {0,1} Interval Normalization.
#'
#' Verifies if the data is normalized in the range {0,1}.
#' If they are not, the normalization is performed and a warning issued.
#'
#' @param data A matrix of data
#'
#' @return The data matrix normalized in the range {0,1}.
#'
#' @examples
#' set.seed(2000)
#'
#' simuData = runif(100, min = 0.5, max=7)
#' sprintf("The minimum value is %.2f and the maximum is %.2f.", min(simuData), max(simuData))
#'
#' simuData = checkRange01(simuData)
#' sprintf("Now the minimum value is %.2f and the maximum is %.2f.", min(simuData), max(simuData))
#'
#' @export
#'
checkRange01 = function(data){
  if (min(data) < 0 || max(data) > 1){
    normalized = (data - min(data)) / (max(data) - min(data))
    return(normalized)
  }
  return(data)
}

#' Verifies if a clustering function is set. If not,
#' uses the spectral clustering as default clustering method.
#'
#' @param clusteringFunction The actual clustering function. NULL if its unknown.
#'
#' @return A clustering function
#'
#' @keywords internal
#'
checkClusteringFunction = function(clusteringFunction){
  if(is.null(clusteringFunction)){
    #uses spectral clustering as default clustering function
    return(spectralClustering)
  }
  return(clusteringFunction)
}

#' ANalysis Of Cluster VAriability
#'
#' The ANOCVA (ANalysis Of Cluster VAriability) is a non-parametric statistical test to compare clusters with applications in functional magnetic resonance imaging data. The ANOCVA allows us to compare the clustering structure of multiple groups simultaneously and also to identify features that contribute to the differential clustering.
#'
#' The test statistic used is the one proposed by Caetano de Jesus (2017).
#'
#' @param dataDist A matrix with multiple matrices of dissimilarites. Given that a subject with N items (e.g. ROIs)
#'        has a matrix of dissimilarities of size NxN, the dataDist matrix should contain the dissimilarity matrix of
#'        all subjects (n) of all populations, resulting in a three-dimensional (nxNxN) matrix.
#' @param id A list in range {1,2,...,n}, where id[i] identifies the population id for the i-th subject.
#' @param replicates The number of bootstrap replicates. The default value is 1000.
#' @param r The optimal number of clusters. If NULL, then it will be estimated by the slope criterion in the interval {2..20}.
#' @param clusteringFunction Determines the clustering function that will be used. The default function is 'spectralClustering'. The clustering function is suposed to return the clustering labels.
#' @param p Slope adjust parameter. Only used if r is unknown.
#' @param maxClust The maximum number of clusters to be tried if estimating optimal number of clusters. The default value is 20.
#' @param criterion The criterion that will be used for estimating the number of clusters (if r is unknown). The options are "slope" or "silhouette". If not defined, "slope" will be used.
#' @param showElapTime Determines whether the total processing time should be displayed. The default value is TRUE.
#'
#' @return ANOCVA p-values
#'
#' @references
#' Fujita A, Takahashi DY, Patriota AG, Sato JR (2014a) A non-parametric statistical test to
#' compare clusters with applications in functional magnetic resonance imaging data. Statistics
#' in Medicine 33: 4949–4962
#'
#' Vidal MC, Sato JR, Balardin JB, Takahashi DY, Fujita A (2017) ANOCVA in R: a software to
#' compare clusters between groups and its application to the study of autism spectrum disorder.
#' Frontiers in Neuroscience 11:1–8
#'
#' Caetano de Jesus DA. (2017) Evaluation of ANOCVA test for cluster comparison through simulations.
#' Master Dissertation. Institute of Mathematics and Statistics, University of São Paulo.
#'
#' @import cluster
#'
#' @examples
#'
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
#' # Defines a k-means function that returns cluster labels directly
#' myKmeans = function(dist, k){
#'   return(kmeans(dist, k, iter.max = 50, nstart = 5)$cluster)
#' }
#'
#' # Number of subjects in each population
#' nsub = 20
#' # Number of items in each subject
#' nitem = 30
#'
#' # Generate simulated data
#' data = array(NA, c(nsub*2, nitem*2, 2))
#' dataDist = array(NA, c(nsub*2, nitem*2, nitem*2))
#' meanx = 2
#' delta = 0.5
#' # Covariance matrix
#' sigma = matrix(c(0.03, 0, 0, 0.03), 2)
#' for (i in seq(nsub*2)){
#'   sub = rbind(mvrnorm(nitem, mu = c(0, 0), Sigma = sigma ),
#'               mvrnorm(nitem, mu = c(meanx,0), Sigma = sigma))
#'   data[i,,] = sub
#'   # If it's a sample of population 2.
#'   if (i > nsub){
#'     data[i,10,1] = data[i,10,1] + delta
#'   }
#'   # Euclidian distance
#'   dataDist[i,,] = as.matrix(dist(data[i,,]))
#' }
#'
#' # Population 1 subject
#' plot(data[5,,], asp = 1, xlab = '', ylab = '', main = 'Population 1 - subject example')
#'
#' # Population 2 subject
#' plot(data[35,,], asp = 1, xlab = '', ylab = '', main = 'Population 2 - subject example')
#'
#' # The first nsub subjects belong to population 1 while the next nsub subjects belong to population 2
#' id = c(rep(1, nsub), rep(2, nsub))
#'
#' \dontrun{
#' # ANOCVA call with different clustering function (myKmeans) and inside estimation of
#' # the number of clusters (r)
#' res1 = anocva(dataDist, id, replicates=500, r = NULL,
#'               clusteringFunction = myKmeans,
#'               p = 1, criterion = "slope")
#' }
#'
#' # Estimate the number of clusters previously by using Spectral Clustering and Slope criterion
#' r = nClustMulti(dataDist, clusteringFunction = spectralClustering, criterion = 'slope')
#'
#' # Calls ANOCVA statistical test
#' res = anocva(dataDist, id, replicates=500, r = r,
#'              clusteringFunction = spectralClustering,
#'              p = 1, criterion = "slope")
#'
#' # DeltaS p-value
#' res$pValueDeltaS
#'
#' # DeltaSq p-values
#' res$pValueDeltaSq
#'
#' # Identifies which items have significant p-values with a significance level of 0.05.
#' which(res$pValueDeltaSq < 0.05)
#'
#' # Identifies which items have significant FDR adjusted p-values (q-values)
#' # with a significance level of 0.05.
#' qValue = p.adjust(res$pValueDeltaSq, "fdr")
#' which(qValue < 0.05)
#'
#' @export
#'
anocva = function(dataDist, id, replicates = 1000, r=NULL, clusteringFunction = NULL, p = 1, maxClust = 20, criterion = c("slope", "silhouette"), showElapTime = TRUE){
  startTime = proc.time()

  N = dim(dataDist)[2]
  k = max(id)

  dataDist = checkRange01(dataDist)
  clusteringFunction = checkClusteringFunction(clusteringFunction)
  r = checkNClust(dataDist, r = r, p = p, maxClust = maxClust, clusteringFunction = clusteringFunction,
                   criterion = criterion)

  boot = sapply(1:(replicates + 1), function(x) anocvaStats(x, dataDist, id, k, N, r, clusteringFunction))
  mDeltaSq = matrix(unlist(boot[2,]), ncol = N, byrow = TRUE)
  vDeltaS = unlist(boot[1,])

  pValueDeltaS = (sum(vDeltaS[2:length(vDeltaS)] >= vDeltaS[1])) / (replicates + 1)

  pValueDeltaSq = array(NA, N)
  for (q in 1:N){
    pValueQi = (sum(mDeltaSq[2:dim(mDeltaSq)[1], q] >= mDeltaSq[1, q])) / (replicates + 1)
    pValueDeltaSq[q] = pValueQi
  }

  finishTime = proc.time()
  if(showElapTime){
    elapTime = (finishTime - startTime)[3]
    print(elapTime)
  }

  return(list("pValueDeltaS" = pValueDeltaS, "pValueDeltaSq" = pValueDeltaSq))
}
