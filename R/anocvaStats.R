#' Calculates ANOVA statistics for original data and bootstrap replicates.
#'
#' @param idx Identifies the bootstrap replicate. 1 means that original data should be used. 2 or more leads to resampling.
#' @param dataDist A matrix with n subjects. Each subject has the size of NxN and represents the distances between the elements of the sample.
#' @param id A list in range {1,2,...,n}, where id[i] identifies the population id for i-th subject.
#' @param k The number of populations.
#' @param N Subjects size.
#' @param r Optimal number of clusters.
#' @param clusteringFunction determines the clustering function that you want to use. The default function is spectralClustering.
#'
#' @return ANOCVA deltaS and deltaSq statistics
#'
#' @keywords internal
#'
anocvaStats <-function(idx, dataDist, id, k, N, r, clusteringFunction){

  Ab = array(0, c(k, N, N))
  Abb = array(0, c(N, N))
  Sj = array(0, c(k, N))
  Sjm = array(0, N)
  S = array(0, N)

  for (j in seq(k)){
    if (idx != 1){
      resample = sample(seq(dim(dataDist)[1]), sum(id == j), replace = TRUE)
    }else{
      resample = which(id == j)
    }

    dataJ = dataDist[resample, , ]
    Ab[j, , ] = colMeans(dataJ)
    Abb = Abb + colSums(dataJ)
  }

  Abb = (Abb / dim(dataDist)[1])

  labels = clusteringFunction(1 - Abb, r)
  S = silhouette(labels, dmatrix = Abb)[, 3]

  for (j in seq(k)){
    Sj[j, ] = silhouette(labels, dmatrix = Ab[j, , ])[, 3]
    Sjm = Sjm + Sj[j, ]
  }
  Sjm = Sjm / k

  deltaS = 0
  for (j in seq(k)){
    deltaS = deltaS + ((S - Sj[j, ]) %*% (S - Sj[j, ]))
  }

  deltaSq = ( S - Sjm) ^ 2

  return(list("deltaS" = deltaS, "deltaSq" = deltaSq))
}
