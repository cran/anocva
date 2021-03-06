% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anocvaStats.R
\name{anocvaStats}
\alias{anocvaStats}
\title{Calculates ANOVA statistics for original data and bootstrap replicates.}
\usage{
anocvaStats(idx, dataDist, id, k, N, r, clusteringFunction)
}
\arguments{
\item{idx}{Identifies the bootstrap replicate. 1 means that original data should be used. 2 or more leads to resampling.}

\item{dataDist}{A matrix with n subjects. Each subject has the size of NxN and represents the distances between the elements of the sample.}

\item{id}{A list in range {1,2,...,n}, where id[i] identifies the population id for i-th subject.}

\item{k}{The number of populations.}

\item{N}{Subjects size.}

\item{r}{Optimal number of clusters.}

\item{clusteringFunction}{determines the clustering function that you want to use. The default function is spectralClustering.}
}
\value{
ANOCVA deltaS and deltaSq statistics
}
\description{
Calculates ANOVA statistics for original data and bootstrap replicates.
}
\references{
Fujita A, Takahashi DY, Patriota AG, Sato JR (2014a) A non-parametric statistical test to
compare clusters with applications in functional magnetic resonance imaging data. Statistics
in Medicine 33: 4949–4962

Caetano de Jesus DA. (2017) Evaluation of ANOCVA test for cluster comparison through simulations.
Master Dissertation. Institute of Mathematics and Statistics, University of São Paulo.
}
\keyword{internal}
