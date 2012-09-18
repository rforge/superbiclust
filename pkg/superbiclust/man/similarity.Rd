\name{similarity}
\alias{similarity}
\title{
Similarity Matrix for bicluster output
}
\description{
computes similarity matrix for the biclustering output based on one of the indices}
\usage{
similarity(x, index = "jaccard")
}

\arguments{
  \item{x}{
BiclustSet object containing row and column indicators of bicluster membership, number of biclusters
}
  \item{index}{
similarity index for the biclusters in output
}
}
\details{
This function operates on BiclustSet object and computes pairwise similarity based on the common elements between biclusters. Several options for similarity index are offered: jaccard (default), kulczynski, sensitivity, specificity  and ochiai indices.
}
\value{
a "similarity object" containing similarity matrix
}
\references{
Later references
}
\author{
Tatsiana Khamiakova
}

