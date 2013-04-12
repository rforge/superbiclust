\name{similarity}
\alias{similarity}
\title{
Similarity Matrix for bicluster output
}
\description{
computes similarity matrix for the biclustering output based on one of the pairwise similarity indices of biclusters in a given bicluster set}
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
This function operates on BiclustSet object and computes pairwise similarity based on the common elements between biclusters. 
Several options for similarity index are offered: jaccard (default), kulczynski, sensitivity, specificity, sorensen  and ochiai indices.
}
\value{
a "similarity object" containing similarity matrix
}
\author{
Tatsiana Khamiakova
}
\seealso{
\code{\linkS4class{BiclustSet}}, \code{\link{HCLtree}}, \code{\link{plotSuper}} 
}
\examples{
#compute sensitivity for BiMAX biclusters
 test <- matrix(rnorm(5000), 100, 50)
 test[11:20,11:20] <- rnorm(100, 3, 0.1)
 test[17:26,21:30] <- rnorm(100, 3, 0.1)
 testBin <- binarize(test,2)
 res <- biclust(x=testBin, method=BCBimax(), minr=4, minc=4, number=10)
 BiMaxBiclustSet <-  BiclustSet(res)
 SensitivityMatr<- similarity(BiMaxBiclustSet,index="sensitivity")
 SensitivityMatr
}
