\name{specificityMat}
\alias{specificityMat}
\title{
Specificity Matrix for bicluster output
}
\description{
Computes specificity matrix for biclusters in two bicluster sets}
\usage{
specificityMat(x, y)
}

\arguments{
  \item{x}{
BiclustSet object containing row and column indicators of bicluster membership, number of biclusters
}
 \item{y}{
BiclustSet object containing row and column indicators of bicluster membership, number of biclusters
}}
\details{
This function operates on BiclustSet objects and computes pairwise similarity based on the common elements between biclusters. 
Give more detail on specificity}
\value{
matrix of pairwise specificities
}
\author{
Tatsiana Khamiakova
}
\seealso{
\code{\link{similarity}},\code{\link{jaccardMat}}, \code{\link{ochiaiMat}}, \code{\link{kulczynskiMat}},
\code{\link{sensitivityMat}},\code{\link{sorensenMat}} 
}
