\name{sorensenMat}
\alias{sorensenMat}
\title{
Sorensen similarity Matrix for bicluster output
}
\description{
Computes Sorensen similarity matrix for biclusters in two bicluster sets}
\usage{
sorensenMat(x, y)
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
Give more detail on sorensenMat similarity score}
\value{
matrix of pairwise Sorensen indices
}
\author{
Tatsiana Khamiakova
}
\seealso{
\code{\link{similarity}},\code{\link{jaccardMat}}, \code{\link{ochiaiMat}}, \code{\link{sensitivityMat}},
\code{\link{specificityMat}},\code{\link{kulczynskiMat}} 
}
