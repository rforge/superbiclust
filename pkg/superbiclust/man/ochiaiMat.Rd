\name{ochiaiMat}
\alias{ochiaiMat}
\title{
Ochiai similarity Matrix for bicluster output
}
\description{
Computes Ochiai similarity matrix for biclusters in two bicluster sets}
\usage{
ochiaiMat(x, y)
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
Give more detail on ochiaiMat similarity score}
\value{
matrix of pairwise Ochiai indices
}
\author{
Tatsiana Khamiakova
}
\seealso{
\code{\link{similarity}},\code{\link{jaccardMat}}, \code{\link{ochiaiMat}}, \code{\link{sensitivityMat}},
\code{\link{specificityMat}},\code{\link{sorensenMat}} 
}
