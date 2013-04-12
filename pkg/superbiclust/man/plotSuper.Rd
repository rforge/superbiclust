\name{plotSuper}
\alias{plotSuper}
\title{
Plot gene profiles within biclusters
}

\description{
Function for plotting gene profiles for compounds within constructed super-bicluster}
\usage{
plotSuper(x, data, BiclustSet)
}

\arguments{
  \item{x}{
	a vector, containing indices of biclusters, to be joined for obtaining the robust bicluster  
}
  \item{data}{
	matrix, dataset, from which the bicluster results are obtained  
}
  \item{BiclustSet}{
	a BiclustSet object containing bicluster output  
}
}

\details{
This function plots gene profiles within a robust bicluster 
}
\value{
biclust object
}
\author{
Tatsiana Khamiakova
}
\seealso{
\code{\link{HCLtree}}, \code{\link{plotSuperAll}}, \code{\link{plotProfilesWithinBicluster}}
}
