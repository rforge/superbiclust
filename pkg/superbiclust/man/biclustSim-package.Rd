\name{superbiclust-package}
\alias{superbiclust-package}
\alias{superbiclust}
\docType{package}
\title{
generating super biclusters form the set of biclusters

}
\description{
The package contains a number of functions for computing similarity matrix of the biclusters obtained by a variety of methods, initialization seeds or various parameter settings.
The matrix is used for the construction of hierarchical tree based on overall similarity, row similarity or column similarity to obtain cut-off points for the similarity metric of choice.
Various statistics are output per bicluster set: a number of a given gene(compound) or gene (compound) set has been present in any bicluster of output or per run.
After the tree is cut, th esuper biclusters are obtained in a form of biclust object, which can further be used for plotting of biclusters.}
\details{
\tabular{ll}{
Package: \tab biclustSim\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2012-08-23\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}
}
\author{
Tatsiana Khamiakova <tatsiana.khamiakova@uhasselt.be>
}
\references{
Shi et al. (2010) A bi-ordering approach to linking gene expression	with clinical annotations in gastric cancer. BMC Bioinformatics. 11. pages 477. 
}
\keyword{biclustering, similarity}
