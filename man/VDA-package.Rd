\name{VDA}
\alias{VDA-package}
\docType{package}
\title{
Multicategory Vertex Discriminant Analysis
}
\description{
This package provides functions to optimize and execute Multicategory Vertex Discriminant Analysis, a method of supervised learning for an outcome with \code{k} predictor categories.  Outcome classification is based on linear discrimination among the vertices of a regular simplex in a \code{k-1}-dimension Euclidean space, where each vertex represents a different category. 

}
\details{
\tabular{ll}{
Package: \tab VDA\cr
Type: \tab Package\cr
Version: \tab 1.02\cr
Date: \tab 2012-02-27\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
}
\author{
Edward Grant, Xia Li, Kenneth Lange, Tong Tong Wu

Maintainer: Edward Grant \email{edward.m.grant@gmail.com}
}
\references{
Lange, K. and Wu, T.T. (2008) An MM Algorithm for Multicategory Vertex Discriminant Analysis. Journal of Computational and Graphical Statistics, Volume 17, No 3, 527-544.
}

\keyword{package}

\examples{
#load dataset from package
data(zoo)

#matrix containing all predictor vectors
x <- zoo[,2:17]

#outcome class vector
y <- zoo[,18]

#run VDA
out <- vda.r(x, y)
}