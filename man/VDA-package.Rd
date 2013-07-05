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
Version: \tab 1.3\cr
Date: \tab 2013-Jul-05\cr
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

#run VDA (ridge penalty)
out <- vda.r(x, y)

#Predict five cases based on VDA 
fivecases <- matrix(0,5,16)
fivecases[1,] <- c(1,0,0,1,0,0,0,1,1,1,0,0,4,0,1,0)
fivecases[2,] <- c(1,0,0,1,0,0,1,1,1,1,0,0,4,1,0,1)
fivecases[3,] <- c(0,1,1,0,1,0,0,0,1,1,0,0,2,1,1,0)
fivecases[4,] <- c(0,0,1,0,0,1,1,1,1,0,0,1,0,1,0,0)
fivecases[5,] <- c(0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0)
predict(out, fivecases)

#run vda.le (lasso and euclidean penalty)
outLE <- vda.le(x, y)

#Predict five cases based on VDA 
fivecases <- matrix(0,5,16)
fivecases[1,] <- c(1,0,0,1,0,0,0,1,1,1,0,0,4,0,1,0)
fivecases[2,] <- c(1,0,0,1,0,0,1,1,1,1,0,0,4,1,0,1)
fivecases[3,] <- c(0,1,1,0,1,0,0,0,1,1,0,0,2,1,1,0)
fivecases[4,] <- c(0,0,1,0,0,1,1,1,1,0,0,1,0,1,0,0)
fivecases[5,] <- c(0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0)
predict(outLE, fivecases)
}