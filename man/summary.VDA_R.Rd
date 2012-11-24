\name{summary.VDA_R}
\alias{summary.VDA_R}
\title{
Summary for a VDA_R object
}
\description{
Takes a fitted VDA_R object produced by \code{\link{VDA_R}} and produces various useful summaries from it. 
}
\usage{
\method{summary}{VDA_R}(object, ...)
}
\arguments{
  \item{object}{
An object of class 'VDA_R', usually the result of a call to \code{\link{VDA_R}}.
}
  \item{\dots}{
Not used. 
}
}
\details{
The function prints the number of cases, the number of classes, and the number of features in \code{object}, of class VDA_R. It also prints the lambda used in the analysis. Additionally, it prints the coefficients and the resulting predictions made by Vertex Discriminant Analysis on the training set and the following training error.
}

\author{
Edward Grant, Xia Li, Kenneth Lange, Tong Tong Wu

Maintainer: Edward Grant \email{egrant@umd.edu}
}
\seealso{
\code{\link{VDA_R}}, \code{\link{print.VDA_R}}
}
\examples{

# load data
data(zoo)

# feature matrix without intercept
x<-zoo[,2:17]

# class vector
y<-zoo[,18]

#run VDA
out<-VDA_R(x, y)

summary(out)
}