\name{summary.vda.le}
\alias{summary.vda.le}
\title{
Summary for a vda.le object
}
\description{
Takes a fitted VDA object produced by \code{\link{vda.le}} and produces various useful summaries from it. 
}
\usage{
\method{summary}{vda.le}(object, ...)
}
\arguments{
  \item{object}{
An object of class 'vda.le', usually the result of a call to \code{\link{vda.le}}.
}
  \item{\dots}{
Not used. 
}
}
\details{
The function prints the number of cases, the number of classes, and the number of features in \code{object}, of class vda.le. It also prints the lambda used in the analysis. Additionally, it prints the coefficients and the resulting predictions made by Vertex Discriminant Analysis on the training set and the following training error.
}

\author{
Edward Grant, Xia Li, Kenneth Lange, Tong Tong Wu

Maintainer: Edward Grant \email{edward.m.grant@gmail.com}
}
\seealso{
\code{\link{vda.le}}, \code{\link{print.vda.le}}
}
\examples{

# load data
data(zoo)

# feature matrix without intercept
x<-zoo[,2:17]

# class vector
y<-zoo[,18]

#run VDA
out<-vda.le(x, y)

summary(out)
}