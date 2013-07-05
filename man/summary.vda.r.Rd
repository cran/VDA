\name{summary.vda.r}
\alias{summary.vda.r}
\alias{summary.vda}
\title{
Summary for a vda.r object
}
\description{
Takes a fitted vda.r object produced by \code{\link{vda.r}} and produces various useful summaries from it. 
}
\usage{
\method{summary}{vda.r}(object, ...)
}
\arguments{
  \item{object}{
An object of class 'vda.r', usually the result of a call to \code{\link{vda.r}}.
}
  \item{\dots}{
Not used. 
}
}
\details{
The function prints the number of cases, the number of classes, and the number of features in \code{object}, of class vda.r. It also prints the lambda used in the analysis. Additionally, it prints the coefficients and the resulting predictions made by Vertex Discriminant Analysis on the training set and the following training error.
}

\author{
Edward Grant, Xia Li, Kenneth Lange, Tong Tong Wu

Maintainer: Edward Grant \email{edward.m.grant@gmail.com}
}
\seealso{
\code{\link{vda.r}}, \code{\link{print.vda.r}}
}
\examples{

# load zoo data
# column 1 is name, columns 2:17 are features, column 18 is class
data(zoo)

# feature matrix without intercept
x<-zoo[,2:17]

# class vector
y<-zoo[,18]

#run VDA
out<-vda.r(x, y)

summary(out)
}