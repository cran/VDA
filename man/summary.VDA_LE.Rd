\name{summary.VDA_LE}
\alias{summary.VDA_LE}
\title{
Summary for a VDA_LE object
}
\description{
Takes a fitted VDA object produced by \code{\link{VDA_LE}} and produces various useful summaries from it. 
}
\usage{
\method{summary}{VDA_LE}(object, ...)
}
\arguments{
  \item{object}{
An object of class 'VDA_LE', usually the result of a call to \code{\link{VDA_LE}}.
}
  \item{\dots}{
Not used. 
}
}
\details{
The function prints the number of cases, the number of classes, and the number of features in \code{object}, of class VDA_LE. It also prints the lambda used in the analysis. Additionally, it prints the coefficients and the resulting predictions made by Vertex Discriminant Analysis on the training set and the following training error.
}

\author{
Edward Grant, Xia Li, Kenneth Lange, Tong Tong Wu

Maintainer: Edward Grant \email{egrant@umd.edu}
}
\seealso{
\code{\link{VDA_LE}}, \code{\link{print.VDA_LE}}
}
\examples{

# load data
data(zoo)

# feature matrix without intercept
x<-zoo[,2:17]

# class vector
y<-zoo[,18]

#run VDA
out<-VDA_LE(x, y)

summary(out)
}