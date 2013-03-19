
\name{print.vda.le}
\alias{print.vda.le}
\title{
Print a vda.le object
}
\description{
The default print method for a vda.le object.
}
\usage{
\method{print}{vda.le}(x, ...)
}
\arguments{
  \item{x}{
Object of class 'vda.le', usually the result of a call to \code{\link{vda.le}}.
}
  \item{\dots}{
Not used. 
}
}
\details{
Prints out the predicted classes for given training data found using Vertex Discriminant Analysis. \code{\link{summary.vda.le}} provides more detailed information about the VDA object \code{x}.
}

\author{
Edward Grant, Xia Li, Kenneth Lange, Tong Tong Wu

Maintainer: Edward Grant \email{edward.m.grant@gmail.com}
}


\references{
Lange, K. and Wu, T.T.  (2008) An MM Algorithm for Multicategory Vertex Discriminant Analysis.  Journal of Computational and Graphical Statistics, Volume 17, No 3, 527-544.
}


\seealso{
\code{\link{vda.le}}, \code{\link{summary.vda.le}}
}

\examples{
# load data
data(zoo)

# feature matrix without intercept
x <- zoo[,2:17]

# class vector
y <- zoo[,18]

#run VDA
out <- vda.le(x, y)

print(out)
}
