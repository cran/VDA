
\name{print.vda.r}
\alias{print.vda.r}
\alias{print.vda}
\title{
Print a vda.r object
}
\description{
The default print method for a vda.r object.
}
\usage{
\method{print}{vda.r}(x, ...)
}
\arguments{
  \item{x}{
Object of class 'vda.r', usually the result of a call to \code{\link{vda.r}}.
}
  \item{\dots}{
Not used. 
}
}
\details{
Prints out the predicted classes for given training data found using Vertex Discriminant Analysis. \code{\link{summary.vda.r}} provides more detailed information about the VDA \code{object}.
}

\author{
Edward Grant, Xia Li, Kenneth Lange, Tong Tong Wu

Maintainer: Edward Grant \email{edward.m.grant@gmail.com}
}


\references{
Lange, K. and Wu, T.T.  (2008) An MM Algorithm for Multicategory Vertex Discriminant Analysis.  Journal of Computational and Graphical Statistics, Volume 17, No 3, 527-544.
}


\seealso{
\code{\link{vda.r}}, \code{\link{summary.vda.r}}
}

\examples{
# load zoo data
# column 1 is name, columns 2:17 are features, column 18 is class
data(zoo)

# feature matrix without intercept
x <- zoo[,2:17]

# class vector
y <- zoo[,18]

#run VDA
out <- vda.r(x, y)

print(out)
}
