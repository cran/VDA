
\name{print.VDA_R}
\alias{print.VDA_R}
\title{
Print a VDA_R object
}
\description{
The default print method for a VDA_R object.
}
\usage{
\method{print}{VDA_R}(x, ...)
}
\arguments{
  \item{x}{
Object of class 'VDA_R', usually the result of a call to \code{\link{VDA_R}}.
}
  \item{\dots}{
Not used. 
}
}
\details{
Prints out the predicted classes for given training data found using Vertex Discriminant Analysis. \code{\link{summary.VDA_R}} provides more detailed information about the VDA \code{object}.
}

\author{
Edward Grant, Xia Li, Kenneth Lange, Tong Tong Wu

Maintainer: Edward Grant \email{egrant@umd.edu}
}


\references{
Lange, K. and Wu, T.T.  (2008) An MM Algorithm for Multicategory Vertex Discriminant Analysis.  Journal of Computational and Graphical Statistics, Volume 17, No 3, 527-544.
}


\seealso{
\code{\link{VDA_R}}, \code{\link{summary.VDA_R}}
}

\examples{
# load data
data(zoo)

# feature matrix without intercept
x <- zoo[,2:17]

# class vector
y <- zoo[,18]

#run VDA
out <- VDA_R(x, y)

print(out)
}
