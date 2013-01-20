
\name{plot.cv.VDA_R}
\alias{plot.cv.VDA_R}
\title{
Plot a  cv.VDA_R object
}
\description{
Plot a the cross validation error across lambda values
}
\usage{
\method{plot}{cv.VDA_R}(x, ...)
}
\arguments{
  \item{x}{
Object of class 'cv.VDA_R', the result of a call to \code{\link{cv.VDA_R}}.
}
  \item{\dots}{
Not used. 
}
}
\details{
Plots the k-fold cross validation testing error for values across a different lambda values. Use \code{\link{cv.VDA_R}} to produce the object of class "cv.VDA_R."
}

\author{
Edward Grant, Xia Li, Kenneth Lange, Tong Tong Wu

Maintainer: Edward Grant \email{egrant@umd.edu}
}


\references{
Lange, K. and Wu, T.T.  (2008) An MM Algorithm for Multicategory Vertex Discriminant Analysis.  Journal of Computational and Graphical Statistics, Volume 17, No 3, 527-544.
}


\seealso{
\code{\link{VDA_R}}, \code{\link{cv.VDA_R}}
}

\examples{
# load data
data(zoo)

# feature matrix without intercept
x <- zoo[,2:17]

# class vector
y <- zoo[,18]

# lambda vector
lam.vec <- (1:10)/10

#run 10 fold cross validation across lambdas
object <- cv.VDA_R(x, y, 10, lam.vec)

#plot
plot(object)
}
