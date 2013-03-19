
\name{plot.cv.vda.r}
\alias{plot.cv.vda.r}
\alias{plot.cv.vda}
\title{
Plot a  cv.vda.r object
}
\description{
Plot a the cross validation error across lambda values
}
\usage{
\method{plot}{cv.vda.r}(x, ...)
}
\arguments{
  \item{x}{
Object of class 'cv.vda.r', the result of a call to \code{\link{cv.vda.r}}.
}
  \item{\dots}{
Not used. 
}
}
\details{
Plots the k-fold cross validation testing error for values across a different lambda values. Use \code{\link{cv.vda.r}} to produce the object of class "cv.vda.r."
}

\author{
Edward Grant, Xia Li, Kenneth Lange, Tong Tong Wu

Maintainer: Edward Grant \email{edward.m.grant@gmail.com}
}


\references{
Lange, K. and Wu, T.T.  (2008) An MM Algorithm for Multicategory Vertex Discriminant Analysis.  Journal of Computational and Graphical Statistics, Volume 17, No 3, 527-544.
}


\seealso{
\code{\link{vda.r}}, \code{\link{cv.vda.r}}
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
object <- cv.vda.r(x, y, 10, lam.vec)

#plot
plot(object)
}
