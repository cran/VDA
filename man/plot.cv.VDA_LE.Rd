
\name{plot.cv.VDA_LE}
\alias{plot.cv.VDA_LE}
\title{
Plot a cv.VDA_LE object
}
\description{
Plot a the cross validation error across lambda values
}
\usage{
\method{plot}{cv.VDA_LE}(x, ...)
}
\arguments{
  \item{x}{
Object of class 'cv.VDA_LE', the result of a call to \code{\link{cv.VDA_LE}}.
}
  \item{\dots}{
Not used. 
}
}

\details{
Plots the k-fold cross validation testing error for values across a different lambda1 and lambda2 values. Use \code{\link{cv.VDA_LE}} to produce the object of class "cv.VDA_LE".
}

\author{
Edward Grant, Xia Li, Kenneth Lange, Tong Tong Wu

Maintainer: Edward Grant \email{egrant@umd.edu}
}


\references{
Lange, K. and Wu, T.T.  (2008) An MM Algorithm for Multicategory Vertex Discriminant Analysis.  Journal of Computational and Graphical Statistics, Volume 17, No 3, 527-544.
}


\seealso{
\code{\link{VDA_LE}}, \code{\link{cv.VDA_LE}}
}

\examples{
# load data
data(zoo)

# feature matrix without intercept
x <- zoo[,2:17]

# class vector
y <- zoo[,18]

# lambda vector
lam.vec.1 <- log(2:4)/100
lam.vec.2 <- (1:3)/100

#run 3 fold cross validation across lambdas
object <- cv.VDA_LE(x, y, 3, lam.vec.1, lam.vec.2)

#plot
plot(object)
}
