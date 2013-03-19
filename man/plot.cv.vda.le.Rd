
\name{plot.cv.vda.le}
\alias{plot.cv.vda.le}
\title{
Plot a cv.vda.le object
}
\description{
Plot a the cross validation error across lambda values
}
\usage{
\method{plot}{cv.vda.le}(x, ...)
}
\arguments{
  \item{x}{
Object of class 'cv.vda.le', the result of a call to \code{\link{cv.vda.le}}.
}
  \item{\dots}{
Not used. 
}
}

\details{
3D plots the k-fold cross validation testing error for values across a different lambda1 and lambda2 values. Use \code{\link{cv.vda.le}} to produce the object of class "cv.vda.le".

When \code{lam.vec.1} or \code{lam.vec.2} is set to 0, a 2D plot will be produced. 
}

\author{
Edward Grant, Xia Li, Kenneth Lange, Tong Tong Wu

Maintainer: Edward Grant \email{edward.m.grant@gmail.com}
}


\references{
Lange, K. and Wu, T.T.  (2008) An MM Algorithm for Multicategory Vertex Discriminant Analysis.  Journal of Computational and Graphical Statistics, Volume 17, No 3, 527-544.
}


\seealso{
\code{\link{vda.le}}, \code{\link{cv.vda.le}}
}

\examples{
# load data
data(zoo)

# feature matrix without intercept
x <- zoo[,2:17]

# class vector
y <- zoo[,18]

# searching for the best pair, using both lasso and euclidean penalizations
cv <- cv.vda.le(x, y, kfold=3, lam.vec.1=(1:5)/100, lam.vec.2=(1:5)/100)
plot(cv)
out <- vda.le(x,y,cv$lam.opt[1],cv$lam.opt[2])

# searching for the best pair, using ONLY lasso penalization
cvlasso <- cv.vda.le(x, y, kfold=3, lam.vec.1=exp(1:10)/1000, lam.vec.2=0)
plot(cvlasso)
cvlasso$lam.opt

# searching for the best pair, using ONLY euclidean penalization (remove #)
#cveuclidian <- cv.vda.le(x, y, kfold=3, lam.vec.1=0, lam.vec.2=exp(1:10)/1000)
#plot(cveuclidian)
#cveuclidian$lam.opt

}
