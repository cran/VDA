
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

# run 10 fold cross validation across lambdas
cv <- cv.vda.r(x, y, 10, lam.vec)

# plot CV results
plot(cv)

# Perform VDA with CV-selected optimal lambda
out <- vda.r(x,y,cv$lam.opt)

# Predict five cases based on VDA
fivecases <- matrix(0,5,16)
fivecases[1,] <- c(1,0,0,1,0,0,0,1,1,1,0,0,4,0,1,0)
fivecases[2,] <- c(1,0,0,1,0,0,1,1,1,1,0,0,4,1,0,1)
fivecases[3,] <- c(0,1,1,0,1,0,0,0,1,1,0,0,2,1,1,0)
fivecases[4,] <- c(0,0,1,0,0,1,1,1,1,0,0,1,0,1,0,0)
fivecases[5,] <- c(0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0)
predict(out, fivecases)
}
