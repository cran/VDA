\name{predict.VDA_R}
\alias{predict.VDA_R}
\title{
Predict a VDA_R object.
}
\description{
The predict function for a VDA_R object.
}
\usage{
\method{predict}{VDA_R}(object, newdata=NULL, ...)
}
\arguments{
  \item{object}{
An object of class 'VDA_R', usually the result of a call to \code{\link{VDA_R}}.
}
  \item{newdata}{
An optional \emph{n} x \emph{p} matrix or data frame containing new data to be classified using VDA. The data must contain the same number of attributes as the training data. If \code{newdata} is omitted, the training data is used.
}
  \item{\dots}{
Not used. 
}
}
\details{
The prediction funtion for Vertex Discriminant Analysis (\code{\link{VDA_R}}). Returns 1 x \emph{n} vector in which each element represents the predicted value for the corresponding case.
}

\author{
Edward Grant, Xia Li, Kenneth Lange, Tong Tong Wu

Maintainer: Edward Grant \email{egrant@umd.edu}
}


\references{
Lange, K. and Wu, T.T. (2008) An MM Algorithm for Multicategory Vertex Discriminant Analysis. Journal of Computational and Graphical Statistics, Volume 17, No 3, 527-544.
}


\seealso{
\code{\link{VDA_R}}, \code{\link{summary.VDA_R}}, \code{\link{print.VDA_R}}
}

\examples{
# load data
data(zoo)

# feature matrix without intercept
x <- zoo[,2:17]

# class vector
y <- zoo[,18]

#run VDA
out <- VDA_R(x,y)

onecase <- matrix(c(0,1,1,0,0,0,1,0,1,1,0,0,2,1,0,1),1,16)
predict(out, onecase)

fivecases <- matrix(0,5,16)
fivecases[1,] <- c(0,0,1,0,0,1,0,1,1,0,0,1,0,1,0,0)
fivecases[2,] <- c(0,0,1,0,1,0,1,0,0,1,0,0,6,0,0,0)
fivecases[3,] <- c(0,0,1,0,0,1,1,0,0,0,0,0,4,0,0,0)
fivecases[4,] <- c(1,0,0,1,0,0,0,1,1,1,0,0,4,1,1,1)
fivecases[5,] <- c(0,0,0,0,0,1,1,1,1,0,1,0,0,1,0,0)
predict(out, fivecases)
}

