\name{cv.VDA_LE}
\alias{cv.VDA_LE}
\alias{cv.VDA_LE.default}
\title{
Choose the optimal pair of lambdas, \eqn{\lambda_1} and \eqn{\lambda_2}
}
\description{
Use k-fold validation to choose the optmial values for the tuning parameters \eqn{\lambda_1} and \eqn{\lambda_2} to be used in Multicategory Vertex Discriminant Analysis (\code{VDA_LE}).
}
\usage{
cv.VDA_LE(x, y, kfold, lam.vec.1, lam.vec.2)
}
\arguments{
  \item{x}{
\emph{n} x \emph{p} matrix or data frame containing the cases for each feature. The rows correspond to cases and the columns to the features. Intercept column is not included in this. 
}
  \item{y}{
\emph{n} x 1 vector representing the outcome variable.  Each element denotes which one of the \emph{k} classes that case belongs to.
}
  \item{kfold}{
The number of folds to use for the k-fold validation for each set of \eqn{lambda_1} and \eqn{lambda_2}
}
  \item{lam.vec.1}{
A vector containing the set of all values of \eqn{\lambda_1}, from which VDA will be conducted. 
}
  \item{lam.vec.2}{
A vector containing the set of all values of \eqn{\lambda_2}, from which VDA will be conducted. VDA_LE is relatively insensitive to lambda values, so it is recommended that a vector of few values is used. The default value is 0.01. 
}
}

\details{
For each pair of \eqn{(\lambda_1,\lambda_2)}, k-fold cross-validation will be conducted and the corresponding average testing error over the k folds will be recorded. \eqn{\lambda_1} represents the parameter for the lasso penalization, while \eqn{\lambda_2} represents the parameter for the group euclidean penalization. The optimal pair is considered the pair of values that give the smallest testing error over the cross validation. 
}

\value{
  \item{kfold}{
  	The number of folds used in k-fold cross validation
  }	
    \item{lam.vec.1}{
  	The user supplied vector of \eqn{\lambda}_1 values
  }	
    \item{lam.vec.2}{
  	The user supplied vector of \eqn{\lambda}_2 values
  }	
  \item{error.cv}{
  A matrix of average testing errors. The rows correspond to \eqn{\lambda}_1 values and the columns correspond to \eqn{\lambda}_2 values.
  }
  \item{lam.opt}{
  The pair of \eqn{\lambda}_1 and \eqn{\lambda}_2 values that return the lowest testing error across k-fold cross validation.
  }
}

\author{
Edward Grant, Xia Li,  Kenneth Lange, Tong Tong Wu

Maintainer: Edward Grant \email{egrant@umd.edu}
}


\references{
Wu, T.T. and Lange, K. (2010) Multicategory Vertex Discriminant Analysis for High-Dimensional Data. Annals of Applied Statistics, Volume 4, No 4, 1698-1721.

Lange, K. and Wu, T.T. (2008) An MM Algorithm for Multicategory Vertex Discriminant Analysis. Journal of Computational and Graphical Statistics, Volume 17, No 3, 527-544.
}


\seealso{
\code{\link{VDA_LE}}. 
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

# searching for the best pair
lambda <- cv.VDA_LE(x, y, kfold=3, lam.vec.1, lam.vec.2)$lam.opt

# run VDA
out <- VDA_LE(x,y,lambda[1], lambda[2])
}

