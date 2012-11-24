\name{cv.VDA_R}
\alias{cv.VDA_R}
\alias{cv.VDA_R.default}
\title{
Choose \eqn{\lambda} using K-fold cross validation
}
\description{
Choose the optimal tuning parameter \eqn{\lambda} for Vertex Discriminant Analyis by using K-fold cross validation.
}
\usage{
cv.VDA_R(x, y, k, lam.vec)
}
\arguments{
  \item{x}{
\emph{n} x \emph{p} matrix or data frame containing the cases for each feature. The rows correspond to cases and the columns to the features. Intercept column is not included in this. 
}
  \item{y}{
\emph{n} x 1 vector representing the outcome variable.  Each element denotes which one of the \emph{k} classes that case belongs to.
}
  \item{k}{
The number of folds to be used in cross-validation.
}
  \item{lam.vec}{
A vector containing the set of all values of \eqn{\lambda}, from which VDA will be conducted. 
}
}

\details{
K-fold cross validation to select optimal \code{lambda} for use in Vertex Disciminant Analysis (\code{VDA_R}). The optimal value is considered the \code{lamda} value that retuns the lowest testing error over the cross validation. If more than one \code{lambda} value give the minumum testing error, the largest \code{lambda} is selected. 

A plot of the cross validation errors can be viewed through \code{\link{plot.cv.VDA_R}}.
}

\value{
  \item{k}{
  The value of \code{K} used for the K-fold cross validation.
  }
  \item{lam.vec}{
  The values of lambda tested.
  }
  \item{mean.error}{
  The mean error corresponding to each lambda across k-folds
  }  
  \item{lam.opt}{
  The determined lambda value among \code{lam.vec} that returns the smallest prediction error. This value is the optimal \code{lambda} value for use in \code{link{VDA_R}}.
  }
  \item{error.cv}{
  The prediction error matrix returned by cross validation method.
  }
}

\author{
Edward Grant, Xia Li, Kenneth Lange, Tong Tong Wu

Maintainer: Edward Grant \email{egrant@umd.edu}
}


\references{
Lange, K. and Wu, T.T. (2008) An MM Algorithm for Multicategory Vertex Discriminant Analysis. Journal of Computational and Graphical Statistics, Volume 17, No 3, 527-544.
}


\seealso{
\code{\link{VDA_R}}

\code{\link{plot.cv.VDA_R}}
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

# searching for the best lambda with 10-fold cross validation
obj <- cv.VDA_R(x, y, 10, lam.vec)

# run VDA
out <- VDA_R(x,y,obj$lam.opt)
}

