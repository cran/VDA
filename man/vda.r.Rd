\name{vda.r}
\alias{vda.r}
\alias{vda.r.default}
\alias{vda}
\title{
Vertex Discriminant Analysis
}
\description{
	Multicategory Vertex Discriminant Analysis (VDA) for classifying an outcome with k possible categories and p features, based on a data set of n cases. The default penalty function is Ridge. Lasso, Euclidean, and a mixture of Lasso and Euclidean are also available. Please refer to \code{\link{vda.le}}
}
\usage{
vda.r(x, y, lambda)
vda(x, y, lambda)
}
\arguments{
  \item{x}{
\emph{n} x \emph{p} matrix or data frame containing the cases for each feature. The rows correspond to cases and the columns to the features. Intercept column is not included in this. 
}
  \item{y}{
\emph{n} x 1 vector representing the outcome variable.  Each element denotes which one of the \emph{k} classes that case belongs to
}
  \item{lambda}{
Tuning constant. The default value is set as \eqn{1/n}.  Can also be found using \code{cv.vda.r}, which uses K-fold cross validation to determine the optimal value.
}
}

\details{
Outcome classification is based on linear discrimination among the vertices of a regular simplex in a \emph{k-1}-dimension Euclidean space, where each vertex represents one of the categories.  Discrimination is phrased as a regression problem involving \eqn{\epsilon-}insensitive residuals and a L2 quadratic ("ridge") penalty on the coefficients of the linear predictors.  The objective function can by minimized by a primal Majorization-Minimization (MM) algorithm that 
\enumerate{
\item relies on quadratic majorization and iteratively re-weighted least squares,
\item is simpler to program than algorithms that pass to the dual of the original optimization problem, and
\item can be accelerated by step doubling.
}
Comparisons on real and simulated data suggest that the MM algorithm for VDA is competitive in statistical accuracy and computational speed with the best currently available algorithms for discriminant analysis, such as linear discriminant analysis (LDA), quadratic discriminant analysis (QDA), \emph{k}-nearest neighbor, one-vs-rest binary support vector machines, multicategory support vector machines, classification and regression tree (CART), and random forest prediction.
}
  
\value{
  \item{feature}{
Feature matrix \code{x} with an intercept vector added as the first column. All entries in the first column should equal 1.
}
  \item{stand.feature}{
  The feature matrix where the all columns are standardized, with the exception of the intercept column which is left unstandardized.
  }
  \item{class}{
  Class vector \code{y}. All elements should be integers between 1 and \code{classes}.
  }
  \item{cases}{
  Number of cases, \emph{n}.
  }
  \item{classes}{
  Number of classes, \emph{k}.
  }
  \item{features}{
  Number of feautres, \emph{p}.
  }
  \item{lambda}{
  Tuning constant \code{lambda} that was used during analysis.
  }
  \item{predicted}{
  Vector of predicted category values based on VDA.
  }
  \item{coefficient}{
  The estimated coefficient matrix where the columns represent the coefficients for each predictor variable corresponding to \code{k-1} outcome categories. The coefficient matrix is used for classifying new cases.
  }
  \item{training_error_rate}{
  The percentage of instances in the training set where the predicted outcome category is not equal to the case's true category.
  }
  \item{call}{
  The matched call
  }
  \item{attr(,"class")}{
  The function results in an object of class "vda.r"
  }
}

\author{
Edward Grant, Xia Li, Kenneth Lange, Tong Tong Wu

Maintainer: Edward Grant \email{edward.m.grant@gmail.com}
}


\references{
Lange, K. and Wu, T.T. (2008) An MM Algorithm for Multicategory Vertex Discriminant Analysis. Journal of Computational and Graphical Statistics, Volume 17, No 3, 527-544.
}

\seealso{
For determining the optimal values for \code{lambda}, refer to \code{\link{cv.vda.r}}.

For high-dimensional setting and conduct variable selection, please refer to \code{\link{vda.le}}.
}

\examples{
# load zoo data
# column 1 is name, columns 2:17 are features, column 18 is class
data(zoo)

#matrix containing all predictor vectors
x <- zoo[,2:17]

#outcome class vector
y <- zoo[,18]

#run VDA
out <- vda.r(x, y)

#Predict five cases based on VDA
fivecases <- matrix(0,5,16)
fivecases[1,] <- c(1,0,0,1,0,0,0,1,1,1,0,0,4,0,1,0)
fivecases[2,] <- c(1,0,0,1,0,0,1,1,1,1,0,0,4,1,0,1)
fivecases[3,] <- c(0,1,1,0,1,0,0,0,1,1,0,0,2,1,1,0)
fivecases[4,] <- c(0,0,1,0,0,1,1,1,1,0,0,1,0,1,0,0)
fivecases[5,] <- c(0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0)
predict(out, fivecases)
}

