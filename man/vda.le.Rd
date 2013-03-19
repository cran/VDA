\name{vda.le}
\alias{vda.le}
\title{
Multicategory Vertex Discriminant Analysis (VDA) For High-Dimensional Data 
}
\description{
The method of vertex discriminant analysis (VDA) is ideally suited to handle multiple categories and an excess of predictors over training cases. \code{vda.le} is an elaboration of VDA that simultaneously conducts classification of \eqn{k} possible categories and variable selection of \emph{p} features, based on a data set of \emph{n} cases.  Variable selection is imposed using \emph{L1} (Lasso) and group Euclidean penalties. To use only Lasso penalization, set \eqn{\lambda2}=0. To use only Euclidean penalization, set \eqn{\lambda1}=0. 

}
\usage{
  vda.le(x, y, lambda1, lambda2)
}
\arguments{
  \item{x}{
    \emph{n} x \emph{p} matrix or data frame containing the cases for each feature.  The rows correspond to cases and the columns to the features. Intercept column is not included in this. 
  }
  \item{y}{
    \emph{n} x 1 vector representing the outcome variable.  Each element denotes which one of the \emph{k} classes that case belongs to
  }
  \item{lambda1}{
  Tuning parameter to control the lasso penalty.  The default value is \eqn{1/n}.  For determining the optimal \code{lambda1}, refer to \code{\link{cv.vda.le}}.
  }
  \item{lambda2}{
  Tuning parameter to control the Euclidean penalty.  The default value is 0.01.  For determining the optimal \code{lambda1}, refer to \code{\link{cv.vda.le}}
  }
}

\details{
\code{vda.le} carries out cyclic coordinate descent in the context of VDA to minimize the loss function.  By adding lasso (\eqn{L1}-norm) and group Euclidean penalties to the VDA loss function, unnecessary predictors are elimnated, adding parsimony and making the model more interpretable.  Lasso penalties are applied to each predictor coefficient separately, while Euclidean penalties couples the coefficients of a single predictor and penalize the group collectively.  If \code{lambda1=0}, then the overall penalty reduces to only group penalties. When \code{lambda2=0}, then the overall penalty reduces to the lasso.  With these penalties in place, cyclic coordinate descent accelerates estimation of all coefficients.}

\value{
  \item{feature}{
    Feature matrix \code{x} with an intercept vector added as the first column.  All entries in the first column should equal \code{1}.
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
    Vector of tuning constants where the first component is \code{lambda1} and the second is \code{lambda2}
    }
    \item{predicted}{
    Vector of predicted category values based on VDA.
    }
    \item{coefficient}{
    The estimated coefficient matrix where the columns represent the coefficients for each predictor variable corresponding to \code{k-1} outcome categories.  The coefficient matrix is used for classifying new cases.
    }
    \item{training_error_rate}{
    The percentage of instances in the training set where the predicted outcome category is not equal to the case's true category.
    }
    \item{nonzeros}{
    Number of feature coefficients retained in the model.  Is equal to \emph{p} - number of features eliminated by penalization.
    }
    \item{selected}{
    An integer vector which represents the attributes that were selected after penalization.
    }
    \item{call}{
    The matched call.
    }
}
    
    \author{
    Edward Grant, Xia Li,  Kenneth Lange, Tong Tong Wu
    
    Maintainer: Edward Grant \email{edward.m.grant@gmail.com}
    }
    
    
    \references{
    Wu, T.T. and Lange, K. (2010) Multicategory Vertex Discriminant Analysis for High-Dimensional Data. Annals of Applied Statistics, Volume 4, No 4, 1698-1721.
    
    Lange, K. and Wu, T.T. (2008) An MM Algorithm for Multicategory Vertex Discriminant Analysis. Journal of Computational and Graphical Statistics, Volume 17, No 3, 527-544.
    }
    
    \seealso{
For determining the optimal values for \code{lambda1} and \code{lambda2}, see \code{\link{cv.vda.le}}

For VDA without variable selection, see \code{\link{vda.r}}.
}
  
\examples{
#load dataset from package
data(zoo)

#matrix containing all predictor vectors
x <- zoo[,2:17]

#outcome class vector
y <- zoo[,18]

#run VDA, Only Lasso Penalization
out<-vda.le(x,y,lambda2=0)

#run VDA, Only Euclidean Penalization
out<-vda.le(x,y,lambda1=0)

#run VDA, Lasso and Euclidean Penalization
out<-vda.le(x,y)
}
    