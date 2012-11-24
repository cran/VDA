VDA_LE.default <-
function (x, y, lambda1=1/length(y), lambda2=0.01)
{
  if ((!is.numeric(lambda1))|lambda1<=0)
    stop ("lambda1 should be a positive number")
  
  if ((!is.numeric(lambda2))|lambda2<=0)
    stop ("lambda2 should be a positive number")
  
  if (length(y)!=nrow(x))
    stop("Dimention doesn't match! 
         Rows of feature matrix X must be the number of cases")
  
  # initialize input
  cases <- length(y)
  classes <- length(unique(y))
  features <- ncol(x)
  # add intercept col
  feature_i <- as.matrix(cbind(rep(1,nrow(x)),x))
  colnames(feature_i)[1] <- "intercept"
  
  return_data <- .Fortran ("VDA_LE", 
                           stand.feature = as.double (feature_i),
                           as.integer (as.vector (y)),
                           as.integer (cases),
                           as.integer (classes),
                           as.integer (features),
                           as.double (lambda1),
                           as.double (lambda2),
                           post_class = as.integer (rep (0,cases)),
                           coefficient = as.double (matrix (0,classes-1,features+1)),
                           training_error_rate = as.double (0),
                           PACKAGE = "VDA")
     
  coefficient<-matrix (return_data$coefficient,classes-1,features+1)
  
  #count nonzero
  selected<-c()
  nonzeros<-0
  for (j in 2:(features+1))
  {
    if (sum(abs(coefficient[,j]))>1E-6) 
    {
      nonzeros<-nonzeros+1
      selected<-c(selected,j-1)
    }
  }
  
  lambda<-cbind(lambda1,lambda2)
  colnames(lambda)<-c('lambda1','lambda2')
  
  
  out <- list (feature = as.data.frame(feature_i),
               stand.feature = matrix (return_data$stand.feature,cases,features+1),
               classvec = class,
               cases = cases,
               classes = classes,
               features = features,
               lambda = lambda,
               predicted = return_data$post_class,
               coefficient = coefficient,
               training_error_rate = return_data$training_error_rate,
               nonzeros=nonzeros,
               selected=selected,
               call=sys.call ())
  class (out) <- "VDA_LE"
  
  return(out)
}
