predict.vda.le <-
function(object, newdata = NULL, ...)
{
  
  if (!inherits(object, "vda.le")) 
    stop("predict.vda.le can only be used to predict from vda.le objects")
  
  if (missing(newdata))
    return(object$predicted)
  else {
    newdata <- as.matrix(newdata) 
    pred <- cbind(rep(1,nrow(newdata)),newdata)%*%t(object$coefficient)
    k <- object$classes
    c <- -(1+sqrt(k))/(k-1)^(3/2)
    d <- sqrt(k/(k-1));
    vertex <- matrix(rep(0,k*(k-1)),nrow=k-1)
    vertex[,1] <- 1/sqrt(k-1)
    for (kk in 2:k){
      vertex[,kk] <- c
      vertex[kk-1,kk] <- c+d;
    }
    #y2 <- diag(k)%*%vertex;
    
    norm <- function(x)sqrt(sum(x^2))
    ddd <- numeric();
    for (kk in 1:k){
      ddd <- cbind(ddd,apply(scale(pred,vertex[,kk],T),1,norm))
    }
    class.pred <- apply(ddd,1,which.min)
    return(class.pred)
  }
}
