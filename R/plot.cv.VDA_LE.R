plot.cv.VDA_LE <-
function(x, ...)
{
 plot.args1<-list(x$lam.vec.1,x$lam.vec.2,x$error.cv, xlab="Lambda1", ylab="Lambda2", zlab="CV Error", col="red")
 
 plot.args2<-list(x$lam.vec.1,x$lam.vec.2,x$error.cv, xlab="Lambda1", ylab="Lambda2", zlab="CV Error", col="blue")
   
  do.call("plot3d", plot.args1) 
  do.call("surface3d", plot.args2)
}
