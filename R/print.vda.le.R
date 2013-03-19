print.vda.le <-
function (x, ...)
{
  cat ("\n Call: \n")
  print (x$call)
  
  cat ("\n Predicted classification: \n")
  print (x$predicted)
  
  cat ("\n Training Error: \n")
  print (x$training_error_rate)
  
  cat ("\n Number of Active Variables: \n")
  print (x$nonzeros)
  
  cat ("\n Selected Variables with Nonzero Coefficients: \n")
  print (names(x$feature)[x$selected+1])
}
