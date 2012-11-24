print.VDA_R <-
function (x, ...)
{
  cat ("\n Call: \n")
  print (x$call)
  
  cat ("\n Predicted classification: \n")
  print (x$predicted)
  
  cat ("\n Training Error: \n")
  print (x$training_error_rate)
}
