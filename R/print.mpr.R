print.mpr <-
function(x, ...){
   cat("Call:\n")
   print(x$call)
   cat("\nModel:\n")
   print(x$model)
   cat("\nCoefficients:\n")
   est <- x$coefficients
   cat("Beta:\n")
   print(est$beta)
   cat("Alpha:\n")
   print(est$alpha)
   if(length(est$tau) > 0){
      cat("Tau:\n")
      print(est$tau)
   }
}
