print.summary.mpr <-
function(x, ...){
   cat("Call:\n")
   print(x$call)
   cat("\nModel:\n")
   print(x$model)
   cat("\nCoefficients:\n")

   if(class(x$overallpmat)=="matrix"){
      printCoefmat(x$coefmat, P.values=TRUE, has.Pvalue=TRUE, signif.legend=FALSE)
      cat("---\nOverall:\n")
      printCoefmat(x$overallpmat, P.values=TRUE, has.Pvalue=TRUE)
   }else{
      printCoefmat(x$coefmat, P.values=TRUE, has.Pvalue=TRUE)
   }
}
