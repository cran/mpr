print.summary.mpr <-
function(x, ...){
   cat("Call:\n")
   print(x$call)
   cat("\nModel:\n")
   print(x$model)
   cat("\nCoefficients:\n")

   if(length(dim(x$overallpmat))==2){
      printCoefmat(x$coefmat, P.values=TRUE, has.Pvalue=TRUE, signif.legend=FALSE)
      cat("---\nOverall:\n")
      printCoefmat(x$overallpmat, P.values=TRUE, has.Pvalue=TRUE, signif.legend=FALSE)
      cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
   }else{
      printCoefmat(x$coefmat, P.values=TRUE, has.Pvalue=TRUE, signif.legend=FALSE)
      cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
   }
}
