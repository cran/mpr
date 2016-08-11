summary.mpr <-
function(object, overall=TRUE, ...){
   vc <- vcov(object)

   if(!is.null(vc)){
      se <- sqrt(diag(vc))
   }else{
      warning("cannot compute standard errors without variance-covariance matrix")
      se <- NA
      overall <- FALSE
   }

   est <- coef(object)
   est <- c(est[[1]], est[[2]], est[[3]])
   zval <- est/se
   pval <- 2*pnorm(abs(zval), lower.tail=FALSE)
   
   coefmat <- cbind(Est=est, SE=se, Z=zval, Pvalue=pval)
   overallpmat <- NA

   if(overall){
      overallpmat <- overallpval(object)
   }

   out <- list(call=object$call, model=object$model, coefmat=coefmat,
               overallpmat=overallpmat)

   class(out) <- "summary.mpr"
   out
}
