extractIC <-
function(object, aic=TRUE){
   if(aic){
      ic <- object$model$aic
   }else{
      ic <- object$model$bic
   }
   ic
}
