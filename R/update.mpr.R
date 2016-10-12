update.mpr <-
function(object, new, comp=1:(object$ncomp), ...){
   ncomp <- object$ncomp
   comp <- unique(comp[comp >= 1 & comp <= ncomp])
   
   if(length(comp) == 0){
      comp <- 1:ncomp
   }

   formula <- object$formula
   rhs <- eval(formula[[3]])

   if(class(new)=="formula"){
      rhs[comp] <- lapply(rhs[comp], update, new)
   }else{
      if(class(new)=="list"){
         j <- 1
         for(i in 1:ncomp){
            rhs[i] <- lapply(rhs[i], update, new[[j]])
            j <- j + 1
         }
      }
   }
   formtemp <- paste(rhs[1:ncomp],collapse=",")
   formtemp <- paste(c(.deparse(formula[[2]]),"~list(",formtemp,")"),collapse="")
   formula <- as.formula(formtemp)  
   
   call <- object$call
   
   call$formula <- formula
   
   extras <- match.call(expand.dots = FALSE)$...

   if(length(extras)){
      existing <- !is.na(match(names(extras), names(call)))
      for(a in names(extras)[existing]){
         call[[a]] <- extras[[a]]
      }
      if(any(!existing)){
         call <- c(as.list(call), extras[!existing])
         call <- as.call(call)
      }
   }
   
   eval(call)
}
