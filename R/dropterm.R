dropterm <-
function(object, lower=~1, comp=1:(object$ncomp), 
                              aic=TRUE,  bestmodel=object, ...){

   ncomp <- object$ncomp
   comp <- unique(comp[comp >= 1 & comp <= ncomp])
   
   if(length(comp) == 0){
      comp <- 1:ncomp
   }

   whichcomp <- paste(comp, collapse="&")

   formula <- object$formula
   rhs <- eval(formula[[3]])

   lower <- lapply(rhs[comp], update, lower)
   lowerterms <- lapply(lapply(lower, terms), attr, "term.labels")
   lowerterms <- Reduce(intersect, lowerterms)
   lower <- formula(paste("~",paste(c(1,lowerterms), collapse="+")))

   dropterms <- lapply(rhs[comp], drop.scope, lower)
   dropterms <- Reduce(intersect, dropterms)
   nterms <- length(dropterms) 

   modeltab <- data.frame(matrix(NA, nterms, dim(object$model)[2]+3))
   names(modeltab) <- c("comp","move","term", names(object$model))

   bestic <- extractIC(bestmodel, aic)

   i <- 1
   for(term in dropterms){
      objupdate <- update(object, as.formula(paste("~ . -", term)), comp, ...)
      modeltab[i,] <- data.frame(comp=whichcomp, move="-", term=term, 
                               objupdate$model, stringsAsFactors=FALSE)
      i <- i + 1

      ic <- extractIC(objupdate, aic)
      
      if(ic < bestic){
         bestmodel <- objupdate
         bestic <- ic
      } 
   }
   list(modeltab=modeltab, bestmodel=bestmodel)
}
