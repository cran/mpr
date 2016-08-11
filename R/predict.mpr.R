predict.mpr <-
function(object, newdata, type=c("survivor", "hazard", "percentile"), tvec, prob=0.5, ...){
   family <- match.arg(object$model$family, names(mprdists))
   famlist <- mprdists[family][[1]]
   
   ncomp <- famlist$ncomp

   est <- coef(object)

   beta  <- est$beta
   alpha  <- est$alpha
   tau  <- est$tau

   formula <- object$formula

   rhs <- eval(formula[[3]])

   formb <- rhs[[1]]
   forma <- rhs[[2]]
   formt <- rhs[[3]]

   xvars <- object$xvars
   xlevels <- object$xlevels
   xfac <- names(object$xlevels)
   newnam <- colnames(newdata)

   mvars <- match(xvars, newnam)
   vna <- is.na(mvars)   

   if(any(vna)){
     errmess <- paste("The following variables not found:",
                 paste(xvars[vna], collapse=", ") )
      stop(errmess)
   }

   nums <- match(setdiff(xvars,xfac), newnam)
   facs <- match(xfac, newnam)

   if(length(nums) > 0){
      for(i in 1:length(nums)){
         newdata[,nums[i]] <- as.numeric(newdata[,nums[i]])
      }
   }
  
   if(length(facs) > 0){
      for(i in 1:length(facs)){
         newdata[,facs[i]] <- as.factor(newdata[,facs[i]])
      }
   }

   blevels <- xlevels[match(attr(terms(formb), "term.labels"), xfac)]
   alevels <- xlevels[match(attr(terms(forma), "term.labels"), xfac)]
   tlevels <- xlevels[match(attr(terms(formt), "term.labels"), xfac)]

   if(!is.null(blevels)){ blevels <- blevels[!is.na(names(blevels))] }
   if(!is.null(alevels)){ alevels <- alevels[!is.na(names(alevels))] }
   if(!is.null(tlevels)){ tlevels <- tlevels[!is.na(names(tlevels))] }

   Xb <- model.matrix(terms(formb), data=newdata, xlev=blevels)
   Xa <- model.matrix(terms(forma), data=newdata, xlev=alevels)
   Xt <- model.matrix(terms(formt), data=newdata, xlev=tlevels)

   parmat <- cbind(Xb%*%beta, Xa%*%alpha, Xt%*%tau)

   type <- match.arg(type)

   switch(type,
      survivor = {
         mprsurv <- Vectorize(famlist$surv, vectorize.args="ti")
         out <- mprsurv(parmat, tvec)
      },
      hazard = {
         mprhaz <- Vectorize(famlist$haz, vectorize.args="ti")
         out <- mprhaz(parmat, tvec)
      },
      percentile = {
         out <- as.matrix(famlist$sim(parmat, 1-prob))
      }, )

   out
}
