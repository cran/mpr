mpr <-
function(formula,data,family="Weibull",init,iterlim=1000,...){

   family <- match.arg(family, names(mprdists))
   famlist <- mprdists[family][[1]]
   ncomp <- famlist$ncomp

   if( !(inherits(formula, "formula") & (length(formula) == 3)) ){
      stop("\"formula\" must be a two sided formula")
   }

   rhs <- formula[[3]]

   if(length(rhs)==1){
      formtemp <- as.formula(paste("~list(~",deparse(rhs),")"))
      rhs <- formtemp[[2]]
   }else{
      if(deparse(rhs[1])!="list()"){
         formtemp <- as.formula(paste("~list(~",deparse(rhs),")"))
         rhs <- formtemp[[2]]
      }
   }

   rhs <- eval(rhs)

   if(!all(unlist(lapply(rhs, inherits, "formula")))){
      stop("RHS of formula must be a list of formula objects")
   }

   if(length(rhs)==1){
      rhs <- rep(rhs, ncomp)
   }

   formb <- rhs[[1]]
   forma <- rhs[[2]]
   formt <- ~0

   if(ncomp == 3 & length(rhs)>=3){
      formt <- rhs[[3]]
   }

   if(ncomp == 3 & length(rhs)==2){
      formt <- ~1
   }

   formtemp <- paste(deparse(formb),deparse(forma),deparse(formt),sep=",")
   formtemp <- paste(c(deparse(formula[[2]]),"~list(",formtemp,")"),collapse="")
   formula <- as.formula(formtemp)

   allterms <- lapply(lapply(list(formb,forma,formt), terms), attr, "term.labels")
   allterms <- Reduce(union, allterms)
  
   if(length(allterms) == 0){
      formall <- "~1"
      formall <- formula(paste(deparse(formula[[2]]), formall))
   }else{
      formall <- paste("~",paste(allterms, collapse="+"), sep="")
      formall <- formula(paste(deparse(formula[[2]]), formall))
   }

   mf <- cl <- match.call()
   m <- match(c("formula", "data"), names(mf), 0)
   mf <- mf[c(1, m)]
   mf$drop.unused.levels <- TRUE
   mf[[1]] <- as.name("model.frame")
   mf$formula <- formall
   mf <- eval(mf, parent.frame())

   resp <- model.response(mf)

   if( !(inherits(resp, "Surv") & attr(resp,"type")=="right") ){
      stop("Response must be a right-censored survival object")
   }

   xlevels <- .getXlevels(terms(formall), mf)
   xvars <- sapply(attr(terms(formall), "variables"), deparse)[-c(1,2)]

   tideltai <- model.extract(mf, "response")[,1:2]

   Xb <- model.matrix(formb,data=mf)
   Xa <- model.matrix(forma,data=mf)
   Xt <- model.matrix(formt,data=mf)

   k.beta <- dim(Xb)[2]
   k.alp <- dim(Xa)[2]
   k.tau <- dim(Xt)[2]
   k <- c(k.beta,k.alp)[1:(ncomp-1)]

   surdat <- cbind(tideltai, Xb, Xa, Xt)

   if( missing(init) ){
      init <- c(0.5, rep(0.01, k.beta - 1),
                -0.2, rep(0.01, k.alp - 1) )
      if(k.tau > 0){
         init <- c(init, -0.2, rep(0.01, k.tau - 1) )
      }
   }else{
      if(any(init == "random")){
         init <- rnorm(k.beta+k.alp+k.tau,0,0.2)
         }
      }

   mprloglike <-  famlist$loglike

   nlmres <- suppressWarnings(.nlm(mprloglike,init,surdat=surdat,k=k,
                                 iterlim=iterlim, ...))
 
   loglike <- -nlmres$min
   npar <- length(init)
   aic <- 2*npar - 2*loglike
   bic <- log(dim(surdat)[1])*npar - 2*loglike                          
   
   code <- nlmres$code
   gradient <- nlmres$gradient
   iter <- nlmres$iterations

   est <- nlmres$estimate
   beta <- est[1:k.beta]
   alpha <- est[(k.beta+1):(k.beta+k.alp)]
   tau <- est[-(1:(k.beta+k.alp))]

   hessian <- nlmres$hessian

   bname <- paste(colnames(Xb),".b",sep="")
   aname <- paste(colnames(Xa),".a",sep="")
   tname <- paste(colnames(Xt),".t",sep="")

   if(length(tau) == 0){
      tname <- character(0)
   }

   names(beta) <- bname
   names(alpha) <- aname
   names(tau) <- tname
   names(gradient) <- c(bname, aname, tname)
   if(!is.null(hessian)){
      rownames(hessian) <- colnames(hessian) <- c(bname, aname, tname)
      vc <- solve(hessian)
   }else{
      vc <- NULL
   }  

   model <- data.frame(family, npar, loglike, aic, bic, code,
                       row.names="", stringsAsFactors=FALSE)

   out <- list(model=model)
   out$coefficients <- list(beta=beta, alpha=alpha, tau=tau)
   out$vcov <- vc
   out$gradient <- gradient
   out$ncomp <- ncomp
   out$formula <- formula
   out$xvars <- xvars
   out$xlevels <- xlevels
   out$call <- cl

   class(out) <- "mpr"
   
   out
}
