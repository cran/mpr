stepmpr <-
function(object, scope=list(lower=~1, upper=~.), comp=1:(object$ncomp),
               direction=c("both", "backward", "forward"),
               joint=TRUE, jointonly=FALSE, aic=TRUE,
               trace=3, ...){

   ncomp <- object$ncomp
   comp <- unique(comp[comp >= 1 & comp <= ncomp])
   
   if(length(comp) == 0){
      comp <- 1:ncomp
   }

   if(class(scope)=="formula"){
     scope <- list(~1, scope)
   }

   lower <- scope[[1]]
   upper <- scope[[2]]

   rhs <- eval(object$formula[[3]])

   lower <- lapply(rhs[comp], update, lower)
   lowerterms <- lapply(lapply(lower, terms), attr, "term.labels")
   lowerterms <- Reduce(intersect, lowerterms)
   lower <- formula(paste("~",paste(c(1,lowerterms), collapse="+")))

   upper <- lapply(rhs[comp], update, upper)
   upperterms <- lapply(lapply(upper, terms), attr, "term.labels")
   upperterms <- Reduce(union, upperterms)
   upper <- formula(paste("~",paste(c(1,upperterms), collapse="+")))


   direction <- match.arg(direction)
   backward <- direction == "both" | direction == "backward"
   forward <- direction == "both" | direction == "forward"

   sortby <- c("aic", "bic")
   sortby <- sortby[c(aic,!aic)]
   prntbrk <- paste(rep("-", 64), collapse="")

   bestmodel <- object
   bestic <- extractIC(bestmodel, aic)
   besticprev <- bestic + 10

   if(trace>0){
      cat(prntbrk,"\n", "Step(0)\n\n", sep="")
      rhs <- eval(bestmodel$formula[[3]])
      cat("Components:\n")
      for(i in 1:ncomp){
         cat(i, .deparse(rhs[[i]]), "\n")
      }
   }
   j <- 1

   while(bestic != besticprev){
      besticprev <- bestic
      bestmodelprev <- bestmodel

      modeltab <- data.frame(comp="*", move="*", term="*", bestmodelprev$model)

      if(backward){
         if(!jointonly){
            for(i in comp){
               d1 <- dropterm(bestmodelprev, lower, i, aic, bestmodel, 
                              hessian=FALSE, ...)
               modeltab <- rbind(modeltab, d1$modeltab)
               bestmodel <- d1$bestmodel
               bestic <- extractIC(bestmodel, aic)
            }
         }
         if((joint | jointonly) & length(comp) > 1){
            d1 <- dropterm(bestmodelprev, lower, comp, aic, bestmodel,
                           hessian=FALSE, ...)
            modeltab <- rbind(modeltab, d1$modeltab)
            bestmodel <- d1$bestmodel
            bestic <- extractIC(bestmodel, aic)
         }
      }

      if(forward){
         if(!jointonly){
            for(i in comp){
               a1 <- addterm(bestmodelprev, upper, i, aic, bestmodel,
                             hessian=FALSE, ...)
               modeltab <- rbind(modeltab, a1$modeltab)
               bestmodel <- a1$bestmodel
               bestic <- extractIC(bestmodel, aic)
            }
         }
         if((joint | jointonly) & length(comp) > 1){
            a1 <- addterm(bestmodelprev, upper, comp, aic, bestmodel,
                          hessian=FALSE, ...)
            modeltab <- rbind(modeltab, a1$modeltab)
            bestmodel <- a1$bestmodel
            bestic <- extractIC(bestmodel, aic)
         }
      }   
   taborder <- order(modeltab[[sortby]])
   modeltab <- modeltab[taborder,]
   row.names(modeltab) <- 1:(dim(modeltab)[1])
   modeltab[c("loglike","aic","bic")] <- round(modeltab[c("loglike","aic","bic")],1)


   if(trace == 2){
      modeltab <- modeltab[1,]
   }

   if(trace>0){
      cat(prntbrk,"\n", "Step(", j, ")\n\n", sep="")
      if(trace >= 2){
         print(modeltab)
         cat("\n")
      }
      rhs <- eval(bestmodel$formula[[3]])
      cat("Components:\n")
      for(i in 1:ncomp){
         cat(i, .deparse(rhs[[i]]), "\n")
      }
   }
   j <- j + 1
   }

   if(trace>0){
      cat(prntbrk,"\n\n")
   }

   bestmodel$call$hessian <- TRUE
   bestmodel <- eval(bestmodel$call)

   bestmodel
}
