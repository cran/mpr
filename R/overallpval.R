overallpval <-
function(object){
   est <- coef(object)
   beta <- est$beta; alpha <- est$alpha; tau <- est$tau

   bnam <- names(beta);   bnam <- substr(bnam, 1, nchar(bnam)-2)
   anam <- names(alpha);  anam <- substr(anam, 1, nchar(anam)-2)
   tnam <- names(tau);    tnam <- substr(tnam, 1, nchar(tnam)-2)
   onam <- Reduce(union, list(bnam,anam,tnam))


   onam <- onam[onam!="(Intercept)"]

   vc <- vcov(object)
   vcnam <- rownames(vc)
   vcnam <- substr(vcnam, 1, nchar(vcnam)-2)

   no <- length(onam)

   chisq <- df <- bat <- rep(NA, length=no)

   coefmat <- NULL

   if(no > 0){
      for(i in 1:no){
         onami <- onam[i]

         vcpos <- vcnam==onami
         vci <- vc[vcpos, vcpos]

         bpos <- bnam==onami
         apos <- anam==onami
         tpos <- tnam==onami

         esti <- rbind(beta[bpos], alpha[apos], tau[tpos])

         chisq[i] <- t(esti) %*% solve(vci) %*% esti
         df[i] <- sum(bpos,apos,tpos)
         bat[i] <- paste(c("b","a","t")[c(any(bpos),any(apos),any(tpos))],collapse="")
      }

      pval <- pchisq(chisq, df, lower.tail=FALSE)
      coefmat <- cbind(Chisq=chisq, Df=df, Pvalue=pval)
      rownames(coefmat) <- paste(onam,bat,sep=".")
   }
   coefmat
}
