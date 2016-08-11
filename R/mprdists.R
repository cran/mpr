
mprdists <- list(

   Weibull = list(
      ncomp = 2,
      sim = function(parmat, u){
         lam <- exp(parmat[,1])
         gam <- exp(parmat[,2])
         (-log(u)/lam)^(1/gam)
      },
      haz = function(parmat, ti){
         lam <- exp(parmat[,1])
         gam <- exp(parmat[,2])
         lam*gam*(ti^(gam-1))
      },
      surv = function(parmat, ti){
         lam <- exp(parmat[,1])
         gam <- exp(parmat[,2])
         exp(-lam*(ti^gam))
      },
      loglike = function(param, surdat, k){
         p <- length(param)
         q <- dim(surdat)[2]
         beta <- param[1:k]
         alpha <- param[(k+1):p]

         ti <- surdat[,1]       
         deltai <- surdat[,2]

         Xb <- as.matrix(surdat[,3:(k+2)])
         Xa <- as.matrix(surdat[,(k+3):q]) 

         xbi <- Xb%*%beta  
         xai <- Xa%*%alpha

         lam <- exp(xbi)
         gam <- exp(xai)

         dldbet <- t( deltai-lam*(ti^gam) )%*%Xb
         dldalp <- t( deltai*(1+gam*log(ti)) - lam*gam*(ti^gam)*log(ti) )%*%Xa
   
         loglike <- -sum(deltai*log(gam*lam*(ti^(gam-1)))-lam*(ti^gam))
         attr(loglike, "gradient") <- -c(dldbet, dldalp)

         loglike
      }
   ),

#######

   WeibullAFT = list(
      ncomp = 2,
      sim = function(parmat, u){
         lam <- exp(parmat[,1])
         gam <- exp(parmat[,2])
         lam <- (lam^gam)
         (-log(u)/lam)^(1/gam)
      },
      haz = function(parmat, ti){
         lam <- exp(parmat[,1])
         gam <- exp(parmat[,2])
         lam <- (lam^gam)
         lam*gam*(ti^(gam-1))
      },
      surv = function(parmat, ti){
         lam <- exp(parmat[,1])
         gam <- exp(parmat[,2])
         lam <- (lam^gam)
         exp(-lam*(ti^gam))
      },
      loglike = function(param, surdat, k){
         p <- length(param)
         q <- dim(surdat)[2]
         beta <- param[1:k]
         alpha <- param[(k+1):p]

         ti <- surdat[,1]       
         deltai <- surdat[,2]

         Xb <- as.matrix(surdat[,3:(k+2)])
         Xa <- as.matrix(surdat[,(k+3):q]) 

         xbi <- Xb%*%beta  
         xai <- Xa%*%alpha

         lam <- exp(xbi)
         gam <- exp(xai)

         lam0 <- lam
         lam <- (lam0^(gam))

         dldbet <- t(gam*(deltai-lam*(ti^gam)))%*%Xb
         dldalp <- t(deltai*(gam*log(lam0*ti)+1)-gam*lam*(ti^gam)*log(lam0*ti))%*%Xa
   
         loglike <- -sum(deltai*log(gam*lam*(ti^(gam-1)))-lam*(ti^gam))
         attr(loglike, "gradient") <- -c(dldbet, dldalp)

         loglike
      }
   ),

#######

   Gompertz = list(
      ncomp = 2,
      sim = function(parmat, u){
         lam <- exp(parmat[,1])
         gam <- parmat[,2]
         ifelse(gam == 0, -log(u)/lam, (1/gam)*log(-(gam/lam)*log(u)+1) )
      },
      haz = function(parmat, ti){
         lam <- exp(parmat[,1])
         gam <- parmat[,2] 
         lam*exp(gam*ti)
      },
      surv = function(parmat, ti){
         lam <- exp(parmat[,1])
         gam <- parmat[,2]
         ifelse(gam == 0, exp(-lam*ti), exp(-(lam/gam)*(exp(gam*ti)-1)) )
      },
      loglike = function(param,surdat,k){
         p <- length(param)
         q <- dim(surdat)[2]
         beta <- param[1:k]
         alpha <- param[(k+1):p]

         ti <- surdat[,1]       
         deltai <- surdat[,2]

         Xb <- as.matrix(surdat[,3:(k+2)])
         Xa <- as.matrix(surdat[,(k+3):q]) 

         xbi <- Xb%*%beta  
         xai <- Xa%*%alpha

         lam <- exp(xbi)
         gam <- xai

         dldbet <- t( deltai - (lam/gam)*(exp(gam*ti)-1) )%*%Xb   
         dldalp <- t( deltai*ti + (lam/(gam^2))*(exp(gam*ti)-1) - (lam/gam)*exp(gam*ti)*ti)%*%Xa
   
         loglike <- -sum(deltai*(log(lam*exp(gam*ti))) -(lam/gam)*(exp(gam*ti)-1) )
         attr(loglike, "gradient") <- -c(dldbet,dldalp)

         loglike
      }
   ),

#######

   Loglogistic = list(
      ncomp = 2,
      sim = function(parmat, u){
         lam <- exp(parmat[,1])
         gam <- exp(parmat[,2])
         (((1/u)-1)/lam)^(1/gam)
      },
      haz = function(parmat, ti){
         lam <- exp(parmat[,1])
         gam <- exp(parmat[,2])
         (lam*gam*(ti^(gam-1)))/(1+lam*(ti^gam))
      },
      surv = function(parmat, ti){
         lam <- exp(parmat[,1])
         gam <- exp(parmat[,2])
         1/(1+lam*(ti^gam))
      },
      loglike = function(param,surdat, k){
         p <- length(param)
         q <- dim(surdat)[2]
         beta <- param[1:k]
         alpha <- param[(k+1):p]

         ti <- surdat[,1]       
         deltai <- surdat[,2]

         Xb <- as.matrix(surdat[,3:(k+2)])
         Xa <- as.matrix(surdat[,(k+3):q]) 

         xbi <- Xb%*%beta  
         xai <- Xa%*%alpha

         lam <- exp(xbi)
         gam <- exp(xai)

         dldbet <- t( deltai - (deltai+1)*(lam*(ti^gam))/(1+lam*(ti^gam)) )%*%Xb
         dldalp <- t( deltai*(1+gam*log(ti)) - (deltai+1)*(lam*gam*(ti^gam)*log(ti))/(1+lam*(ti^gam)) )%*%Xa
   
         loglike <- -sum( deltai*log((lam*gam*(ti^(gam-1)))/(1+lam*(ti^gam))) - log(1+lam*(ti^gam)) )
         attr(loglike, "gradient") <- -c(dldbet,dldalp)

         loglike
      }
   ),

#######

   TDL = list(
      ncomp = 2,
      sim = function(parmat, u){
         lam <- parmat[,1]
         gam <- parmat[,2]
         ifelse(gam == 0, -log(u)/(exp(lam)/(1+exp(lam))),
                         (1/gam)*(log((u^(-gam))*(1+exp(lam))-1)-lam) )
      },
      haz = function(parmat, ti){
         lam <- parmat[,1]
         gam <- parmat[,2]
         exp(gam*ti+lam)/(1+exp(gam*ti+lam))
      },
      surv = function(parmat, ti){
         lam <- parmat[,1]
         gam <- parmat[,2]
         ifelse(gam == 0, exp(-(exp(lam)/(1+exp(lam)))*ti),
                           ((1+exp(gam*ti+lam))/(1+exp(lam)))^(-1/gam) )
      },
      loglike = function(param,surdat,k){
         p <- length(param)
         q <- dim(surdat)[2]
         beta <- param[1:k]
         alpha <- param[(k+1):p]

         ti <- surdat[,1]       
         deltai <- surdat[,2]

         Xb <- as.matrix(surdat[,3:(k+2)])
         Xa <- as.matrix(surdat[,(k+3):q]) 

         xbi <- Xb%*%beta  
         xai <- Xa%*%alpha 

         lam <- xbi
         gam <- xai

         wi <- exp(lam)
         zi <- exp(gam*ti+lam)

         dldbet <- t( deltai - (deltai+(1/gam))*(zi/(1+zi)) + (1/gam)*(wi/(1+wi)) )%*%Xb   
         dldalp <- t( deltai*ti - (deltai+(1/gam))*(ti*zi/(1+zi)) + (1/(gam^2))*log((1+zi)/(1+wi)) )%*%Xa
   
         loglike <- -sum(deltai*log(zi/(1+zi))-(1/gam)*log((1+zi)/(1+wi)))
         attr(loglike, "gradient") <- -c(dldbet, dldalp)

         loglike
      }
   ),

#######

   Burr = list(
      ncomp = 3,
      sim = function(parmat, u){
         lam <- exp(parmat[,1])
         gam <- exp(parmat[,2])
         rho <- exp(parmat[,3])
         (((u^(-rho))-1)/(lam*rho))^(1/gam)
      },
      haz = function(parmat, ti){
         lam <- exp(parmat[,1])
         gam <- exp(parmat[,2])
         rho <- exp(parmat[,3])
         (lam*gam*(ti^(gam-1)))/(1+lam*rho*(ti^gam))
      },
      surv = function(parmat, ti){
         lam <- exp(parmat[,1])
         gam <- exp(parmat[,2])
         rho <- exp(parmat[,3])
         (1+lam*rho*(ti^gam))^(-1/rho)
      },
      loglike = function(param,surdat,k){
         p <- length(param)  
         q <- dim(surdat)[2]
         k1 <- k[1]
         k2 <- k[2]
         beta <- param[1:k1]
         alpha <- param[(k1+1):(k1+k2)]
         tau <- param[(k1+k2+1):p]

         ti <- surdat[,1]       
         deltai <- surdat[,2]

         Xb <- as.matrix(surdat[,3:(k1+2)]) 
         Xa <- as.matrix(surdat[,(k1+3):(k1+k2+2)])
         Xt <- as.matrix(surdat[,(k1+k2+3):q])

         xbi <- Xb%*%beta
         xai <- Xa%*%alpha
         xti <- Xt%*%tau  

         lam <- exp(xbi)
         gam <- exp(xai)
         rho <- exp(xti)

         dldbet <- t( deltai - (deltai+1/rho)*(lam*rho*(ti^gam))/(1+lam*rho*(ti^gam)) )%*%Xb
         dldalp <- t( deltai*(1+gam*log(ti)) - (deltai+1/rho)*(lam*gam*rho*(ti^gam)*log(ti))/(1+lam*rho*(ti^gam)) )%*%Xa
         dldtau <- t( (1/rho)*log(1+lam*rho*(ti^gam)) - (deltai+1/rho)*(lam*rho*(ti^gam))/(1+lam*rho*(ti^gam)) )%*%Xt
   
         loglike <- -sum( deltai*log((lam*gam*(ti^(gam-1)))/(1+lam*rho*(ti^gam))) - (1/rho)*log(1+lam*rho*(ti^gam)) )
         attr(loglike, "gradient") <- -c(dldbet, dldalp, dldtau)

         loglike
      }
   ),

#######

   PGW = list(
      ncomp = 3,
      sim = function(parmat, u){
         lam <- exp(parmat[,1])
         gam <- exp(parmat[,2])
         rho <- exp(parmat[,3])
         (((1-log(u))/lam)^(1/rho)-1)^(1/gam)
      },
      haz = function(parmat, ti){
         lam <- exp(parmat[,1])
         gam <- exp(parmat[,2])
         rho <- exp(parmat[,3])
         lam*gam*rho*(ti^(gam-1))*((1+ti^gam)^(rho-1))
      },
      surv = function(parmat, ti){
         lam <- exp(parmat[,1])
         gam <- exp(parmat[,2])
         rho <- exp(parmat[,3])
         exp( -lam*(((1+ti^gam)^rho)-1) )
      },
      loglike = function(param, surdat, k){
         p <- length(param)  
         q <- dim(surdat)[2]
         k1 <- k[1]
         k2 <- k[2]
         beta <- param[1:k1]
         alpha <- param[(k1+1):(k1+k2)]
         tau <- param[(k1+k2+1):p]

         ti <- surdat[,1]       
         deltai <- surdat[,2]

         Xb <- as.matrix(surdat[,3:(k1+2)]) 
         Xa <- as.matrix(surdat[,(k1+3):(k1+k2+2)])
         Xt <- as.matrix(surdat[,(k1+k2+3):q])

         xbi <- Xb%*%beta
         xai <- Xa%*%alpha
         xti <- Xt%*%tau  

         lam <- exp(xbi)
         gam <- exp(xai)
         rho <- exp(xti)

         wi <- 1 + ti^gam

         dldbet <- t( deltai - lam*(wi^rho - 1) )%*%Xb
         dldalp <- t( deltai*(1 + gam*log(ti) + gam*(rho-1)*(ti^gam)*(log(ti))/wi) - lam*gam*rho*(wi^(rho-1))*(ti^gam)*log(ti) )%*%Xa
         dldtau <- t( deltai*(1 + rho*log(wi)) - lam*rho*(wi^rho)*log(wi) )%*%Xt
   
         loglike <- -sum( deltai*log(lam*gam*rho*(ti^(gam-1))*(wi^(rho-1)) ) - lam*((wi^rho)-1) )
         attr(loglike, "gradient") <- -c(dldbet, dldalp, dldtau)

         loglike
      }
   )
)
