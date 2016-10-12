.deparse <- function(x){
   paste(trimws(deparse(x, width.cutoff=500)),collapse="")
}

.nlm <-
function(..., hessian=TRUE){
   nlm(..., hessian=hessian)
}
