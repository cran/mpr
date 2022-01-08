\name{extractIC}
\alias{extractIC}
\title{
Extracting Information Criteria for MPR models.
}
\description{
This function extracts the Akaike Information Criterion (AIC) or Bayesian Information Criterion (BIC)
from a fitted \code{\link{mpr}} model.
}
\usage{
extractIC(object, aic = TRUE)
}
\arguments{
  \item{object}{
  fitted \code{\link{mpr}} model.
}
  \item{aic}{
  logical. If \code{TRUE}, AIC is extracted. If \code{FALSE}, BIC is extracted.
}
}
\value{
The AIC or BIC value.
}
\author{
Kevin Burke.
}
\seealso{
\code{\link{mpr}}, \code{\link{addterm}}, \code{\link{dropterm}}.
}
\examples{
# Veterans' administration lung cancer data
veteran <- survival::veteran
head(veteran)

mod1 <- mpr(Surv(time, status) ~ list(~ trt, ~ trt), data=veteran, family="Weibull")
mod1
extractIC(mod1) # AIC
extractIC(mod1, aic=FALSE) # BIC


}

\keyword{internal}

