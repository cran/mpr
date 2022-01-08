\name{vcov.mpr}
\alias{vcov.mpr}
\title{
Extracting variance-covariance matrix for MPR models.
}
\description{
This function extracts the variance-covariance matrix for the estimated parameters from a fitted
\code{\link{mpr}} model.
}
\usage{
\method{vcov}{mpr}(object, ...)
}
\arguments{
  \item{object}{
  fitted \code{\link{mpr}} model.
}
  \item{\dots}{
  further arguments passed to or from other methods.
}
}
\value{
The variance-covariance matrix.
}
\author{
Kevin Burke.
}
\seealso{
\code{\link{mpr}}.
}
\examples{
# Veterans' administration lung cancer data
veteran <- survival::veteran
head(veteran)

mod1 <- mpr(Surv(time, status) ~ list(~ trt, ~ trt), data=veteran, family="Weibull")
vcov(mod1)
mod1$vcov

}

\keyword{internal}


