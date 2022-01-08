\name{predict.mpr}
\alias{predict.mpr}
\title{
Predict method for Multi-Parameter Regression (MPR) Fits
}
\description{
Survival predictions based on \code{mpr} objects.
}
\usage{
\method{predict}{mpr}(object, newdata, type = c("survivor", "hazard", "percentile"),
        tvec, prob = 0.5, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  an object of class \dQuote{\code{mpr}} which is the result of a call to \code{\link{mpr}}.
}
  \item{newdata}{
  \code{data.frame} in which to look for variables with which to predict. 
}
  \item{type}{
  type of prediction which may be a survivor function, hazard function or percentile value.
}
  \item{tvec}{
  vector of times at which the predicted survivor or hazard function will be evaluated. Only required
  if \code{type} is \code{"survivor"} or \code{"hazard"}.
}
  \item{prob}{
  numeric value between 0 and 1 (i.e., probability) indicating the percentile to be predicted. By
  default \code{prob = 0.5} which corresponds to the median value. Only required if \code{type} is 
  \code{"percentile"}.
}
  \item{\dots}{
  further arguments passed to or from other methods.
}
}

\value{
A matrix of predictions whose rows correspond to the rows of \code{newdata}. When \code{type}
is \code{"survivor"} or \code{"hazard"}, this matrix of predictions has columns
corresponding to \code{tvec}. However, when \code{type} is \code{"percentile"}, the matrix only
has one column.
}

\author{
Kevin Burke.
}

\seealso{
\code{\link{mpr}}, \code{\link{summary.mpr}}
}
\examples{
library(survival)

# Veterans' administration lung cancer data
veteran <- survival::veteran
head(veteran)

# Weibull MPR treatment model
mod1 <- mpr(Surv(time, status) ~ list(~ trt, ~ trt), data=veteran,
            family="Weibull")

# predicted survivor function evaluated at four times
predict(mod1, newdata=data.frame(trt=c(1,2)), type="survivor",
        tvec=c(25, 50, 100, 150))

# predicted percentiles
predict(mod1, newdata=data.frame(trt=c(1,2)), type="percentile", prob=0.5)
predict(mod1, newdata=data.frame(trt=c(1,2)), type="percentile", prob=0.1)

# comparing predicted survivor functions to Kaplan-Meier curves
KM <- survfit(Surv(time, status) ~ trt, data=veteran)
plot(KM, col=1:2)
tvec <- seq(0, max(KM$time), length=100)
Stpred <- predict(mod1, newdata=data.frame(trt=c(1,2)), type="survivor",
                  tvec=tvec)
lines(tvec, Stpred[1,])
lines(tvec, Stpred[2,], col=2)
}