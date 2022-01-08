
\name{Distributions}


\alias{Distributions} \alias{distributions} \alias{Distribution} \alias{distribution}
\alias{mprdists}
\alias{family}

\alias{Weibull} \alias{WeibullAFT} \alias{Gompertz} \alias{Loglogistic} \alias{TDL} \alias{Burr}
\alias{PGW}


\title{Distributions in the \code{mpr} Package}

\description{
Information on the distributions currently available within the \code{mpr} package.
}


\details{
When fitting a Multi-Parameter Regression (MPR) model to data, the underlying distribution is selected
using the \dQuote{\code{family}} argument in the \code{\link{mpr}} function.

Currently the \code{mpr} package includes distributions which have upto three parameters:
\describe{
  \item{1}{
  \eqn{\lambda}. This is a \emph{scale} parameter which controls the overall magnitude of the hazard
  function and is typically the \dQuote{interest} parameter in standard Single-Parameter Regression
  (SPR) models. The Multi-Parameter Regression (MPR) framework is more general and considers all
  parameters to be of interest. 
  }
  \item{2}{
  \eqn{\gamma}. This is a \emph{shape} parameter which controls the time evolution of the hazard.
  }
  \item{3}{
  \eqn{\rho}. This is an additional \emph{shape} parameter which controls the time evolution of the
  hazard (available within the Burr and PGW distributions).
  }
}
The MPR framework allows these parameters to depend on covariates as follows:
 \deqn{g_1(\lambda) = x^T \beta}{g1(\lambda) = x' \beta}
 \deqn{g_2(\gamma) = z^T \alpha}{g2(\gamma) = z' \alpha}
 \deqn{g_3(\rho) = w^T \tau}{g3(\rho) = w' \tau}
where \eqn{g_1(.)}{g1(.)}, \eqn{g_2(.)}{g2(.)} and \eqn{g_3(.)}{g3(.)} are appropriate link functions 
(log-link for positive parameters and identity-link for unconstrained parameters), \eqn{x}, \eqn{z} and 
\eqn{w} are covariate vectors, which may or may not contain covariates in common, and \eqn{\beta},
\eqn{\alpha} and \eqn{\tau} are the corresponding vectors of regression coefficients.

The distributions currently available are described below in terms of their hazard functions:

\tabular{lclclcl}{

\code{family} \tab\verb{  }||\verb{  }\tab  \bold{Hazard \eqn{h(t)}} \tab\verb{  }|\verb{  }\tab
\bold{Parameters} \tab\verb{  }|\verb{  }\tab \bold{Note}\cr
\tab\verb{  }||\verb{  }\tab \tab\verb{  }|\verb{  }\tab \tab\verb{  }|\verb{  }\tab \cr

\code{Weibull} \tab\verb{  }||\verb{  }\tab
\eqn{
 \lambda \gamma t^{\gamma - 1} }{
 \lambda \gamma t^(\gamma - 1) } \tab\verb{  }|\verb{  }\tab
\eqn{\lambda > 0}, \eqn{\gamma > 0} \tab\verb{  }|\verb{  }\tab SPR\eqn{(\lambda)} \eqn{=} PH \cr
\tab\verb{  }||\verb{  }\tab \tab\verb{  }|\verb{  }\tab \tab\verb{  }|\verb{  }\tab \cr

\code{WeibullAFT} \tab\verb{  }||\verb{  }\tab
\eqn{
 \lambda \gamma (\lambda t)^{\gamma - 1} }{
 \lambda \gamma (\lambda t)^(\gamma - 1) } \tab\verb{  }|\verb{  }\tab
\eqn{\lambda > 0}, \eqn{\gamma > 0} \tab\verb{  }|\verb{  }\tab  SPR\eqn{(\lambda)} \eqn{=} AFT  \cr
\tab\verb{  }||\verb{  }\tab \tab\verb{  }|\verb{  }\tab \tab\verb{  }|\verb{  }\tab \cr

\code{Gompertz} \tab\verb{  }||\verb{  }\tab
\eqn{
 \lambda \exp(\gamma t) }{
 \lambda exp(\gamma t) } \tab\verb{  }|\verb{  }\tab
\eqn{\lambda > 0}, \eqn{\gamma \in (-\infty, \infty)} \tab\verb{  }|\verb{  }\tab  SPR\eqn{(\lambda)} \eqn{=} PH \cr
\tab\verb{  }||\verb{  }\tab \tab\verb{  }|\verb{  }\tab \tab\verb{  }|\verb{  }\tab \cr

\code{Loglogistic} \tab\verb{  }||\verb{  }\tab
\eqn{
 \frac{\lambda \gamma t^{\gamma - 1}}{1 + \lambda t^\gamma} }{
 \lambda \gamma t^(\gamma - 1) / (1 + \lambda t^\gamma) } \tab\verb{  }|\verb{  }\tab
\eqn{\lambda > 0}, \eqn{\gamma > 0} \tab\verb{  }|\verb{  }\tab  SPR\eqn{(\lambda)} \eqn{=} PO  \cr
\tab\verb{  }||\verb{  }\tab \tab\verb{  }|\verb{  }\tab \tab\verb{  }|\verb{  }\tab \cr

\code{TDL} \tab\verb{  }||\verb{  }\tab
\eqn{
 \frac{\exp(\gamma t + \lambda)}{1+\exp(\gamma t + \lambda)} }{
 exp(\gamma t + \lambda) / (1+exp(\gamma t + \lambda)) } \tab\verb{  }|\verb{  }\tab
\eqn{\lambda \in (-\infty, \infty)}, \eqn{\gamma \in (-\infty, \infty)} \tab\verb{  }|\verb{  }\tab ---  \cr
\tab\verb{  }||\verb{  }\tab \tab\verb{  }|\verb{  }\tab \tab\verb{  }|\verb{  }\tab \cr

\code{Burr} \tab\verb{  }||\verb{  }\tab
\eqn{
 \frac{\lambda \gamma t^{\gamma - 1}}{1 + \lambda \rho t^\gamma} }{
 \lambda \gamma t^(\gamma - 1) / (1 + \lambda \rho t^\gamma) } \tab\verb{  }|\verb{  }\tab
\eqn{\lambda > 0}, \eqn{\gamma > 0}, \eqn{\rho > 0} \tab\verb{  }|\verb{  }\tab  ---  \cr
\tab\verb{  }||\verb{  }\tab \tab\verb{  }|\verb{  }\tab \tab\verb{  }|\verb{  }\tab \cr

\code{PGW} \tab\verb{  }||\verb{  }\tab
\eqn{
 \lambda \gamma \rho t^{\gamma-1} (1+t^\gamma)^{\rho-1} }{
 \lambda \gamma \rho t^(\gamma-1) (1+t^\gamma)^(\rho-1) } \tab\verb{  }|\verb{  }\tab
\eqn{\lambda > 0}, \eqn{\gamma > 0}, \eqn{\rho > 0} \tab\verb{  }|\verb{  }\tab  SPR\eqn{(\lambda)} \eqn{=} PH 
}
The acronymns which appear in the table above are:
\describe{
  \item{SPR\eqn{(\lambda)}}{
  a Single-Parameter Regression (SPR) model where covariates enter through the scale parameter,
  \eqn{\lambda}. For example, in the row corresponding to the \code{Weibull} model, \dQuote{SPR\eqn{(\lambda)}
  \eqn{=} PH} means that the Weibull SPR\eqn{(\lambda)} model is a PH model. Thus, this standard parametric
  PH model is generalised via the Weibull MPR model.
  }
  \item{PH}{
  proportional hazards.
  }
  \item{AFT}{
  accelerated failure time.
  }
  \item{PO}{
  proportional odds.
  }
  \item{TDL}{
  time-dependent logistic.
  }
  \item{PGW}{
  power generalised Weibull.
  }
}
}

\author{
Kevin Burke.
}

\seealso{
\code{\link{mpr}}
}

\examples{
# Veterans' administration lung cancer data
veteran <- survival::veteran
head(veteran)

# Weibull MPR treatment model
mpr(Surv(time, status) ~ list(~ trt, ~ trt), data=veteran, family="Weibull")

# Burr MPR treatment model
mpr(Surv(time, status) ~ list(~ trt, ~ trt, ~ trt), data=veteran,
    family="Burr")

}