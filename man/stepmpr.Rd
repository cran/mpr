\name{stepmpr}
\alias{stepmpr}

\title{
Stepwise Selection Procedure for Multi-Parameter Regression (MPR) Models
}
\description{
Applies a stepwise selection procedure to an object of class \dQuote{\code{mpr}} to find the best model
in the sense of AIC (or BIC).
}
\usage{
stepmpr(object, scope = list(lower = ~ 1, upper = ~ .),
        comp = 1:(object$ncomp), direction = c("both", "backward", "forward"),
        joint = TRUE, jointonly = FALSE, aic = TRUE, trace = 3, ...)
}
\arguments{
  \item{object}{
  an object of class \dQuote{\code{mpr}} which is the result of a call to \code{\link{mpr}}.
}
  \item{scope}{
  either a single formula defining the \code{upper} (maximal) model or a list containing two formulae - 
  the \code{lower} (minimal) and \code{upper} (maximal) models respectively. See \dQuote{Details} for
  further information.
}
  \item{comp}{
  a numeric vector indicating the regression component(s) to which the selection procedure should be
  applied. Note that \dQuote{\code{1}} \eqn{ = \lambda}, \dQuote{\code{2}} \eqn{ = \gamma} and
  \dQuote{\code{3}} \eqn{ = \rho}. For more information on the various components, see \code{\link{mpr}}
  and \code{\link{distributions}}. 
}
  \item{direction}{
  the mode of stepwise search, which can be one of \code{"both"}, \code{"backward"}, or \code{"forward"},
  with a default of \code{"both"}.
}
  \item{joint}{
  logical. If \code{TRUE}, the selection procedure carries out \emph{joint} component (i.e., 
  simultaneous) steps in addition to \emph{individual} component steps. If \code{FALSE}, only individual
  component steps are carried out. See \dQuote{Details} for more information.
}
  \item{jointonly}{
  logical. If \code{TRUE}, the selection procedure only carries out \emph{joint} component steps.  
}
  \item{aic}{
  logical. If \code{TRUE}, AIC is used as the basis for determining the best model among those
  considered. If \code{FALSE}, BIC is used.
}
  \item{trace}{
  if positive, information is printed during the running of \code{stepmpr}. Larger values may give more
  detailed information. 
}
  \item{\dots}{
  additional arguments to be passed to internal methods.
}
}
\details{
The function \code{stepmpr} uses repeated calls to \code{\link{addterm}} and
\code{\link{dropterm}} and is based on the idea that variable selection should be applied to
each component individually \emph{and} to all components jointly (when \code{joint = TRUE}). As an example,
consider the case where forward selection (\code{direction = "forward"}) will be carried out in components
1 and 2 individually (\code{comp = 1:2}) \emph{and} jointly (\code{joint = TRUE}). At a given iteration of
the algorithm, the following single-term additions are then carried out:
\describe{
  \item{individual step \code{1}:}{
  each term currently absent from component 1 will be considered.
}
  \item{individual step \code{2}:}{
  each term currently absent from component 2 will be considered.
}
  \item{joint step \code{1&2}:}{
  each term currently absent from both components 1 and 2 will be considered.
}
}
The reason for the joint step is to account for the possibility that a covariate may only appear significant
when it is present simultaneously in \emph{both} regression components. This situation can arise as the
variance-covariance matrix for the estimated regression coefficients is typically not block diagonal with
respect to the regression components and, in particular, coefficients for the same covariate in
different components are typically highly correlated. Of course, the \code{stepmpr} function has the
flexibility to carry individual steps only, joint steps only or individual steps in a particular component
only as the end user prefers. See \dQuote{Examples} below.


The set of models searched is determined by the \code{scope} argument which is either a single
\code{upper} formula or a list whose elements are \code{lower} and \code{upper} formulae. The
\code{upper} formula must \emph{contain} each formula corresponding to the components under
consideration (as indicated by \code{comp}) whereas the \code{lower} formula must \emph{be contained
within} each of these formulae. For more information on the use of \code{lower} and \code{upper}, see
\code{\link{addterm}}.

If \code{scope} is missing, the \code{lower} model is simply the null model (i.e., a model with no 
covariates) and the \code{upper} model is formed using terms from the initial model; specifically, 
all terms from the regession components under consideration (as indicated by \code{comp}) are used
in the \code{upper} formula.
}
\value{
A final \dQuote{best} \code{mpr} model as selected using the stepwise procedure.
}

\author{
Kevin Burke.
}

\seealso{
\code{\link{mpr}}, \code{\link{dropterm}}, \code{\link{addterm}}, \code{\link{update.mpr}}.
}
\examples{
# Veterans' administration lung cancer data
veteran <- survival::veteran
head(veteran)

#######
mod0 <- mpr(Surv(time, status) ~ 1, data=veteran)
mod0 # family = "Weibull" by default

# the "upper" model formula (by default the lower will be ~ 1)
scope <- ~ trt + celltype

stepmpr(mod0, scope)
stepmpr(mod0, scope, direction="forward", aic=FALSE)

# individual steps only
stepmpr(mod0, scope, joint=FALSE)

# joint steps only
stepmpr(mod0, scope, jointonly=TRUE)

# component 1 only (and, hence, only individual steps)
stepmpr(mod0, scope, comp=1)

#######
mod1 <- mpr(Surv(time, status) ~ trt + celltype, data=veteran)
mod1

stepmpr(mod1)
stepmpr(mod1, scope = ~ .^2)

# "lower" model formula forces trt to stay in
stepmpr(mod1, scope = list(~trt, ~.))

}

