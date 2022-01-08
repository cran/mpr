\name{update.mpr}
\alias{update.mpr}
\title{
Update and Re-fit a Multi-Parameter Regression (MPR) Model Call
}
\description{
Updates the right-hand side of the \code{formula} and re-fits the \code{mpr} model.
}
\usage{
\method{update}{mpr}(object, new, comp = 1:(object$ncomp), ...)
}

\arguments{
  \item{object}{
  an object of class \dQuote{\code{mpr}} which is the result of a call to \code{\link{mpr}}.
}
  \item{new}{
  either a one-sided \code{formula} (in which case the \code{comp} argument is also required) or a
  \code{list} of one-sided \code{formula} objects whose length is equal
  to the number of regression components in the \code{mpr} object, e.g., the Weibull model has two components. 
}
  \item{comp}{
  a numeric vector indicating the regression component(s) to be updated (only needed when \code{new} is
  a one-sided \code{formula}) where \dQuote{\code{1}} \eqn{ = \lambda}, \dQuote{\code{2}}
  \eqn{ = \gamma} and \dQuote{\code{3}} \eqn{ = \rho}.  For more information on the various components,
  see \code{\link{mpr}} and \code{\link{distributions}}.
}
  \item{\dots}{
  additional arguments to be passed to the updated \code{mpr} call.
}
}
\details{
There are two ways in whcih the \code{update.mpr} function can be used. The first specificies which
component(s) will be updated
(via the \code{comp} argument) along with the update to be applied via \code{new} which must then be a one-sided
\code{formula}. The second approach specifies both the components in question and the updates to be applied
through \code{new} which is a \code{list} of one-sided \code{formula} objects (in this case \code{comp}
is ignored). See \dQuote{Examples} below.

In the \code{new} formula (or list of formulae) . means \dQuote{what is already there}.
}
\value{
The fitted, updated \code{mpr} object.
}

\author{
Kevin Burke.
}

\seealso{
\code{\link{mpr}}, \code{\link{addterm}}, \code{\link{dropterm}}, \code{\link{stepmpr}}.
}
\examples{
# Veterans' administration lung cancer data
veteran <- survival::veteran
head(veteran)

# Weibull MPR treatment model
mod1 <- mpr(Surv(time, status) ~ list(~ trt, ~ trt), data=veteran,
            family="Weibull")

# remove trt from first component
update(mod1, ~ . - trt, comp=1)
update(mod1, list(~ . - trt, ~ .))

# remove trt from both components
update(mod1, ~ . - trt, comp=1:2)
update(mod1, list(~ . - trt, ~ . - trt))

# add celltype to second component
update(mod1, ~ . + celltype, comp=2)
update(mod1, list(~ . , ~ . + celltype))

# simultaneously remove trt from first component and add celltype to second
# component. This is only possible using the approach where "new" is a list. 
update(mod1, list(~ . - trt, ~ . + celltype))

# can also update other things, e.g. "family"
update(mod1, ~ ., family="Gompertz")
update(mod1, ~ . + celltype, family="Loglogistic")

mod2 <- update(mod1, ~ ., family="Burr") # change to Burr model
mod2
update(mod2, ~ . + celltype, comp=2:3) # add celltype to components 2 and 3
}

