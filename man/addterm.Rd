\name{addterm}
\alias{addterm}
\alias{dropterm}

\title{
All Possible Single-Term Additions to / Deletions from a Multi-Parameter Regression (MPR) Model
}
\description{
Identifies all models which arise via single-term additions to / deletions from a
component (or, simultaneously, multiple components) of the supplied \code{mpr} model, fits all such
models and summarises these models in a table.
}
\usage{
addterm(object, upper = ~ ., comp = 1:(object$ncomp),
        aic = TRUE, bestmodel = object, ...)

dropterm(object, lower = ~ 1, comp = 1:(object$ncomp),
         aic = TRUE, bestmodel = object, ...)
}

\arguments{
  \item{object}{
  an object of class \dQuote{\code{mpr}} which is the result of a call to \code{\link{mpr}}.
}
  \item{upper}{
  a one-sided \code{formula} (used in \code{addterm}) specifying a maximal model which must
  include the current one.
}
  \item{lower}{
  a one-sided \code{formula} (used in \code{dropterm}) specifying a minimal model which must
  be within the current one.
}
 \item{comp}{
  a numeric value (or vector) indicating the regression component (or components) where (simultaneous) 
  additions / deletions occur. Note that \dQuote{\code{1}} \eqn{ = \lambda}, \dQuote{\code{2}}
  \eqn{ = \gamma} and \dQuote{\code{3}} \eqn{ = \rho}.  For more information on the various components, 
  see \code{\link{mpr}} and \code{\link{distributions}}.
}
  \item{aic}{
  logical. If \code{TRUE}, AIC is used as the basis for determining the best model among those
  considered. If \code{FALSE}, BIC is used.
}
  \item{bestmodel}{
  an initial best model which, by default, is the supplied current model. This argument is used within the
  \code{\link{stepmpr}} function but is unlikely to be used directly by the end user.
}
  \item{\dots}{
  additional arguments to be passed to internal methods.
}
}
\details{
 The hierarchy is respected when considering terms to be added or dropped, e.g., all main effects
 contained in a second-order interaction must remain.

 When using \code{addterm}, the terms in the \code{upper} model formula must be a \emph{superset}
 of the terms for each regression component indicated by \code{comp}. For example, if component
 1 is \code{~ a + b + c} and component 2 is \code{~ a + b} (and terms are to be added to both
 simultaneously, i.e., \code{comp=1:2}),
 then \code{upper = ~ a + b + c + d} is acceptable and means that the variable
 \code{d} will be added simultaneously to components 1 and 2 (this can be written more compactly as
 \code{upper = ~ . + d}). On the other hand, \code{~ a + b + d} is not acceptable since its terms do
 not form a superset of the terms in component 1 (however, this would be acceptable if we were only
 considering component 2, i.e., if \code{comp=2}).

 When using \code{dropterm}, the terms in the \code{lower} model formula must be a \emph{subset}
 of the terms for each regression component indicated by \code{comp}. Again, if component 1 is
 \code{~ a + b + c} and component 2 is \code{~ a + b} (and terms are to be dropped from both
 simultaneously, i.e., \code{comp=1:2}), then \code{lower = ~ a} is
 acceptable and means that the variable \code{b} will be dropped simultaneously from components 1 and 2 
 (this can be written more compactly as \code{lower = ~ . - b}). On the other hand, \code{~ c} is not
 acceptable since its terms do not form a subset of the terms in component 2 (however, this would be
 acceptable if we were only considering component 1, i.e., if \code{comp=1}).

 To summarise the above two paragraphs, the \code{upper} formula must \emph{contain} each formula
 corresponding to the components under consideration whereas the \code{lower} formula must \emph{be
 contained within} each of these formulae.

 %In the \code{upper} or \code{lower} model formulae . means \dQUote{what is already there}.
 %Specifically, in the \code{upper} formula, it is the \emph{union} of terms from the components
 %selected via \code{comp} whereas, in the \code{lower} formula, it is the \emph{intersection} of such
 %terms.
}
\value{
A \code{list} containing the following components:
  \item{modeltab}{
  a table containing information about each of the fitted models. This information comes from the
  \dQuote{\code{model}} element of each of the \code{mpr} objects - see \code{\link{mpr}} for details.
  }
  \item{bestmodel}{
  the model with the lowest AIC (or BIC if \code{aic = FALSE}) among the fitted models including the
  initial \code{bestmodel} passed to the \code{addterm} / \code{dropterm} function.
  }
}

\author{
Kevin Burke.
}

\seealso{
\code{\link{mpr}}, \code{\link{stepmpr}}, \code{\link{update.mpr}}.
}
\examples{
# Veterans' administration lung cancer data
veteran <- survival::veteran
head(veteran)

# null model
mod1 <- mpr(Surv(time, status) ~ list(~ 1, ~ 1), data=veteran)
mod1 # family = "Weibull" by default

# consider adding trt and celltype to component 1 
addterm(mod1, ~ trt + celltype, comp=1)

# consider adding trt and celltype to components 1 and 2 simultaneously
addterm(mod1, ~ trt + celltype, comp=1:2)$modeltab

# further examples
mod2 <- mpr(Surv(time, status) ~ list(~ trt + celltype, ~ trt + karno),
            data=veteran)
dropterm(mod2, ~ 1, comp=1:2)$modeltab
dropterm(mod2, ~ 1, comp=1)$modeltab
dropterm(mod2, ~ 1, comp=2)$modeltab

# does nothing since celltype is only in component 1
dropterm(mod2, ~ . - celltype, comp=1:2)$modeltab 

# removes celltype from component 1
dropterm(mod2, ~ . - celltype, comp=1)$modeltab
}