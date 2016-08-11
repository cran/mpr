\name{summary.mpr}
\alias{summary.mpr}
\alias{print.summary.mpr}
\alias{overallpval}

\title{
Summarising Multi-Parameter Regression (MPR) Fits
}
\description{
\code{summary} method for class \dQuote{\code{mpr}}
}
\usage{
\method{summary}{mpr}(object, overall = TRUE, ...)
}
\arguments{
  \item{object}{
  an object of class \dQuote{\code{mpr}} which is the result of a call to \code{\link{mpr}}.
}
  \item{overall}{
  logical. If \code{TRUE}, p-values testing the overall effect of a covariate are shown. See
  \dQuote{Details} for more information.
}
  \item{\dots}{
  further arguments passed to or from other methods.
}
}
\details{
The function \code{print.summary.lm} produces a typical table of coefficients, standard errors and
p-values along with \dQuote{significance stars}. In addition, a table of overall p-values are shown.

Multi-Parameter Regression (MPR) models are defined by allowing mutliple distributional parameters to
depend on covariates. The regression components are:
 \deqn{g_1(\lambda) = x^T \beta}{g1(\lambda) = x' \beta}
 \deqn{g_2(\gamma) = z^T \alpha}{g2(\gamma) = z' \alpha}
 \deqn{g_3(\rho) = w^T \tau}{g3(\rho) = w' \tau}
and the table of coefficients displayed by \code{print.summary.lm} follows this ordering.
Furthermore, the names of the coefficients in the table are proceeded by \dQuote{\code{.b}} for
\eqn{\beta} coefficients, \dQuote{\code{.a}} for \eqn{\alpha} coefficients and \dQuote{\code{.t}} for
\eqn{\tau} coefficients to avoid ambiguity.

Let us assume that a covariate \eqn{c}, say, appears in both the \eqn{\lambda} and \eqn{\gamma}
regression components. The standard table of coefficients provides p-values corresponding to the following
null hypotheses:
\deqn{H_0: \beta_c = 0}{H0: \beta_c = 0}
\deqn{H_0: \alpha_c = 0}{H0: \alpha_c = 0} 
where \eqn{\beta_c} and \eqn{\alpha_c} are the regression coefficients of \eqn{c} (one for each of the
two components in which \eqn{c} appears). However, in the context of MPR models, it may be of interest
to test the hypothesis that the \bold{overall} effect of \eqn{c} is zero, i.e., that its \eqn{\beta}
\emph{and} \eqn{\alpha} effects are jointly zero:
\deqn{H_0: \beta_c = \alpha_c = 0}{H0: \beta_c = \alpha_c = 0}
Thus, if \code{overall=TRUE}, \code{print.summary.lm} displays a table of such \dQuote{overall p-values}.

}
\value{
The function \code{summary.mpr} returns a \code{list} containing the following components: 
  \item{call}{
  the matched call from the \code{mpr} object.
  }
  \item{model}{
  a \code{data.frame} containing useful information about the fitted model. This is the same as the
  \dQuote{\code{model}} element of the \code{mpr} object - see \code{\link{mpr}} for details.
  }
  \item{coefmat}{
  a typical coefficient matrix whose columns are the estimated regression coefficients,
  standard errors and p-values.
  }
  \item{overallpmat}{
  a matrix containing the overall p-values as described above in \dQuote{Details}.
  }
}

\author{
Kevin Burke.
}

\seealso{
\code{\link{mpr}}, \code{\link{predict.mpr}}.
}
\examples{
# Veterans' administration lung cancer data
data(veteran, package="survival")
head(veteran)

# Weibull MPR treatment model (family = "Weibull" by default)
mod1 <- mpr(Surv(time, status) ~ list(~ trt, ~ trt), data=veteran)

summary(mod1)
}

