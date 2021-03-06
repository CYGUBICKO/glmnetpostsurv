% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/postsurvplots.R
\name{plot.varImp}
\alias{plot.varImp}
\title{Generic method for plotting variable importance of various models}
\usage{
\method{plot}{varImp}(x, ...)
}
\arguments{
\item{x}{a \code{\link[glmnetsurv]{glmnetsurv}} object.}

\item{...}{for future implementations}
}
\description{
Generic method for plotting variable importance of various models
}
\examples{

data(veteran, package="survival")
# glmnet
gfit1 <- glmnetsurv(Surv(time, status) ~ factor(trt) + karno + diagtime + age + prior
	, data = veteran
	, lambda = 0.02
	, alpha = 0.8
)
imp1 <- varImp(gfit1)
plot(imp1)
imp2 <- varImp(gfit1, show_sign = TRUE, scale = TRUE)
plot(imp2)

}
\seealso{
\code{\link[glmnetsurv]{varImp.glmnetsurv}}
}
