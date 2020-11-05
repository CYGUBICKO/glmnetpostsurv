#' Compute survival curve and cumulative hazard from a glmnet model
#'
#' Compute the predicted survivor and cumulative hazard function for a penalized Cox proportional hazard model.
#'
#' @aliases glmnetsurvfit
#'
#' @details
#' \code{glmnetsurvfit} and \code{glmnetbasehaz} functions produce survival curves and estimated cumulative hazard, respectively, for the fitted \code{\link[glmnet]{glmnet}} model. They both return the estimated survival probability and the estimated cumulative hazard, which are both Breslow estimate.
#'
#' The \code{glmnetbasehaz} is an alias for \code{glmnetsurvfit} which simply computed the predicted survival estimates (baseline).
#'
#' If the \code{newX} argument is missing, the "average" survival or cumulative hazard estimates are produced with the predictor values equal to means of the data set. See \code{\link[survival]{survfit.coxph}} for warning against this. If the \code{newX} is specified, then the returned object will contain a matrix of both survival and cumulative hazard estimates with each column for each row in the \code{newX}.
#'
#' @param fit fitted \code{\link[glmnet]{glmnet}} object
#' @param Y \code{link[survival]{Surv}} object used in the \code{fit}.
#' @param X input matrix used in the \code{fit}.
#' @param newX a data frame containing the variables appearing on the right hand side of \code{\link[glmnet]{glmnet}} formula.
#' @param s value of \code{lambda} at which predictions are required.
#' @param wt option weight applied.
#' @param ... for future implementations
#'
#' @return \code{glmnetsurvfit} and \code{glmnetbasehaz} return S3 objects of class \code{\link[glmnetpostsurv]{glmnetsurvfit.glmnet}} and \code{\link[glmnetpostsurv]{glmnetbasehaz.glmnet}}, respectively:
#' \item{n}{number of observations used in the fit.}
#' \item{events}{total number of events of interest in the fit.}
#' \item{time}{time points defined by the risk set.}
#' \item{n.risk}{the number of individuals at risk at time \code{t}.}
#' \item{n.event}{the number of events that occur at time \code{t}.}
#' \item{n.censor}{the number of subjects who exit the risk set, without an event, at time \code{t}.}
#' \item{surv}{a vector or a matrix of estimated survival function.}
#' \item{cumhaz, hazard}{a vector or a matrix of estimated cumulative hazard.}
#' \item{xmeans}{column means for X.}
#' \item{s}{lambda value.}
#' \item{call}{the call that produced the object.}
#'
#' @seealso
#' \code{\link[glmnet]{glmnet}}, \code{\link[glmnetpostsurv]{plot.glmnetsurvfit}}
#'
#' @rdname glmnetsurvfit.glmnet
#'
#' @examples
#' data(veteran, package="survival")
#' lam <- 0
#' alp <- 1
#' sobj <- with(veteran, Surv(time, status))
#' X <- model.matrix(~factor(trt) + karno + diagtime + age + prior, data = veteran)[,-1]
#' gmodel <- glmnet(X, sobj, 'cox', alpha = alp, lambda = lam)
#'
#' # Survival estimate
#' gsurv <- glmnetsurvfit(fit = gmodel, Y = sobj, X = X, s = lam)
#' 
#' # Baseline survival estimate
#' gbsurv <- glmnetbasehaz(gmodel, centered = FALSE)
#'
#' @import glmnet
#' @export 

glmnetsurvfit.glmnet <- function(fit, Y, X, newX, s, wt = NULL, ...) {
	afit <- glmnetHazard(fit, Y, X, s, wt)
	chaz <- afit$chaz
	surv.est <- exp(-chaz)
	xmeans <- apply(X, 2, mean)
	if (!missing(newX)) {
		beta.hat <- drop(predict(fit, type = "coefficients", s = s))
		X.centered <- X - rep(xmeans, each = NROW(X))
		lp <- unname(drop(X.centered %*% beta.hat))
		surv.est <- t(sapply(surv.est, function(x) x^exp(lp)))
		chaz <- -log(surv.est)
	}
	out <- list(n = afit$n
		, events = sum(afit$n.event)
		, time = afit$time
		, n.risk = afit$n.risk
		, n.event = afit$n.event
		, n.censor = afit$n.censor
		, surv = surv.est
		, cumhaz = chaz
		, xmeans = xmeans
		, s = s
	)
	out$call <- match.call()
	class(out) <- "glmnetsurvfit"
	out
}

#' Compute baseline survival and cumulative hazard
#'
#' @aliases glmnetbasehaz
#'
#' @param centered if \code{TRUE} (default), return data from a predicted survival function at the mean values of the predictors, if \code{FALSE} returns prediction for all predictors equal to zero (baseline hazard).
#'
#' @rdname glmnetsurvfit.glmnet
#' @import glmnet
#' @export
#'

glmnetbasehaz.glmnet <- function(fit, centered = TRUE){
	sfit <- glmnetsurvfit.glmnet(fit)
	## Expected cummulative hazard rate sum of hazard s.t y(t)<=t
	chaz <- sfit$cumhaz
	surv.est <- exp(-chaz)
	## Compute the cumhaz with the mean of the covariates otherwise set
	## all covariates to 0 (above)
	if (!centered) {
		beta.hat <- drop(predict(fit, type = "coefficients", s = sfit$s))
		## Centered estimates
		X.mean <- sfit$xmeans
		offset <- drop(X.mean %*% beta.hat)
		chaz <- chaz * exp(-offset)
		surv.est <- exp(-chaz)
	}
	out <- list(time = sfit$time, hazard = chaz, surv = surv.est)
	class(out) <- c("glmnetsurvfit", "glmnetbasehaz")
	out
}

