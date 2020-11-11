#' Compute survival curve and cumulative hazard from a glmnet model through glmnetsurv
#'
#' Compute the predicted survivor and cumulative hazard function for a penalized Cox proportional hazard model.
#'
#' @aliases glmnetsurvfit
#'
#' @details
#' \code{glmnetsurvfit} and \code{glmnetbasehaz} functions produce survival curves and estimated cumulative hazard, respectively, for a \code{\link[glmnet]{glmnet}} model fitted through \code{\link[glmnetsurv]{glmnetsurv}}. They both return the estimated survival probability and the estimated cumulative hazard, which are both Breslow estimate.
#'
#' The \code{glmnetbasehaz} is an alias for \code{glmnetsurvfit} which simply computed the predicted survival estimates (baseline).
#'
#' If the \code{newdata} argument is missing, the "average" survival or cumulative hazard estimates are produced with the predictor values equal to means of the data set. See \code{\link[survival]{survfit.coxph}} for warning against this. If the \code{newdata} is specified, then the returned object will contain a matrix of both survival and cumulative hazard estimates with each column for each row in the \code{newdata}.
#'
#' @param fit fitted \code{\link[glmnetsurv]{glmnetsurv}} object
#' @param newdata a matrix containing the variables appearing on the right hand side of the model formula of \code{\link[glmnetsurv]{glmnetsurv}} model.
#' @param ... for future implementations
#'
#' @return \code{glmnetsurvfit} and \code{glmnetbasehaz} return S3 objects of class \code{\link[glmnetsurv]{glmnetsurvfit.glmnetsurv}} and \code{\link[glmnetsurv]{glmnetbasehaz.glmnetsurv}}, respectively:
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
#' \code{\link[glmnet]{glmnet}}, \code{\link[glmnetsurv]{plot.glmnetsurvfit}}
#'
#' @rdname glmnetsurvfit.glmnetsurv
#'
#' @examples
#' data(veteran, package="survival")
#' lam <- 0
#' alp <- 1
#' gmodel <- glmnetsurv(Surv(time, status) ~ factor(trt) + karno + diagtime + age + prior
#'		, data = veteran
#'		, alpha = alp
#'		, lambda = lam
#' )
#'
#' # Survival estimate
#' gsurv <- glmnetsurvfit(fit = gmodel)
#' 
#' # Baseline survival estimate
#' gbsurv <- glmnetbasehaz(gmodel, centered = FALSE)
#'
#' @import glmnet
#' @export 

glmnetsurvfit.glmnetsurv <- function(fit, newdata, ...) {
	mfit <- fit$fit
	if(!inherits(mfit, "coxnet"))stop("The object should be a cox model. Use glmnetsurv to fit the model first.")
	s <- fit$s
	if (length(s)>1)stop("Refit the glmnetsurv model with a single lambda (optimal). See ?glmnetsurvcv")
	afit <- glmnetHazard(fit)
	chaz <- afit$chaz
	surv.est <- exp(-chaz)
	if (!missing(newdata)) {
		beta.hat <- as.vector(predict(mfit, type = "coefficients", s = s))
		new_form <- delete.response(fit$terms)
		m <- model.frame(new_form, data = newdata, xlev = fit$xlevels
			, na.action = fit$na.action, drop.unused.levels = TRUE
		)
		newX <- model.matrix(new_form, m, contr = fit$contrasts, xlev = fit$xlevels)
		xnames <- colnames(newX)
		assign <- setNames(attr(newX, "assign"), xnames)[-1]
		xnames <- names(assign)
		newX <- newX[ , xnames, drop=FALSE]
		xmeans <- apply(newX, 2, mean)
		X.centered <- newX - rep(xmeans, each = NROW(newX))
		lp <- as.vector(X.centered %*% beta.hat)
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
#' @rdname glmnetsurvfit.glmnetsurv
#' @import glmnet
#' @export
#'

glmnetbasehaz.glmnetsurv <- function(fit, centered = TRUE){
	sfit <- glmnetsurvfit.glmnetsurv(fit = fit)
	## Expected cummulative hazard rate sum of hazard s.t y(t)<=t
	chaz <- sfit$cumhaz
	surv.est <- exp(-chaz)
	## Compute the cumhaz with the mean of the covariates otherwise set
	## all covariates to 0 (above)
	if (!centered) {
		beta.hat <- as.vector(predict(fit$fit, type = "coefficients", s = fit$s))
		## Centered estimates
		X.mean <- apply(fit$X, 2, mean)
		offset <- as.vector(X.mean %*% beta.hat)
		chaz <- chaz * exp(-offset)
		surv.est <- exp(-chaz)
	}
	out <- list(time = sfit$time, hazard = chaz, surv = surv.est)
	class(out) <- c("glmnetsurvfit", "glmnetbasehaz")
	out
}


#' Predict survival probabilities at various time points
#'
#' The function extracts the survival probability predictions from a \code{glmnet} model.
#'
#' @aliases predictSurvProb
#'
#' @param object fitted \code{\link[glmnetsurv]{glmnetsurv}}.
#' @param newdata a matrix containing the variables appearing in model \code{\link[glmnetsurv]{glmnetsurv}} formula.
#' @param times a vector of times in the range of the response, at which to return the survival probabilities.
#' @param ... for future implementations.
#'
#' @return a matrix of probabilities with as many rows as the rows of the \code{newdata} and as many columns as number of time points (\code{times}). 
#'
#' @examples
#'
#' data(veteran, package="survival")
#' # Penalized
#' lam <- 0.02
#' alp <- 1
#' gfit1 <- glmnetsurv(Surv(time, status) ~ factor(trt) + karno + diagtime + age + prior
#'		, data = veteran
#'		, alpha = alp
#'		, lambda = lam
#' )
#' p1 <- predictSurvProb.glmnetsurv(gfit1, newdata = veteran[1:80,], time = 10)
#'
#' # Unpenalized model
#' lam2 <- 0
#' alp2 <- 1
#' gfit2 <- glmnetsurv(Surv(time, status) ~ factor(trt) + karno + diagtime + age + prior
#'		, data = veteran
#'		, alpha = alp2
#'		, lambda = lam2
#' )
#' p2 <- predictSurvProb.glmnetsurv(gfit2, newdata = veteran[1:80,], times = 10)
#'
#' plot(p1, p2, xlim=c(0,1), ylim=c(0,1)
#' 	, xlab = "Penalized predicted survival chance at 10"
#' 	, ylab="Unpenalized predicted survival chance at 10"
#' )
#'
#' @importFrom prodlim sindex
#' @export

predictSurvProb.glmnetsurv <- function(object, newdata, times, ...){
	N <- NROW(newdata)
	sfit <- glmnetsurvfit(object, newdata)
	S <- t(sfit$surv)
	Time <- sfit$time
	if(N == 1) S <- matrix(S, nrow = 1)
	p <- S[, prodlim::sindex(Time, times), drop = FALSE]
	p
}

#' Extract predictions from glmnet model
#'
#' Extract event probabilities from the fitted model.
#'
#' @aliases predictRisk
#'
#' @details
#' For survival outcome, the function predicts the risk, \eqn{1 - S(t|x)}, where \eqn{S(t|x)} is the survival chance of an individual characterized by \eqn{x}.
#'
#' @param object fitted \code{\link[glmnetsurv]{glmnetsurv}}.
#' @param newdata a matrix containing the variables appearing in model \code{\link[glmnetsurv]{glmnetsurv}} formula.
#' @param times a vector of times in the range of the response, at which to return the survival probabilities.
#' @param ... for future implementations.
#'
#' @return a matrix of probabilities with as many rows as the rows of the \code{newdata} and as many columns as number of time points (\code{times}).
#'
#' @examples
#'
#' data(veteran, package="survival")
#' # Penalized
#' lam <- 0.02
#' alp <- 1
#' gfit1 <- glmnetsurv(Surv(time, status) ~ factor(trt) + karno + diagtime + age + prior
#'		, data = veteran
#'		, alpha = alp
#'		, lambda = lam
#' )
#' r1 <- predictRisk.glmnetsurv(gfit1, newdata = veteran[1:80,], times = 10)
#'
#' # Unpenalized model
#' lam2 <- 0
#' alp2 <- 1
#' gfit2 <- glmnetsurv(Surv(time, status) ~ factor(trt) + karno + diagtime + age + prior
#'		, data = veteran
#'		, alpha = alp2
#'		, lambda = lam2
#' )
#' r2 <- predictRisk.glmnetsurv(gfit2, newdata = veteran[1:80,], times = 10)
#' plot(r1, r2, xlim=c(0,1), ylim=c(0,1)
#' 	, xlab = "Penalized predicted survival chance at 10"
#' 	, ylab="Unpenalized predicted survival chance at 10"
#' )
#'
#' @export

predictRisk.glmnetsurv <- function(object, newdata, times, ...){
	p <- 1 - predictSurvProb.glmnetsurv(object, newdata, times)
	p
}

