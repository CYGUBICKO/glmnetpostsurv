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
	
	y <- fit$y
	x <- model.matrix(fit)
	risk <- predict(fit, type = "risk")
	afit <- glmnetHazard(y = y, x = x, risk = risk)
	chaz <- afit$chaz
	surv.est <- exp(-chaz)
	if (!missing(newdata)) {
		lp <- predict(fit, newdata = newdata, type = "lp")
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
	chaz <- sfit$cumhaz
	surv.est <- exp(-chaz)
	if (!centered) {
		beta.hat <- as.vector(fit$coefficients) 
		## Centered estimates
		X.mean <- fit$means
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
	p <-  cbind(1, S)[, 1 + prodlim::sindex(Time, times),drop = FALSE]
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
#' r1 <- predictRisk(gfit1, newdata = veteran[1:80,], times = 10)
#'
#' # Unpenalized model
#' lam2 <- 0
#' alp2 <- 1
#' gfit2 <- glmnetsurv(Surv(time, status) ~ factor(trt) + karno + diagtime + age + prior
#'		, data = veteran
#'		, alpha = alp2
#'		, lambda = lam2
#' )
#' r2 <- predictRisk(gfit2, newdata = veteran[1:80,], times = 10)
#' plot(r1, r2, xlim=c(0,1), ylim=c(0,1)
#' 	, xlab = "Penalized predicted survival chance at 10"
#' 	, ylab="Unpenalized predicted survival chance at 10"
#' )
#'
#' @importFrom riskRegression predictRisk
#' @export predictRisk
#' @export

predictRisk.glmnetsurv <- function(object, newdata, times, ...){
	p <- 1 - predictSurvProb.glmnetsurv(object, newdata, times)
	p
}

#' Prediction for glmnetsurv model
#'
#' Compute fitted values and model terms for the glmnetsurv model.
#'
#' @details
#' The computation of these predictions similar to those in \code{\link[survival]{predict.coxph}}. Our current implementation does not incorporate stratification.
#'
#' @param object fitted \code{\link[glmnetsurv]{glmnetsurv}} object
#' @param ... for future implementations.
#' @param newdata optional data frame containing the variables appearing on the right hand side of \code{\link[glmnetsurv]{glmnetsurv}} formula. If absent, the predictions are for the data frame used in the original fit.
#' @param type the type of predicted value. Either linear predictor (\code{"lp"}), the risk score (\code{"risk"} equivalently \code{exp(lp)}) and the terms of linear predictor (\code{"terms"}).
#' @param terms if \code{type = "terms"}, this argument can be used to specify which terms to be return. Default is all.
#' @param na.action defines the missing value action for the \code{newdata}. If \code{newdata} is absent, then the behavior of missing is dictated by the \code{na.action} option of the original fit.
#'
#' @return a vector of predictions, depending on the \code{type}.
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
#' predict(gfit1)
#'
#' @export

predict.glmnetsurv <- function(object, ..., newdata = NULL
	, type = c("lp", "risk", "terms"), terms = names(object$assign)
	, na.action = na.pass){
	
	type <- match.arg(type)
	xmeans <- object$means
	beta.hat <- object$coefficients
	beta.hat <- beta.hat[rownames(beta.hat),,drop = TRUE]
	all_terms <- terms
	new_form <- terms(object)
	if (is.null(newdata)){
		newX <- model.matrix(object)
		assign <- object$assign
		xnames <- names(assign)
		newX <- newX[ , xnames, drop=FALSE]
		xmeans <- xmeans[xnames]
		beta.hat <- beta.hat[xnames]
		newX.centered <- newX - rep(xmeans, each = NROW(newX))
		lp <- as.vector(drop(newX.centered %*% beta.hat))
#		lp <- object$linear.predictors
	} else {
		x_form <- delete.response(new_form)
		m <- model.frame(x_form, data = newdata, xlev = object$xlevels
			, na.action = na.action, drop.unused.levels = TRUE
		)
		newX <- model.matrix(object, m, contr = object$contrasts, xlev = object$xlevels)
		xnames <- colnames(newX)
		assign <- setNames(attr(newX, "assign"), xnames)
		xnames <- names(assign)
		newX <- newX[ , xnames, drop=FALSE]
		xmeans <- xmeans[xnames]
		beta.hat <- beta.hat[xnames]
		newX.centered <- newX - rep(xmeans, each = NROW(newX))
		lp <- as.vector(drop(newX.centered %*% beta.hat))
	}
	## Terms
	if (type == "terms"){
		term_list <- list()
		tvals <- unique(assign)
		for (i in seq_along(tvals)){
			w <- assign == tvals[i]
			term_list[[i]] <- newX.centered[, w, drop = FALSE] %*% beta.hat[w]
		}
		terms_df <- do.call("cbind", term_list)
		colnames(terms_df) <- all_terms
		if(!missing(terms)){	terms_df <- terms_df[, terms, drop = FALSE]}
		terms_df <- terms_df[, beta.hat!=0, drop = FALSE]
	}
	out <- switch(type
		, lp = lp
		, risk = exp(lp)
		, terms = terms_df
	)
	return(out)
}

#' Compute the concordance statistic for the pcoxtime model
#'
#' The function computes the agreement between the observed response and the predictor.
#'
#' @aliases concordScore
#'
#' @details
#' Computes Harrel's C index for predictions for \code{\link[glmnetsurv]{glmnetsurv}}, \code{\link[pcoxtime]{pcoxtime}}, \code{\link[survival]{coxph}}, etc, object and takes into account censoring. See \code{\link[survival]{survConcordance}}.
#'
#' @param fit fitted \code{\link[glmnetsurv]{glmnetsurv}}, \code{\link[pcoxtime]{pcoxtime}}, \code{\link[survival]{coxph}}, etc.
#' @param newdata optional data frame containing the variables appearing on the right hand side of \code{\link[glmnetsurv]{glmnetsurv}} formula.
#' @param stats logical. If \code{TRUE} all the related concordance statistics are returned.
#'
#' @return an object containing the concordance, followed by the number of pairs that agree, disagree, are tied, and are not comparable.
#'
#' @examples
#'
#' data(veteran, package="survival")
#' # Penalized
#' lam <- 0.1
#' alp <- 0.5
#' pfit1 <- glmnetsurv(Surv(time, status) ~ factor(trt) + karno + diagtime + age + prior
#'		, data = veteran
#'		, lambda = lam
#'		, alpha = alp
#'	)
#' c1 <- concordScore(pfit1)
#' c1
#'
#' # Unpenalized
#' lam <- 0
#' alp <- 1
#' pfit2 <- glmnetsurv(Surv(time, status) ~ factor(trt) + karno + diagtime + age + prior
#'		, data = veteran
#'		, lambda = lam
#'		, alpha = alp
#'	)
#' c2 <- concordScore(pfit2)
#' c2
#'
#' @export

concordScore <- function(fit, newdata = NULL, stats = FALSE){
	if (is.null(newdata)) {
		risk <- predict(fit, type = "risk")
		y <- fit$y
	} else {
		risk <- predict(fit, newdata = newdata, type = "risk")
		y <- model.extract(model.frame(fit$terms, data = newdata), "response")
	}

	conindex <- survival::survConcordance(y ~ risk)
	if (!stats){
		conindex <- conindex$concordance
	}
	return(conindex)
}

#' Permutation variable importance
#'
#' Computes the relative importance based on random permutation of focal variable for various survival models.
#'
#' @details 
#' Given predictors \code{x_1, x_2, ..., x_n} used to predict the survival outcome, \code{y}. Suppose, for example, \code{x_1} has low predictive power for the response. Then, if we randomly permute the observed values for \code{x_1}, then the prediction for \code{y} will not change much. Conversely, if any of the predictors highly predicts the response, the permutation of that specific predictor will lead to a considerable change in the predictive measure of the model. In this case, we conclude that this predictor is important. In our implementation, Harrel's concordance index is used to measure the prediction accuracy.
#' @param model fitted \code{\link[glmnetsurv]{glmnetsurv}}, \code{\link[pcoxtime]{pcoxtime}}, \code{\link[survival]{coxph}}, etc.
#' @param newdata optional data frame containing the variables appearing on the right hand side of \code{\link[glmnetsurv]{glmnetsurv}} formula.
#' @param nrep number of replicates for permulations
#'
#' @return a named vector of variable scores
#'
#' @examples
#'
#' data(veteran, package="survival")
#' # Penalized
#' lam <- 0.1
#' alp <- 0.5
#' pfit1 <- glmnetsurv(Surv(time, status) ~ factor(trt) + karno + diagtime + age + prior
#'		, data = veteran
#'		, lambda = lam
#'		, alpha = alp
#'	)
#' imp <- permuteImp(pfit1, newdata = veteran, nrep = 50)
#' imp
#'
#' @export

permuteImp <- function(model, newdata, nrep = 50){
	# Overall score
	overall_c <- concordScore(model, newdata = newdata, stats = FALSE)
	Terms <- terms(model)
	yvars <- formula(Terms)[[2]]
	y <- with(newdata, eval(yvars))
	xvars <- all.vars(formula(delete.response(Terms)))
	N <- NROW(newdata)
	newdata <- newdata[, xvars, drop = FALSE]
	vi <- sapply(xvars, function(x){
		permute_df <- newdata[rep(seq(N), nrep), ]
		permute_var <- as.vector(replicate(nrep, sample(newdata[,x], N, replace = FALSE)))
		index <- rep(1:nrep, each = N)
		permute_df[, x] <- permute_var
		risk <- predict(model, newdata = permute_df, type = "risk")
		perm_c <- tapply(risk, index, function(r){
			survConcordance(y~r)$concordance
		})
		mean((overall_c - perm_c)/overall_c)
	})
	return(vi)
}

#' Compute Breslow estimates of survival functions
#'
#' @details
#' This function computes the Breslow's survival functions. It requires a model object with the \code{y}, a \code{Surv} object, estimated model coefficients, \code{coefficients} and the covariate mean values (estimate from model.matrix), \code{means}.
#'
#' @param fit fitted model objects with the objects described in details section.
#' @param centered if \code{TRUE} (default), return data from a predicted survival function at the mean values of the predictors, if \code{FALSE} returns prediction for all predictors equal to zero (baseline hazard).
#' @export
breslow <- function(fit, centered = FALSE){
	beta.hat <- fit$coefficients
	xmeans <- fit$means
	relhaz_bar <- as.vector(exp(drop(xmeans %*% beta.hat)))
	y <- fit$y
	X <- model.matrix(fit)
	lp <- as.vector(drop(X %*% beta.hat))
#	lp <- fit$linear.predictors
	rset <- robustHazard(y = y, lp = lp)
	time <- rset$time
	risk <- rset$risk
	events <- rset$events
	hazard <- events/risk
	time_order <- order(time)
	time <- time[time_order]
	chaz <- unname(cumsum(hazard[time_order]))
	surv.est <- exp(-chaz)
	if (centered){
		surv.est <- surv.est^relhaz_bar
		chaz <- -log(surv.est)
	}
	indexkeep <- sindex(sort(time), unique(sort(time)))
	out <- list(n = rset$n
		, events = sum(rset$n.event)
		, time = unique(sort(time[indexkeep]))
		, n.risk = rset$n.risk
		, n.event = rset$n.event
		, n.censor = rset$n.censor
		, surv = surv.est[indexkeep]
		, cumhaz = chaz[indexkeep]
	)
	out
}
