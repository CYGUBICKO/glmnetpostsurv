#' Compute predicted hazard for glmnet 
#' 
#' This code is borrowed from internal function agsurv from survival package. 
#'
#' @param fit fitted \code{\link[glmnet]{glmnet}}.
#' @param Y \code{link[survival]{Surv}} object used in the \code{fit}.
#' @param X input matrix used in the \code{fit}.
#' @param s value of \code{lambda} at which predictions are required.
#' @param wt option weight applied.
#' @return A list of S3 objects. 
#' \item{n}{number of observations used in the fit.}
#' \item{events}{total number of events of interest in the fit.}
#' \item{time}{time points defined by the risk set.}
#' \item{n.risk}{the number of individuals at risk at time \code{t}.}
#' \item{n.event}{the number of events that occur at time \code{t}.}
#' \item{n.censor}{the number of subjects who exit the risk set, without an event, at time \code{t}.}
#' \item{surv}{a vector or a matrix of estimated survival function.}
#' \item{chaz, hazard}{a vector or a matrix of estimated cumulative hazard.}
#' @keywords internal

glmnetHazard <- function(fit, Y, X, s, wt = NULL){
	if (is.null(wt)) wt <- rep(1, NROW(Y))
	beta.hat <- drop(predict(fit, type = "coefficients", s = s))
	xmeans <- apply(X, 2, mean)
	X.centered <- X - rep(xmeans, each = NROW(X))
	oldlp <- unname(drop(X.centered %*% beta.hat))
	oldrisk <- exp(oldlp)
	status <- Y[, ncol(Y)]
	dtime <- Y[, ncol(Y) - 1]
	death <- (status == 1)
	time <- sort(unique(dtime))
	nevent <- as.vector(rowsum(wt * death, dtime))
	ncens <- as.vector(rowsum(wt * (!death), dtime))
	wrisk <- wt * oldrisk
	rcumsum <- function(x) rev(cumsum(rev(x)))
	nrisk <- rcumsum(rowsum(wrisk, dtime))
	irisk <- rcumsum(rowsum(wt, dtime))
	if (NCOL(Y) != 2){
		delta <- min(diff(time))/2
		etime <- c(sort(unique(Y[, 1])), max(Y[, 1]) + delta)
		indx <- approx(etime, 1:length(etime), time, method = "constant", rule = 2, f = 1)$y
		esum <- rcumsum(rowsum(wrisk, Y[, 1]))
		nrisk <- nrisk - c(esum, 0)[indx]
		irisk <- irisk - c(rcumsum(rowsum(wt, Y[, 1])), 0)[indx]
	}
	haz <- nevent/nrisk
	result <- list(n = NROW(Y), time = time, n.event = nevent
		, n.risk = irisk, n.censor = ncens, hazard = haz
		, chaz = cumsum(haz)
	)
	return(result)
}

