#' Compute predicted hazard for glmnet 
#' 
#' This code is borrowed from internal function agsurv from survival package. 
#'
#' @param object fitted \code{\link[glmnetsurv]{glmnetsurv}} object.
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

glmnetHazard <- function(y, x = NULL, wt = rep(1, NROW(y)), risk=NULL, survtype=NULL, vartype=NULL){
	status <- y[, ncol(y)]
	dtime <- y[, ncol(y) - 1]
	death <- (status == 1)
	time <- sort(unique(dtime))
	nevent <- as.vector(rowsum(wt * death, dtime))
	ncens <- as.vector(rowsum(wt * (!death), dtime))
	wrisk <- wt * risk
	rcumsum <- function(x) rev(cumsum(rev(x)))
	nrisk <- rcumsum(rowsum(wrisk, dtime))
	irisk <- rcumsum(rowsum(wt, dtime))
	if (NCOL(y) != 2){
		delta <- min(diff(time))/2
		etime <- c(sort(unique(y[, 1])), max(y[, 1]) + delta)
		indx <- approx(etime, 1:length(etime), time, method = "constant", rule = 2, f = 1)$y
		esum <- rcumsum(rowsum(wrisk, y[, 1]))
		nrisk <- nrisk - c(esum, 0)[indx]
		irisk <- irisk - c(rcumsum(rowsum(wt, y[, 1])), 0)[indx]
	}
	haz <- nevent/nrisk
	result <- list(n = NROW(y), time = time, n.event = nevent
		, n.risk = irisk, n.censor = ncens, hazard = haz
		, chaz = cumsum(haz)
	)
	return(result)
}


#' Breslow estimator for the individuals at risk
#' 
#' Computes the number of indiduals ar risk given the linear predictor
#'
#' @param y \code{\link[survival]{Surv}} object.
#' @param lp a vector of linear predictor. Computed from the \code{X} matrix (predictors) and the estimated coefficients (\code{beta}).
#' @return A list of S3 objects. 
#' \item{n}{number of observations used in the fit.}
#' \item{time}{time points defined by the risk set.}
#' \item{n.event}{the number of events that occur at time \code{t}.}
#' \item{n.risk}{the number of individuals at risk at time \code{t}.}
#' \item{n.censor}{the number of subjects who exit the risk set, without an event, at time \code{t}.}
#'
#' @keywords internal
robustHazard <- function(y, lp) {

	## Number of observations
	N <- NROW(y)
	## Single or double time time vars as per Surv
	p <- NCOL(y)

	## Relative hazard
	relhaz <- exp(lp)
	
	## Initialize risk score for each patient
	risk <- numeric(N)

	## Individuals at risk
	n.risk <- numeric(N)

	if (p == 2) {
		endtime <- y[,1]
		events <- y[,2]
		for (i in 1:N){
			cond <- endtime[i] >= endtime
			indicator <- ifelse(cond, 1, 0)
			risk <- risk + (indicator * relhaz[[i]])
			n.risk <- n.risk + indicator
		}
	} else {
		starttime <- y[,1]
		endtime <- y[,2]
		events <- y[,3]
		for (i in 1:N){
			cond <- (endtime[i] >= endtime) & (starttime[i] < endtime)
			indicator <- ifelse(cond, 1, 0)
			risk <- risk + (indicator * relhaz[[i]])
			n.risk <- n.risk + indicator
		}
	}
	death <- (events == 1)
	nevent <- as.vector(rowsum(1*death, endtime))
	ncens <- as.vector(rowsum(1*!death, endtime))
	n.risk <- sort(unique(drop(n.risk)), decreasing = TRUE)
	return(list(n = N, time = endtime, risk = risk, n.risk = n.risk, n.event = nevent, n.censor = ncens, events = events))
}



