#' Fit and perform post fitting procedures to glmnet survival models
#' 
#' Extends functionality of survival models in \code{\link[glmnet]{glmnet}} to compute survival curves and other calibrations.
#' @details
#' This functions offers a user friendly formular-data interface for fitting survival models using \code{\link[glmnet]{glmnet}}. It performs both the cross-validation and model fitting, depending on which option is specified in \code{fittype}. Any additional \code{\link[glmnet]{glmnet}} or \code{\link[glmnet]{cv.glmnet}} arguments can be specified in \code{...}.
#'
#' @param formula Object of class formula describing 
#' the model. The response and terms are specified 
#' similar to \code{\link[survival]{Surv}} function.
#' @param data optional data frame containing 
#' variables specified in the formula.
#' @param fittype whether to perform cross-validation \code{(fittype = "cv")} or fit a \code{glmnet} model \code{(fittype = "fit")}.
#' @param family currently, only \code{glmnet} \code{"cox"} family (survival model) is allowed.
#' @param alpha the elasticnet mixing parameter, see \code{\link[glmnet]{glmnet}}.
#' @param lambda optional user supplied lambda sequence, \code{\link[glmnet]{glmnet}}.
#' @param contrasts.arg an optional list. See 
#' the contrasts.arg of
#' \code{[stats]{model.matrix.default}}.
#' @param xlevs a named list of character vectors 
#' giving the full set of levels to be assumed 
#' for each factor. See \code{[stats]{model.frame}}.
#' @param na.action a function which indicates 
#' what should happen when the data contain NAs. 
#' See \code{[stats]{model.frame}}.
#' @param ... any of the options in \code{\link[glmnet]{glmnet}} or \code{\link[glmnet]{cv.glmnet}}.
#'
#' @seealso
#' \code{\link[survival]{Surv}}, \code{\link[glmnet]{glmnet}}, \code{\link[glmnet]{cv.glmnet}}
#'
#' @export
#' @importFrom stats .getXlevels aggregate approx as.formula coef coefficients delete.response model.extract model.frame model.matrix na.omit na.pass predict setNames terms
#' @docType package
#' @name glmnetsurv
#' @import glmnet
#' @examples
#'
#' # 
#' data(veteran, package="survival")
#' ## Fit unpenalized Cox using glmnet
#' lam <- 0 # Should fit unpenalized Cox model
#' gfit1 <- glmnetsurv(Surv(time, status) ~ factor(trt) + karno + diagtime + age + prior
#'		, data = veteran
#'		, lambda = lam
#'		, alpha = 1
#' 	, fittype = "fit"
#'	)
#' print(gfit1)
#'
#' # Perform cross-validation
#' gfit2 <- glmnetsurv(Surv(time, status) ~ factor(trt) + karno + diagtime + age + prior
#'		, data = veteran
#'		, lambda = NULL
#'		, alpha = 1
#' 	, fittype = "cv"
#'	)
#' plot(gfit2)

glmnetsurv <- function(formula = formula(data), data = sys.parent()
	, fittype = c("fit", "cv"), family = "cox", alpha = 1
	, lambda = NULL, contrasts.arg = NULL, xlevs = NULL
	, na.action = na.omit, ...){
	if(family != "cox")stop("Only cox family allowed currently!")
	fittype <- match.arg(fittype)
	call <- match.call()
	m <- match.call(expand.dots = FALSE)
	temp <- c("", "formula", "data")
	m <- m[match(temp, names(m), nomatch = 0)]
	Terms <- terms(formula, data = data)
	m$formula <- Terms
	m$na.action <- na.action
	m[[1]] <- as.name("model.frame")
	m <- eval(m, sys.parent())
	Terms2 <- terms(m)
	y <- model.extract(m, "response")
	term.labels <- attr(Terms, "term.labels")
	xlevels <- .getXlevels(Terms, m)
	
	if(!inherits(y, "Surv")) stop("formula: must be a survival formula. ?survival")
	X <- model.matrix(Terms, m, contr = contrasts.arg, xlev = xlevs)
	contrasts <- attr(X, "contrasts")
	xnames <- colnames(X)
	assign <- setNames(attr(X, "assign"), xnames)[-1]
	X <- X[,-1, drop = FALSE]
	glmnet_args <- list(x = X, y = y, family = family, alpha = alpha, lambda = lambda)
	new_args <- list(...)
	if (length(new_args))glmnet_args[names(new_args)] = new_args
	
	if (fittype=="fit"){
		fit <- do.call("glmnet", glmnet_args)
	} else {
		fit <- do.call("cv.glmnet", glmnet_args)
	}
	fit$call <- match.call()
	result <- list(fit = fit, X = X, y = y, s = fit$lambda
		, contrasts = contrasts, na.action = na.action
		, xlevels = xlevels, terms = Terms2
	)
	result$call <- match.call()
	class(result) <- "glmnetsurv"
  	return(result)
}


