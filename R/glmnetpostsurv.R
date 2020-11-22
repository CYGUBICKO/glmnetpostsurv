#' Fit and perform post fitting procedures to glmnet survival models
#' 
#' Extends functionality of survival models in \code{\link[glmnet]{glmnet}} to compute survival curves and other calibrations.
#' @details
#' This functions offers a user friendly formular-data interface for fitting survival models using \code{\link[glmnet]{glmnet}}. Any additional \code{\link[glmnet]{glmnet}} arguments can be specified in \code{...}.
#'
#' @param formula Object of class formula describing 
#' the model. The response and terms are specified 
#' similar to \code{\link[survival]{Surv}} function.
#' @param data optional data frame containing 
#' variables specified in the formula.
#' @param family currently, only \code{glmnet} \code{"cox"} family (survival model) is allowed.
#' @param alpha the elasticnet mixing parameter, see \code{\link[glmnet]{glmnet}}.
#' @param lambda optional user supplied lambda sequence, see \code{\link[glmnet]{glmnet}}. It is recommended that you supply a sequence of \code{lambdas.optimal} from \code{\link[glmnetsurv]{glmnetsurvcv}} object.
#'	@param s a single value of lambda over which predictions or extractions are made. Ideally, this value should be obtained from  \code{\link[glmnetsurv]{glmnetsurvcv}} if not known. This can be \code{NULL} (or not specified) if a single value of \code{lambda} is specified, otherwise required if \code{lambda = NULL} or if \code{lambda} is a vector.
#' @param contrasts.arg an optional list. See 
#' the contrasts.arg of
#' \code{[stats]{model.matrix.default}}.
#' @param xlevs a named list of character vectors 
#' giving the full set of levels to be assumed 
#' for each factor. See \code{[stats]{model.frame}}.
#' @param na.action a function which indicates 
#' what should happen when the data contain NAs. 
#' See \code{[stats]{model.frame}}.
#' @param ... any of the options in \code{\link[glmnet]{glmnet}}.
#'
#' @seealso
#' \code{\link[glmnetsurv]{glmnetsurvcv}}, \code{\link[survival]{Surv}}, \code{\link[glmnet]{glmnet}}, \code{\link[glmnet]{cv.glmnet}}
#'
#' @return A list of \code{glmnetsurv} objects:
#' \item{fit}{fitted \code{\link[glmnet]{glmnet}} model object}
#' \item{X}{model matrix of model terms.}
#' \item{y}{Surv object defining the event times and event status.}
#' \item{s}{lambda used}
#'
#' @export
#' @importFrom stats .getXlevels aggregate approx as.formula coef coefficients delete.response model.extract model.frame model.matrix na.omit na.pass predict setNames terms reorder formula
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
#'	)
#' print(gfit1)
#'
#' # Perform cross-validation
#' gfit2 <- glmnetsurv(Surv(time, status) ~ factor(trt) + karno + diagtime + age + prior
#'		, data = veteran
#'		, lambda = NULL
#'		, alpha = 1
#' 	, s = 0.002
#'	)
#' plot(gfit2)

glmnetsurv <- function(formula = formula(data), data = sys.parent()
	, family = "cox", alpha = 1, lambda = NULL, s = NULL
	, contrasts.arg = NULL, xlevs = NULL, na.action = na.omit, ...){
	if(family != "cox")stop("Only cox family allowed currently!")
	if ((is.null(lambda)|length(lambda)>1) & is.null(s))stop("s is required for predictions.")
	if (any(s<0) | length(s)>1)stop("s is a non-negative single value.")
	
	if (is.null(s))s <- lambda
	sobj <- glmnetsurvdata(formula, data, contrasts.arg, xlevs, na.action)
	X <- sobj$X
	y <- sobj$y
	contrasts <- sobj$contrasts
	na.action <- sobj$na.action
	xlevels <- sobj$xlevels
	Terms2 <- sobj$terms
	glmnet_args <- list(x = X, y = y, family = family, alpha = alpha, lambda = lambda)
	new_args <- list(...)
	
	if (length(new_args))glmnet_args[names(new_args)] <- new_args
	fit <- do.call("glmnet", glmnet_args)
	
	fit$call <- match.call()
	result <- list(fit = fit, X = X, y = y, s = s
		, contrasts = contrasts, na.action = na.action
		, xlevels = xlevels, terms = Terms2
	)
	result$call <- match.call()
	class(result) <- "glmnetsurv"
  	return(result)
}

#' Cross-validation for glmnet survival models via glmnetsurv
#'
#' Performs \code{k}-fold cross-validation for \code{\link[glmnet]{glmnet}} via \code{\link[glmnetsurv]{glmnetsurvcv}}, plots
#' solution path plots, and returns optimal value of lambda
#' (and optimal alpha if more than one is given).
#'
#' @details
#' Performs cross-validation as illustrated in \code{\link[glmnet]{cv.glmnet}} but has additional capability to support more than one \code{alpha}.
#'
#' If more than one \code{alpha} is specified, say code{(0.2, 0.5, 1)}, the \code{glmnetsurvcv} will search for optimal values for alpha with respect to the corresponding lambda values. In this case, optimal alpha and lambda sequence will be returned, i.e., the \code{(alpha, lambda)} pair that corresponds to the lowest predicted cross-validated error (likelihood deviance).
#'
#'
#' @param formula Object of class formula describing 
#' the model. The response and terms are specified 
#' similar to \code{\link[survival]{Surv}} function.
#' @param data optional data frame containing 
#' variables specified in the formula.
#' @param family currently, only \code{glmnet} \code{"cox"} family (survival model) is allowed.
#' @param alpha elasticnet mixing parameter, with
#' \code{0 <= alpha <= 1}. If a vector of  \code{alpha} is supplied, cross-validation will be performed for each of the \code{alpha} and optimal value returned. The default is \code{1}.
#' @param lambda optional user supplied lambda sequence, \code{\link[glmnet]{cv.glmnet}}.
#' @param nfolds number of folds. Default is \code{10}.
#' @param foldid an optional sequence of values between \code{1} and {nfolds} specifying what fold each observation is in. This is important when comparing performance across models. If specified, \code{nfolds} can be missing.
#' @param refit logical. Whether to return solution path based on optimal lambda and alpha picked by the model. Default is \code{refit = TRUE}.
#' @param contrasts.arg an optional list. See 
#' the contrasts.arg of
#' \code{[stats]{model.matrix.default}}.
#' @param xlevs a named list of character vectors 
#' giving the full set of levels to be assumed 
#' for each factor. See \code{[stats]{model.frame}}.
#' @param na.action a function which indicates 
#' what should happen when the data contain NAs. 
#' See \code{[stats]{model.frame}}.
#' @param ... any of the options in \code{\link[glmnet]{cv.glmnet}}.
#'
#' @return An S3 object of class \code{\link[glmnetsurv]{glmnetsurvcv}}:
#' \item{lambda.min}{the value of lambda that gives minimum cross-validated error.}
#' \item{lambda.1se}{largest value of lambda such that error is within \code{1} standard error of the minimum.}
#' \item{alpha.optimal}{optimal alpha corresponding to \code{lambda.min}.}
#' \item{lambdas.optimal}{the sequence of lambdas containing \code{lambda.min}.}
#' \item{foldids}{the fold assignment used.}
#' \item{dfs}{list of data frames containing mean cross-validated error summaries and estimated coefficients in each fold.}
#' \item{fit}{if \code{refit = TRUE}, summaries corresponding to the optimal \code{alpha} and \code{lambdas}. This is used to plot solution path}.
#'
#' @seealso
#' \code{\link[glmnetsurv]{plot.glmnetsurvcv}}, \code{\link[glmnetsurv]{glmnetsurvcv}}, \code{\link[glmnet]{cv.glmnet}}
#'
#' @examples
#'
#' data(veteran, package="survival")
#' cv1 <- glmnetsurvcv(Surv(time, status) ~ factor(trt) + karno + diagtime + age + prior
#'		, data = veteran
#'		, alpha = 1
#'		, refit = FALSE
#'	)
#' print(cv1)
#'
#' # Train model using optimal alpha and lambda
#' fit1 <- glmnetsurv(Surv(time, status) ~ factor(trt) + karno + diagtime + age + prior
#'		, data = veteran
#'		, alpha = cv1$alpha.optimal
#'		, lambda = cv1$lambda.min
#'	)
#' print(fit1)
#'
#' @export

glmnetsurvcv <- function(formula = formula(data), data = sys.parent(), family = "cox"
	, alpha = 1, lambda = NULL, nfolds = 10, foldid = NULL, refit = TRUE
	, contrasts.arg = NULL, xlevs = NULL, na.action = na.omit, ...){
  	
	if(family != "cox")stop("Only cox family allowed currently!")
	sobj <- glmnetsurvdata(formula, data, contrasts.arg, xlevs, na.action)
	X <- sobj$X
	y <- sobj$y
	glmnet_args <- list(x = X, y = y, family = family, alpha = alpha, lambda = lambda)
  new_args <- list(...)
	if (length(new_args))glmnet_args[names(new_args)] <- new_args
	glmnet_args$nfolds <- nfolds
	glmnet_args$foldid <- foldid
	
	if(is.null(foldid)) foldid <- sample(rep(seq(nfolds), length.out = length(y)))
	glmnet_args$foldid <- foldid
	fit1 <- lapply(alpha, function(al){
		glmnet_args$alpha <- al
		mod <- do.call("cv.glmnet", glmnet_args)
		min_lambdas_df <- data.frame(lambda.min = mod$lambda.min
			, lambda.1se = mod$lambda.1se
			, cv.min = mod$cvm[which(mod$lambda==mod$lambda.min)]
			, alpha = al
		)
		tune_df <- data.frame(lambda = mod$lambda
			, alpha = al
			, cvm = mod$cvm
			, cvsd = mod$cvsd
			, cvlo = mod$cvlo
			, cvup = mod$cvup
		)
		beta_df <- lapply(mod$lambda, function(lam){
			bb <-	coef(mod, s = lam)
			p <- length(as.vector(bb))
			df <- data.frame(term = rownames(bb)
				, fold = rep(1, p)
				, estimate = as.vector(bb)
				, lambda = rep(lam, p)
				, alpha = rep(al, p)
				, l1_norm = rep(sum(abs(as.vector(bb))), p)
				, nzero = rep(mod$nzero[mod$lambda==lam], p)
			)
		})
		beta_df <- do.call("rbind", beta_df)
		res <- list(min_lambdas_df = min_lambdas_df
			, tune_df = tune_df, beta_df = beta_df
		)
		return(res)
	})
	fit1 <- do.call("rbind", fit1)
	min_metrics_df <- do.call("rbind", fit1[, "min_lambdas_df"])
	rownames(min_metrics_df) <- NULL
	cvm_df <- do.call("rbind", fit1[, "tune_df"])
	rownames(cvm_df) <- NULL
	beta_df <- do.call("rbind", fit1[, "beta_df"])
	rownames(beta_df) <- NULL

	min_lambdas <- min_metrics_df[which.min(min_metrics_df$cv.min), ]

	### Min lambda: optimal lambda
	lambda.min <- min_lambdas$lambda.min
	### 1 std error
	lambda.1se <- min_lambdas$lambda.1se
	### Optimal alpha
	alpha.optimal <- min_lambdas$alpha
	lambdas.optimal <- cvm_df$lambda[cvm_df$alpha==alpha.optimal]
	
	if(refit){
		glmnet_args$lambda <- lambdas.optimal
		glmnet_args$alpha <- alpha.optimal
		mod <- do.call("glmnet", glmnet_args)
		beta_est <- lapply(mod$lambda, function(lam){
			bb <- coef(mod, s = lam)
			p <- length(as.vector(bb))
			df <- data.frame(term = rownames(bb)
				, estimate = as.vector(bb)
				, lambda = rep(lam, p)
				, alpha = rep(alpha.optimal, p)
				, l1_norm = rep(sum(abs(as.vector(bb))), p)
				, nzero = rep(sum(bb!=0), p)
			)
		})
		beta_refit_df <- do.call("rbind", beta_est)
	} else {
		beta_refit_df <- NULL
	}

	out <- list(lambda.min = lambda.min, lambda.1se = lambda.1se
		, alpha.optimal = alpha.optimal, lambdas.optimal = lambdas.optimal
		, foldids = foldid, dfs = list(min_metrics_df = min_metrics_df
		, cvm_df = cvm_df, beta = beta_df), fit = list(beta = beta_refit_df)
	)
	out$call <- match.call()
	class(out) <- "glmnetsurvcv"
	return(out)

}

#' Prepare data for glmnet model
#'
#' @param formula Object of class formula describing 
#' the model. The response and terms are specified 
#' similar to \code{\link[survival]{Surv}} function.
#' @param data optional data frame containing 
#' variables specified in the formula.
#' @param contrasts.arg an optional list. See 
#' the contrasts.arg of
#' \code{[stats]{model.matrix.default}}.
#' @param xlevs a named list of character vectors 
#' giving the full set of levels to be assumed 
#' for each factor. See \code{[stats]{model.frame}}.
#' @param na.action a function which indicates 
#' what should happen when the data contain NAs. 
#' See \code{[stats]{model.frame}}.
#'
#' @keywords internal

glmnetsurvdata <- function(formula = formula(data), data = sys.parent()
	, contrasts.arg = NULL, xlevs = NULL, na.action = na.omit){
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
	result <- list(X = X, y = y, contrasts = contrasts
		, na.action = na.action, xlevels = xlevels, terms = Terms2
	)
  	return(result)
}


#' Compute variable importance of various survival models object
#'
#' @aliases varImp
#'
#' @details
#' Absolute value of the coefficients (parameters) corresponding the tuned model are used \code{type = param}. Otherwise, variable level importance is computed using permutation. See  \code{\link[glmnetsurv]{permuteImp}}.
#'
#' @param object fitted \code{\link[glmnetsurv]{glmnetsurv}}, \code{\link[pcoxtime]{pcoxtime}}, \code{\link[survival]{coxph}}, etc, object.
#' @param type if \code{type = "param"} absolute value of estimated coefficients are used. If \code{type = "variable"} variable level importance is computed using permutation.
#' @param scale if \code{TRUE} the scores are divided by the absolute sum of the coefficients.
#' @param newdata optional data frame containing the variables appearing on the right hand side of \code{\link[glmnetsurv]{glmnetsurv}} formula. Required if \code{type = "variable"}
#' @param nrep number of replicates for permulations. Only if \code{type = "variable"}.
#' @param ... not implemented. 
#'
#' @seealso
#' \code{\link[glmnetsurv]{plotImp}}
#'
#' @examples
#'
#' data(veteran, package="survival")
#' # glmnet
#' gfit1 <- glmnetsurv(Surv(time, status) ~ factor(trt) + karno + diagtime + age + prior
#'		, data = veteran
#'		, lambda = 0.02
#'		, alpha = 0.8
#'	)
#' imp1 <- varImp(gfit1, type = "param")
#' print(imp1)
#' imp2 <- varImp(gfit1, type = "variable", newdata = veteran)
#' print(imp2)
#'
#' @export

varImp <- function(object, type = c("param", "variable"), scale = TRUE, newdata, nrep = 20, ...){
	type <- match.arg(type)
	if (type=="param"){
		if (inherits(object, "glmnetsurv")){
			s <- object$s
			beta <- predict(object$fit, s = s, type = "coef")
			if(is.list(beta)) {
				out <- do.call("cbind", lapply(beta, function(x) x[,1]))
				out <- as.data.frame(out)
			} else out <- data.frame(Overall = beta[,1])	
		} else {
			beta <- coef(object)
			out <- data.frame(Overall = beta)
		}
		out <- out[rownames(out) != "(Intercept)",,drop = FALSE]
	} else {
		out <- data.frame(Overall = permuteImp(object, newdata, nrep))
	}
	out$sign <- sign(out$Overall)
	out$Overall <- abs(out$Overall)
	if (scale){
		out$Overall <- out$Overall/sum(out$Overall, na.rm = TRUE)
	}
	return(out)
}
