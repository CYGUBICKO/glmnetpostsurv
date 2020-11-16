#' Print method fot  glmnetsurv object
#'
#' This function prints a summary of the glmnet object from  glmnetsurv.
#'
#' @details The call that produced \code{\link[glmnetsurv]{glmnetsurv}} is printed, followed by \code{\link[glmnet]{glmnet}} objects. 
#'
#' @param x fitted \code{\link[glmnetsurv]{glmnetsurv}} model object 
#' @param ... for future implementations
#'
#' @return see \code{\link[glmnet]{glmnet}}.
#'
#' @method print glmnetsurv
#' @export
#' @export print.glmnetsurv

print.glmnetsurv <- function(x, ...){
	x <- x$fit
	print(x, ...)
}

#' Plot glmnet objects through glmnetsurv
#'
#' @details
#' Internally extract \code{\link[glmnet]{glmnet}} and call plot method
#'
#' @param x \code{\link[glmnetsurv]{glmnetsurv}} object.
#' @param ... additional \code{\link[glmnet]{plot.glmnet}} plots methods arguments.
#'
#' @method plot glmnetsurv
#' @export
#' @export plot.glmnetsurv

plot.glmnetsurv <- function(x, ...){
	x <- x$fit
	plot_args <- list(x = x)
  	new_args <- list(...)
  	if(length(new_args))plot_args[names(new_args)] = new_args
  	do.call("plot", plot_args)
}

#' Print cross-validated glmnetsurvcv object
#'
#' Print the summary of the result of cross-validation for a glmnetsurv object.
#'
#' @details
#' A summary of optimal lambda and alpha for training glmnetsurv model.
#'
#' @param x \code{\link[glmnetsurv]{glmnetsurvcv}} object
#' @param ... for future implementations
#'
#' @method print glmnetsurvcv
#' @export
#' @export print.glmnetsurvcv
print.glmnetsurvcv <- function(x, ...){
	cat("Call:\n")
	print(x$call)
	cat("\nOptimal parameter values\n")
	out <- data.frame(cbind(lambda.min = x$lambda.min, lambda.1se = x$lambda.1se, alpha.optimal = x$alpha.optimal))
	print(out, row.names = FALSE)
	cat("\n")
}

#' Print a short summary of survival function
#'
#' Print the number of observations and number of events.
#'
#' @details Provide a summary of \code{\link[glmnetsurv]{glmnetsurvfit.glmnetsurv}} object.
#'
#' @param x the result of a call to the \code{\link[glmnetsurv]{glmnetsurvfit.glmnetsurv}} function.
#' @param ... for future implementations
#'
#' @return The call to the \code{\link[glmnetsurv]{glmnetsurvfit.glmnetsurv}} and the summary of the survival function.
#'
#' @method print glmnetsurvfit
#' @export
#' @export print.glmnetsurvfit
print.glmnetsurvfit <- function(x, ...){
	if (!inherits(x, "glmnetbasehaz")){
		cat("Call:\n")
		print(x$call)
		out <- data.frame(cbind(n = x$n, events = sum(x$events)))
		print(out, row.names = FALSE)
		cat("\n")
	}
}

#' @export
glmnetsurvfit <- function(fit, newdata, ...) UseMethod("glmnetsurvfit")

#' @export
glmnetbasehaz <- function(fit, centered = TRUE) UseMethod("glmnetbasehaz")

#' @export 
varImp <- function(object, show_sign = FALSE, scale = FALSE, ...) UseMethod("varImp")

