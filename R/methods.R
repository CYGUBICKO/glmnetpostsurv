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

#' @export
glmnetsurvfit <- function(fit, newdata, ...) UseMethod("glmnetsurvfit")

#' @export
glmnetbasehaz <- function(fit, centered = TRUE) UseMethod("glmnetbasehaz")

