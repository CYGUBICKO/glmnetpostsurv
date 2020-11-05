#' @export
glmnetsurvfit <- function(fit, Y, X, newX, s, wt = NULL,...) UseMethod("glmnetsurvfit")

#' @export
glmnetbasehaz <- function(fit, Y, X, s, wt = NULL, centered = TRUE) UseMethod("glmnetbasehaz")

