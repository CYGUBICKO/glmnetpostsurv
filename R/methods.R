#' @export
glmnetsurvfit <- function(fit, Y, X, newX, s, wt = NULL,...) UseMethod("glmnetsurvfit")

#' @export
glmnetbasehaz <- function(fit, centered = TRUE) UseMethod("glmnetbasehaz")

