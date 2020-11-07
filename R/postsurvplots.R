#' Plot survival and cumulative hazard curves
#'
#' Plot estimated survival and cumulative  hazard curves for \code{glmnet} model.
#'\code{\link[glmnetpostsurv]{glmnetsurvfit.glmnet}}
#' @details
#' Depending on the specification in \code{\link[glmnetpostsurv]{glmnetsurvfit.glmnet}}, this function plots either average or individual survival or cumulative hazard curves. The plot is a \code{\link[ggplot2]{ggplot}} object, hence can be be customized further, see example below.
#'
#' @param x a \code{\link[glmnetpostsurv]{glmnetsurvfit.glmnet}} or \code{\link[glmnetpostsurv]{glmnetbasehaz.glmnet}} object.
#' @param ... for future implementations
#' @param type type of curve to generate. Either \code{type = "surv"} for survival curves or \code{type = "cumhaz"} for cumulative hazard curve.
#' @param lsize line size for the curves. Default is \code{0.3}.
#' @param compare logical. Whether to return plot with labels to add additional \code{geom} object for comparison. Default is \code{FALSE}.
#'
#' @return a \code{\link[ggplot2]{ggplot}} object.
#'
#' @examples
#'
#' library(ggplot2)
#' data(veteran, package="survival")
#'
#' lam <- 0
#' alp <- 1
#' sobj <- with(veteran, Surv(time, status))
#' X <- model.matrix(~factor(trt) + karno + diagtime + age + prior, data = veteran)[,-1]
#' gmodel <- glmnet(X, sobj, 'cox', alpha = alp, lambda = lam)
#'
#' # Survival estimate
#' gsurv <- glmnetsurvfit(fit = gmodel, Y = sobj, X = X, s = lam)
#' 
#' # Plot survival curves
#' plot(gsurv)
#'
#' # Baseline survival estimate
#' gbsurv <- glmnetbasehaz(gmodel, Y = sobj, X = X, s = lam, centered = FALSE)
#' plot(gbsurv)
#'
#' # Compare overall and baseline cumulative hazard
#' p1 <- plot(gsurv, type = "cumhaz", compare = TRUE)
#' df2 <- data.frame(time = gbsurv$time, cumhaz = gbsurv$hazard)
#' p2 <- (p1
#'		+ geom_line(data = df2, aes(x = time, y = cumhaz, group = 1, col = "baseline"))
#'		+ scale_colour_manual(name = "C. hazard"
#'			, values = c("#E41A1C", "#000000")
#'			, labels = c("baseline", "overall")
#'		)
#' )
#' print(p2)
#'
#' @import ggplot2
#' @export

plot.glmnetsurvfit <- function(x, ..., type = c("surv", "cumhaz"), lsize = 0.3, compare = FALSE) {
	type <- match.arg(type)
	if (inherits(x, "glmnetbasehaz")){
		cumhaz <- x$hazard
	} else {
		cumhaz <- x$cumhaz
	}
	surv <- x$surv
	time <- x$time
	plot_df <- data.frame(id = 1, time = time, surv = surv, cumhaz = cumhaz)
	if (NCOL(surv) > 1){
		nindivs <- NCOL(surv)
		individ <- as.factor(rep(1:nindivs, each = length(time)))
		surv <- as.vector(surv)
		cumhaz <- as.vector(cumhaz)
		time <- rep(time, nindivs)
		plot_df <- data.frame(id = individ, time = time, surv = surv, cumhaz = cumhaz)
	}

	id <- NULL
	p0 <- (ggplot(plot_df, aes(x = time, group = id), colour = "grey")
		+ labs(x = "Time")
		+ theme_bw()
		+ theme(panel.spacing = grid::unit(0,"lines"), legend.position = "bottom")
	)

	if (type == "surv"){
		p1 <- p0 + geom_line(aes(y = surv), size = lsize) + labs(y = "Survival prob.")
		if (compare){
			p1 <- p0 + geom_line(aes(y = surv, col = "glmnet"), size = lsize) + labs(y = "Survival prob.")
		}
	} else {
		p1 <- p0 + geom_line(aes(y = cumhaz), size = lsize) + labs(y = "Cummualtive hazard")
		if (compare){
			p1 <- p0 + geom_line(aes(y = cumhaz, col = "glmnet"), size = lsize) + labs(y = "Cummualtive hazard")
		}
	}
	return(p1)
}


#' Prediction performance
#'
#' Plots predictive performance of \code{glmnet} in survival analysis in comparison to other models. It uses risk scoring from \code{\link[riskRegression]{Score}}. This extension allows \code{glmnet} to support performance measure scoring by R package \code{pec}. See examples.
#'
#' @details
#' Implements plot method for \code{\link[riskRegression]{Score}} for time-dependent Brier score, AUC and ROC. However, currently, no support for time-dependent covariate models.
#'
#' @param x \code{\link[riskRegression]{Score}} object. See examples.
#' @param ... for future implementations.
#' @param type metric to return. Choices are \code{"roc", "auc", "brier"}.
#' @param pos spacing between the lines.
#'
#' @return a \code{\link[ggplot2]{ggplot}} object.
#'
#' @examples
#'
#' data(veteran, package="survival")
#' # glmnet
#' lam <- 0.02
#' alp <- 0.8
#' sobj <- with(veteran, Surv(time, status))
#' X <- model.matrix(~factor(trt) + karno + diagtime + age + prior, data = veteran)[,-1]
#' gfit1 <- glmnet(X, sobj, 'cox', alpha = alp, lambda = lam)
#'
#' # coxph
#' cfit1 <- coxph(Surv(time, status) ~ factor(trt) + karno + diagtime + age + prior
#'		, data = veteran
#'		, method = "breslow"
#'		, x = TRUE
#'		, y = TRUE
#'	)
#'
#' # Evaluate model performance at 90, 180, 365 time points
#' score_obj <- Score(list("coxph" = cfit1, "pcox" = gfit1)
#' 	, Surv(time, status) ~ 1
#' 	, data = veteran
#' 	, plots = "roc"
#'		, metrics = c("auc", "brier")
#' 	, B = 10
#'		, times = c(90, 180, 365)
#'		, newdata = X
#'		, Y = sobj
#'		, X = X
#' 	, s = lam
#' )
#'
#' # Plot AUC
#' plot(score_obj, type = "auc")
#' # Plot ROC
#' plot(score_obj, type = "roc")
#' # Plot brier
#' plot(score_obj, type = "brier")
#'
#' # Prediction error using pec package
#'\dontrun{
#' 	if (require("pec")) {
#'			pec_fit <- pec(list("coxph" = cfit1, "pcox" = gfit1)
#'				, Surv(time, status) ~ 1
#' 			, data = veteran
#' 			, splitMethod = "Boot632plus"
#'				, keep.matrix = TRUE
#'			)
#'			plot(pec_fit)
#' 	}
#'}
#'
#' @export

plot.Score <- function(x, ..., type = c("roc", "auc", "brier"), pos = 0.3){
	if (!inherits(x, "Score"))
		stop("Object should be score. See ?riskRegression::Score")
	type <- match.arg(type)
	if (type == "roc"){
		df <- x$ROC$plotframe
		FPR <- TPR <- model <- AUC <- lower <- upper <- Brier <- NULL
		p1 <- (ggplot(df, aes(x = FPR, y = TPR, color = as.factor(times)))
			+ geom_line(size = 1)
			+ geom_abline(size = 1, colour = "grey")
			+ facet_wrap(~model)
			+ labs(x = "1-Specificity", y = "Sensitivity", colour = "Time")
			+ scale_colour_viridis_d(option = "inferno")
			+ theme(legend.position = "right")
		)
	} else if (type == "auc"){
		df <- x$AUC$score
		p1 <- (ggplot(df, aes(x = model, y = AUC, colour = as.factor(times)))
			+ geom_point(position = position_dodge(pos))
			+ geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(pos))
			+ scale_colour_viridis_d(option = "inferno")
			+ labs(x = "AUC", y = "Model", colour = "Time")
			+ theme(legend.position = "right")
		)
	} else {
		df <- x$Brier$score
		p1 <- (ggplot(df, aes(x = model, y = Brier, colour = as.factor(times)))
			+ geom_point(position = position_dodge(pos))
			+ geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(pos))
			+ scale_colour_viridis_d(option = "inferno")
			+ labs(x = "Brier", y = "Model", colour = "Time")
			+ theme(legend.position = "right")
		)
	}
	return(p1)
}

