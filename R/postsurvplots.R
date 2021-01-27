#' Plot survival and cumulative hazard curves
#'
#' Plot estimated survival and cumulative  hazard curves for \code{glmnet} model.
#'\code{\link[glmnetsurv]{glmnetsurvfit.glmnetsurv}}
#' @details
#' Depending on the specification in \code{\link[glmnetsurv]{glmnetsurvfit.glmnetsurv}}, this function plots either average or individual survival or cumulative hazard curves. The plot is a \code{\link[ggplot2]{ggplot}} object, hence can be be customized further, see example below.
#'
#' @param x a \code{\link[glmnetsurv]{glmnetsurvfit.glmnetsurv}} or \code{\link[glmnetsurv]{glmnetbasehaz.glmnetsurv}} object.
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
#' gmodel <- glmnetsurv(Surv(time, status) ~ factor(trt) + karno + diagtime + age + prior
#'		, data = veteran
#'		, lambda = 0
#'		, alpha = 1
#'	)
#'
#' # Survival estimate
#' gsurv <- glmnetsurvfit(fit = gmodel)
#' 
#' # Plot survival curves
#' plot(gsurv)
#'
#' # Baseline survival estimate
#' gbsurv <- glmnetbasehaz(gmodel, centered = FALSE)
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
		p1 <- p0 + geom_step(aes(y = surv), size = lsize) + labs(y = "Survival prob.")
		if (compare){
			p1 <- p0 + geom_step(aes(y = surv, col = "glmnet"), size = lsize) + labs(y = "Survival prob.")
		}
	} else {
		p1 <- p0 + geom_step(aes(y = cumhaz), size = lsize) + labs(y = "Cummualtive hazard")
		if (compare){
			p1 <- p0 + geom_step(aes(y = cumhaz, col = "glmnet"), size = lsize) + labs(y = "Cummualtive hazard")
		}
	}
	return(p1)
}


#' Plot solution path for glmnetsurvcv
#'
#' Plots the cross-validation curve, and upper and lower standard deviation curves, as a function of the optimal lambdas. Also, plots the solution path as a function of optimal lambdas (or randomly picked fold, if \code{refit = FALSE}) or \code{l1}-norm.
#'
#' @details
#' To plot solution path corresponding to optimal alpha and lambda, set \code{refit = TRUE} in \code{\link[glmnetsurv]{glmnetsurvcv}}. The plot is a \code{\link[ggplot2]{ggplot}} object, hence can be be customized further.
#'
#' @param x fitted \code{\link[glmnetsurv]{glmnetsurvcv}} object.
#' @param ... for future implementations
#' @param type which plot to return. \code{type = "cve"} (default) return a cross-validation curve and \code{type = "fit"} returns coefficient profiles (solution path). See details.
#' @param xvar only if \code{type = "fit"}. Plot coefficients a function of either lambda (\code{xvar = "lambda"}) or l1-norm (\code{xvar = "l1"}).
#' @param show_nzero logical. Whether to show number of nonzero coefficients on the plot. Default is \code{show_nzero = FALSE}. Still experimental for \code{type = "cve"}.
#' @param seed random number generator. Important if \code{refit = FALSE} in \code{\link[glmnetsurv]{glmnetsurvcv}}.
#'
#' @return a \code{\link[ggplot2]{ggplot}} object.
#'
#' @examples
#'
#' library(ggplot2)
#' data(veteran, package="survival")
#' # Using a vector of alphas = (0.8, 1)
#' cv1 <- glmnetsurvcv(Surv(time, status) ~ factor(trt) + karno + diagtime + age + prior
#' 	, data = veteran
#'		, alpha = c(0.8, 1)
#'		, refit = TRUE
#'	)
#' # Plot cross-validation curves
#' plot(cv1, type = "cve")
#'
#' # Plot
#' plot(cv1, type = "fit")
#'
#' @import ggplot2
#' @export

plot.glmnetsurvcv <- function(x, ..., type = c("cve", "fit"), xvar = c("lambda", "l1"), show_nzero = FALSE, seed = 1234){

	type <- match.arg(type)
	set.seed(seed)
	lambda <- cvm <- cvlo <- cvup <- NULL
	lambda.min <- lambda.1se <- NULL
	estimate <- term <- l1_norm <- NULL
	xvar_breaks <- nzero <- NULL
	sec_axisLabs <- function(var){
		lamb_tmp_df <- beta_df[, c(var, "nzero")]
		lamb_tmp_df <- lamb_tmp_df[!duplicated(lamb_tmp_df[[var]]), ]
		lamb_tmp_df <- lamb_tmp_df[order(-lamb_tmp_df$nzero),]
		if (var=="lambda"){
			lambdavals <- log(lamb_tmp_df[[var]])
		} else {
			lambdavals <- lamb_tmp_df[[var]]
		}
		nzero_breaks <- base::pretty(lambdavals, 5)
		clossest_lambda_breaks <- unlist(lapply(nzero_breaks, function(x)which.min(abs(lambdavals - x))))
		closest_lambda_to_breaks <- lambdavals[clossest_lambda_breaks]
		nzero_labels <- lamb_tmp_df[lambdavals %in% closest_lambda_to_breaks, ]
		return(data.frame(xvar_breaks = nzero_labels[[var]], nzero = nzero_labels$nzero))
	}
	if (type == "cve") {
		cvm_df <- x$dfs$cvm_df
		min_df <- x$dfs$min_metrics_df
		cvm_df$optimal <- ifelse(cvm_df$alpha==x$alpha.optimal, "  (Optimal)", "")
		min_df$optimal <- ifelse(min_df$alpha==x$alpha.optimal, "  (Optimal)", "")
		cvm_df$alpha <- as.factor(cvm_df$alpha)
		cvm_df$alpha_labels <- paste0("alpha== ", cvm_df$alpha, cvm_df$optimal)
		min_df$alpha_labels <- paste0("alpha== ", min_df$alpha, min_df$optimal)
		cvm_plot <- (ggplot(cvm_df, aes(x = log(lambda), y = cvm))
			+ geom_point(colour = "red", size = 0.2)
			+ geom_errorbar(aes(ymin = cvlo, ymax = cvup)
				, width = 0.01
				, colour = "red"
				, alpha = 0.4
			)
			+ facet_wrap(~alpha_labels, labeller = label_parsed, scales = "free_x")
			+ geom_vline(data = min_df, aes(xintercept = log(lambda.min)), lty = 2, size = 0.2)
			+ geom_vline(data = min_df, aes(xintercept = log(lambda.1se)), lty = 2, size = 0.2)
			+ labs(x = expression(log(lambda)), y = "Partial Likelihood Deviance")
		)
		if (show_nzero){
			beta_df <- x$dfs$beta
			lamb_tmp_df <- beta_df[, c("lambda", "alpha", "nzero")]
			lamb_tmp_df <- lamb_tmp_df[!duplicated(lamb_tmp_df[c("lambda","alpha", "nzero")]), ]
			labels_df <- lapply(split(lamb_tmp_df, lamb_tmp_df$alpha), function(dd){
				df <- sec_axisLabs("lambda")
				df$alpha <- unique(dd$alpha)
				df$optimal <- ifelse(df$alpha==x$alpha.optimal, "  (Optimal)", "")
				df$alpha_labels <- paste0("alpha== ", df$alpha,	df$optimal)
				return(df)
			})
			labels_df <- do.call("rbind", labels_df)
			cvm_plot <- (cvm_plot
				+ geom_text(data = labels_df, aes(y = Inf, x = log(xvar_breaks), label = nzero)
					, hjust = 0, vjust = 1
				)
			)
		}
		return(cvm_plot)
	} else {

		xvar <- match.arg(xvar)
		beta_df <- x$fit$beta
		facet <- FALSE
		if (is.null(beta_df)){
			facet <- TRUE
			beta_df = x$dfs$beta
			rand_fold <- sample(unique(beta_df$fold), 1)
			beta_df <- beta_df[beta_df$fold==rand_fold, ]
		}

		base_plot <- (ggplot(beta_df, aes(y = estimate, group = term, colour = term))
			+ scale_colour_viridis_d(option = "inferno")
			+ labs(y = "Coefficeint estimate", colour = "Predictor")
			+ theme(legend.position = "none")
		)

		if (xvar == "lambda"){
			coef_plot <- (base_plot + geom_line(aes(x = log(lambda)))
				+ scale_x_continuous(sec.axis = sec_axis(~.
						, breaks = log(sec_axisLabs("lambda")$xvar_breaks)
						, labels = sec_axisLabs("lambda")$nzero
					)
				)
				+ labs(x = expression(log(lambda)))
			)
			if (!facet){
				coef_plot <- (coef_plot
					+ geom_vline(xintercept = log(x$lambda.min), lty = 2, size = 0.2)
					+ geom_vline(xintercept = log(x$lambda.1se), lty = 2, size = 0.2)
				)
			}
		} else {
			coef_plot <- (base_plot
				+ geom_line(aes(x = l1_norm))
				+ scale_x_continuous(sec.axis = sec_axis(~.
							, breaks = sec_axisLabs("l1_norm")$xvar_breaks
							, labels = sec_axisLabs("l1_norm")$nzero
						)
					)
				+ labs(x = "L1 Norm")
			)
		}
		if (facet) {
			coef_plot <- coef_plot + facet_wrap(~alpha, scales = "free_x")
			message("These are CV coefficient plots. \nSet refit = TRUE to plot estimates based on whole dataset")
		}
		return(coef_plot)
	}
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
#' gfit1 <- glmnetsurv(Surv(time, status) ~ factor(trt) + karno + diagtime + age + prior
#'		, data = veteran
#'		, lambda = 0.02
#'		, alpha = 0.8
#'	)
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
#' score_obj <- Score(list("coxph" = cfit1, "glmnet" = gfit1)
#' 	, Surv(time, status) ~ 1
#' 	, data = veteran
#' 	, plots = "roc"
#'		, metrics = c("auc", "brier")
#' 	, B = 10
#'		, times = c(90, 180, 365)
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
#'			pec_fit <- pec(list("coxph" = cfit1, "glmnet" = gfit1)
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
	times <- NULL
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

#' Generic method for plotting variable importance of various models
#'
#' @param x a \code{\link[glmnetsurv]{glmnetsurv}} object. 
#' @param drop_zero if \code{TRUE} only nonzero estimates are shown.
#' @param label_size the size of text labels 
#' @param ... for future implementations
#'
#' @seealso
#' \code{\link[glmnetsurv]{varImp}}
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
#' plotImp(imp1)
#' imp2 <- varImp(gfit1, type = "variable", newdata = veteran)
#' plotImp(imp2)
#'
#' @export

plotImp <- function(x, drop_zero = FALSE, label_size = 3, ...){
	x$terms <- rownames(x)
	oldsign <- any(x$sign<0)
	x$sign <- ifelse(x$sign==1, "+", ifelse(x$sign==-1, "-", "0"))
	x <- x[order(x$Overall), ]
	if (drop_zero){
		x <- x[x$Overall!=0, ]
	}
	Overall <- NULL
	p1 <- (ggplot(x, aes(x = reorder(terms, Overall), y = Overall, group = 1, colour = sign))
		+ geom_point()
		+ geom_segment(aes(x = terms, xend = terms
				, y = 0, yend = Overall
				, colour = sign
			)
		)
		+ geom_text(aes(y = 0, label = paste0(terms, "  "))
			, size = label_size
			, vjust = 0.2
			, hjust = 1
			, show.legend = FALSE
		)
		+ scale_colour_manual(values = c("+" = "black", "-" = "red", "0" = "grey")
			, limits = c("+", "0", "-")
		)
		+ guides(fill = FALSE)
		+ labs(x = "", y = "Importance")
		+ coord_flip(clip = "off", expand = TRUE)
		+ theme_minimal()
		+ theme(plot.title = element_text(hjust = 0.5, face = "bold", colour = "grey")
			, axis.ticks.y = element_blank()
			, axis.text.y = element_blank()
			, panel.grid.major.y = element_blank()
		)
	)
	if (!oldsign) {
		p1 <- p1 + guides(colour = FALSE)
	}
	return(p1)
}

