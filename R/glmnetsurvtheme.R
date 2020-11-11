#' Set theme for glmnetsurv plots
#' 
#' Sets a theme for glmnetsurv and other ggplot objects
#'
#' @examples
#' library(ggplot2)
#' gsurvtheme()
#' data(veteran, package="survival")
#' lam <- 0.02
#' alp <- 1
#' gfit <- glmnetsurv(Surv(time, status) ~ factor(trt) + karno + diagtime + age + prior
#'		, data = veteran
#'		, lambda = 0
#'		, alpha = 1
#'	)
#'
#' # Plot survival curves
#' gsurv <- glmnetsurvfit(gfit)
#' plot(gsurv)
#' @import ggplot2
#' @export

gsurvtheme <- function(){
   theme_set(theme_bw() +
      theme(panel.spacing = grid::unit(0,"lines")
      	, plot.title = element_text(hjust = 0.5)
			, legend.position = "bottom"
			, axis.ticks.y = element_blank()
			, axis.text.x = element_text(size = 12)
			, axis.text.y = element_text(size = 12)
			, axis.title.x = element_text(size = 12)
			, axis.title.y = element_text(size = 12)
			, legend.title = element_text(size = 13, hjust = 0.5)
			, legend.text = element_text(size = 13)
			, panel.grid.major = element_blank()
			, legend.key.size = unit(0.8, "cm")
			, legend.key = element_rect(fill = "white")
			, panel.spacing.y = unit(0.3, "lines")
			, panel.spacing.x = unit(1, "lines")
			, strip.background = element_blank()
			, panel.border = element_rect(colour = "grey"
				, fill = NA
				, size = 0.8
			)
			, strip.text.x = element_text(size = 11
				, colour = "black"
				, face = "bold"
			)
      )
   )
}


