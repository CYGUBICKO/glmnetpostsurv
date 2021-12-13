library(shellpipes)
library(prodlim)
library(glmnetsurv)
library(splines)

loadEnvironments()

glmnetsurv::gsurvtheme()

# Evaluate model performance at 90, 180, 365 time points
score_obj <- Score(list("coxph" = coxph_mod, "glmnet" = glmnet_mod)
	, formula=Hist(entry=tstart,time=tstop,event=status)~1
	, data = df[sample(1:nrow(df),100),]
	, plots = "roc"
	, metrics = c("auc", "brier")
	, B = 100
	, times = c(90, 180, 365)
)

# Plot AUC
plot(score_obj, type = "auc")
# Plot ROC
plot(score_obj, type = "roc")
# Plot brier
plot(score_obj, type = "brier")
