library(shellpipes)
library(glmnetsurv)
library(splines)

loadEnvironments()

lam <- glmnetcv_mod$lambda.min
lams <- glmnetcv_mod$lambdas.optimal
alp <- glmnetcv_mod$alpha.optimal
glmnet_mod <- glmnetsurv(mod_form
	, data=df
	, lambda=lams
	, s=lam
	, alpha=alp
)
print(glmnet_mod)
plot(glmnet_mod, xvar="lambda", label=TRUE)

saveVars(df
	, mod_form
	, glmnet_mod
)
