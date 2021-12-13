library(shellpipes)
library(glmnetsurv)
library(splines)

loadEnvironments()

glmnetcv_mod <- glmnetsurvcv(mod_form
	, data=df
)
print(glmnetcv_mod)

saveVars(df
	, mod_form
	, glmnetcv_mod
)
