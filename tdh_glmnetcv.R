library(shellpipes)
library(glmnetsurv)
library(splines)

loadEnvironments()

glmnetcv_mod <- glmnetsurvcv(mod_form
	, data=train_df
)
print(glmnetcv_mod)

saveVars(train_df
	, test_df
	, mod_form
	, glmnetcv_mod
)
