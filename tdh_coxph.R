library(shellpipes)
library(survival)
library(splines)

loadEnvironments()

coxph_mod <- coxph(mod_form
	, data=train_df
	, method="breslow"
	, x=TRUE
)
print(coxph_mod)

saveVars(train_df
	, test_df
	, mod_form
	, coxph_mod
)
