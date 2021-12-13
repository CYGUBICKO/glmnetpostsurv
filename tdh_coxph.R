library(shellpipes)
library(survival)
library(splines)

loadEnvironments()

coxph_mod <- coxph(mod_form
	, data=df
	, method="breslow"
	, x=TRUE
)
print(coxph_mod)

saveVars(df
	, mod_form
	, coxph_mod
)
