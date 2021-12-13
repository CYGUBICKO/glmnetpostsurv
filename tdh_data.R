library(shellpipes)
library(survival)
library(splines)

loadEnvironments()

df <- survival::cgd

mod_form <- as.formula(Surv(tstart, tstop, status) ~ treat + sex +
  ns(age,3) + height + weight + inherit + steroids + propylac + hos.cat
)

saveVars(df
	, mod_form
)
