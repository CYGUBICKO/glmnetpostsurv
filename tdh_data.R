library(shellpipes)
library(dplyr)
library(survival)
library(splines)

loadEnvironments()

tdhPartition <- function(df, id, prop) {
	require(dplyr)
	index <- (df
		%>% group_by_at(id)
		%>% mutate(..index=runif(1,0,1)<prop)
		%>% pull(..index)
	)
	train_df <- df[index,]
	test_df <- df[!index,]
	out <- list(train_df=train_df, test_df=test_df)
}

df <- (survival::cgd
	%>% tdhPartition(id="id", prop=0.60)
)
train_df <- df$train_df
test_df <- df$test_df


mod_form <- as.formula(Surv(tstart, tstop, status) ~ treat + sex +
  ns(age,3) + height + weight + inherit + steroids + propylac + hos.cat
)

saveVars(train_df
	, test_df
	, mod_form
)
