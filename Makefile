## This is posthoc function for survival analysis using glmnet

current: target
-include target.mk

# -include makestuff/perl.def

vim_session:
	bash -cl "vmt"

######################################################################

Sources += $(wildcard *.md *.R *.rmd)
Sources += $(wildcard R/*.R)
Sources += $(wildcard man/*.Rd) NAMESPACE DESCRIPTION

automatic_makeR = defined

######################################################################

glmnetpostsurv.Rout: R/glmnetpostsurv.R

glmnetHazard.Rout: R/glmnetHazard.R

glmnetsurvfit.Rout: R/glmnetsurvfit.R

methods.Rout: R/methods.R

######################################################################

### Makestuff

Sources += Makefile


Ignore += makestuff
msrepo = https://github.com/dushoff
Makefile: makestuff/Makefile
makestuff/Makefile:
	git clone $(msrepo)/makestuff
	ls $@

-include makestuff/os.mk
-include makestuff/git.mk
-include makestuff/visual.mk
-include makestuff/projdir.mk
-include makestuff/makeR.mk
-include makestuff/pandoc.mk

