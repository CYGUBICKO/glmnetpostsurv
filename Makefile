## This is posthoc function for survival analysis using glmnet

current: target
-include target.mk

# -include makestuff/perl.def

vim_session:
	bash -cl "vmt"

autopipeR = defined

######################################################################

Sources += $(wildcard *.md *.R *.rmd)
Sources += $(wildcard R/*.R)
Sources += $(wildcard man/*.Rd) NAMESPACE DESCRIPTION

automatic_makeR = defined

######################################################################

glmnetpostsurv.Rout: R/glmnetpostsurv.R

glmnetHazard.Rout: R/glmnetHazard.R

glmnetsurvfit.Rout: R/glmnetsurvfit.R

postsurvplots.Rout: R/postsurvplots.R

methods.Rout: R/methods.R

glmnetsurvtheme.Rout: R/glmnetsurvtheme.R

pkgsExport.Rout: R/pkgsExport.R

######################################################################

## Time-dependent covariate with glmnet and coxph
tdh_data.Rout: tdh_data.R

## coxph
tdh_coxph.Rout: tdh_coxph.R tdh_data.rda

## glmnet
tdh_glmnetcv.Rout: tdh_glmnetcv.R tdh_data.rda
tdh_glmnetmod.Rout: tdh_glmnetmod.R tdh_glmnetcv.rda

## Compare measures
tdh_compare.Rout: tdh_compare.R tdh_coxph.rda tdh_glmnetmod.rda

######################################################################

## Package installation and checks
Ignore += glmnetsurv_1*

build-package:
	R CMD build .

install-tarball:
	R CMD INSTALL glmnetsurv_1*

check-package:
	echo "devtools::check('.')" | R --slave

update-doc:
	echo "devtools::document('.')" | R --slave

install:
	$(MAKE) update-doc build-package install-tarball

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

-include makestuff/texi.mk
-include makestuff/pipeR.mk

-include makestuff/git.mk
-include makestuff/visual.mk

