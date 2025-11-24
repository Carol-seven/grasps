objects := $(wildcard R/*.R) DESCRIPTION
version := $(shell egrep "^Version:" DESCRIPTION | awk '{print $$NF}')
pkg := $(shell egrep "^Package:" DESCRIPTION | awk '{print $$NF}')
tar := $(pkg)_$(version).tar.gz
tinytest := $(wildcard tests/testthat/*.R)
checkLog := $(pkg).Rcheck/00check.log
qmd := $(wildcard vignettes/*.qmd)
rmd := $(wildcard vignettes/*.Rmd)
vignettes := $(patsubst %.qmd,%.html,$(qmd)) $(patsubst %.Rmd,%.html,$(rmd))


.PHONY: check
check: $(checkLog)


.PHONY: build
build: $(tar)


.PHONY: install
install: $(tar)
	R CMD INSTALL $(tar)


.PHONY: preview
preview: $(vignettes)


.PHONY: pkgdown
pkgdown:
	Rscript -e "library(methods); pkgdown::build_site();"


$(tar): $(objects) $(qmd) $(rmd)
	@$(RM) -rf src/RcppExports.cpp R/RcppExports.R
	@Rscript -e "library(methods);" \
	-e "Rcpp::compileAttributes()" \
	-e "devtools::document();";
	R CMD build --compact-vignettes=gs+qpdf .


$(checkLog): $(tar) $(tinytest)
	R CMD check $(tar)


.PHONY: check-as-cran
check-as-cran: $(tar)
	R CMD check --as-cran $(tar)


.PHONY: check-revdep
check-revdep: $(tar)
	@mkdir -p revdep
	@rm -rf revdep/*.Rcheck
	@cp $(tar) revdep
	R CMD BATCH --no-save --no-restore misc/revdep_check.R &


vignettes/%.html: vignettes/%.Rmd
	Rscript -e "library(methods); rmarkdown::render('$?')"


vignettes/%.html: vignettes/%.qmd
	quarto render $< --to html --output $@


.PHONY: readme
readme: README.md
README.md: README.Rmd
	@Rscript -e "rmarkdown::render('$<')"


.PHONY: TAGS
TAGS:
	Rscript -e "utils::rtags(path = 'R', ofile = 'TAGS')"
	gtags

.PHONY: clean
clean:
	@$(RM) -rf *~ */*~ *.Rhistroy *.tar.gz src/*.so src/*.o \
	*.Rcheck/ *.Rout .\#* *_cache
