all: unmarked_presentation.html unmarked_presentation.R

unmarked_presentation.html: unmarked_presentation.Rmd
	Rscript -e "rmarkdown::render('unmarked_presentation.Rmd')"

unmarked_presentation.Rmd: unmarked_presentation.ipynb
	Rscript -e "rmarkdown:::convert_ipynb('unmarked_presentation.ipynb')"

unmarked_presentation.R: unmarked_presentation.Rmd
	Rscript -e "knitr::purl('unmarked_presentation.Rmd', documentation=2)"
