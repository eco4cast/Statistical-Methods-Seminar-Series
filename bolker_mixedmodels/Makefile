
notes.slides.html: notes.Rmd
	. ./mkslides

notes.html: notes.Rmd
	Rscript -e "rmarkdown::render('notes.Rmd')"

full_notes.html: notes.Rmd
	export SHOW_RESULTS=TRUE; Rscript -e "rmarkdown::render('notes.Rmd', output_file = 'full_notes.html')"
