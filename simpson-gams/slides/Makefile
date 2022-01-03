all: slides purl copy

slides: gams.Rmd slides.css
	Rscript -e 'library("rmarkdown"); render("gams.Rmd")'

purl: gams.Rmd
	Rscript -e "knitr::purl(\"gams.Rmd\")"

copy: gams.html slides.css macros.js
	cp -R -u gams_files gams.html macros.js slides.css libs resources ~/work/web/jekyll/blog/slides/gams-webinar-efi-esa-2022/
