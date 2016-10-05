logo.png: logo.pdf
	convert -density 300 logo.pdf -quality 90 logo.png

logo.pdf: logo.tex
	pdflatex logo.tex

clean:
	rm -vf logo.aux logo.log logo.png logo.pdf

