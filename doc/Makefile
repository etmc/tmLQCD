.PHONY: main.pdf all clean

all: main.pdf

main.pdf: main.tex
	latexmk -f -bibtex -pdf -pdflatex="pdflatex -interaction=nonstopmode" -use-make main.tex

clean:
	latexmk -CA
	rm -f main.spl main.tdo main.pdf *converted-to.pdf

distclean:
	latexmk -CA
	rm -f main.bbl main.spl main.tdo *converted-to.pdf
