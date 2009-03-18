all: c-code

c-code:
	gcc -lgsl -lgslcblas -Wall avl.c ht.c thermo-community.c -o tc

pdf:
	pdflatex effdQ.tex
	bibtex effdQ
	pdflatex effdQ.tex
	pdflatex effdQ.tex
	open effdQ.pdf

clean:
	rm -f *.aux *.bbl *.blg *.log tc 2> /dev/null
