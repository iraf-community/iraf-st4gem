doc=stdata_install

# batchmode means don't ask questions about errors; I don't know the answers anyway
all:
	pdflatex -interaction=batchmode ${doc}.tex

clean:
	rm -f *.dvi *.aux *.idx *.log *.out *.toc .log *.ilg *.ind
	rm -f Index.tex

distclean:      clean
	rm -f *.pdf

