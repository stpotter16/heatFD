docfiles = proj2-model.tex

if HAVE_PDFLATEX
docfiles += proj2-model.pdf
doc_DATA = proj2-model.pdf
#MANNAME = proj2-model
#MANTEXSRC = $(MANNAME).tex
#MANAUX = $(MANNAME).aux
#MANPDF = $(MANNAME).pdf
#CLEAFILES = $(MANPDF) $(MANNAME).log $(MANNAME).out $(MANNAME).toc
$(MANPDF): $(MANTEXSRC)
	$(PDFLATEX) $<
endif
