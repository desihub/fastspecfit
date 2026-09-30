TEX = emline-model continuum-modeling fastbayes
PDFS = $(TEX:%=%.pdf)

LATEXMK = latexmk

all: $(PDFS)

%.pdf: %.tex
	$(LATEXMK) -pdf -quiet $<
	$(LATEXMK) -c $*

clean:
	$(LATEXMK) -C
	rm -f $(PDFS)

.PHONY: all clean
