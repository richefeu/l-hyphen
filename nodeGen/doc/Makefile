# Tutoriel nodegen : compilation LaTeX (latexmk + pdflatex)
#   make          tutorial_nodegen.pdf
#   make figures  régénère les figures (make_figures.sh, nécessite inkscape)
#   make clean    fichiers auxiliaires de LaTeX ; make distclean : aussi le PDF
LATEXMK ?= latexmk
DOC      = tutorial_nodegen

all: $(DOC).pdf

$(DOC).pdf: $(DOC).tex validation_results.tex $(wildcard figures/*)
	$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error $(DOC).tex

figures:
	./make_figures.sh

clean:
	$(LATEXMK) -c $(DOC).tex

distclean:
	$(LATEXMK) -C $(DOC).tex

.PHONY: all figures clean distclean
