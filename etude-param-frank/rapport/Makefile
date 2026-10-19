# Compte rendu de etude-param-frank : compilation LaTeX (latexmk + pdflatex)
#   make          rapport.pdf (figures prises ici, dans ../cases/ et ../echantillon/)
#   make figures  régénère les figures à partir des résultats de ../cases/ (étude lancée au préalable)
#   make clean    fichiers auxiliaires de LaTeX ; make distclean : aussi le PDF
LATEXMK ?= latexmk
PYTHON  ?= python3
DOC      = rapport

all: $(DOC).pdf

$(DOC).pdf: $(DOC).tex $(wildcard *.png ../cases/*.png ../echantillon/*.png)
	$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error $(DOC).tex

figures:
	cd ../cases && for s in plot_maps plot_curves plot_snapshots plot_breaks plot_crackwatch; do \
		$(PYTHON) $$s.py || exit 1; done
	$(PYTHON) figures_rapport.py

clean:
	$(LATEXMK) -c $(DOC).tex

distclean:
	$(LATEXMK) -C $(DOC).tex

.PHONY: all figures clean distclean
