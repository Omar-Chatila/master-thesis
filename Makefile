.PHONY: all clean watch designs

all:
	latexmk main.tex

watch:
	latexmk -pvc main.tex

designs:
	mkdir -p build/designs
	latexmk -jobname=design-habitat-grid designs/habitat-grid.tex
	latexmk -jobname=design-field-notes designs/field-notes.tex
	latexmk -jobname=design-technical-atlas designs/technical-atlas.tex
	latexmk -jobname=design-quiet-academic designs/quiet-academic.tex
	latexmk -jobname=design-image-watermark designs/image-watermark.tex
	latexmk -jobname=design-tengwar-book designs/tengwar-book.tex
	mv build/design-*.pdf build/designs/

clean:
	latexmk -C main.tex
