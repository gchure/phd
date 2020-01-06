PY=python
PANDOC=pandoc

INPUTDIR=$(CURDIR)
OUTPUTDIR=../dst
STYLEDIR=$(INPUTDIR)/styles
BIBFILE=$(INPUTDIR)/references.bib

pdf:
	pandoc $(INPUTDIR)/metadata.yaml \
	$(INPUTDIR)/chapter*/*.md \
	$(INPUTDIR)/ref_heading.tex \
	-o $(OUTPUTDIR)/thesis.pdf \
	--default-image-extension=.pdf \
	--template=$(STYLEDIR)/template.tex \
	--filter pandoc-xnos \
	--filter pandoc-citeproc \
	--lua-filter=$(INPUTDIR)/short-captions.lua \
	--top-level-division chapter \
	--resource-path='.:chapter_01/figs:chapter_02/figs' \

	

