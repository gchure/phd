PY=python
PANDOC=pandoc

INPUTDIR=$(CURDIR)/src
OUTPUTDIR=$(CURDIR)dst
STYLEDIR=$(INPUTDIR)/styles
BIBFILE=$(INPUTDIR)/references.bib

pdf:
	pandoc $(INPUTDIR)/metadata.yaml \
	$(INPUTDIR)/chapter*/*.md \
	$(INPUTDIR)/ref_heading.md \
	-o $(OUTPUTDIR)/thesis.pdf \
	--default-image-extension=.pdf \
	--template=$(STYLEDIR)/template.tex \
	--filter pandoc-xnos \
	--filter pandoc-citeproc \
	--lua-filter=$(INPUTDIR)/frontmatter/short-captions.lua \
	--top-level-division chapter \
	--resource-path='.:chapter_01/figs:chapter_02/figs:chapter_06/figs:' \

html:
	bundle exec jekyll build --destination site;\
	sh copyfigs.sh ;\
	echo "done"


	

