PY=python
PANDOC=pandoc



pdf:
	cd src && \
	pandoc metadata.yaml \
	chapter*/*.md \
	ref_heading.md \
	-o ../dst/thesis.pdf \
	--default-image-extension=.pdf \
	--template=styles/template.tex \
	--filter pandoc-xnos \
	--filter pandoc-citeproc \
	--filter pandoc-crossref \
	--lua-filter=frontmatter/short-captions.lua \
	--top-level-division chapter \
	--resource-path='.:chapter_01/figs:chapter_02/figs:chapter_06/figs:' \
	--bibliography='references.bib' \
	&& cd - \

html:	
	JEKYLL_ENV=production bundle exec jekyll build --destination docs;\
	sh copyfigs.sh ;\


	

