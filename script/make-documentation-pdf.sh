#!/bin/bash -v
# run from project root directory; NOT THIS DIRECTORY

doc_mdpp="doc/manual.mdpp"
doc_md="doc/manual.md"
doc_md2="doc/manual2.md"
doc_tpl="doc/manual-template.tex"
pdf_engine="lualatex"
dst_doc_pdf="doc/manual.pdf"

markdown-pp "$doc_mdpp" -o "$doc_md"

# remove lines: %... (needed for title in manpage)
# replace lines: ?... with %...
cat "$doc_md" | sed '/^%/d' | sed 's/^?\(.*\)/%\1/g' > "$doc_md2"
# generate pdf
pandoc --template="$doc_tpl" "$doc_md2" --toc --pdf-engine=$pdf_engine -o "$dst_doc_pdf"
# remove temporary
rm "$doc_md" "$doc_md2"

echo "mkdoc $dst_doc_pdf"