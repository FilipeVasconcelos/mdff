.SUFFIXES: .tex .aux .bib

filename=mdff

pspdfd:	psd
	ps2pdf -sPAPERSIZE=a4 ${filename}.ps

pspdf:	ps 
	ps2pdf -sPAPERSIZE=a4 ${filename}.ps

pdfd: 
	pdflatex -draftmode ${filename}
	bibtex ${filename} ||true
	makeindex ${filename}.idx
	makeindex ${filename}.gls 
	pdflatex -draftmode ${filename}
	pdflatex ${filename}
	
pdf: 
	pdflatex ${filename}
	bibtex ${filename} ||true
	makeindex ${filename}.idx
	makeindex ${filename}.gls 
	pdflatex ${filename}
	pdflatex ${filename}

pdf-print: ps
	ps2pdf -sPAPERSIZE=a4 -dColorConversionStrategy=/LeaveColorUnchanged -dPDFSETTINGS=/printer ${filename}.ps

text: html
	html2text -width 100 -style pretty ${filename}/${filename}.html | sed -n '/./,$$p' | head -n-2 >${filename}.txt

html:
	@#latex2html -split +0 -info "" -no_navigation ${filename}
	htlatex ${filename}

ps:	dvi
	dvips -t a4 ${filename}.dvi

psd: dvid
	dvips -t a4 ${filename}.dvi

dvi:
	latex ${filename}
	bibtex ${filename} ||true
	latex ${filename}
	latex ${filename}

dvid:
	latex -draftmode ${filename}
	bibtex ${filename} 
	makeindex ${filename}
	makeglossaries ${filename}
	latex -draftmode ${filename}
	latex -draftmode ${filename}

read:  pdf
	xpdf ${filename}.pdf 

aread: pdf
	acroread ${filename}.pdf

gvread: ps
	gv ${filename}.ps 

clean:
	rm -f *.log *.aux *.out *.bbl *.blg *.mtc* *.toc *.maf *.idx *.gls *.ilg *.glo *.glg
