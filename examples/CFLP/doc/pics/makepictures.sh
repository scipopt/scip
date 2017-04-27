latex $1.tex
dvips $1.dvi
rm $1.dvi
rm $1.aux
rm $1.log
mv $1.ps $1_bla.eps
eps2eps $1_bla.eps $1.eps
rm $1_bla.eps
epstopdf $1.eps
echo I make pictures, photographic pictures
