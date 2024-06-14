pdflatex -shell-escape -interaction=nonstopmode figure.tex

# rename files
mv figure-figure0.pdf lift.pdf
mv figure-figure1.pdf drag.pdf

# clean up
rm figure.pdf figure.log figure.auxlock figure.aux figure.dvi
rm figure-figure0.* figure-figure1.*
