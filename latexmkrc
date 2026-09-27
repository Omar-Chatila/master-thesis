$pdf_mode = 1;                    # pdfLaTeX (widely available and reproducible)
$pdflatex = 'pdflatex --synctex=1 --interaction=nonstopmode --halt-on-error --file-line-error %O %S';
$lualatex = 'lualatex --synctex=1 --interaction=nonstopmode --halt-on-error --file-line-error %O %S';
$xelatex = 'xelatex --synctex=1 --interaction=nonstopmode --halt-on-error --file-line-error %O %S';
$bibtex = 'biber %O %B';
$out_dir = 'build';
$aux_dir = 'build';
$max_repeat = 5;

# Keep the source tree clean while retaining the PDF and SyncTeX output.
$clean_ext = 'acn acr alg glg glo gls ist loa lol nav snm vrb';
