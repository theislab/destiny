destiny.pdf: man/*.Rd
	RD2PDF_INPUTENC=inputenx R_RD4PDF=ae,hyper R CMD Rd2pdf --force --batch --no-preview --encoding=UTF-8 --output=$@ .
