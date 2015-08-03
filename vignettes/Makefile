all: destiny.pdf

destiny.pdf: destiny.tex bibliography.bib
	rubber --pdf destiny.tex

destiny.tex: destiny.ipynb destiny.tplx mystyle_processors.py
	ipython nbconvert --Exporter.preprocessors="['mystyle_processors.MyLatexPreprocessor']" --to=latex --template=destiny.tplx destiny.ipynb

clean:
	rubber --clean destiny.tex
	rm -rf __pycache__

notebook:
	ipython notebook destiny.ipynb

watch:
	inotifyexec --filter 'destiny.ipynb|destiny.tplx|mystyle_processors.py' make

clean-all: clean
	rm destiny.tex
	rm -rf destiny_files
