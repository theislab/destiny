from IPython.nbconvert.preprocessors.base import Preprocessor
from pygments.formatters import LatexFormatter

class MyLatexPreprocessor(Preprocessor):
	def preprocess(self, nb, resources):
		resources['latex']['pygments_definitions'] = LatexFormatter(style='trac').get_style_defs()
		return nb, resources 
