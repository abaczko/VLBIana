from fractions import Fraction
def set_size(width, fraction=1,subplots=(1,1)):
	""" Set aesthetic figure dimensions to avoid scaling in latex.
	Taken from https://jwalton.info/Embed-Publication-Matplotlib-Latex/
	
	Parameters
	----------
	width: float or string
	        Width in pts, or string of predined document type
	fraction: float,optional
	        Fraction of the width which you wish the figure to occupy
	subplots: array-like, optional
	The nu,ber of rows and columns of subplots
	
	Returns
	-------
	fig_dim: tuple
	        Dimensions of figure in inches
	"""
	if width.find('_')!=-1:
		w				= width.split('_')
		width		= w[0]
		fraction= float(w[1])
	if width =='aanda':
		width_pt = 256.0748
	elif width =='aanda*':
		width_pt = 523.5307
	elif width == 'beamer':
		width_pt = 720
	elif width == 'screen':
		width_pt = 600
	else:
		width_pt = width
	# Width of figure
	fig_width_pt = width_pt * fraction
	
	# Convert from pt to inches
	inches_per_pt = 1 / 72.27
	
	# Golden ratio to set aesthetic figure height
	golden_ratio = (5 ** 0.5 - 1) / 2.
	
	# Figure width in inches
	fig_width_in = fig_width_pt * inches_per_pt
	# Figure height in inches
	fig_height_in = fig_width_in * golden_ratio* (subplots[0] / subplots[1])
	
	return fig_width_in, fig_height_in
