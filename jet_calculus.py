from astropy.io import fits
from colormaps import cmaps
from matplotlib.colors import LogNorm
from astropy.io import fits
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from scipy.ndimage import fourier_shift
from itertools import groupby
import scipy.constants as const

#for drawing ellipses/circles
import math
from skimage.draw import (circle,circle_perimeter,ellipse,ellipse_perimeter)
def keyfunc(s):
	    return [int(''.join(g)) if k else ''.join(g) for k, g in groupby(s, str.isdigit)]

def read_fits(file):
	with fits.open(file) as hdulist:
		comps   = hdulist[1].data
		header  = hdulist[0].header
		img     = hdulist[0].data
		img			= img.reshape(img.shape[2],img.shape[3])
		#img			= img[0][0]
	return comps,header,img

def read_header(file):
	with fits.open(file) as hdulist:
		header  = hdulist[0].header
	return header

def px_to_mas(data,head):
	#convert pixel to mas
	xx,yy=data
	mas=3.6e6
	return round((xx - head['crpix1'])*head['cdelt2']*mas,2),round((yy - head['crpix2'])*head['cdelt2']*mas,2)

def PXPERBEAM(b_maj,b_min,px_inc):
	#print (str(b_maj))
	#print (str(px_inc))
	beam_area = np.pi/(4*np.log(2))*b_min*b_maj
	PXPERBEAM = beam_area/(px_inc**2)
	return PXPERBEAM
#
def axis_along_pa(major,minor,ellang,pa):
	'''
	Evaluates the axis of the beam along a given angle in polar coordinates.
	'''
	t =  (pa - ellang)*np.pi/180
	return (major*minor)/np.sqrt((minor*np.cos(t))**2+(major*np.sin(t))**2)
#
def mod_peak_flux(S,amaj,amin,bmaj,bmin):
	'''
	A way to derive the Model component peak flux
	'''
	return S/np.sqrt((amaj**2+bmaj**2)*(amin**2+bmin**2))*(bmaj*bmin)

def derive_positional_error(beam,compsize,snr_c,snr_m):
	'''
	resolution depend on beam size and component size as well as clean snr and model component snr.
	snr_c: S/N of map peak to post-fit rms noise
	snr_m: S/N of component peak flux density to post-fit rms noise
	'''
	return np.sqrt(np.square(beam/2./snr_c)+np.square(compsize/2./snr_m))
#
def derive_comp_width_error(W,snr):
	'''
	W: component FWHM
	snr: S/N of component peak flux density to post-fit rms noise
	'''
	return W/snr
#
def derive_tb(S,f,z,major,ratio=1,**kwargs):
	'''
	Derives brightness temperature (see Kovalev et al. 2005.
	S			: Flux density of component/Region [Jy]
	f			: Observing Frequency [GHz]
	major : Major axis ofcomponent/region[mas]
	minor : Minor axis of component/region[mas]
	'''
	Const = 2*np.log(2)*np.square(const.c)*1e-26/(const.pi*const.k*np.square(np.pi/180/3.6e6))
	tb = Const*S*(1+z)/(major**2*ratio*np.square(f*1e9))
	if 'Serr' in kwargs.keys():
		Serr = kwargs.get('Serr')
		if 'majorerr' in kwargs.keys():
			majorerr = kwargs.get('majorerr')
		else:
			raise Exception('Please give a value for majorerr as well as Serr to derive error for tb.\n')
		dtb = Const*(1+z)/(ratio*np.square(f*1e9))*np.sqrt(np.square(Serr/np.square(major))+np.square(2*S*majorerr/major**3))
		print('derived error as well.\n')
		return (tb,dtb)
	else:
		return tb

def draw_circ(img,x,y,radius):
	rr,cc = circle_perimeter(x,y,radius)
	img[rr,cc] = 1
	return img

def draw_ell(img,x,y,rmin,rmaj,posang):
	#returns rr,cc: index px that belong to the ellipse perimeter.
	rr,cc =ellipse_perimeter(x,y,rmin,rmaj,orientation=posang) #(x,y) center coord., posang major axis orientation clockwise [rad]
	img[rr,cc] = 1
	return img
#

def apply_shift(img,shift):
	offset_image = fourier_shift(np.fft.fftn(img), shift)
	imgalign = np.fft.ifftn(offset_image)
	img2 = imgalign.real

	return img2

