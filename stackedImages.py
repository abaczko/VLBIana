import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from pylab import *
import os,sys
import ehtim as eh
from glob import glob
from astropy.io import fits
from itertools import groupby
from astropy.table import Table
from jet_calculus import *
from astropy.time import Time
from skimage.restoration import estimate_sigma
import argparse

workDir		= os.getcwd()+'/'
plotDir = workDir+'plots/'
if not os.path.isdir(plotDir):
  os.makedirs(plotDir)


def all_equal(iterable):
    g = groupby(iterable)
    return next(g, True) and not next(g, False)

def read_shift(shiftFile):
	'''
	Returns shift in RA and DEC in mas
	'''
	with fits.open(shiftFile) as hdul:
		data=hdul[1].data
	return data.dx,data.dy,data.date
	

def apply_shift(img,shift):
	offset_image = fourier_shift(np.fft.fftn(img), shift)
	imgalign = np.fft.ifftn(offset_image)
	img2 = imgalign.real

	return img2

def estimate_noise(img):
#	if type(img) == str:
#		img = cv2.imread(img)
	return estimate_sigma(img, multichannel=False,average_sigmas=False)

class stackedImage():

	def __init__(self,imgFiles,shiftFile=None):
		self.header		= [read_header(m) for m in imgFiles]
		self.ehtFiles	= [eh.image.load_fits(m, aipscc=True) for m in imgFiles]
		self.beam			= [[h['BMAJ']*np.pi/180,h['BMIN']*np.pi/180,h['BPA']*np.pi/180] for h in self.header]
		self.fovx			= np.array([m.fovx() for m in self.ehtFiles])
		self.fovy			= np.array([m.fovy() for m in self.ehtFiles])
		self.naxis1		= np.array([h['NAXIS1'] for h in self.header])
		self.naxis2		= np.array([h['NAXIS2'] for h in self.header])
		self.px_inc		= np.array([h['cdelt2'] for h in self.header])
		self.noise_difmap	= np.array([h['noise'] for h in self.header])
		self.ppb			= [PXPERBEAM(b[0],b[1],pxi*np.pi/180) for b,pxi in zip(self.beam,self.px_inc)]
		self.mean_beam= np.array([np.sqrt(b[0]*b[1]) for b in self.beam])
#		self.stacked_beam	= np.sum(self.mean_beam)/len(self.mean_beam)
		self.stacked_beam = max(self.mean_beam)
		self.stacked_naxis= ma.amax([self.naxis1,self.naxis2])
		self.stacked_fov	= ma.amax([self.fovx,self.fovy])
		self.stacked_px_inc = self.stacked_fov/self.stacked_naxis*180/np.pi
		self.stacked_ppb	= PXPERBEAM(self.stacked_beam,self.stacked_beam,self.stacked_px_inc*np.pi/180)
		self.date_im	= None
		self.stackedImg = None
		self.blurImg	= None
		self.stacked_noise_diff = None
		#check whether all images have the same field of view. If not, regrid them to the smalles fov.
		if len(np.unique(np.concatenate([self.fovx,self.fovy])))>1:
			maxfov = [max(x,y) for x,y in zip(self.fovx,self.fovy)]
			self.ehtFiles = [f.pad(mf,mf) for f,mf in zip(self.ehtFiles,maxfov)]
			sys.stdout.write('Regrid images as not same FOV for all\n')
			self.ehtFiles = [f.regrid_image(self.stacked_fov,self.stacked_naxis,interp='linear') for f in self.ehtFiles]
	
		# blur with a circular gauss
		blurFiles		= [im.blur_circ(self.stacked_beam) for im in self.ehtFiles]
		self.blurImg			= np.array([im.imarr(pol='I') for im in blurFiles])
		if shiftFile is not None:
			shift_ra,shift_dec,shift_date = read_shift(shiftFile)
			self.date_im= [np.round(Time(h['DATE-OBS']).to_value('decimalyear'),4) for h in self.header]
			date_shift	= [np.round(Time(d, format='decimalyear').value,4) for d in shift_date]
			date_intersection = list(set(self.date_im).intersection(date_shift))
			if date_intersection.sort() == self.date_im.sort():
				sys.stdout.write('Shifts for all epochs are there. Will continue\n')
				shift_ra	= -np.array([a for a,b in zip(shift_ra,date_shift) if b in date_intersection])
				shift_dec = -np.array([a for a,b in zip(shift_dec,date_shift) if b in date_intersection])
				shift_date= np.array([a for a in date_shift if a in date_intersection])
			else:
				sys.stdout.write('Dates/Number of shifts are not equal to the given image observing dates. Please check.\n')
				return
############# CHECK if shift is estimated correctly! Do I need self.stacked+px_inc?
######################
			print(shift_ra)
			shift_ra = shift_ra/self.stacked_px_inc/3.6e6
			shift_dec= shift_dec/self.stacked_px_inc/3.6e6
			print(shift_ra)
			self.blurImg = np.array([apply_shift(f,[sd,sr]) for f,sd,sr in zip(self.blurImg,shift_dec,shift_ra)])
	
		self.stackedImg	= self.blurImg.sum(axis=0)/len(self.blurImg)*self.stacked_ppb
	
		print(np.sum(self.noise_difmap)/len(self.noise_difmap))	
	def write_out(self,outfile='stackedImage.fits'):
		hdr	= fits.Header(self.header[0])
		hdr['RMS']	= (estimate_noise(self.stackedImg),'estimated Gaussian noise')
		hdr['NOISE'] = (np.sum(self.noise_difmap)/len(self.noise_difmap),'averaged mean of map noises')
		hdr['CDELT1']	= (-self.stacked_px_inc,'Pixel increment (degrees)')
		hdr['CDELT2']	= (self.stacked_px_inc,'Pixel increment (degrees)')
		hdr['BMAJ']		= (self.stacked_beam*180/np.pi,'Stacked beam major (degrees)')
		hdr['BMIN']		= (self.stacked_beam*180/np.pi,'Stacked beam minor (degrees)')
		hdr['BPA']		= 0.0
		hdr['naxis1'] = (self.stacked_naxis,'Stacked naxis 1)')
		hdr['naxis2'] = (self.stacked_naxis,'Stacked naxis 2')
		hdr['datamin'] = (ma.min(self.stackedImg),'Min data of stacked')
		hdr['datamax'] = (ma.max(self.stackedImg),'Max data of stacked')
		hdr['crpix1'] = (self.stacked_naxis/2.,'Reference pixel')
		hdr['crpix2'] = (self.stacked_naxis/2.,'Reference pixel')

		fits.writeto(outfile, np.flipud(self.stackedImg), hdr, overwrite=True)

	def plot_stacked(self,img=None):
		if img==None:
			img = self.stackedImg
		#plt.ioff()
		fig1 = plt.figure()
		plt.imshow(img)
		plt.colorbar()
		plt.show()


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="If run as '__main__': Produce stacked image of given file list. If shifts for each image need to be applyed give path to a file containing the shifts.")
	parser.add_argument('f',help='File list, provide path to all files as xx* or aa,bb,cc...',nargs='+',type=argparse.FileType('rb'))
	parser.add_argument('-s','--shift_file',help='Path to an ascii table containing the shifts for RA and DEC',default=None,nargs=1,type=argparse.FileType('r'))
	parser.add_argument('-VK','--vlba_K',help='Only a shortcut if VLBA K values of Anne are loaded',default=None,type=str)

	args=parser.parse_args()

	if args.vlba_K is not None:
		imgFiles= glob('BR*.fits')
		imgFiles = sorted(imgFiles)
		if args.vlba_K == 'VLBAK':
			args.shift_file = 'shifts_k_q.fits'
			sys.stdout.write('VLBAK given to function. Hence, shift is applied.\n')
		elif args.vlba_K == 'VLBAQ':
			sys.stdout.write('VLBAQ given to function. Hence, no shift is applied.\n')
		sys.stdout.write('vlba_K given to function. Assume it is VLBA K or Q data and data has to be excluded\n')
		exclude =['BR099A','BR099B','BR099H','BR119B','BR130E','BR130F','BR130G','BR120A','BR120F','BR120G']
		imf=[]
		for im in imgFiles:
			ims = im.split('_')[0]
			if ims in exclude:
				print('remove {}'.format(im.split('_')[0]))
			else:
				imf.append(im)
		imgFiles = imf
	else:
		imgFiles = args.f
	if args.shift_file is not None:
		sys.stdout.write('Will shift data.\n')
		stacked_image = stackedImage(imgFiles,shiftFile=args.shift_file)
	else:
		stacked_image = stackedImage(imgFiles)
	stacked_image.write_out()
	stacked_image.plot_stacked()
