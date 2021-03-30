#!/usr/bin/env python
import string,math,sys,fileinput,glob,os,time
from numpy import *
from scipy import *
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from pylab import *
import os

##from skimage.feature import register_translation
from skimage.registration import phase_cross_correlation
from scipy.ndimage import fourier_shift
from matplotlib.patches import Circle,Ellipse
from skimage.draw import circle_perimeter,ellipse_perimeter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from VLBIana.modules.plot_functions import *
#import matplotlib.colors as colors
#font = {'family' : 'normal',
#				'weight' : 'normal',
#				'size'	 : 5}
#matplotlib.rc('font', **font)
import ehtim as eh

from VLBIana.modules.jet_calculus import *

#plt.rcParams.update({'font.size':7})

workDir		= os.getcwd()+'/'
plotDir = workDir+'plots/'
if not os.path.isdir(plotDir):
  os.makedirs(plotDir)

def apply_shift(img,shift):
	offset_image = fourier_shift(np.fft.fftn(img), shift)
	imgalign = np.fft.ifftn(offset_image)
	img2 = imgalign.real

	return img2

def align(img1,img2,inc1,inc2,mask1=False,mask2=False):

	#align image1 with image2 
	##shift, error, diffphase = register_translation((img1), (img2),100)
	if mask1 is False:
		shift,error,diffphase = phase_cross_correlation((img1),(img2),upsample_factor=100)
	else:
		shift = phase_cross_correlation((img1),(img2),upsample_factor=100,reference_mask=mask1,moving_mask=mask2)

	if mask1 is False:
		print ('register images new shift (y,x): [{} : {}] mas'.format(shift[0]*inc1,shift[1]*inc2))
		print ('register images new shift (y,x): {} px +- {}'.format(shift,error))

	else:
		print ('register images new shift (y,x): [{} : {}] mas'.format(shift[0]*inc1,shift[1]*inc2))
		print ('register images new shift (y,x): {} px'.format(shift))

	
	#shift img2 to found position
	img2 = apply_shift(img2,shift)
	
	#return aligned image and the computed shift
	if mask1 is False:
		return {'align':img2, 'shift':shift, 'error':error, 'diffphase': diffphase}
	else:
		return {'align':img2, 'shift':shift}
#==============================================

###run the script
'''
How to run the scipt:

def plot_aligned_maps(maps,flux_cut=False,radius=False,npix_x=False,npix_y=False,e_maj=False,e_min=False,e_pa=False,masked_shift=False):

def plot_aligned_maps(eht_maps,maps,masked_shift=False, **kwargs):	
	read in images
	files		= eht_maps
'''
def plot_aligned_maps(eht_maps,maps,masked_shift=False, **kwargs):
	'''read in images
	'''
	files   = eht_maps
	header	= [read_header(m) for m in maps]
	maps_beam=[[h['BMAJ']*np.pi/180,h['BMIN']*np.pi/180,h['BPA']*np.pi/180] for h in header]
	maps_ps = [dm.psize/eh.RADPERUAS/1e3 for dm in files]
	maps_fov = [dm.fovx()/eh.RADPERUAS/1e3 for dm in files]
	px_inc_ra=[h['cdelt1'] for h in header]
	px_inc_dec=[h['cdelt2'] for h in header]
	ppb=[PXPERBEAM(bp[0],bp[1],pxi*np.pi/180) for bp,pxi in zip(maps_beam,px_inc_dec)]
	ppb_r=PXPERBEAM(maps_beam[0][0],maps_beam[0][1],px_inc_dec[1]*np.pi/180)
	noise1=header[0]['NOISE'] 
	noise2=header[1]['NOISE']
	noise1_r = noise1/ppb[0]*ppb_r
	noise2_r = noise2/ppb[1]*ppb_r

	if 'sigma' in kwargs.keys():
		sigma = kwargs.get('sigma',False)
	else:
		sigma=1

	# regrid and blur the clean maps
	#check for same fov and pixel size
	if np.logical_and(maps_ps[0] != maps_ps[1], maps_fov[0]!=maps_fov[1]):
	
		beam	= maps_beam[0]
		fov		= files[1].fovx()
		naxis = header[1]['NAXIS1']
		file1regrid			=	files[0].regrid_image(fov, naxis, interp='linear')
		file1regridblur	=	file1regrid.blur_gauss(beam,frac=1)
		file2regrid			= files[1].regrid_image(fov, naxis, interp='linear')
		file2regridblur	= file2regrid.blur_gauss(beam, frac=1)
	
	file1rb = file1regridblur.imarr(pol='I').copy()
	file2rb = file2regridblur.imarr(pol='I').copy()

	mask1=np.zeros_like(file1rb, dtype=np.bool)
	mask2=np.zeros_like(file2rb, dtype=np.bool)
	# cut out inner, optically thick part of the image
	if 'npix_x' in kwargs.keys():
		npix_x = kwargs.get('npix_x',False)
		npix_x = kwargs.get('npix_x',False)
		px_min_x = int(naxis/2-npix_x)
		px_max_x = int(naxis/2+npix_x)
		px_min_y = int(naxis/2-npix_y)
		px_max_y = int(naxis/2+npix_y)

		px_range_x = np.arange(px_min_x,px_max_x+1,1)
		px_range_y = np.arange(px_min_y,px_max_y+1,1)

		index=np.meshgrid(px_range_y,px_range_x)
		mask1[tuple(index)] = True
		mask2[tuple(index)] = True

	if 'cut_left' in kwargs.keys():
		cut_left = kwargs.get('cut_left',False)
		px_max = int(naxis/2+cut_left)
		px_range_x = np.arange(0,px_max,1)
		mask1[:,px_range_x] = True
		mask2[:,px_range_x] = True

	if 'cut_right' in kwargs.keys():
		cut_right = kwargs.get('cut_right',False)
		px_max = int(naxis/2-cut_right)
		px_range_x = np.arange(px_max,naxis,1)
		mask1[:,px_range_x] = True
		mask2[:,px_range_x] = True

	if 'radius' in kwargs.keys():
		radius = kwargs.get('radius',False)
		rr,cc = circle(int(len(file1rb)/2),int(len(file1rb)/2),radius)
		mask1[rr,cc] = True
		mask2[rr,cc] = True

	if 'e_maj' in kwargs.keys():
		e_maj = kwargs.get('e_maj',False)
		e_min = kwargs.get('e_min',False)
		e_pa	= kwargs.get('e_pa',False)
		e_xoffset = kwargs.get('e_xoffset',False)
		if e_xoffset!=False:
			x,y = int(len(file1rb)/2)+e_xoffset,int(len(file1rb)/2)
		else:
			x,y = int(len(file1rb)/2),int(len(file1rb)/2)
		if e_pa!=False:
			e_pa = 90-e_pa
		else:
			e_pa = maps_beam[0][2]
		rr,cc =ellipse(y,x,e_maj,e_min,rotation=e_pa*np.pi/180)
		mask1[rr,cc] = True
		mask2[rr,cc] = True

	if 'flux_cut' in kwargs.keys():
		flux_cut = kwargs.get('flux_cut',False)
		mask1[file1>flux_cut*ma.amax(file1rb)] = True
		mask2[file2>flux_cut*ma.amax(file2rb)] = True

#	file1rbm = np.ma.masked_array(file1rb.copy(),mask1)
#	file2rbm = np.ma.masked_array(file2rb.copy(),mask2)
	file1rbm = file1rb.copy()
	file2rbm = file2rb.copy()
	file1rbm[mask1]=0
	file2rbm[mask2]=0


	#align image
	if masked_shift:
		file2_shift = align(file1rb,file2rb,px_inc_dec[1]*3.6e6,px_inc_ra[1]*3.6e6,mask1=~mask1,mask2=~mask2)
	else:
#		file2_shift=align(file1copy.imarr(pol='I'),file2copy.imarr(pol='I'),px_inc_dec[1]*3.6e6,px_inc_ra[1]*3.6e6)
		file2_shift=align(file1rbm,file2rbm,px_inc_dec[1]*3.6e6,px_inc_ra[1]*3.6e6)


######################
	#shift file2regridblur to found position
	#file2regridblur_shift = apply_shift(np.flipud(file2regridblur.imarr(pol='I')),file2_shift['shift'])
#	file2rb_shift = apply_shift(file2rb,file2_shift['shift'])
#	file2rb_shift *= ppb_r
###########################

	file1_plt = files[0].regrid_image(files[0].fovx(),header[0]['NAXIS1']).blur_gauss(maps_beam[0],frac=1).imarr(pol='I')
	file2_plt = files[1].regrid_image(files[1].fovx(),header[1]['NAXIS1']).blur_gauss(maps_beam[1],frac=1).imarr(pol='I')
	file1_plt		= file1_plt*ppb[0]
	file2_plt		= file2_plt*ppb[1]
	file1rb_plt	= file1rb*ppb_r
	file2rb_plt	= file2rb*ppb_r
	file1rbm_plt	= file1rbm*ppb_r
	file2rbm_plt	= file2rbm*ppb_r
	file2rb_shift_plt = apply_shift(file2rb,file2_shift['shift'])* ppb_r

	ra=file2regridblur.fovx()/eh.RADPERUAS/1e3/2
	dec=file2regridblur.fovy()/eh.RADPERUAS/1e3/4
	ra_min=-ra
	ra_max=ra
	dec_min=-dec
	dec_max=dec
	scale1	= px_inc_dec[0]*3.6e6
	scale2	= px_inc_dec[1]*3.6e6

	freq1 = np.round(header[0]['CRVAL3']*1e-9,1)
	freq2 = np.round(header[1]['CRVAL3']*1e-9,1)

	x1=np.linspace(-naxis*0.5*scale1,(naxis*0.5-1)*scale1,naxis)
	y1=np.linspace(naxis*0.5*scale1,-(naxis*0.5-1)*scale1,naxis)
	x2=np.linspace(-naxis*0.5*scale2,(naxis*0.5-1)*scale2,naxis)
	y2=np.linspace(naxis*0.5*scale2,-(naxis*0.5-1)*scale2,naxis)
	extent1 = np.max(x1), np.min(x1), np.min(y1), np.max(y1)
	extent2 = np.max(x2), np.min(x2), np.min(y2), np.max(y2)
	
	f,ax = plt.subplots(2,2,figsize=(set_size('aanda*',subplots=(2,2))),gridspec_kw={'hspace': 0.3, 'wspace':0.3})

	axe_ratio='scaled'
	
	ax[0,0].set_title('{} GHz original'.format(freq1))
	ax[0,1].set_title('{} GHz original'.format(freq2))
	ax[1,0].set_title('{} GHz regrid $+$ blur'.format(freq1))
	ax[1,1].set_title('{} GHz regrid $+$ blur'.format(freq2))

	level0	= min([noise1,noise2,noise1_r,noise2_r])*sigma
	lev=[]
	for i in range(0,10):
		lev.append(level0*2**i)

	level1r	= noise1_r*sigma
	lev1_r=[]
	for i in range(0,10):
		lev1_r.append(level1r*2**i)
	level2r	= noise2_r*sigma
	lev2_r=[]
	for i in range(0,10):
		lev2_r.append(level2r*2**i)

	imax = max([ma.amax(ii) for ii in [file1_plt,file2_plt,file1rb_plt,file2rb_plt]])
	im = ax[0,0].imshow(file1_plt,cmap=colormap,norm=mpl.colors.SymLogNorm(linthresh=level0,linscale=0.5,vmin=level0,vmax=0.5*imax),extent=extent1)
	im = ax[0,1].imshow(file2_plt,cmap=colormap,norm=mpl.colors.SymLogNorm(linthresh=level0,linscale=0.5,vmin=level0,vmax=0.5*imax),extent=extent2)
	im = ax[1,0].imshow(file1rb_plt,cmap=colormap,norm=mpl.colors.SymLogNorm(linthresh=level0,linscale=0.5,vmin=level0,vmax=0.5*imax),extent=extent2)
	im = ax[1,1].imshow(file2rb_plt,cmap=colormap,norm=mpl.colors.SymLogNorm(linthresh=level0,linscale=0.5,vmin=level0,vmax=0.5*imax),extent=extent2)

	#privImshow(file1_plt,noise1,extent1,ax[0,0])
	#privImshow(file2_plt,noise1,extent2,ax[0,1])
	#privImshow(file1rb_plt,noise1_r,extent2,ax[1,0])
	#privImshow(file2rb_plt,noise1_r,extent2,ax[1,1])

	#divider = make_axes_locatable(ax)
	#cax = divider.append_axes('right', size='5%', pad='3%')
	cbar = f.colorbar(im,ax=ax[:])
	cbar.set_label('Flux Density [Jy/beam]')

	for aa in ax.flat:
		aa.set(xlabel='RA [mas]', ylabel='Relative DEC [mas]')
		aa.xaxis.set_minor_locator(AutoMinorLocator())
		aa.yaxis.set_minor_locator(AutoMinorLocator())
		aa.axis(axe_ratio)
		aa.set_xlim(ra_max,ra_min)
		aa.set_ylim(dec_min,dec_max)
	f.savefig(plotDir+'{:d}GHz_convolved_with_{:d}GHz.pdf'.format(int(freq2),int(freq1)),bbox_inches='tight')
	
######################################################	
	
#	plt.cla()
	f,ax = plt.subplots(2,2,figsize=(set_size('aanda*',subplots=(2,2))),gridspec_kw={'hspace': 0.3, 'wspace':0.3})

	
	ax[0,0].set_title('{}GHz regrid $+$ blur'.format(freq1))
	ax[0,1].set_title('{}GHz regrid $+$ blur'.format(freq2))
	ax[1,0].set_title('{0}GHz/ {1}GHz not shifted'.format(freq1,freq2))
	ax[1,1].set_title('{0}GHz/ {1}GHz shifted'.format(freq1,freq2))

	imax = max([ma.amax(ii) for ii in [file1rbm_plt,file2rbm_plt]])
	im = ax[0,0].imshow(file1rbm_plt,cmap=colormap,norm=mpl.colors.SymLogNorm(linthresh=level0,linscale=0.5,vmin=level0,vmax=0.5*imax),extent=extent2)
	im = ax[0,1].imshow(file2rbm_plt,cmap=colormap,norm=mpl.colors.SymLogNorm(linthresh=level0,linscale=0.5,vmin=level0,vmax=0.5*imax),extent=extent2)
	#divider = make_axes_locatable(ax[0,1])
	#cax = divider.append_axes('right', size='5%', pad=0.05)
	#cbar = f.colorbar(col, use_gridspec=True,cax=cax)
	cbar = f.colorbar(im,ax=ax[0:1])
	cbar.set_label('Flux Density [Jy/beam]')
	#
	cntr1=ax[1,0].contour(np.flipud(file1rb_plt),linewidths=0.5,levels=lev1_r,colors='grey',extent=extent2,alpha=1)
	cntr2=ax[1,0].contour(np.flipud(file2rb_plt),linewidths=0.5,levels=lev2_r,colors='darkblue',extent=extent2,alpha=0.6)
	h1,_ = cntr1.legend_elements()
	h2,_ = cntr2.legend_elements()
	ax[1,0].legend([h1[0],h2[0]],['{}GHz'.format(freq1),'{}GHz'.format(freq2)],loc='upper right')
	#
	cntr1=ax[1,1].contour(np.flipud(file1rb_plt),linewidths=0.5,levels=lev1_r,colors='grey',extent=extent2)
	cntr2=ax[1,1].contour(np.flipud(file2rb_shift_plt),linewidths=0.5,levels=lev2_r,colors='darkblue',extent=extent2,alpha=0.6)
	h1,_ = cntr1.legend_elements()
	h2,_ = cntr2.legend_elements()
	ax[1,1].legend([h1[0],h2[0]],['{}GHz'.format(freq1),'{}GHz'.format(freq2)],loc='upper right')

	for aa in ax.flat:
		aa.set(xlabel='RA [mas]', ylabel='Relative Declination [mas]')
		aa.axis(axe_ratio)
		aa.xaxis.set_minor_locator(AutoMinorLocator())
		aa.yaxis.set_minor_locator(AutoMinorLocator())
		aa.set_xlim(ra_max,ra_min)
		aa.set_ylim(dec_min,dec_max)

	f.savefig(plotDir+'shifted_maps_{:d}GHz_{:d}GHz.pdf'.format(int(freq1),int(freq2)),bbox_inches='tight')
#########
	plt.cla()
	shift_export=file2_shift['shift'].copy()
	sys.stdout.write('final shift: {}'.format(shift_export))
	shift_export[0]*=px_inc_dec[1]*3.6e6
	shift_export[1]*=px_inc_ra[1]*3.6e6
	if masked_shift: 
	#error_export[1]*=px_inc_ra[1]*3.6e6
		sys.stdout.write('shift in mas: {}'.format(shift_export))
	else:
		error_export=file2_shift['error'].copy()
		error_export*=px_inc_dec[1]*3.6e6
		sys.stdout.write('shift in mas: {}\pm{}'.format(shift_export,error_export))
#########################
	# plot spix map
	file1rb_plt = np.flipud(file1rb_plt)
	file2rb_shift_plt = np.flipud(file2rb_shift_plt)
	spix1 = file1rb_plt*(file1rb_plt > noise1*sigma) #replaces indices where condition is not met with 0
	spix2 = file2rb_shift_plt*(file2rb_shift_plt > noise2*sigma)
	spix1[spix1==0] = noise1*sigma
	spix2[spix2==0] = noise2*sigma

	frequ1 = header[0]['CRVAL3']*1e-9
	frequ2 = header[1]['CRVAL3']*1e-9

	a = np.log10(spix2/spix1)/np.log10(freq2/freq1)

	spix_vmin,spix_vmax=-3,5
	sys.stdout.write('\nSpectral index max(alpha)={} - min(alpha)={}\nCutoff {}<alpha<{}\n'.format(ma.amax(a),ma.amin(a),spix_vmin,spix_vmax))
	a[a<spix_vmin] = spix_vmin
	a[a>spix_vmax] = spix_vmax
	a[spix2==noise2*sigma] =spix_vmin

	level10	= noise1*sigma
	lev1=[]
	level20	= noise2*sigma
	lev2=[]
	
	for i in range(0,10):
		lev1.append(level10*2**i)
		lev2.append(level20*2**i)

	f,ax = plt.subplots(figsize=(set_size('aanda')))
	
	cset = ax.contour(spix1,linewidths=[0.5],levels=lev1_r,colors=['grey'], extent=extent2,origin='lower',alpha=0.7)
	im = ax.imshow(a,cmap='hot_r',origin='lower',extent= extent2,vmin=spix_vmin,vmax=spix_vmax)
	divider = make_axes_locatable(ax)
	cax = divider.append_axes('right', size='5%', pad=0.05)
	cbar = f.colorbar(im, use_gridspec=True,cax=cax)
	cbar.set_label(r'Spectral index $\alpha$')
	h1,_ = cset.legend_elements()
	ax.legend([h1[0]],['{}GHz'.format(freq1)],loc='upper right')

	ax.axis('scaled')
	ax.set_xlabel('RA [mas]')
	ax.set_ylabel('Relative Dec [mas]')
	ax.set_xlim(ra_max,ra_min)
	ax.set_ylim(dec_min,dec_max)
	ax.minorticks_on()
	#ax.tick_params('both', length=8, width=2, which='major')
	#f.tick_params(axis='both',which='both',direction='in', labelsize=13)

#	cb = plt.colorbar(im,cax=cax,cmap='jet')
#	cb.ax.tick_params(labelsize=13, width=2)
#	cb.set_label(r'$\alpha$',fontsize=15)

	plt.savefig(plotDir+'spectral_index_between_{:d}_{:d}.pdf'.format(int(freq1),int(freq2)),bbox_inches='tight')
	plt.close('all')
############################
	if masked_shift:
		return {'file1':file1regridblur,'file2':file2regridblur,'shift':shift_export,'increment_dec':px_inc_dec[1]*3.6e6,'increment_ra':px_inc_ra[1]*3.6e6}
	else:
		return {'file1':file1regridblur,'file2':file2regridblur,'shift':shift_export,'increment_dec':px_inc_dec[1]*3.6e6,'increment_ra':px_inc_ra[1]*3.6e6,'error':error_export,'diffphase':file2_shift['diffphase']}
