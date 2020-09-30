#!/usr/bin/env python
import string,math,sys,fileinput,glob,os,time
from numpy import *
from scipy import *
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from pylab import *

from skimage.feature import register_translation
from skimage.feature.register_translation import _upsampled_dft
from scipy.ndimage import fourier_shift
from matplotlib.patches import Circle,Ellipse
from skimage.draw import circle_perimeter,ellipse_perimeter
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

#font = {'family' : 'normal',
#				'weight' : 'normal',
#				'size'	 : 5}
#matplotlib.rc('font', **font)
import ehtim as eh

from jet_calculus import *

plt.rcParams.update({'font.size':7})

def apply_shift(img,shift):
	offset_image = fourier_shift(np.fft.fftn(img), shift)
	imgalign = np.fft.ifftn(offset_image)
	img2 = imgalign.real

	return img2

def align(img1,img2,inc1,inc2):

	#align image1 with image2 
	shift, error, diffphase = register_translation((img1), (img2),100)
	print ('register images new shift (y,x): [{} : {}] mas'.format(shift[0]*inc1,shift[1]*inc2))
	print ('register images new shift (y,x): {} px +- {}'.format(shift,error))

	
	#shift img2 to found position
	img2 = apply_shift(img2,shift)
	
	#return aligned image and the computed shift
	return {'align':img2, 'shift':shift}
#==============================================

###run the script
'''
How to run the scipt:

def plot_aligned_maps(maps,flux_cut=False,radius=False,npix_x=False,npix_y=False,e_maj=False,e_min=False,e_pa=False):

'''
def plot_aligned_maps(eht_images,maps, **kwargs):	
	'''read in images'''
	files		= eht_images
	header	= [read_header(m) for m in maps]
	maps_beam=[[h['BMAJ']*np.pi/180,h['BMIN']*np.pi/180,h['BPA']*np.pi/180] for h in header]
	maps_ps = [dm.psize/eh.RADPERUAS/1e3 for dm in files]
	maps_fov = [dm.fovx()/eh.RADPERUAS/1e3 for dm in files]
	px_inc_ra=[h['cdelt1'] for h in header]
	px_inc_dec=[h['cdelt2'] for h in header]
	ppb=[PXPERBEAM(bp[0],bp[1],pxi*np.pi/180) for bp,pxi in zip(maps_beam,px_inc_dec)]
	

	#check for same fov and pixel size
	if np.logical_and(maps_ps[0] != maps_ps[1], maps_fov[0]!=maps_fov[1]):
	
		beam	= maps_beam[0]
		fov		= files[1].fovx()
		naxis = header[1]['NAXIS1']
		file2regrid			= files[1].regrid_image(fov, naxis, interp='linear')
		file2regridblur	= file2regrid.blur_gauss(beam, frac=1)
		file1regrid			=	files[0].regrid_image(fov, naxis, interp='linear')
		file1regridblur	=	file1regrid.blur_gauss(beam,frac=1)

	ppb_r=PXPERBEAM(maps_beam[0][0],maps_beam[0][1],px_inc_dec[1]*np.pi/180)
	file1 = file1regridblur.imarr(pol='I').copy()
	file2 = file2regridblur.imarr(pol='I').copy()

	file1=np.flipud(file1)
	file2=np.flipud(file2)

	file1_plt=files[0].blur_gauss(maps_beam[0],frac=1).imarr(pol='I')
	file2_plt=files[1].blur_gauss(maps_beam[1],frac=1).imarr(pol='I')
	
	file1_plt = np.flipud(file1_plt)*ppb[0]
	file2_plt = np.flipud(file2_plt)*ppb[1]

	sys.stdout.write('max file1={}\n'.format(np.amax(file1_plt)))
	sys.stdout.write('max file2={}\n'.format(np.amax(file2_plt)))

	file1regridblur_plt		= np.flipud(file1regridblur.imarr(pol='I'))*ppb_r
	file2regridblur_plt		= np.flipud(file2regridblur.imarr(pol='I'))*ppb_r
	sys.stdout.write('max file1_regrid={}\n'.format(np.amax(file1regridblur_plt)))
	sys.stdout.write('max file2_regrid={}\n'.format(np.amax(file2regridblur_plt)))

	noise1=header[0]['NOISE'] # in the beginning I divided by ppb[0]
	noise2=header[1]['NOISE']
	noise1_r = noise1/ppb[0]*ppb_r
	noise2_r = noise2/ppb[1]*ppb_r

	sys.stdout.write ('noise1={0}\nnoise2={1}\n'.format(noise1,noise2))
	sys.stdout.write('noise1_regrid={0}\nnoise2_regrid={1}\n'.format(noise1_r,noise2_r))
	if 'sigma' in kwargs.keys():
		sigma = kwargs.get('sigma',False)
	else:
		sigma=1
	level10	= noise1*sigma
#	lev1		= [-level10]
	lev1=[]
	level20	= noise2*sigma
#	lev2		= [-level20]
	lev2=[]
	
	for i in range(0,10):
		lev1.append(level10*2**i)
		lev2.append(level20*2**i)

	level10_r	= noise1_r*sigma
#	lev1_r		= [-level10_r]
	lev1_r=[]
	level20_r	= noise2_r*sigma
#	lev2_r		= [-level20_r]
	lev2_r=[]
	
	for i in range(0,10):
		lev1_r.append(level10_r*2**i)
		lev2_r.append(level20_r*2**i)
	
	# cut out inner, optically thick part of the image
	if 'npix_x' in kwargs.keys():
		npix_x = kwargs.get('npix_x',False)
		npix_y = kwargs.get('npix_y',False)
		px_min_x = int(naxis/2-npix_x)
		px_max_x = int(naxis/2+npix_x)
		px_min_y = int(naxis/2-npix_y)
		px_max_y = int(naxis/2+npix_y)

		px_range_x = np.arange(px_min_x,px_max_x+1,1)
		px_range_y = np.arange(px_min_y,px_max_y+1,1)

		index=np.meshgrid(px_range_y,px_range_x)
		file1[tuple(index)]=0
		file2[tuple(index)]=0
	
	if 'cut_left' in kwargs.keys():
		cut_left = kwargs.get('cut_left',False)
		px_max = int(naxis/2+cut_left)
		px_range_x = np.arange(0,px_max,1)
		file1[:,px_range_x]=0
		file2[:,px_range_x]=0

	if 'cut_right' in kwargs.keys():
		cut_right = kwargs.get('cut_right',False)
		px_max = int(naxis/2-cut_right)
		px_range_x = np.arange(px_max,naxis,1)
		file1[:,px_range_x]= 0
		file2[:,px_range_x]= 0

	if 'radius' in kwargs.keys():
		radius = kwargs.get('radius',False)
		rr,cc = circle(int(len(file1)/2),int(len(file1)/2),radius)
		file1[rr,cc]= 0
		file2[rr,cc]= 0

	if 'e_maj' in kwargs.keys():
		e_maj = kwargs.get('e_maj',False)
		e_min = kwargs.get('e_min',False)
		e_pa	= kwargs.get('e_pa',False)
		e_xoffset = kwargs.get('e_xoffset',False)
		if e_xoffset!=False:
			x,y = int(len(file1)/2)+e_xoffset,int(len(file1)/2)
		else:
			x,y = int(len(file1)/2),int(len(file1)/2)
		if e_pa!=False:
			e_pa = 90-e_pa
		else:
			e_pa = maps_beam[0][2]
		rr,cc =ellipse(y,x,e_maj,e_min,rotation=e_pa*np.pi/180)
		file1[rr,cc]= 0
		file2[rr,cc]= 0

	if 'flux_cut' in kwargs.keys():
		flux_cut = kwargs.get('flux_cut',False)
		file2[file2>flux_cut*ma.amax(file2)]=0
		file1[file1>flux_cut*ma.amax(file1)]=0

	file1[file1<=0] = 0
	file2[file2<=0] = 0

	file1copy = file1regridblur.copy()
	file2copy = file2regridblur.copy()

	file1copy.imvec=file1.flatten()
	file2copy.imvec=file2.flatten()

	#align image
	file2_shift=align(file1copy.imarr(pol='I'),file2copy.imarr(pol='I'),px_inc_dec[1]*3.6e6,px_inc_ra[1]*3.6e6)

	file1=file1*ppb_r
	file2=file2*ppb_r
	file2_shift['align']=file2_shift['align']*ppb_r

######################
	#shift file2regridblur to found position
	file2regridblur_shift = apply_shift(np.flipud(file2regridblur.imarr(pol='I')),file2_shift['shift'])
	file2regridblur_shift = file2regridblur_shift*ppb_r
###########################

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
	#image colormap
	im_color = 'inferno' #string for matplotlib colormap
	widths = [float(file1.shape[1]) for i in range(2)]
	heights = [dec_max/float(ra_max)*widths[i] for i in range(2)]

	fig_width = 10.
	fig_height = fig_width * sum(heights) / sum(widths)
	ratios = [fh/heights[0] for fh in heights]
	f,ax = plt.subplots(2,2,figsize=(fig_width,fig_height))

	axe_ratio='scaled'

	ax[0,0].set_title('{} GHz original'.format(freq1))
	col = ax[0,0].imshow(file1_plt,cmap=im_color,norm=colors.SymLogNorm(linthresh=level10,linscale=0.5,vmin=lev1[0],vmax=0.5*ma.amax(file1_plt)),extent=extent1,origin='lower')
	divider = make_axes_locatable(ax[0,0])
	cax = divider.append_axes('right', size='5%', pad=0.05)
	cbar = f.colorbar(col, use_gridspec=True,cax=cax)
	cbar.set_label('Flux Density [Jy/beam]')

	ax[0,1].set_title('{} GHz original'.format(freq2))
	col = ax[0,1].imshow(file2_plt,cmap=im_color,norm=colors.SymLogNorm(linthresh=level20,linscale=0.5,vmin=lev1[0],vmax=0.5*ma.amax(file2_plt)),extent=extent2,origin='lower')
	divider = make_axes_locatable(ax[0,1])
	cax = divider.append_axes('right', size='5%', pad=0.05)
	cbar = f.colorbar(col, use_gridspec=True,cax=cax)
	cbar.set_label('Flux Density [Jy/beam]')
		
	ax[1,0].set_title('{} GHz regrid & blur'.format(freq1))
	col = ax[1,0].imshow(file1regridblur_plt,cmap=im_color,norm=colors.SymLogNorm(linthresh=level10_r,linscale=0.5,vmin=lev1_r[0],vmax=0.5*ma.amax(file1regridblur_plt)),extent=extent2,origin='lower')
	divider = make_axes_locatable(ax[1,0])
	cax = divider.append_axes('right', size='5%', pad=0.05)
	cbar = f.colorbar(col, use_gridspec=True,cax=cax)
	cbar.set_label('Flux Density [Jy/beam]')

	ax[1,1].set_title('{} GHz regrid & blur'.format(freq2))
	col = ax[1,1].imshow(file2regridblur_plt,cmap=im_color,norm=colors.SymLogNorm(linthresh=level20_r,linscale=0.5,vmin=lev2_r[0],vmax=0.5*ma.amax(file2regridblur_plt)),extent=extent2,origin='lower')
	divider = make_axes_locatable(ax[1,1])
	cax = divider.append_axes('right', size='5%', pad=0.05)
	cbar = f.colorbar(col, use_gridspec=True,cax=cax)
	cbar.set_label('Flux Density [Jy/beam]')

	for aa in ax.flat:
		aa.set(xlabel='Right Ascension  [mas]', ylabel='Relative DEC [mas]')
		aa.axis(axe_ratio)
		aa.set_xlim(ra_max,ra_min)
		aa.set_ylim(dec_min,dec_max)

	f.subplots_adjust(hspace=0.3,wspace=0.3)
	f.savefig('{0}GHz_convolved_with_{1}GHz.pdf'.format(freq2,freq1),bbox_inches='tight')
	
######################################################	
	
	plt.cla()
	widths = [float(file1.shape[1]) for i in range(2)]
	heights = [dec_max/float(ra_max)*widths[i] for i in range(2)]

	fig_width = 10.
	fig_height = fig_width * sum(heights) / sum(widths)
	ratios = [fh/heights[0] for fh in heights]
	f,ax = plt.subplots(2,2,figsize=(fig_width,fig_height))

	ax[0,0].set_title('{}GHz regrid & blur'.format(freq1))
	col = ax[0,0].imshow(file1,cmap=im_color,norm=colors.SymLogNorm(linthresh=level10_r,linscale=0.5,vmin=ma.amin(file1),vmax=ma.amax(file1)),extent=extent2,origin='lower')
	divider = make_axes_locatable(ax[0,0])
	cax = divider.append_axes('right', size='5%', pad=0.05)
	cbar = f.colorbar(col, use_gridspec=True,cax=cax)
	cbar.set_label('Flux Density [Jy/beam]')
	#
	ax[0,1].set_title('{}GHz regrid & blur'.format(freq2))
	col=ax[0,1].imshow(file2,cmap=im_color,norm=colors.SymLogNorm(linthresh=level20_r,linscale=0.5,vmin=ma.amin(file2),vmax=ma.amax(file2)),extent=extent2,origin='lower')
	divider = make_axes_locatable(ax[0,1])
	cax = divider.append_axes('right', size='5%', pad=0.05)
	cbar = f.colorbar(col, use_gridspec=True,cax=cax)
	cbar.set_label('Flux Density [Jy/beam]')
	#
	ax[1,0].set_title('{0}GHz/ {1}GHz not shifted'.format(freq1,freq2))
	cntr1=ax[1,0].contour(file1regridblur_plt,linewidths=0.5,levels=lev1_r,colors='grey',extent=extent2,alpha=1)
	cntr2=ax[1,0].contour(file2regridblur_plt,linewidths=0.5,levels=lev2_r,colors='darkblue',extent=extent2,alpha=0.6)
	h1,_ = cntr1.legend_elements()
	h2,_ = cntr2.legend_elements()
	ax[1,0].legend([h1[0],h2[0]],['{}GHz'.format(freq1),'{}GHz'.format(freq2)],loc='upper right')
	#
	ax[1,1].set_title('{0}GHz/ {1}GHz shifted'.format(freq1,freq2))
	file2subtr=file2regridblur_shift-file2_shift['align']
	cntr1=ax[1,1].contour(file1regridblur_plt,linewidths=0.5,levels=lev1_r,colors='grey',extent=extent2)
	cntr2=ax[1,1].contour(file2subtr,linewidths=0.5,levels=lev2_r,colors='darkblue',extent=extent2,alpha=0.6)
	cntr3=ax[1,1].contour(file2_shift['align'],linewidths=0.5,levels=lev2_r,colors='firebrick',extent=extent2,alpha=0.6)
	h1,_ = cntr1.legend_elements()
	h2,_ = cntr2.legend_elements()
	h3,_ = cntr3.legend_elements()
	ax[1,1].legend([h1[0],h2[0],h3[0]],['{}GHz'.format(freq1),'{}GHz excluded'.format(freq2),'{}GHz included'.format(freq2)],loc='upper right')

	for aa in ax.flat:
		aa.set(xlabel='Right Ascension [mas]', ylabel='Relative Declination [mas]')
		aa.axis(axe_ratio)
		aa.set_xlim(ra_max,ra_min)
		aa.set_ylim(dec_min,dec_max)

	f.subplots_adjust(hspace=0.3,wspace=0.3)

	f.savefig('shifted_maps_{0}GHz_{1}GHz.pdf'.format(freq1,freq2),bbox_inches='tight')
	plt.cla()
	shift_export=file2_shift['shift'].copy()
	sys.stdout.write('final shift: {}'.format(shift_export))
	shift_export[0]*=px_inc_dec[1]*3.6e6
	shift_export[1]*=px_inc_ra[1]*3.6e6
	sys.stdout.write('shift in mas: {}'.format(shift_export))
#########################
	# plot spix map
	spix1 = file1regridblur_plt*(file1regridblur_plt > noise1*sigma) #replaces indices where condition is not met with 0
	spix2 = file2regridblur_shift*(file2regridblur_shift > noise2*sigma)
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

	f,ax = plt.subplots()
	
	cset = ax.contour(spix1,linewidths=[0.5],levels=lev1_r,colors=['grey'], extent=extent2,origin='lower',alpha=0.7)
	im = ax.imshow(a,cmap='hot_r',origin='bottom',extent= extent2,vmin=spix_vmin,vmax=spix_vmax)
	divider = make_axes_locatable(ax)
	cax = divider.append_axes('right', size='5%', pad=0.05)
	cbar = f.colorbar(im, use_gridspec=True,cax=cax)
	cbar.set_label(r'Spectral index $\alpha$', fontsize=12)
	h1,_ = cset.legend_elements()
	ax.legend([h1[0]],['{}GHz'.format(freq1)],loc='upper right',fontsize=11)

	ax.axis('scaled')
	ax.set_xlabel('Right Ascension [mas]',fontsize=12)
	ax.set_ylabel('Relative Declination [mas]',fontsize=12)
	ax.set_xlim(ra_max,ra_min)
	ax.set_ylim(dec_min,dec_max)
	ax.minorticks_on()
	#ax.tick_params('both', length=8, width=2, which='major')
	#f.tick_params(axis='both',which='both',direction='in', labelsize=13)

#	cb = plt.colorbar(im,cax=cax,cmap='jet')
#	cb.ax.tick_params(labelsize=13, width=2)
#	cb.set_label(r'$\alpha$',fontsize=15)

	plt.savefig('spectral_index_between_{:.1f}_{:.1f}.pdf'.format(freq1,freq2),bbox_inches='tight')
	plt.close('all')
############################
	return {'file1':file1regridblur,'file2':file2regridblur,'shift':shift_export,'increment_dec':px_inc_dec[1]*3.6e6,'increment_ra':px_inc_ra[1]*3.6e6}

'''do some ploting here for example compare pixel by pixel difference
for original images and the aligned ones'''
