#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import VLBIana.modules.fit_functions as ff 
from VLBIana.modules.plot_functions import *
from VLBIana.modules.jet_calculus import *

plt.ion()
#plt.style.use('talkstyle')
class CleanMap(object):
	'''
	Useage example:
	Map = CleanMap('maps/NGC1052_C.fits')
	Map.plotMap(plot_mod='models/NGC1052_C_model.fits')
	'''
	def __init__(self,mapFile,ccomp=False):
		self.map = mapFile
		self.head = None
		self.comp = None
		self.cmap = None
		self.cmaph = dict()
		self.modFile = None
		self.modh = dict()
		if ccomp:
			self.comp,self.head,self.cmap = read_fits(self.map)
		else:
			with fits.open(mapFile) as hdulist:
				self.head = hdulist[0].header
				img = hdulist[0].data
				if img.ndim == 4:
					self.cmap = img.reshape(img.shape[2],img.shape[3])
				else:
					self.cmap = img

		self.cmaph['noise']	= self.head['NOISE']
		self.cmaph['px_inc']= self.head['CDELT2']
		self.cmaph['naxis']	= self.head['NAXIS1']
		self.cmaph['beam']	= [self.head['BMAJ']*3.6e6,self.head['BMIN']*3.6e6,self.head['BPA']]
		self.cmaph['freq']	= np.around(self.head['crval3']*1e-9,1)
		self.cmaph['ppb']		= [PXPERBEAM(self.cmaph['beam'][0]*np.pi/180,self.cmaph['beam'][1]*np.pi/180,self.cmaph['px_inc']*3.6e6*np.pi/180)]
		self.cmaph['fov']		= self.cmaph['px_inc']*self.cmaph['naxis']*3.6e6

	def modelFile(self, modelFile):
		with fits.open(modelFile) as f:
			self.modFile = f[1].data
		self.modh['x']			= self.modFile['DELTAX']*3.6e6
		self.modh['y']			= self.modFile['DELTAY']*3.6e6
		self.modh['maj']		= self.modFile['MAJOR AX']*3.6e6
		self.modh['min']		= self.modFile['MINOR AX']*3.6e6
		self.modh['posa']		= self.modFile['POSANGLE']
		
	def plotMap(self,sigma=3,fig_size='screen',ra_lim=False,dec_lim=False,saveFile=False,plot_mod=False,plot_cntr=True,plot_color=False,cntr_color=False,model_color=False):

		if not cntr_color:
			if plot_cntr:
				cntr_color = 'grey'
			else:
				cntr_color = 'black'
		####################
		# setting all parameters for plotting a clean image
		#####################
		if not ra_lim:
			ra  = self.cmaph['fov']/3.
			dec = self.cmaph['fov']/5.
		else:
			ra	= ra_lim
			dec	= dec_lim

		scale	= -self.cmaph['px_inc']*3.6e6
		sigma	= sigma
		level0= self.cmaph['noise']*sigma
		lev=[]
		for i in range(0,10):
			lev.append(level0*2**i)
		xx	= np.linspace(-self.cmaph['naxis']*0.5*scale,(self.cmaph['naxis']*0.5-1)*scale,self.cmaph['naxis'])
		yy	= np.linspace(self.cmaph['naxis']*0.5*scale,-(self.cmaph['naxis']*0.5-1)*scale,self.cmaph['naxis'])
		vmax=0.5*ma.amax(self.cmap)
		norm = mpl.colors.SymLogNorm(linthresh=level0,linscale=0.5,vmin=level0,vmax=vmax,base=np.e)

		##################
		# Plotting
		###################
		f,ax = plt.subplots()
		ax.axis('scaled')
		ax.set_xlim(-ra,ra)
		ax.set_ylim(-dec,dec)
		ax.invert_xaxis()
		plotBeam(self.cmaph['beam'][1],self.cmaph['beam'][0],self.cmaph['beam'][2],ra,dec,ax)
		if plot_cntr:
			cntr=ax.contour(xx,yy,self.cmap,linewidths=0.5,levels=lev,colors=cntr_color,alpha=1)
		if plot_color:
			extent = np.max(xx),np.min(xx),np.min(yy),np.max(yy)
			im = ax.imshow(self.cmap,cmap=colormap,extent=extent,origin='lower')
			im.set_norm(norm)

		#######
		# plot model if wanted 
		#######
		if plot_mod:
			if not model_color:
				model_color='red'
			self.modelFile(plot_mod)
			Mx = self.modh['x']
			My = self.modh['y']
			Mposa = self.modh['posa']
			Mmaj	= self.modh['maj']
			Mmin	= self.modh['min']
			for j,xx in enumerate(Mx):
				e_comp = Ellipse([Mx[j],My[j]],Mmaj[j],Mmin[j],-Mposa[j]+90, color=model_color, zorder=2, fill=False,lw=0.5)
				ax.add_artist(e_comp)
				maj1_x = Mx[j]-np.sin(-np.pi/180*Mposa[j])*Mmaj[j]*0.5
				maj1_y = My[j]+np.cos(-np.pi/180*Mposa[j])*Mmaj[j]*0.5
				maj2_x = Mx[j]+np.sin(-np.pi/180*Mposa[j])*Mmaj[j]*0.5
				maj2_y = My[j]-np.cos(-np.pi/180*Mposa[j])*Mmaj[j]*0.5
	
				min1_x = Mx[j]-np.sin(-np.pi/180*(Mposa[j]+90))*Mmin[j]*0.5
				min1_y = My[j]+np.cos(-np.pi/180*(Mposa[j]+90))*Mmin[j]*0.5
				min2_x = Mx[j]+np.sin(-np.pi/180*(Mposa[j]+90))*Mmin[j]*0.5
				min2_y = My[j]-np.cos(-np.pi/180*(Mposa[j]+90))*Mmin[j]*0.5
				if maj1_y==maj2_y:
					ax.plot(maj1_x,maj1_y, color = model_color,marker='+', markersize=5,lw=0.5)
				else:
					ax.plot([maj1_x,maj2_x],[maj1_y,maj2_y], color = model_color, lw = 1)
					ax.plot([min1_x,min2_x],[min1_y,min2_y], color = model_color, lw = 1)

		# set axis, labels, etc.
		ax.set(xlabel='RA [mas]', ylabel='Relative DEC [mas]')
		ax.minorticks_on()
		ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
		ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

		if saveFile:
			figsize=set_size(fig_size)
			set_corrected_size(f,figsize)
			plt.savefig(saveFile,bbox_inches='tight')
		else:
			plt.show()
		return f,ax
