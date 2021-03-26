#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

from jet_calculus import *

class CleanMap(object):
	def __init__(self,mapFile,ccomp=False):
		self.map = mapFile
		self.head = None
		self.comp = None
		self.cmap = None
		self.cmaph = dict()
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

	def readHead(self):
		self.cmaph['noise']	= self.head['NOISE']
		self.cmaph['px_inc']	= self.head['CDELT2']
		self.cmaph['naxis']	= self.head['NAXIS1']
		self.cmaph['beam']		= [self.head['BMAJ']*3.6e6,self.head['BMIN']*3.6e6,self.head['BPA']]
		self.cmaph['freq']		= np.around(self.head['crval3']*1e-9,1)
		self.cmaph['ppb']		= [PXPERBEAM(self.cmaph['beam'][0]*np.pi/180,self.cmaph['beam'][1]*np.pi/180,self.cmaph['px_inc']*3.6e6*np.pi/180)]
		self.cmaph['fov']		= self.cmaph['px_inc']*self.cmaph['naxis']*3.6e6
		return self.cmaph


