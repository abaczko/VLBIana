#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
To plot some diagnostic plots and fit a Powerlaw or borkenPowerlaw to the ridgelines
'''
from ridgeline.plot_ridgeline import *
from ridgeline.fit_width import *
from ridgeline.scattering import *

mapFile		= ['xx','xx']
ridgeLine = ['xx','xx']
errorFile	= ['xx','xx']
logFile		= ['xx','xx']

label =['freq1','fre2']

shiftFile = 'xx.txt' # file containing shifts to align both ridgelines
saveFile = 'Ridglines' #files will start with that string
theta	= 24	#rotation angle used when estimating the ridge-line
inc		= 80	# inclination of the source to automatically derive the de-projected distances

ridgeline_plotter(mapFile,ridgeLine,shiftFile,saveFile,label,theta,logFile,plot_rl=True,incl=inc)
ridgeline_plotter(mapFile,ridgeLine,shiftFile,saveFile,label,theta,logFile,plot_width=True,plot_deconvolved=True,errorbars=True,incl=inc)
ridgeline_plotter(mapFile,ridgeLine,shiftFile,saveFile,label,theta,logFile,plot_width=True,plot_reslimit=False,errorbars=True,incl=inc)
ridgeline_plotter(mapFile,ridgeLine,shiftFile,saveFile,label,theta,logFile,plot_map=True,incl=inc,plot_stacked=True)
rl_fit(mapFile,ridgeLine,shiftFile,saveFile,label,theta,logFile,fit='Powerlaw',plot_res=True,incl=inc,plot_flux=True)
