#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from astropy.table import Table
from astropy.io import ascii

import sys

def read_slfit_log(logFile,outfile=False,writeTable=False):
	'''
	Read logfiles written by aips task SLFIT and saves results to an astropy table
	'''
	lines=[]
	isl = [int(logFile.split('.')[0].split('_')[-1])]
	if not outfile and writeTable:
		outfile = logfile.split('.')[0]
	with open(logFile,'r') as f:
		for line in f:
			if line.find('LOGFILE') == -1:
				lines.append(line)

	i=0
	peak,peak_err,peak_unit,pos,pos_err,pos_unit,fwhm,fwhm_err,fwhm_unit=[],[],[],[],[],[],[],[],[]
	for line in lines:
		if line.find('Results in physical units')!=-1:
			peak.append(float(lines[i+2].split()[3]))
			peak_err.append(float(lines[i+2].split()[4]))
			peak_unit.append(lines[i+2].split()[5])
			pos.append(float(lines[i+3].split()[3]))
			pos_err.append(float(lines[i+3].split()[4]))
			pos_unit.append(lines[i+3].split()[5])
			fwhm.append(float(lines[i+4].split()[3]))
			fwhm_err.append(float(lines[i+4].split()[4]))
			fwhm_unit.append(lines[i+4].split()[5])
			if peak_unit[-1] == 'MicroJY/BEAM':
				peak_err[-1]*=1e-3
				peak[-1]*=1e-3
		#	elif peak_unit[-1] == 'MillJy/BEAM':
		#		peak_err[-1]*=1e-3
		#		peak[-1]*=1e-3

		i+=1
	slfit = Table([isl,peak,peak_err,pos,pos_err,fwhm,fwhm_err,peak_unit,pos_unit,fwhm_unit],names=['isl','peak','peak_err','pos','pos_err','fwhm','fwhm_err','peak_unit','pos_unit','fwhm_unit'])
	if writeTable:
		ascii.write(slfit, outfile+'.dat', overwrite=True)
	return slfit

def log_to_table(logFiles,outfile=False,writeTable=True):
	'''
	Read multiple logfiles from SLFIT and save the info to an astropy table
	'''
	if not outfile and writeTable:
		rot = logFiles[0].split('/')[-2].split('_')[-1]
		outfile = '_'.join(logFiles[0].split('/')[-1].split('.')[0].split('_')[:-1])+'_'+rot

	SLfits = Table(names=['isl','peak','peak_err','pos','pos_err','fwhm','fwhm_err','peak_unit','pos_unit','fwhm_unit'],dtype=['i4','f4','f4','f4','f4','f4','f4','U25','U25','U25'])

	for lf in logFiles:
		slfit = read_slfit_log(lf)
		SLfits.add_row(slfit[0])

	indx = SLfits['isl'].argsort()
	SLfits = SLfits[indx]
	if writeTable:
		ascii.write(SLfits, outfile+'.dat', overwrite=True)

	return SLfits
