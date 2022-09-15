#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from astropy.table import Table
import sys,warnings,statistics,collections
#
from VLBIana.modules.cleanMap import *

def readShift(shiftFile):
    relative_shifts = Table.read(shiftFile, format='ascii')
    return relative_shifts['RA'], relative_shifts['DEC']

class RidgeLine(object):
    def __init__(self,ridgeFile,logFile,mapFile,shift=None,incl=90):
        self.ridgeline= ridgeFile
        self.log            = logFile
        self.map            = mapFile
        self.cmap           = None
        self.peak_err   = []
        self.pos_err    = []
        self.fwhm_err   = []
        self.log_isl    = []
        self.RL             = None
        self.RLbinned = None
        self.shift      = shift
        self.i              = incl*np.pi/180

        print('Applying an inclination of the jets of i={} rad'.format(self.i))

    def readMap(self):
        Map = CleanMap(self.map)
        self.cmap = Map.cmaph
        return Map

    def readLog(self):
        '''
        To get the errors on fit parameters from SLFIT
        '''
        logFile = Table.read(self.log,format='ascii')
        self.peak_err = logFile['peak_err']
        self.fwhm_err = logFile['fwhm_err']
        self.pos_err    = logFile['pos_err']
        self.log_isl    = logFile['isl']
        return self.peak_err,self.fwhm_err,self.pos_err#,self.log_isl

    def ReadErrors(self,errFile):
        width_errors = Table.read(errFile,format='ascii',data_start=1)
        return width_errors

    def readRidgeline(self,theta=0,widthErr=False,add_image_err=10.):
        '''
        Read the ridgeline files and return a table with all values
        give theta in rad
        '''
        logFile     = self.readLog()
        self.RL = Table.read(self.ridgeline,format='ascii',data_start=1,names=['isl','m','n','Xpos','Ypos','Dist','Peak','FWHM'])
        self.readMap()
        # convert properties to the required units and derive uncertainties

        self.RL['FWHM'] = np.abs(self.RL['FWHM'])*1e3

        self.RL['fomalont'] = self.RL['FWHM']/(self.RL['Peak']/self.cmap['noise'])
        # change numbers again
        self.RL['Peak'] *= 1e3
        self.RL['Xpos'] *= 1e3
        self.RL['Ypos'] *= 1e3
        self.RL['Freq'] = [self.cmap['freq']]*len(self.RL['Xpos'])
        self.RL['YposErr']      = self.pos_err
        self.RL['FitErr']           = self.fwhm_err
        self.RL['PeakFitErr']   = self.peak_err
        if widthErr:
            sys.stdout.write('Use stddev calculated from changing the rotation angle.\n')
            self.RL['FWHMErr']  = self.ReadErrors(widthErr)['RL_stddev']
        else:
            self.RL['FWHMErr'] = self.RL['FitErr']
        # removing some values
        mask        = np.logical_and.reduce((self.RL['Peak']> 5*self.cmap['noise'],self.RL['FWHM']>self.cmap['beam'][0],self.RL['FWHM']>self.RL['FWHMErr']))
        nn1 = len(self.RL['Xpos'])
        self.RL = self.RL[mask]
        nn2 = len(self.RL['Xpos'])

        if nn1 != nn2:
            print('WARNING: '+str(nn1-nn2)+' Entries have been deleted (either flux is below 3*noise or fwhm>beam')
        # continue filling the table
        self.RL['Dist']             = np.sign(self.RL['Xpos'])*np.sqrt(self.RL['Xpos']**2+self.RL['Ypos']**2)/np.sin(self.i)
        self.RL['DistErr']      = np.abs(self.RL['Ypos']/np.sqrt(self.RL['Xpos']**2+self.RL['Ypos']**2)*self.RL['YposErr'])
        self.RL['FWHMDeconvolved']      = np.sqrt(self.RL['FWHM']**2-self.cmap['beam'][0]**2)
        if add_image_err:
            self.RL['FWHMDeconvolvedErr']   = ((self.RL['FWHM']/self.RL['FWHMDeconvolved'])*self.RL['FWHMErr'])+self.cmap['beam'][0]/add_image_err
        else:
            self.RL['FWHMDeconvolvedErr']   = ((self.RL['FWHM']/self.RL['FWHMDeconvolved'])*self.RL['FWHMErr'])

        self.RL['RA']           = self.RL['Xpos']*np.cos(theta)-self.RL['Ypos']*np.sin(theta)
        self.RL['Dec']      = self.RL['Ypos']*np.cos(theta)+self.RL['Xpos']*np.sin(theta)

        self.RL['RAErr']    = np.abs(self.RL['YposErr']*np.sin(theta))
        self.RL['DecErr']   = np.abs(self.RL['YposErr']*np.cos(theta))
        if self.shift is not None:
            print('shifting data\n')
            self.RL['RAshift']  = self.RL['RA']+self.shift[0]
            self.RL['Decshift'] = self.RL['Dec']-self.shift[1]
            self.RL['RAshiftErr'] = self.RL['RAErr']
            self.RL['DecshiftErr'] = self.RL['DecErr']
            self.RL['Distshift']= np.sign(self.RL['RAshift'])*np.sqrt(self.RL['RAshift']**2+self.RL['Decshift']**2)
            self.RL['DistshiftErr']= self.RL['DistErr']
        return self.RL

    def binRidgeLine(self,rlbin=1/4.):
        ''' To bin the ridgeline into distance bins equally to rlbin x beam '''
        bins = np.arange(min(self.RL['Dist']),max(self.RL['Dist']),self.cmap['beam'][0]*rlbin)
        self.RLbinned = Table(dtype=self.RL.dtype)
        for i,dd in enumerate(bins,start=1):
            rl = self.RL[np.logical_and(self.RL['Dist'] >= bins[i-1], self.RL['Dist'] < bins[i])]
            if len(rl['RA'])>1:
                self.RLbinned.add_row((i,1,1,np.mean(rl['Xpos']),
                    np.mean(rl['Ypos']),
                    np.mean(rl['Dist']),
                    np.mean(rl['Peak']),
                    np.median(rl['FWHM']),
                    np.median(rl['FWHM'])/(np.mean(rl['Peak'])/self.cmap['noise']), #fomalont
                    np.mean(rl['Freq']),
                    statistics.stdev(rl['Ypos']),
                    statistics.stdev(rl['FWHM']),
                    statistics.stdev(rl['Peak']),
                    (np.median(rl['FWHM'])/(np.mean(rl['Peak'])/self.cmap['noise']))+(self.cmap['beam'][0]/5.), #fwhm error
                    statistics.stdev(rl['Dist']),
                    np.median(rl['FWHMDeconvolved']),
                    self.cmap['beam'][0]/10.,#(np.median(rl['FWHM'])/np.median(rl['FWHMDeconvolved']))*np.median(rl['FWHM'])/(np.mean(rl['Peak'])/self.cmap['noise']),#statistics.stdev(rl['FWHMDeconvolved']), #error deconvolved
                    np.mean(rl['RA']),np.mean(rl['Dec']),
                    statistics.stdev(rl['RA']),
                    statistics.stdev(rl['Dec']),
                    np.mean(rl['RAshift']),
                    np.mean(rl['Decshift']),
                    statistics.stdev(rl['RAshift']),
                    statistics.stdev(rl['Decshift']),
                    np.mean(rl['Distshift']),
                    statistics.stdev(rl['Distshift'])))
            elif len(rl['RA'])==1:
                #print('only one value per bin, add original line')
                self.RLbinned.add_row(rl[0])
            else:
                print('No value in this bin')
                pass
            if i == len(bins)-1:
                break
        return self.RLbinned
