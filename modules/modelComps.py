#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from astropy.table import Table,vstack
from astropy.io import fits
from astropy.time import Time
from glob import glob
from jet_calculus import *
import matplotlib.pyplot as plt

plt.ion()
cm = plt.get_cmap('gist_rainbow')

def trim_axs(axs, N):
  """little helper to massage the axs list to have correct length..."""
  axs = axs.flat
  for ax in axs[N:]:
    ax.remove()
  return axs[:N]
def keyfunc(s):
  return [int(''.join(g)) if k else ''.join(g) for k, g in groupby(s, str.isdigit)]

class modelComp(object):
    '''
    To read and plot modelfit components
    If you need the initial moFiles array to be sorted: try:
    modFiles = sorted(modFiles,key-keyfunc)
    It should sort as you typically expect, e.g. 0,1,2,3 or a,b,c

    Instruction:
    modelFile= glob('model/*.fits')
    modF = modelComp(modelFile)
    modF.plot_epoch_xy()

    For sorting by id:
    modF.sort_by_id()

    Then you can use modF.model_sorted for plotting by id
    e.g. as is done in modF.plot_comp_xy(xax='date',yax='FLUX')
    '''
    def __init__(self,modFiles,cleanFiles=False,shift=False):
        self.modFiles =modFiles
        self.model  =   dict()
        self.cchead = dict()
        self.shift  = None

        for modFile in modFiles:
            with fits.open(modFile) as f:
                h = f[0].header
                key = h['OBJECT']+'.'+str(h['CRVAL3']/1e9)+'.'+h['DATE-OBS']
                self.model[key] = dict()
                self.model[key]['data_table'] = f[1].data

            self.model[key]['head'] = h
            self.model[key]['data'] = Table(self.model[key]['data_table'])
            self.model[key]['date'] = Time(h['DATE-OBS'])
            self.model[key]['date'].format = 'mjd'
            self.model[key]['date'] = self.model[key]['date'].value

            self.model[key]['freq'] = h['CRVAL3']/1e9
            self.model[key]['source'] = h['OBJECT']
            self.model[key]['beam'] = [h['BMAJ'],h['BMIN'],h['BPA']]
            self.model[key]['noise'] = h['NOISE']
            self.model[key]['px_inc'] = h['CDELT2']
            self.model[key]['data']['DELTAX'] *= 3.6e6
            self.model[key]['data']['DELTAY'] *= 3.6e6
            self.model[key]['data']['DIST'] = np.sign(self.model[key]['data']['DELTAX'])*np.sqrt(self.model[key]['data']['DELTAX']**2 + self.model[key]['data']['DELTAY']**2 )
            self.model[key]['data']['id'] = np.arange(1,len(self.model[key]['data'])+1,1)

        if cleanFiles:
            for cleanFile in cleanFiles:
                with fits.open(cleanFile) as f:
                    h = f[0].header
                key = h['OBJECT']+'.'+str(h['CRVAL3']/1e9)+'.'+h['DATE-OBS']
                self.cchead[key] = dict()
                self.cchead[key]['head'] = h
                self.cchead[key]['beam'] = [h['BMAJ'],h['BMIN'],h['BPA']]
                self.cchead[key]['noise'] =h['NOISE']
                self.cchead[key]['px_inc'] = h['CDELT2']

        if shift:
            self.shift = Table.read(shiftFile, format='ascii')
            self.model[key]['data_shifted'] = self.model[key]['data'].copy()
            self.model[key]['data_shifted']['DELTAX'] += self.shift['RA']
            self.model[key]['data_shifted']['DELTAY'] += self.shift['DEC']
            self.model[key]['data_shifted']['DIST'] += np.sqrt(self.shift['RA']**2 + self.shift['DEC']**2)

    def sort_by_id(self):
        ids = np.unique(vstack([self.model[key]['data']['id'] for key in self.model])['id'])

        self.model_sorted = dict()
        for ID in ids:
            date = []
            for i,key in enumerate(self.model):
                if i ==0:
                    self.model_sorted[ID] = Table(dtype=self.model[key]['data'])
                row = self.model[key]['data'][self.model[key]['data']['id']==ID]
                if len(row) > 0:
                    self.model_sorted[ID].add_row(row[0])
                    date.append(self.model[key]['date'])
            self.model_sorted[ID]['date']= date

###
    def plot_comp_xy(self,xax='DIST',yax='FLUX',out=False):
        NUM_COLORS = len(self.model_sorted)
        linecolors = [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)]
        symbols = ['o','v','*','P','s','D'] * np.int(np.ceil(NUM_COLORS/6))

        fig,ax = plt.subplots(figsize=(12,8))
        for i,comp in enumerate(self.model_sorted):
            color = linecolors[i]
            symbol = symbols[i]
            xx = self.model_sorted[comp][xax]
            yy = self.model_sorted[comp][yax]

            ax.scatter(xx,yy,marker=symbol, color = color, label = comp)

        ax.set(xlabel=xax, ylabel=yax)
        fig.legend(loc=7, fontsize='small')
        fig.tight_layout()
        fig.subplots_adjust(right=0.9, top=0.98)

        if out:
            if type(out)==bool:
                outf = self.model[comp]['source']+'_components_all.pdf'
            fig.savefig(outf,bbox_inches='tight')
###

    def plot_epoch_xy(self,xax='DIST',yax='FLUX',plot_all=True,out=True,ylog=False,xlog=False):
        if plot_all:
            nn = len(self.model)
            nx = 4
            ny = int(np.ceil(nn/nx))
        else:
            nx = 1
            ny = 1

        fig,axs = plt.subplots(ny,nx, figsize=(12,8))
        axs   = trim_axs(axs,len(self.model))

        if plot_all:
            for ax,modF in zip(axs,self.model):
                xx = self.model[modF]['data'][xax]
                yy = self.model[modF]['data'][yax]
                if ylog:
                    ax.set_yscale('log')
                if xlog:
                    ax.set_xscale('log')

                ax.scatter(xx,yy)
                ax.minorticks_on()
                ax.set(xlabel=xax,ylabel=yax)
                ax.annotate('{} - {:.1f} GHz - {}'.format(self.model[modF]['source'],self.model[modF]['freq'],self.model[modF]['date']),xy=(0.03,0.9),xycoords='axes fraction',size=8)

            plt.tight_layout(pad=0.2,w_pad=0.2,h_pad=0.2)

            if out:
                if type(out)==bool:
                    outf = self.model[modF]['source']+'_components.pdf'
                fig.savefig(outf,bbox_inches='tight')
            plt.cla()


    def plot_all(self,xax='date',yax='DIST',out=False):
        fig,ax = plt.subplots(figsize=(12,8))
        dates = []
        ymin = 0
        ymax = 0
        for modF in self.model:
            dates.append(self.model[modF]['date'].value)
            if np.min(self.model[modF]['data'][yax]) < ymax :
                ymin = np.min(self.model[modF]['data'][yax])
            if np.max(self.model[modF]['data'][yax]) > ymax :
                ymax = np.max(self.model[modF]['data'][yax])
        xmin = np.min(dates) - 0.5
        xmax = np.max(dates) + 0.5
        ax.axis([xmin,xmax,ymax,ymin])
        print([xmin,xmax,ymin,ymax])
        print(dates)
        for i,modF in enumerate(self.model):
            ax.scatter([dates[i]]*len(self.model[modF]['data'][yax]),self.model[modF]['data'][yax])
        ax.minorticks_on()
        if xax=='date':
            ax.set(xlabel='Date [MJD]',ylabel='Distance [mas]')
        else:
            ax.set(xlabel=xax,ylabel='Distance [mas]')
        plt.tight_layout(pad=0.4,w_pad=0.5,h_pad=1.0)
        if out:
            if type(out)==bool:
                outf = self.model[modF]['source']+'_all.pdf'
            fig.savefig(outf,bbox_inches='tight')
        plt.cla()
