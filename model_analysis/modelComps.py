#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from astropy.table import Table,vstack
from astropy.io import fits,ascii
from astropy.time import Time
import scipy as sp
from glob import glob
import sys,os
from itertools import cycle,chain,compress
from VLBIana.modules.jet_calculus import *
from VLBIana.modules.plot_functions import *
import VLBIana.modules.fit_functions as ff
import ehtim as eh
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
plt.ion()

def trim_axs(axs, N):
  """little helper to massage the axs list to have correct length..."""
  axs = axs.flat
  for ax in axs[N:]:
    ax.remove()

  return axs[:N]
def keyfunc(s):
  return [int(''.join(g)) if k else ''.join(g) for k, g in groupby(s, str.isdigit)]

def apply_shift(img,shift):
    offset_image = fourier_shift(np.fft.fftn(img), shift)
    imgalign = np.fft.ifftn(offset_image)
    img2 = imgalign.real

    return img2

class modelComp(object):
    '''
    To read and plot modelfit components
    If you need the initial moFiles array to be sorted: try:
    modFiles = sorted(modFiles,keyfunc)
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
    def __init__(self,modFiles,cleanFiles=None,shift=None,colormap=None,z=None):
        sys.stdout.write('initiate modelComp class.\n')
        self.modFiles =modFiles
        self.model  =   dict()
        self.model_sorted = None
        self.cchead = dict()
        self.ccmap  = dict()
        self.ccmap_shifted = dict()
        self.keys   = []
        self.ids    = []
        self.eps    = [chr(value) for value in range(ord('a'),ord('a')+len(self.modFiles))]
        self.z      = z

        if self.z==None:
            sys.stdout.write('Please provide a redshift if Tb should be calculated, giving the parameter of z when calling the modelComp class.\n')
        else:
            sys.stdout.write('Redshift provided, hence Tb will be calculated.\n')

        if colormap:
            self.colormap = colormap
        else:
            self.colormap = 'inferno'
        self.id_colors = []
        self.symbols = []

        for ii,modFile in enumerate(modFiles):
            with fits.open(modFile) as f:
                h = f[0].header
                key = h['OBJECT']+'.'+str(h['CRVAL3']/1e9)+'.'+h['DATE-OBS']
                self.keys.append(key)
                self.model[key] = dict()
                self.model[key]['data_table'] = f[1].data

            self.model[key]['head'] = h
            self.model[key]['data'] = Table(self.model[key]['data_table'])
            self.model[key]['date'] = Time(h['DATE-OBS'])
            self.model[key]['date_obs'] =  self.model[key]['date'].value.split()[0]
            self.model[key]['date'].format = 'mjd'
            self.model[key]['date'] = self.model[key]['date'].value

            self.model[key]['freq'] = np.around(h['CRVAL3']/1e9,2)
            self.model[key]['source'] = h['OBJECT']
            self.model[key]['beam'] = [h['BMAJ'],h['BMIN'],h['BPA']]
            self.model[key]['noise'] = h['NOISE']
            self.model[key]['px_inc'] = h['CDELT2']
            self.model[key]['data']['DELTAX'] *= 3.6e6
            self.model[key]['data']['DELTAY'] *= 3.6e6
            self.model[key]['data']['MAJOR AX'] *= 3.6e6
            self.model[key]['data']['MINOR AX'] *= 3.6e6
            self.model[key]['data']['ratio'] = self.model[key]['data']['MINOR AX'] / self.model[key]['data']['MAJOR AX']
            self.model[key]['data']['DIST'] = np.sign(self.model[key]['data']['DELTAX'])*np.sqrt(self.model[key]['data']['DELTAX']**2 + self.model[key]['data']['DELTAY']**2 )
            self.model[key]['data'].add_column(np.arange(1,len(self.model[key]['data'])+1,1).astype(str)
, name='id', index=0)
            nid = len(self.model[key]['data'])
            self.model[key]['data']['ep'] = nid*[self.eps[ii]]
            self.model[key]['data']['freq'] = nid*[self.model[key]['freq']]
            self.model[key]['data']['cid'] = mpl.colors.to_hex('red')
#include tb here
            if self.z:
                self.model[key]['data']['tb'] = derive_tb(self.model[key]['data']['FLUX'],self.model[key]['freq'],self.z,self.model[key]['data']['MAJOR AX'],self.model[key]['data']['ratio'])
                for ii,value in enumerate(self.model[key]['data']['tb']):
                    if value>1e13:
                        self.model[key]['data']['tb'] [ii] = np.nan
                self.model[key]['data']['logtb'] = np.log10(self.model[key]['data']['tb'])

    #generate a list of the all components
        [self.ids.extend(self.model[key]['data']['id']) for key in self.keys]
        self.ids = np.unique(self.ids)
        self.ids = np.array(sorted(self.ids, key=keyfunc))
        # generate a standard color and symbol for each id
        NUM_COLORS = len(self.ids)
        cmap = mpl.colormaps.get_cmap(self.colormap)
        colors = cmap(np.linspace(0,0.85,NUM_COLORS))

        #cmap = cm.get_cmap(colormap,NUM_COLORS)
        for i in range(len(colors)):
            self.id_colors.append(mpl.colors.rgb2hex(colors[i]))
        self.id_colors = np.array(self.id_colors)
        self.symbols = np.array(['o','v','*','P','s','D'] * int(np.ceil(NUM_COLORS/6)))

        for key in self.keys:
            for ii,model in enumerate(self.model[key]['data']):
                mask = self.ids==model['id']
                model['cid'] = self.id_colors[mask][0]

        if cleanFiles:
            sys.stdout.write('Load clean fits files to class.\n')
            for i,cleanFile in enumerate(cleanFiles):
                with fits.open(cleanFile) as f:
                    h = f[0].header
                    img = f[0].data
                    if img.ndim == 4:
                        img = img.reshape(img.shape[2],img.shape[3])

                key = self.keys[i]
                self.cchead[key] = dict()
                self.cchead[key]['head']    = h
                self.cchead[key]['beam']    = [h['BMAJ']*3.6e6,h['BMIN']*3.6e6,h['BPA']]
                self.cchead[key]['noise']   = h['NOISE']
                self.cchead[key]['px_inc']  = h['CDELT2']
                self.cchead[key]['naxis']   = h['NAXIS1']
                self.cchead[key]['freq']    = np.around(h['CRVAL3']/1e9,2)
                self.cchead[key]['fov']     = self.cchead[key]['px_inc']*self.cchead[key]['naxis']*3.6e6 
                self.ccmap[key] = img
 
        if shift:
            sys.stdout.write('shifting model data according to values in shiftFIle.\n')
            shiftFile= shift
            self.shift = Table.read(shiftFile, format='ascii')
            for i,shift in enumerate(self.shift):
                key =self.keys[i]
                self.model[key]['data_shifted'] = self.model[key]['data'].copy()
                self.model[key]['data_shifted']['DELTAX'] += shift['RA']
                self.model[key]['data_shifted']['DELTAY'] += shift['DEC']
                self.model[key]['data_shifted']['DIST'] += np.sqrt(shift['RA']**2 + shift['DEC']**2)
                inc=self.cchead[key]['px_inc']*3.6e6
                self.ccmap_shifted[key] = apply_shift(self.ccmap[key],[shift['DEC']/inc,-shift['RA']/inc])

    def update_ids(self):
        sys.stdout.write('Updatend component id list.\n')
        self.ids = []
        [self.ids.extend(self.model[key]['data']['id']) for key in self.keys]
        self.ids = np.unique(self.ids)
        self.ids = np.array(sorted(self.ids, key=keyfunc))

    def update_cm(self,colormap=False):
        sys.stdout.write('Updating colormap.\n')
        if colormap==False:
            colormap=self.colormap
        else:
            self.colormap = colormap
        self.update_ids()
        NUM_COLORS = len(self.ids)
        colors = cmap(np.linspace(0,0.85,NUM_COLORS))
        self.id_colors=[]
        for i in range(len(colors)):
            self.id_colors.append(mpl.colors.rgb2hex(colors[i]))
        self.id_colors = np.array(self.id_colors)
        self.symbols = np.array(['o','v','*','P','s','D'] * int(np.ceil(NUM_COLORS/6)))
        for key in self.keys:
            for ii,model in enumerate(self.model[key]['data']):
                mask = self.ids==model['id']
                model['cid'] = self.id_colors[mask][0]

    def sort_by_id(self):
        sys.stdout.write('Creating class.model_sorted with components sorted by id.\n')
        self.update_ids()

        self.model_sorted = dict()
        for ID in self.ids:
            date = []
            for i,key in enumerate(self.model):
                if i ==0:
                    self.model_sorted[ID] = Table(dtype=self.model[key]['data'])
                row = self.model[key]['data'][self.model[key]['data']['id']==ID]
                if len(row) > 0:
                    self.model_sorted[ID].add_row(row[0])
                    date.append(self.model[key]['date'])
            self.model_sorted[ID]['date']= date
        self.update_cm()

    def change_id(self,old_ids,new_ids):
        '''
        Usage:
            old_ids = ['a_2','b_5','c_6'] is a list of the component ids to be changed. Consisting of the epoch id (ep) and the component id in that epoch (id)
            new_ids = ['a2','a2','b'] a list of names you want to assign to the components. The id of the component in this epoch will then be changed to this value
        '''
        sys.stdout.write('Change model ids according to imput.\n')
        for ii,oids in enumerate(old_ids):
            epoch,cid = oids.split('_')
            key=self.keys[np.where(np.array(self.eps)==epoch)[0][0]]
            jj=np.where(self.model[key]['data']['id']==cid)[0][0]
            self.model[key]['data']['id'][jj]=new_ids[ii]
            if self.model[key]['data_shifted']:
                self.model[key]['data_shifted']['id'][jj]= new_ids[ii]
        self.update_cm()
        self.sort_by_id()


###
    def plot_comp_xy(self,xax='DIST',yax='FLUX',out=False,comps=False,line=False,ccolor=False):
        sys.stdout.write('Plot component {} vs {}\n'.format(xax,yax))
        if not self.model_sorted:
            self.sort_by_id()
        if comps:
            model_sorted = {key: self.model_sorted[key] for key in comps}
        else:
            model_sorted = self.model_sorted.copy()

        if ccolor:
            ccolors= np.array(len(self.ids)*[ccolor])
        else:
            ccolors = self.id_colors

        fig,ax = plt.subplots(figsize=(12,8))
        if xax == 'freq':
            xlabel='Freq [GHz]'
        if yax == 'FLUX':
            ylabel='Flux [Jy]'
        elif yax == 'tb':
            ylabel =r'$\log(T_\mathrm{B})$ [K]'
            #ax.set_yscale('log')
            ax.loglog()
        elif yax =='DIST':
            ylabel = 'Distance [mas]'
        for i,comp in enumerate(model_sorted):
            color = ccolors[self.ids==comp][0]
            symbol = self.symbols[i]
            xx = self.model_sorted[comp][xax]
            yy = self.model_sorted[comp][yax]

            if line:
                ax.plot(xx,yy,marker=symbol, color = color, label = comp,lw=1)
            else:
                ax.scatter(xx,yy,marker=symbol, color = color, label = comp)

        #ax.set(xlabel=xax, ylabel=yax)
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        keys_list = sorted(by_label, key=keyfunc)
        sorted_by_label = {key: by_label[key] for key in keys_list}
        #pos = ax.get_position()
        #ax.set_position([pos.x0, pos.y0, pos.width, pos.height * 0.85])
        #ax.set_position([pos.x0, pos.y0, pos.width * 0.9, pos.height])
        fig.legend(sorted_by_label.values(), sorted_by_label.keys(),loc='upper left', bbox_to_anchor=(0.9,0.98),ncol=1)
        ax.set(xlabel=xlabel, ylabel=ylabel)
        ax.minorticks_on()
        if not yax=='tb':
            ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

       # fig.legend(loc=7, fontsize='small')
        fig.tight_layout()
        fig.subplots_adjust(right=0.9, top=0.98)

        if out:
            if type(out)==bool:
                outf = self.model[self.keys[0]]['source']+'_{}_vs_{}.pdf'.format(xax,yax)
            fig.savefig(outf,bbox_inches='tight')
            sys.stdout.write('Plot has been written to {}\n'.format(outf))
        else:
            plt.show()

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
            fig.subplots_adjust(right=0.95, top=0.98)

            if out:
                if type(out)==bool:
                    outf = self.model[modF]['source']+'_components.pdf'
                fig.savefig(outf,bbox_inches='tight')
                sys.stdout.write('Plot has been written to {}\n'.format(outf))
            else:
              plt.show()
            plt.cla()

    def plot_evolution_map(self,sigma=4,fig_size='screen',ra=False,dec=False,saveFile=False,plot_color=False,cntr_color=False,cntr_lw=False,ccolor=False,out=True,plot_id=True,plot_below=False,shifted=False,plot_cmp_arrow=False,plot_label=True,plot_ontop=False):
        sys.stdout.write('Plot modelcomponents over clean maps\n')
        nn = len(self.model)
        nx = 1
        ny = int(np.ceil(nn/nx))

        if not cntr_color:
                cntr_color = 'grey'
        if not cntr_lw:
            cntr_lw=1

        if ccolor:
            sys.stdout.write('Set user component color\n')
            ccolors= len(self.ids)*[ccolor]
        else:
            ccolors = self.id_colors
        if dec:
            if type(dec)==list:
                if len(dec)>2:
                    if len(dec) > len(self.keys):
                        sys.stdout.write('please give only 1 parameter for each fits-image ra and dec\n')
                        return
            else:
                ra  = np.mean([self.cclean[key]['fov']/3. for key in self.keys()])
                dec = [self.cclean[key]['fov']/5. for key in self.keys()]

        Dra = [-ra,ra]
        shifty = [0]
        nk = len(self.keys)
        for kk,keys in enumerate(self.keys):
            if kk>0:
                shifty.append(shifty[kk-1]-2*dec[kk-1])
        if plot_ontop:
            Ddec = [-dec[0],dec[0]]
        else:
            Ddec= [shifty[-1]-dec[0]/2.,dec[0]]

        xsize=6
        if plot_ontop:
            ysize=dec[0]*xsize/ra
        else:
            ysize = np.sum(dec)*xsize/ra
        fig,ax = plt.subplots(1,1,figsize=(xsize,ysize), constrained_layout=True)
        ax.set_ylabel('Freq [GHz]',fontsize=16)
        ax.set_xlabel('RA [mas]',fontsize=16)

        lines   = []
        labels  = []
        ####################
        # setting all parameters for plotting a clean image
        #####################
        for kk,key in enumerate(self.keys):
            modelh= self.model[key]
            clean = self.cchead[key]
            if shifted:
                model = self.model[key]['data_shifted'].copy()
            else:
                model = self.model[key]['data'].copy()
            if shifted:
                ccmap = self.ccmap_shifted[key].copy()
            else:
                ccmap =  self.ccmap[key].copy()
            #shift the maps and models
            inc=clean['px_inc']*3.6e6
            if not plot_ontop:
                if key!=self.keys[0]:
                    model['DELTAY'] += shifty[kk]
                    model['DIST'] += shifty[kk]

            if plot_ontop:
                sdec =dec[kk]
            else:
                sdec = shifty[kk]
            # create loglevels
            scale   = -clean['px_inc']*3.6e6
            sigma   = sigma
            level0  = clean['noise']*sigma
            lev=[]
            for i in range(0,10):
                lev.append(level0*2**i)
            xx  = np.linspace(-clean['naxis']*0.5*scale,(clean['naxis']*0.5-1)*scale,clean['naxis'])
            yy  = np.linspace(clean['naxis']*0.5*scale,-(clean['naxis']*0.5-1)*scale,clean['naxis'])
            vmax=0.5*ma.amax(ccmap)
            norm = mpl.colors.SymLogNorm(linthresh=level0,linscale=0.5,vmin=level0,vmax=vmax,base=np.e)

            ##################
            # Plotting
            ###################
            ax.axis('scaled')
            ax.set_xlim(Dra)
            ax.set_ylim(Ddec)
            ax.minorticks_on()
            ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

            #ax.set_aspect(1)
            #ax.set_adjustable("datalim")
            ax.invert_xaxis()
            plotBeam(clean['beam'][0],clean['beam'][1],clean['beam'][2],ra,sdec-dec[kk],ax)

            if plot_ontop:
                colorsc =['red','blue','green','black']
                cntr_color = colorsc[kk]
                extent = np.max(xx),np.min(xx),np.min(yy),np.max(yy)
            else:
                extent = np.max(xx),np.min(xx),np.min(yy)+shifty[kk],np.max(yy)+shifty[kk]
            cntr = ax.contour(ccmap,linewidths=cntr_lw,levels=lev,colors=cntr_color,alpha=1,extent=extent)
            if plot_color:
                im = ax.imshow(ccmap,cmap=self.colormap,extent=extent,origin='lower', interpolation='gaussian')
                im.set_norm(norm)

            Mx  = model['DELTAX']
            My  = model['DELTAY']
            Mposa   = model['POSANGLE']
            Mmaj    = model['MAJOR AX']
            Mmin    = model['MINOR AX']
            for j,xx in enumerate(Mx):
                epid=model['id'][j]
                ccolor = ccolors[self.ids==epid][0]
                e_comp = Ellipse([Mx[j],My[j]],Mmaj[j],Mmin[j],-Mposa[j]+90, color=ccolor, zorder=2, fill=False,lw=0.5)
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
                    ax.plot(maj1_x,maj1_y, color = ccolor,marker='+', markersize=5,lw=1,label=model['id'][j])
                else:
                    ax.plot([maj1_x,maj2_x],[maj1_y,maj2_y], color = ccolor, lw = 1,label=model['id'][j])
                    ax.plot([min1_x,min2_x],[min1_y,min2_y], color = ccolor, lw = 1)
                if plot_cmp_arrow:
                    ax.annotate('{}'.format(model['id'][j]),xy=(Mx[j],My[j]), xycoords='data', xytext=(-5, (-1)**j*40), textcoords='offset points',arrowprops=dict(arrowstyle="->", connectionstyle="arc3"), color=ccolor,fontsize=10)

                if key==self.keys[-1]:
                    ax.axvline(x = maj1_x, color = ccolor)

            # set axis, labels, etc.
#            scalebar = AnchoredSizeBar(ax.transData, 1, '1 mas', 'lower right',frameon=False,size_vertical=0.04)
#            ax.add_artist(scalebar)
#            ax.get_yaxis().set_visible(False)
#            if key != self.keys[-1]:
#                ax.spines['bottom'].set_visible(False)
#            if key != self.keys[0]:
#                ax.spines['top'].set_visible(False)

           # ax.annotate('{:.1f} GHz'.format(modelh['freq']),xy=(-0.5*ra,sdec-0.1*sdec),xycoords='data',size=14)
           # ax.annotate('{}'.format(modelh['date_obs']),xy=(-0.5*ra,sdec-(0.1*sdec)-0.3),xycoords='data',size=14)

# The following lines sort the labels and handles alphabetically 
            if plot_label:
                handles, Labels = ax.get_legend_handles_labels()
                labels.extend(Labels)
                lines.extend(handles)

        if plot_label:
            by_label = dict(zip(labels, lines))
            keys_list = sorted(by_label, key=keyfunc)
            sorted_by_label = {key: by_label[key] for key in keys_list}
            pos = ax.get_position()
            ax.set_position([pos.x0, pos.y0, pos.width, pos.height * 0.85])
            #ax.set_position([pos.x0, pos.y0, pos.width * 0.9, pos.height])
            fig.legend(sorted_by_label.values(), sorted_by_label.keys(),loc='upper left', bbox_to_anchor=(0.12,0.98),ncol=3)
# sorting done, legend plottet
        tlabel = ['{}'.format(self.model[key]['freq']) for key in self.keys]
        plt.yticks(shifty,tlabel,fontsize=14)
        plt.tight_layout()#pad=0.2,w_pad=0.2,h_pad=0.2)
        fig.subplots_adjust(right=0.95, top=0.98)

        if out:
            if type(out)==bool:
                outf = modelh['source']+'_map+model_evolution.pdf'
            if shifted:
                outf = modelh['source']+'_map+model_evolution_shifted.pdf'
            fig.savefig(outf,bbox_inches='tight')
            sys.stdout.write('Plot has been written to {}\n'.format(outf))
        else:
            plt.show()
        plt.clf()
        return fig,ax



    def overplot_model(self,sigma=3,fig_size='screen',ra=False,dec=False,saveFile=False,plot_color=False,cntr_color=False,cntr_lw=False,ccolor=False,out=True,plot_id=True,shifted=False,plot_cmp_arrow=False,plot_label=True):
        sys.stdout.write('Plot modelcomponents over clean maps\n')
#        if plot_all:
        nn = len(self.model)
        if nn>6:
            nx = 4
        else:
            nx = 2
        ny = int(np.ceil(nn/nx))

        if not cntr_color:
            cntr_color = 'grey'
        if not cntr_lw:
            cntr_lw=1

        if ccolor:
            sys.stdout.write('Set user component color\n')
            ccolors= len(self.ids)*[ccolor]
        else:
            ccolors = self.id_colors
        load_gs = False
        if not dec:
            ra  = np.mean([self.cclean[key]['fov']/3. for key in self.keys()])
            dec = np.mean([self.cclean[key]['fov']/5. for key in self.keys()])
        RA = ra
        DEC = dec
        if type(dec)==list:
            if len(dec)>2:
                if len(dec) > len(self.keys):
                    sys.stdout.write('please give only 1 parameter for each fits-image ra and dec\n')
                    return
                load_gs = True
                gs_heights= dec
            elif len(dec)==2:
                Dra  = ra
                Ddec = dec
        else:
            Dra = [-ra,ra]
            Ddec= [-dec,dec]

        xsize =6
        xxs = xsize*nx
        yys=np.mean(dec)*xxs/np.mean(ra)
        yys*=ny

        fig,axs = plt.subplots(ny,nx, figsize=(xxs,yys))
        axs   = trim_axs(axs,len(self.model))


        lines   = []
        labels  = []
        kk=0
        ####################
        # setting all parameters for plotting a clean image
        #####################
        for ax,key in zip(axs,self.keys):
            modelh= self.model[key]
            clean = self.cchead[key]
            if shifted:
                model = self.model[key]['data_shifted']
            else:
                model = self.model[key]['data']
            if shifted:
                ccmap = self.ccmap_shifted[key]
            else:
                ccmap =  self.ccmap[key]
            if load_gs:
                Dra = [-ra[kk],ra[kk]]
                Ddec = [-dec[kk],dec[kk]]
                RA = ra[kk]
                DEC = dec[kk]

            scale   = -clean['px_inc']*3.6e6
            sigma   = sigma
            level0  = clean['noise']*sigma
            lev=[]
            for i in range(0,10):
                lev.append(level0*2**i)
            xx  = np.linspace(-clean['naxis']*0.5*scale,(clean['naxis']*0.5-1)*scale,clean['naxis'])
            yy  = np.linspace(clean['naxis']*0.5*scale,-(clean['naxis']*0.5-1)*scale,clean['naxis'])
            extent = np.max(xx),np.min(xx),np.min(yy),np.max(yy)
            vmax=0.5*ma.amax(ccmap)
            norm = mpl.colors.SymLogNorm(linthresh=level0,linscale=0.5,vmin=level0,vmax=vmax,base=np.e)

            ##################
            # Plotting
            ###################
            #f,ax = plt.subplots()
            ax.axis('scaled')
            ax.set_xlim(Dra)
            ax.set_ylim(Ddec)
            ax.set_aspect(1)
            ax.set_adjustable("datalim")
            ax.invert_xaxis()
            plotBeam(clean['beam'][0],clean['beam'][1],clean['beam'][2],RA,-DEC,ax)

           # cntr=ax.contour(xx,yy,ccmap,linewidths=cntr_lw,levels=lev,colors=cntr_color,alpha=1)
            cntr = ax.contour(ccmap,linewidths=cntr_lw,levels=lev,colors=cntr_color,alpha=1,extent=extent)
            if plot_color:
                im = ax.imshow(ccmap,cmap=self.colormap,extent=extent,origin='lower', interpolation='gaussian')
                im.set_norm(norm)

            Mx  = model['DELTAX']
            My  = model['DELTAY']
            Mposa   = model['POSANGLE']
            Mmaj    = model['MAJOR AX']
            Mmin    = model['MINOR AX']
            for j,xx in enumerate(Mx):
                epid=model['id'][j]
                ccolor = ccolors[self.ids==epid][0]
                e_comp = Ellipse([Mx[j],My[j]],Mmaj[j],Mmin[j],-Mposa[j]+90, color=ccolor, zorder=2, fill=False,lw=0.5)
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
                    ax.plot(maj1_x,maj1_y, color = ccolor,marker='+', markersize=5,lw=1,label=model['id'][j])
                else:
                    ax.plot([maj1_x,maj2_x],[maj1_y,maj2_y], color = ccolor, lw = 1,label=model['id'][j])
                    ax.plot([min1_x,min2_x],[min1_y,min2_y], color = ccolor, lw = 1)
                if plot_cmp_arrow:
                    ax.annotate('{}'.format(model['id'][j]),xy=(Mx[j],My[j]), xycoords='data', xytext=(-5, (-1)**j*40), textcoords='offset points',arrowprops=dict(arrowstyle="->", connectionstyle="arc3"), color=ccolor,fontsize=8)


            # set axis, labels, etc.
            ax.set(xlabel='RA [mas]', ylabel='DEC [mas]')
            ax.minorticks_on()
            ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            ax.set_title('{} - {:.1f} GHz - {}'.format(modelh['source'],modelh['freq'],modelh['date_obs']),size=12)
# The following lines sort the labels and handles alphabetically 
            if plot_label:
                handles, Labels = ax.get_legend_handles_labels()
                labels.extend(Labels)
                lines.extend(handles)
            kk+=1

        if plot_label:
            by_label = dict(zip(labels, lines))
            keys_list = sorted(by_label, key=keyfunc)
            sorted_by_label = {key: by_label[key] for key in keys_list}
            pos = ax.get_position()
            ax.set_position([pos.x0, pos.y0, pos.width, pos.height * 0.85])
            #ax.set_position([pos.x0, pos.y0, pos.width * 0.9, pos.height])
            fig.legend(sorted_by_label.values(), sorted_by_label.keys(),loc='upper left', bbox_to_anchor=(1,0.98),ncol=2)
            #secx = ax.secondary_xaxis('left',)
            #secx.set_xlabel(args['Freq'])
# sorting done, legend plottet

        plt.tight_layout()#pad=0.2,w_pad=0.2,h_pad=0.2)
        fig.subplots_adjust(right=0.95, top=0.98)

        if out:
            if type(out)==bool:
                outf = modelh['source']+'_Model_overplot.pdf'
            if shifted:
                outf = modelh['source']+'_Model_overplot_shifted.pdf'
            fig.savefig(outf,bbox_inches='tight')
            sys.stdout.write('Plot has been written to {}\n'.format(outf))
        else:
            plt.show()
        plt.clf()
        return fig,ax

    def plot_all(self,xax='date',yax='DIST',out=False):
        fig,ax = plt.subplots(figsize=(12,8))
        dates = []
        ymin = 0
        ymax = 0
        for modF in self.model:
            dates.append(self.model[modF]['date'])
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
           sys.stdout.write('Plot has been written to {}\n'.format(outf))
        else:
            plt.show()
        plt.cla()

    def fit_comp_spectrum(self,additional_map=False,plot_areas=False,plot_all_components=False,comps=False,ccolor=False,out=True):
        sys.stdout.write('Fit component spectrum\n')
        if not self.model_sorted:
            self.sort_by_id()
        if comps:
            model_sorted = {key: self.model_sorted[key] for key in comps}
        else:
            model_sorted = self.model_sorted.copy()

        if ccolor:
            ccolors= np.array(len(self.ids)*[ccolor])
        else:
            ccolors = self.id_colors

        NCOMP = len(model_sorted)
        if NCOMP<5:
            nx=2
        if NCOMP<13:
            nx = 3
        else:
            nx = 4
        ny= int(np.ceil(NCOMP/nx))
        print(ny)
        xsize =15
        xxs = xsize
        yys= xsize*ny/nx
        #yys =12
        yys*=ny

      #  figsize = set_size('aanda*',subplots=(ny,nx))
      #  fig,axs = plt.subplots(ny,nx)
        axe_ratio = 'scaled'
        #fig,axs=plt.subplots(nrows=ny,ncols=nx,figsize=set_size('aanda*',subplots=(ny,nx)),sharex=True,sharey=True)
   #     fig,axs=plt.subplots(ny,nx,figsize=set_size('aanda*',subplots=(ny,nx)))#,sharex=True,sharey=True)
        fig,axs = plt.subplots(ny,nx, figsize=(xxs,xxs+2))

        axs = trim_axs(axs,NCOMP)
        #axs.set(xlabel=xlabel, ylabel=ylabel)
        fig.supxlabel('Frequency [GHz]')
        fig.supylabel('Flux density [Jy]')

# Doing the fit
        sn_p,sn_sd,sn_ch2,an_out,pl_p,pl_sd,pl_ch2,pl_out=[],[],[],[],[],[],[],[]
        fit=[]
        CompSN,CompPL,num,Sm,athin,athick,alpha,chi2SN,chi2PL,numE,SmE,athinE,athickE,alphaE=[],[],[],[],[],[],[],[],[],[],[],[],[],[]
        for i,comp in enumerate(model_sorted):
            print('Fit power law to Comp '+comp)
            print(comp)
            cflux = self.model_sorted[comp]['FLUX']
            cfluxerr = 0.15*cflux.copy()
            cfreq = self.model_sorted[comp]['freq']
            cid = self.model_sorted[comp]['id']

            pl_x0 = np.array([np.min(cflux),-1])

            beta,sd_beta,chi2,pl_out = ff.odr_fit(ff.powerlaw,[cfreq,cflux,cfluxerr],pl_x0,verbose=1)
            #ff.odr_fit(ff.powerlaw,[cc['Freq'],cc['Flux'],cc['Fluxerr']],pl_x0,verbose=1)
            pl_p.append(beta)
            pl_sd.append(sd_beta)
            pl_ch2.append(chi2)
            #fit Snu
            print('Fit SSA to Comp '+cid[0])
            sn_x0 = np.array([100,np.max(cflux),3,-1])
            if cid[0]=='B3':
                sn_x0[0]=20
            beta,sd_beta,chi2,sn_out = ff.odr_fit(ff.Snu,[cfreq,cflux,cfluxerr],sn_x0,verbose=1)
            sn_p.append(beta)
            sn_sd.append(sd_beta)
            sn_ch2.append(chi2)

            if np.logical_and(sn_ch2[i]>pl_ch2[i],pl_out.info<5):
                sys.stdout.write('Power law fits better\n')
                CompPL.append(cid[0])
                alpha.append(pl_p[i][1])
                alphaE.append(pl_sd[i][1])
                chi2PL.append(pl_ch2[i])
                fit.append('PL')
            elif np.logical_and(pl_ch2[i]>sn_ch2[i],sn_out.info<5):
                sys.stdout.write('ssa spectrum fits better\n')
                CompSN.append(cid[0])
                num.append(sn_p[i][0])
                Sm.append(sn_p[i][1])
                athin.append(sn_p[i][3])
                athick.append(sn_p[i][2])
                chi2SN.append(sn_ch2[i])
                numE.append(sn_sd[i][0])
                SmE.append(sn_sd[i][1])
                athinE.append(sn_sd[i][3])
                athickE.append(sn_sd[i][2])
                fit.append('SN')
            else:
                sys.stdout.write('NO FIT WORKED, use power law\n')
                CompPL.append(cid[0])
                alpha.append(pl_p[i][1])
                alphaE.append(pl_sd[i][1])
                chi2PL.append(pl_ch2[i])
                fit.append('PL')
        namesSN=['Comp','num','numE','Sm','SmE','athin','athinE','athick','athickE','Chi2']
        SN = Table([CompSN,np.round(num,2),numE,np.round(Sm,2),SmE,np.round(athin,2),athinE,np.round(athick,2),athickE,chi2SN],names=namesSN)
        print(SN)
        namesPL=['Comp','alpha','alphaE','Chi2']
        PL = Table([CompPL,np.round(alpha,2),alphaE,chi2PL],names=namesPL)
        print(PL)
        savefile ='component_spectrum_fit'
        ascii.write(SN,savefile+'_SN.tex',formats=dict.fromkeys(namesSN[1:], '%.2f'),overwrite=True,format='latex')
        ascii.write(PL,savefile+'_PL.tex',formats=dict.fromkeys(namesPL[1:], '%.2f'),overwrite=True,format='latex')
        #exponent = math.floor(np.log10(np.min(vstack([d for d in self.model_sorted[model_sorted]])['Flux'])))
        exponent = -2
        ymin=float('1e{}'.format(exponent))
        props = dict(boxstyle='round',fc='w',alpha=0.5)

        #doing the plotting
        i=0
        for ax,comp in zip(axs,model_sorted):
            cflux = self.model_sorted[comp]['FLUX']
            cfluxerr = 0.15*cflux.copy()
            cfreq = self.model_sorted[comp]['freq']
            cid = self.model_sorted[comp]['id']
            color = ccolors[self.ids==comp][0]
            symbol = self.symbols[i]

          #  ax.axis('scaled')
         #   ax.set_aspect(1)
         #   ax.set_adjustable("datalim")

            ax.set_xscale('log')
            ax.set_yscale('log')
          #  ax.axis(axe_ratio)
            ax.minorticks_on()
            ax.errorbar(cfreq,cflux,cfluxerr,marker='.',ms=3,ls='none',color='red',elinewidth=0.8,label='Comp {}'.format(cid[0]))
            xr=np.arange(1,300,0.01)
            if fit[i]=='PL':
                textstr = '\n'.join((
                    r'$\alpha={:.2f}\pm{:.2f}$'.format(pl_p[i][1],pl_sd[i][1]),
                ))
                ax.annotate(textstr, xy=(0.05,0.1),xycoords='axes fraction',fontsize=6,bbox=props)
                ax.plot(xr,ff.powerlaw(pl_p[i],xr),'k',lw=0.5)
            elif fit[i]=='SN':
                textstr = '\n'.join((
                    r'$\nu_m={:.2f}$'.format(sn_p[i][0]),
                    r'$S_m={:.2f}$'.format(sn_p[i][1]),
                    r'$\alpha_{{thin}}={:.2f}$'.format(sn_p[i][3]),
                    r'$\alpha_{{thick}}={:.2f}$'.format(sn_p[i][2]),
                    r'$\chi_\mathrm{{red}}^2={:.2f}$'.format(sn_ch2[i])
                ))
        
                ax.annotate(textstr, xy=(0.05,0.1),xycoords='axes fraction',fontsize=5,bbox=props)
                ax.plot(xr,ff.Snu(sn_p[i],xr),'k',lw=0.5)
        
            ax.set_xlim(1,300)
            ax.set_ylim(ymin,5)
        
            ax.set_yticks([0.001,0.01,0.1,1])
            ax.set_xticks([1,10,100])
            ax.tick_params(axis="x",direction="in")
            ax.tick_params(axis="y",direction="in")
            ax.legend(loc=1)
            i+=1
#        handles, labels = ax.get_legend_handles_labels()
#        by_label = dict(zip(labels, handles))
#        keys_list = sorted(by_label, key=keyfunc)
#        sorted_by_label = {key: by_label[key] for key in keys_list}
#        fig.legend(sorted_by_label.values(), sorted_by_label.keys(),loc='upper left', bbox_to_anchor=(0.9,0.98),ncol=1)
        fig.subplots_adjust(right=0.9, top=0.98, left=0.08,bottom=0.06)
        #set_corrected_size(fig,figsize)
        if out:
            if type(out)==bool:
                outf = self.model[self.keys[0]]['source']+'Component_spectrum_fit.pdf'
            fig.savefig(outf,bbox_inches='tight')
            sys.stdout.write('Plot has been written to {}\n'.format(outf))
        else:
            plt.show()


    def write_tex_table(self):
        '''Write out a tex table with all modelfit parameters.
        '''
        sys.stdout.write('Writing model component parameters to tex-table.\n')
        i=0
        if os.path.exists('Allmodel.tex'):
            os.remove('Allmodel.tex')
            sys.stdout.write('\nA file Allmodel.tex already exists, will be deleted and created.\n')
        else:
            sys.stdout.write('\nA file Allmodel.tex will be created.\n')

        for i,key in enumerate(self.model):

            ff = ['%s','%1.3f','%2.2f','%2.2f','%2.2f','%2.1f','%2.2f','%2.2f']
            kk = ['id','FLUX','DELTAX','DELTAY','MAJOR AX','ratio','DIST','logtb']
            formats=dict(zip(kk,ff))
            with open('Allmodel.tex',mode='a') as f:
                f.write('\\midrule\\multicolumn{8}{c}{%2.1f\\GHz - Band}\\\\\\midrule\n'%self.model[key]['freq'])
#                if i==0:
                self.model[key]['data'].write(f,format='ascii.no_header', include_names=kk,delimiter='&', formats=formats,fill_values=[(ascii.masked, '  --  ')],overwrite=False)

 #               else:
  #                  dd.write(f,format='ascii.no_header', include_names=kk,delimiter='&', formats=formats,fill_values=[(ascii.masked, '  --  ')],overwrite=False)

            #i+=1

        with open('Allmodel.tex',mode='r') as f:
            file_lines=[''.join([x.strip(),'\\\\', '\n']) for x in f.readlines()]

        freq_list = [str(self.model[key]['freq']) for key in self.keys]
        with open('Allmodel.tex',mode='w') as f:
            f.write('\\begin{table}[!hbt]\n')
            f.write('\\centering\n')
            f.write('\\small\n')
            f.write('\\setlength{\extrarowheight}{1pt}\n')
            f.write('\\setlength{\\tabcolsep}{3pt}\n')
            f.write('\\caption{Parameters of Gaussian model fitting components for '+','.join(freq_list)+'\,GHz observation for '+self.model[self.keys[0]]['source']+'.\label{tab:vlbamodels}\n')
            f.write('\\begin{adjustbox}{width=\linewidth}\n')
            f.write('\\begin{threeparttable}\n')
            f.write('\\begin{tabular}{lSSSSSSc}\\toprule\n')
            f.write('{ID} & {Flux density} & {RA} & {DEC} & {Major \\tnote{1}} & {Ratio\\tnote{2}} & {Distance} & {$\log\,T_\mathrm{b}$} \\\\\n')
            f.write(' & {$[$Jy$]$} & {$[$mas$]$} & {$[$mas$]$} & {$[$mas$]$} & & {$[$mas$]$} & {$[$K$]$} \\\\\n')
#            f.write('\\midrule\n')
            f.writelines(file_lines)
            f.write('\\bottomrule\n')
            f.write('\\end{tabular}\n')
            f.write('\\begin{TableNotes}\n')
            f.write('\\footnotesize\n')
            f.write('\\item[1] FWHM major axis of restoring beam\n')
            f.write('\\item[2] Ratio of minor to major axis\n')
            f.write('\\end{TableNotes}\n')
            f.write('\\end{threeparttable}\n')
            f.write('\\end{adjustbox}\n')
            f.write('\\end{table}\n')
        sys.stdout.write('\nAllmodel.tex has been written.\n')
