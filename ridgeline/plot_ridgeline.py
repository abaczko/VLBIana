#!/usr/bin/env python
# -*- coding: utf-8 -*-

from astropy.table import Table
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition, mark_inset)
import os,sys
from glob import glob
from VLBIana.modules.plot_functions import *
from VLBIana.modules.cleanMap import *
from VLBIana.modules.jet_calculus import *
from VLBIana.ridgeline.ridgeLine import *

def ridgeline_plotter(mapF,ridgeL,saveFile,lab,theta,logF,modFile=False,**kwargs):

	args = {'xr'				: False,
					'yr'				: False,
					'errorbars'	: False,
					'plot_rl'		: False,
					'plot_width': False,
					'plot_deconvolved'	: False,
					'plot_beam'	: False,
					'plot_map'	: False,
					'plot_reslimit' : False,
					'ra_lim'	: False,
					'dec_lim'	: False,
					'asize'		: 8,
					'add_rl'	: False,
					'fig_size':'aanda*',
					'fig_extension': 'pdf',
					'incl' : 90,
					'binRidgeLine' : False,
					'errorFile':False,
					'multiL' : False,
					'plot_stacked' : False,
					'shiftFile': False,
					'overplot_rl': True,
					'sigma' : [4]}
	args.update(kwargs)

	if type(mapF)==str:
		mapF		= [mapF]
		ridgeL	= [ridgeL]
		logF		= [logF]
		modFile	= [modFile]
		lab			= [lab]
	else:
		pass
	mapFile		= mapF.copy()
	ridgeLine	= ridgeL.copy()
	logFile		= logF.copy()
	modFile		= modFile
	label			= lab.copy()

	print('files \n{}\n will be plotted'.format(mapFile))


	if args['add_rl']:
		mapFile.append(args['add_rl'][0])
		ridgeLine.append(args['add_rl'][1])
		logFile.append(args['add_rl'][2])
		label.append(args['add_rl'][3])
		if args['errorFile']:
			args['errorFile'].append(args['add_rl'][4])
		n = len(ridgeLine)
		colors = cmap(np.linspace(0,0.95,n))
		colors = colors[len(colors)-n:]
		colors[n-1]=[0.8,0.5,0.5,1]
		mpl.rcParams['axes.prop_cycle'] = cycler(color=colors)

	header				= [read_header(m) for m in mapFile]
	px_inc				= [h['CDELT2'] for h in header]
	naxis					= [h['NAXIS1'] for h in header]
	maps_beam			= [[h['BMAJ']*3.6e6,h['BMIN']*3.6e6,h['BPA']] for h in header]
	ppb     = [PXPERBEAM(bp[0]*np.pi/180,bp[1]*np.pi/180,pxi*3.6e6*np.pi/180) for bp,pxi in zip(maps_beam,px_inc)]
	noise   = [h['NOISE'] for h in header]
	freqs   = [np.around(h['crval3']*1e-9,1) for h in header]
	fov			= [pxi*na*3.6e6 for pxi,na in zip(px_inc,naxis)]

	flux_uncertainty = 0.15
	theta = theta*np.pi/180
	if args['shiftFile']:
		relative_shifts = Table.read(args['shiftFile'], format='ascii')
		sys.stdout.write('The following shifts have been read and will be used for aligning the ridgelines.\n')
		sys.stdout.write('{}'.format(relative_shifts))
		shiftRA,shiftDEC = readShift(args['shiftFile'])
		DistStr = 'Distshift'
		RAStr = 'RAshift'
		DecStr = 'Decshift'
		Ridge = [RidgeLine(RL,log,cmap,shift=(shiftra,shiftdec),incl=args['incl']) for RL,log,cmap,shiftra,shiftdec in zip(ridgeLine,logFile,mapFile,shiftRA,shiftDEC)]
	else:
		sys.stdout.write('No shifts will be applied to ridge-lines.\n')
		Ridge = [RidgeLine(RL,log,cmap,incl=args['incl']) for RL,log,cmap in zip(ridgeLine,logFile,mapFile)]
		DistStr = 'Dist'
		RAStr = 'RA'
		DecStr = 'Dec'

	if args['errorFile']:
		sys.stdout.write('Use errorfile')
		ridgelines = [RL.readRidgeline(theta,widthErr=eF) for RL,eF in zip(Ridge,args['errorFile'])]
	else:
		ridgelines = [RL.readRidgeline(theta) for RL in Ridge]

	if args['binRidgeLine']:
		ridgelinesbinned = [RL.binRidgeLine() for RL in Ridge]
		ridgelines = ridgelinesbinned
	i=0
	Jet,CJet = [],[]
	for rl in ridgelines:
	#	rl['Freq']			= freqs[i]*len(rl['Xpos'])
		Jet.append(rl[np.sign(rl[DistStr])==1])
		CJet.append(rl[np.sign(rl[DistStr])==-1])
		i+=1
	
	if args['plot_reslimit']:
		mmod,mheader,mimg=[],[],[]
		for ff in modFile:
			cc,hh,ii = read_fits(ff)
			mmod.append(cc)
			mheader.append(hh)
			mimg.append(ii)
		noise_mmap  = [h['NOISE'] for h in mheader]

		mod_maj   = [m['MAJOR AX']*3.6e6 for m in mmod]
		mod_flux  = [m['FLUX'][np.where(mmaj!=0)] for m,mmaj in zip(mmod,mod_maj)]
		mod_min   = [m['MINOR AX'][np.where(mmaj!=0)]*3.6e6 for m,mmaj in zip(mmod,mod_maj)]
		mod_beam  = [[h['BMAJ']*3.6e6,h['BMIN']*3.6e6,h['BPA']] for h in mheader]
		if args['shiftFile']:
			mod_x			= [m['DELTAX'][np.where(mmaj!=0)]*3.6e6+rs for m,rs,mmaj in zip(mmod,relative_shifts['RA'],mod_maj)]
			mod_y			= [m['DELTAY'][np.where(mmaj!=0)]*3.6e6+rs for m,rs,mmaj in zip(mmod,relative_shifts['DEC'],mod_maj)]
		else:
			mod_x			= [m['DELTAX'][np.where(mmaj!=0)]*3.6e6 for m,mmaj in zip(mmod,mod_maj)]
			mod_y			= [m['DELTAY'][np.where(mmaj!=0)]*3.6e6 for m,mmaj in zip(mmod,mod_maj)]

		mod_r			= [np.sign(mx)*np.sqrt(mx**2+my**2) for mx,my in zip(mod_x,mod_y)]
		ind				= [m.argsort() for m in mod_r]
		mod_maj		= [m[np.where(m!=0)] for m in mod_maj]
		mod_maj		= [m[i] for m,i in zip(mod_maj,ind)]
		mod_flux  =	[m[i] for m,i in zip(mod_flux,ind)]
		mod_min   = [m[i] for m,i in zip(mod_min,ind)] 
		mod_x			= [m[i] for m,i in zip(mod_x,ind)]		
		mod_y			= [m[i] for m,i in zip(mod_y,ind)]		
		mod_r			= [m[i] for m,i in zip(mod_r,ind)]		


		mod_peakflux = [mod_peak_flux(S,amaj,amin,bb[0],bb[1]) for S,amaj,amin,bb in zip(mod_flux,mod_maj,mod_min,mod_beam)]
		snr_mmap = [pp/(pb*nn) for pp,nn,pb in zip(mod_flux,noise_mmap,ppb)]
		reslim = [resLim(amaj,amin,snr,0) for amaj,amin,snr in zip(mod_maj,mod_min,snr_mmap)] 

	xData	=	[rl['Xpos'].copy() for rl in ridgelines]
	yData	=	[rl['Ypos'].copy() for rl in ridgelines]
	Dist	= [rl[DistStr].copy() for rl in ridgelines]
	FWHM	= [rl['FWHM'].copy() for rl in ridgelines]
	flux	= [rl['Peak'].copy() for rl in ridgelines]
	FWHM_deconvolved = [rl['FWHMDeconvolved'].copy() for rl in ridgelines]
	fwhm_error_deconvolved = [rl['FWHMDeconvolvedErr'].copy() for rl in ridgelines]
	RA		= [rl[RAStr].copy() for rl in ridgelines]
	Dec		= [rl[DecStr].copy() for rl in ridgelines]
	RA_error	= [rl['RAErr'].copy() for rl in ridgelines]
	Dec_error	= [rl['DecErr'].copy() for rl in ridgelines]

	saveFile+='_shift'

#
	fwhm_error	=  [rl['fomalont'].copy() for rl in ridgelines]
#
	flux_error	= [f*0.15 for f in flux]
#	

	n=len(ridgelines)
	if args['xr']:
		xmin	= args['xr'][0]
		xmax	= args['xr'][1]
	else:
		xmin = np.min([np.min(x) for x in RA])
		xmax = np.max([np.max(x) for x in RA])
	if args['yr']:
		ymin = args['yr'][0]
		ymax = args['yr'][1]
	else:
		ymin = np.min([np.min(y) for y in Dec])
		ymax = np.max([np.max(y) for y in Dec])
	
	##########################################
	###### Produce the plots #################
	##########################################
	if args['plot_map']:
		RAm		= [rl[RAStr].copy() for rl in ridgelines]
		Decm		= [rl[DecStr].copy() for rl in ridgelines]

		nn = len(mapFile)
		if nn>1:
			if args['multiL']:
				sys.stdout.write('Will use image dimension specific for multi-lamda data set.\n')
				ra=[75,35,35,9,9,9]
				dec=[rra/2. for rra in ra]
			else:
				if not args['ra_lim']:
					sys.stdout.write('Will calculate image dimensions. Give parameters ra_lim and dec_lim if other dimensions are required.\n')
					ra  = [FOV/3 for FOV in fov]
					dec = [FOV/5 for FOV in fov]
				else:
					ra	= args['ra_lim']
					dec	= args['dec_lim']
	
		else:
			if not args['ra_lim']:
				ra  = [FOV/3 for FOV in fov]
				dec = [FOV/5 for FOV in fov]
			else:
				ra	= [args['ra_lim']]
				dec	= [args['dec_lim']]
		
		ra_min=[-rr for rr in ra]
		ra_max=ra
		dec_min=[-dd for dd in dec]
		dec_max=dec

		scale=[-pxi*3.6e6 for pxi in px_inc]
		if len(args['sigma'])==1:
			sigma=len(px_inc)*args['sigma']
		else:
			sigma=args['sigma']
		level0=[n*s for n,s in zip(noise,sigma)]
		lev=[]
		for l in level0:
			ll=[]
			for i in range(0,10):
				ll.append(l*2**i)
			lev.append(ll)
		xx=[np.linspace(-nn*0.5*ss,(nn*0.5-1)*ss,nn) for nn,ss in zip(naxis,scale)]
		yy=[np.linspace(nn*0.5*ss,-(nn*0.5-1)*ss,nn) for nn,ss in zip(naxis,scale)]
		if nn>1:
			if nn % 2:
				xs = 1
				ys = nn
				if args['fig_size']=='aanda*':
					args['fig_size']='aanda'

			else:
				xs = 2
				ys = int(nn/2)
		else:
			xs = 1
			ys = 1

		if np.logical_or(xs>1,ys>1):
			figsize=set_size(args['fig_size'],subplots=(ys,xs))
			#figsize = set_scaled_size(args['fig_size'],subplots=(ys,xs))
		else:
			figsize=set_size(args['fig_size'])

		f,ax = plt.subplots(ys,xs)
		axe_ratio='scaled'
		k=0
		l=0
		i=0
		for img in mapFile:
			Map = CleanMap(img)
			cleanmap = Map.cmap	
			vmin=level0[i]
			vmax=0.5*ma.amax(cleanmap)
			norm = mpl.colors.SymLogNorm(linthresh=vmin,linscale=0.5,vmin=vmin,vmax=vmax)
			if np.logical_and(xs>1,ys>1):	
				ax[k,l].axis(axe_ratio)
				ax[k,l].set_xlim(ra_min[i],ra_max[i])
				ax[k,l].set_ylim(dec_min[i],dec_max[i])
				ax[k,l].invert_xaxis()
				ax[k,l].annotate(label[i], xy=(0.65,0.9),xycoords='axes fraction',color='black',fontsize=args['asize'],bbox=bbox_props)
		
				plotBeam(maps_beam[i][1],maps_beam[i][0],maps_beam[i][2],ra_max[i],dec_min[i],ax[k,l])
				cntr=ax[k,l].contour(xx[i],yy[i],cleanmap,linewidths=0.5,levels=lev[i],colors='grey',alpha=1) #,extent=extent[i],origin='lower',alpha=1)
				cntr.set_norm(norm)
				if args['overplot_rl']:
					ax[k,l].scatter(RAm[i],Decm[i],c='red',s=1,marker='.')
					ax[k,l].errorbar(RAm[i],Decm[i],yerr=Dec_error[i],xerr=RA_error[i],fmt='none',c='red',alpha=0.5,elinewidth=0.4)
			elif np.logical_or(xs>1,ys>1):
				ax[k].axis(axe_ratio)
				ax[k].set_xlim(ra_min[i],ra_max[i])
				ax[k].set_ylim(dec_min[i],dec_max[i])
				ax[k].invert_xaxis()
				plotBeam(maps_beam[i][1],maps_beam[i][0],maps_beam[i][2],ra_max[i],dec_min[i],ax[k])
				cntr=ax[k].contour(xx[i],yy[i],cleanmap,linewidths=0.5,levels=lev[i],colors='grey',alpha=1) #,extent=extent[i],origin='lower',alpha=1)
				cntr.set_norm(norm)
				if args['overplot_rl']:
					ax[k].scatter(RAm[i],Decm[i],c='red',s=1,marker='.')
					ax[k].errorbar(RAm[i],Decm[i],yerr=Dec_error[i],xerr=RA_error[i],fmt='none',c='red',alpha=0.5,elinewidth=0.4)

				if args['plot_stacked']:
					ax[k].annotate('VLBA {0:.1f} GHz stacked'.format(freqs[i]), xy=(0.5,0.8),xycoords='axes fraction',color='black',fontsize=args['asize'],bbox=bbox_props)
				else:
					ax[k].annotate(label[i], xy=(0.5,0.8),xycoords='axes fraction',color='black',fontsize=args['asize'],bbox=bbox_props)
			else:	
				ax.axis(axe_ratio)
				ax.set_xlim(ra_min[i],ra_max[i])
				ax.set_ylim(dec_min[i],dec_max[i])
				ax.invert_xaxis()
				plotBeam(maps_beam[i][1],maps_beam[i][0],maps_beam[i][2],ra_max[i],dec_min[i],ax)
				cntr=ax.contour(xx[i],yy[i],cleanmap,linewidths=1,levels=lev[i],colors='grey',alpha=1) #,extent=extent[i],origin='lower',alpha=1)
				cntr.set_norm(norm)
				if args['overplot_rl']:
					ax.scatter(RAm[i],Decm[i],c='red',s=2,marker='.')
					ax.errorbar(RAm[i],Decm[i],yerr=Dec_error[i],xerr=RA_error[i],fmt='none',c='red',alpha=0.5,elinewidth=0.4)
				ax.annotate(label[i], xy=(0.7,0.9),xycoords='axes fraction',color='black',fontsize=args['asize']+3,bbox=bbox_props)

			i+=1
			if k<ys-1:
				k+=1
			elif l<xs-1:
				l+=1
				k=0
		if np.logical_or(xs>1,ys>1):	
			#ax[0,0].set(ylabel='Relative DEC [mas]')
			#ax[1,0].set(ylabel='Relative DEC [mas]')
			#ax[2,0].set(xlabel='RA [mas]', ylabel='Relative DEC [mas]')
			#ax[2,1].set(xlabel='RA [mas]')

			for aa in ax.flat:
				aa.set(xlabel='RA [mas]', ylabel='Relative DEC [mas]')
				aa.minorticks_on()
				aa.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
				aa.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
		else:
				ax.set(xlabel='RA [mas]', ylabel='Relative DEC [mas]')
				ax.minorticks_on()
				ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
				ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
		if args['overplot_rl']:
			saveFile+='_overplot'
		set_corrected_size(f,figsize)
	
###########################
	if args['plot_rl']:
		######
		if args['fig_size']=='aanda*':
			args['fig_size']='aanda'
		f,ax=plt.subplots()
		figsize=set_size(args['fig_size'])
		axe_ratio='scaled'
		ax.axis(axe_ratio)
		plt.xlabel('RA [mas]')
		plt.ylabel('Relative DEC [mas]')

		i=0
		for rl,m in zip(ridgelines,markers):
			ax.scatter(RA[i],Dec[i],s=1,marker='.',label=label[i],alpha=0.6)
			ax.errorbar(RA[i],Dec[i],yerr=Dec_error[i],xerr=RA_error[i],fmt='none',alpha=0.3,label='_nolegend_',elinewidth=0.4)

			i+=1
			
		plt.xlim(xmax*1.04,xmin*1.1)
		plt.ylim(ymin*1.04,ymax*1.5)
		handles, labels = ax.get_legend_handles_labels()
		ax.legend(handles,labels,loc=0,markerscale=2)
		ax.minorticks_on()
		set_corrected_size(f,figsize)

#############################################		
	elif args['plot_width']:
		ymin = 1e-1
		ymax = max([max(f) for f in FWHM])#20
		f,ax=plt.subplots(2,1,sharex=True,gridspec_kw={'hspace': 0})
		figsize=set_size(args['fig_size'],subplots=(2,1))
		axesWidthPlot(ax[1],ylabel='Peak Flux density [mJy/beam]',xlabel='Distance along jet [mas]',xscale='linear')

		i=0
		ax[0].set_xlim(xmax*1.1,xmin*1.1)
		if args['plot_deconvolved']:
			axesWidthPlot(ax[0],ylabel='De-convolved Jet width [mas]',xlabel='Distance along jet [mas]',secxax='Distance along jet [$R_\mathrm{S}$]',secyax='De-convolved Jet width [$R_\mathrm{S}$]',xscale='linear')
			for rl,m in zip(ridgeLine,markers):
				ax[0].scatter(Dist[i],FWHM_deconvolved[i],s=2,marker=m,label='{}'.format(label[i]))
				if args['errorbars']:
					ax[0].errorbar(Dist[i],FWHM_deconvolved[i],yerr=fwhm_error_deconvolved[i],fmt='none',elinewidth=0.2,errorevery=1,label='_nolegend_',alpha=0.4)
				i+=1
			ax[0].set_ylim(1e-3,ymax*1.1)
			saveFile+='_jet_width+flux_deconvolved'
		else:
			axesWidthPlot(ax[0],ylabel='Jet width [mas]',xlabel='Distance along jet [mas]',secxax='Distance along jet [$R_\mathrm{S}$]',secyax='Jet width [$R_\mathrm{S}$]',xscale='linear',yscale='linear')
			for rl,m in zip(ridgeLine,markers):
				ax[0].scatter(Dist[i],FWHM[i],s=2,marker=m,label='{}'.format(label[i]))
				if args['errorbars']:
					ax[0].errorbar(Dist[i],FWHM[i],yerr=fwhm_error[i],fmt='none',elinewidth=0.4,errorevery=1,label='_nolegend_',alpha=0.4)
				if args['plot_beam']:
					ax[0].plot(np.linspace(xmax-0.1*xmax,xmin),np.linspace(maps_beam[i][0],maps_beam[i][0]),linestyle=':',alpha=0.6,label='_nolegend_')
				else:
					print('no beam plotted')
				if args['plot_reslimit']:
					ax[0].plot(mod_r[i],reslim[i],linestyle=':',marker='o',markersize=2,lw=0.5,alpha=0.6,label='_nolegend_')
				else:
					print('no resolution limit plotted')
				ax[0].set_ylim(0,ymax)
				i+=1
			saveFile+='_jet_width+flux'
		
			####################################################################
		# Plot Amplitude 
		
		ymin = 1e-1
		ymax = 1e3
		i=0
	
		ax[1].xaxis.set_minor_locator(AutoMinorLocator())
		for rl,m in zip(ridgeLine,markers):
			ax[1].scatter(Dist[i],flux[i],s=2,marker=m,label='{}'.format(label[i]))
			if args['errorbars']:
				ax[1].errorbar(Dist[i],flux[i],yerr=flux_error[i],fmt='none',elinewidth=0.4,errorevery=1,label='_nolegend_',alpha=0.4)
			i+=1
		
		ax[1].set_xlim(xmax*1.1,xmin*1.1)
		ax[1].set_ylim(ymin*1.1,ymax*1.1)

		handles, labels = ax[0].get_legend_handles_labels()
		ax[0].legend(handles,labels,loc=0,markerscale=2)
		handles, labels = ax[1].get_legend_handles_labels()
		ax[1].legend(handles,labels,loc=0,markerscale=2)

	set_corrected_size(f,figsize)
	saveFile = saveFile+'.'+args['fig_extension']
	plt.savefig(saveFile,bbox_inches='tight')
	sys.stdout.write('Saved plot to file: {}\n'.format(saveFile))
############################################################

if __name__=="__main__":
	''' Maybe print a description on how to use this module'''
