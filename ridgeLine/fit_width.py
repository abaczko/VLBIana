#!/usr/bin/env python
# -*- coding: utf-8 -*-

import VLBIana.modules.fit_functions as ff 
from VLBIana.modules.plot_functions import *
from astropy.io import ascii
import os
from glob import glob
from ridgeLine import *
from VLBIana.modulesjet_calculus import *

#plt.style.use('talkstyle')
def rl_fit(mapFile,ridgeLine,shiftFile,saveFile,label,theta,logFile,fit='Powerlaw',asize=7,plot_hist=False,write_fit_info=True,plot_res=True,plot_ang=True,plot_flux=False,fig_size='aanda*',add_rl=False,fig_extension='pdf',incl=90,binRidgeLine=False,fit_both_jets_together=True,errorFile=False,write_tex=False):
	''' To fit the width of a jet 
	'''
	if type(ridgeLine)==str:
		mapFile		= [mapFile]
		ridgeLine	= [ridgeLine]
		label			= [label]
		logFile		= [logFile]
		if errorFile:
			errorFile		= [errorFile]

	else:
		pass
	#	print('files \n{}\n will be fitted'.format(ridgeLine))

	if add_rl:
		mapFile.append(add_rl[0])
		ridgeLine.append(add_rl[1])
		logFile.append(add_rl[2])
		label.append(add_rl[3])
		if errorFile:
			errorFile.append(add_rl[4])
		n = len(ridgeLine)
		colors = cmap(np.linspace(0,0.95,n))
		colors = colors[len(colors)-n:]
		colors[n-1]=[0.8,0.5,0.5,1]
		mpl.rcParams['axes.prop_cycle'] = cycler(color=colors)

# get number of subplots
	nsub=0
	if plot_flux:
		nsub+=1
		nf = nsub
	if plot_ang:
		nsub+=1
		na = nsub
	if plot_res:
		nsub+=1
		nr=nsub
	if plot_hist:
		nsub+=1
		nh=nsub
	
	nsub+=1

	flux_uncertainty = 0.15
	shiftRA,shiftDEC = readShift(shiftFile)
	theta = theta*np.pi/180

	Ridge = [RidgeLine(RL,log,cmap,shift=(shiftra,shiftdec),incl=incl) for RL,log,cmap,shiftra,shiftdec in zip(ridgeLine,logFile,mapFile,shiftRA,shiftDEC)]
	if errorFile:
		sys.stdout.write('Use errorfile')
		ridgelines = [RL.readRidgeline(theta,widthErr=eF) for RL,eF in zip(Ridge,errorFile)]
	else:
		ridgelines = [RL.readRidgeline(theta) for RL in Ridge]
	if add_rl:
		ridgelines[-1]= ridgelines[-1][np.abs(ridgelines[-1]['Distshift'])>=0.2]
	if binRidgeLine:
		print('use binned ridge line')
		ridgelines = [RL.binRidgeLine() for RL in Ridge]
	i=0
	Jet,CJet = [],[]
	for rl in ridgelines:
		Jet.append(rl[np.sign(rl['Distshift'])==1])
		CJet.append(rl[np.sign(rl['Distshift'])==-1])
		i+=1
	
	fwhmE= 'FWHMDeconvolvedErr'
	fwhm	= 'FWHMDeconvolved'
#
	DistJ,DistCJ,fwhmJ,fwhmCJ,fwhmEJ,fwhmECJ,fluxJ,fluxCJ=[],[],[],[],[],[],[],[]
	i=0
	for J,CJ in zip(Jet,CJet):
		DistJ.append(J['Distshift'].copy())
		fwhmJ.append(J[fwhm].copy())
		fwhmEJ.append(J[fwhmE].copy())
		fluxJ.append(J['Peak'].copy())
		DistCJ.append(np.abs(CJ['Distshift'].copy()))
		fwhmCJ.append(CJ[fwhm].copy())
		fwhmECJ.append(CJ[fwhmE].copy())
		fluxCJ.append(CJ['Peak'].copy())

		i+=1
#
########################################
#### make fits ###
	X	= [rl['Distshift'] for rl in ridgelines] 
	Y	= [rl[fwhm] for rl in ridgelines] 
	DY= [rl[fwhmE] for rl in ridgelines] 

	if fit=='Powerlaw':
		saveFile+='_powerlaw_'
		if add_rl:
			betaJ,sd_betaJ,chi2J,out = ff.fit_pl(DistJ[:-1],fwhmJ[:-1],fwhmEJ[:-1])
			betaCJ,sd_betaCJ,chi2CJ,out = ff.fit_pl(DistCJ[:-1],fwhmCJ[:-1],fwhmECJ[:-1])
			beta,sd_beta,chi2,out = ff.fit_pl(X[:-1],Y[:-1],DY[:-1])
		else:
			betaJ,sd_betaJ,chi2J,out = ff.fit_pl(DistJ,fwhmJ,fwhmEJ)
			betaCJ,sd_betaCJ,chi2CJ,out = ff.fit_pl(DistCJ,fwhmCJ,fwhmECJ)
			beta,sd_beta,chi2,out = ff.fit_pl(X,Y,DY)
		ascii.write([['Eastern Jet','Western Jet','Both Jets'],[betaJ[1],betaCJ[1],beta[1]],[sd_betaJ[1],sd_betaCJ[1],sd_beta[1]],[chi2J,chi2CJ,chi2]],saveFile+'fit_Jet.tex',names=['Fit to','a','da','chi2'],formats={'a':'%.2f','da':'%.2f','chi2':'%.2f'}, overwrite=True,format='latex')
	elif fit=='brokenPowerlaw':
		saveFile+='_broken_powerlaw_'
		if add_rl:
			betaJ,sd_betaJ,chi2J,out = ff.fit_bpl(DistJ[:-1],fwhmJ[:-1],fwhmEJ[:-1],x0=[0.3,0,1,3]) 
			betaCJ,sd_betaCJ,chi2CJ,out = ff.fit_bpl(DistCJ[:-1],fwhmCJ[:-1],fwhmECJ[:-1],x0=[0.1,0,1,3])
		else:
			betaJ,sd_betaJ,chi2J,out = ff.fit_bpl(DistJ,fwhmJ,fwhmEJ,x0=[0.3,0,1,2]) 
			betaCJ,sd_betaCJ,chi2CJ,out = ff.fit_bpl(DistCJ,fwhmCJ,fwhmECJ,x0=[0.3,0,1,2])
		beta,sd_beta,chi2,out = ff.fit_bpl(X,Y,DY,x0=[0.1,0,0.5,2])
		ascii.write([['Eastern Jet','Western Jet','Both Jets'],[betaJ[0],betaCJ[0],beta[0]],[sd_betaJ[0],sd_betaCJ[0],sd_beta[0]],[betaJ[1],betaCJ[1],beta[1]],[sd_betaJ[1],sd_betaCJ[1],sd_beta[1]],[betaJ[2],betaCJ[2],beta[2]],[sd_betaJ[2],sd_betaCJ[2],sd_beta[2]],[betaJ[3],betaCJ[3],beta[3]],[sd_betaJ[3],sd_betaCJ[3],sd_beta[3]],[chi2J,chi2CJ,chi2]],saveFile+'fit_Jet.tex',names=['Fit to','w0','dw0','ku','dku','kd','dkd','bp','dbp','chi2'],formats={'w0':'%.2f','dw0':'%.2f','ku':'%.2f','dku':'%.2f','kd':'%.2f','dkd':'%.2f','bp':'%.2f','dbp':'%.2f','chi2':'%.2f'}, overwrite=True,format='latex')


###########################################
### for plotting ###
	ymin = 1e-3
	ymax = max([max(f) for f in Y])*3
	xmin	= min (np.abs(np.concatenate(X)))
	xmax = max (np.abs(np.concatenate(X)))*1.1
	xr=np.arange(xmin,xmax,0.01)
	yrmin = -1
	yrmax = 1

	f,ax=plt.subplots(nsub,2,figsize=(set_size(fig_size,subplots=(nsub,2))),sharex='col',sharey='row',gridspec_kw={'hspace': 0, 'wspace': 0})
	axesWidthPlot(ax[0,0],secxax=r'Distance from core [$R_\mathrm{S}$]',secyax=r'Jet width \small{[$R_\mathrm{S}$]}')
	axesWidthPlot(ax[0,1],secxax=r'Distance from core [$R_\mathrm{S}$]',secyax=r'Jet width \small{[$R_\mathrm{S}$]}',ylabel=False)
	for axs in ax.flatten():
		axs.set_xscale('log')
		axs.set_yscale('log')

	if plot_flux:
		axesWidthPlot(ax[nf,0],ylabel=r'Peak Flux density \small{[$\mathrm{\frac{mJy}{beam}}$]}')

	if plot_ang:
		axesWidthPlot(ax[na,0],ylabel=r'$\phi_{app}$ [deg]',yscale='linear')
		axesWidthPlot(ax[na,1],yscale='linear',ylabel=False)
		
	if plot_res:
		axesWidthPlot(ax[nr,1],ylabel=False,yscale='linear')
		axesWidthPlot(ax[nr,0],ylabel='Residuals',yscale='linear')

	if plot_hist:
		plot_res = True
		add_subplot_unshare(ax[nh,0])
		add_subplot_unshare(ax[nh,1])
#	if plot_res:
#		axesWidthPlot(ax[1,0],ylabel=r'Peak Flux density \small{[$\mathrm{\frac{mJy}{beam}}$]}', xlabel=False)
#	else:
#		axesWidthPlot(ax[1,0],ylabel=r'Peak Flux density \small{[$\mathrm{\frac{mJy}{beam}}$]}')
#		axesWidthPlot(ax[1,1],ylabel=False)

#	ax[1,1].set_xscale('log')
#	ax[1,1].set_yscale('log')


	ax[0,0].annotate('Eastern Jet', xy=(0.7,0.9),xycoords='axes fraction',size=asize+1)
	ax[0,1].annotate('Western Jet', xy=(0.7,0.9),xycoords='axes fraction',size=asize+1)


########### Now do the plotting ######
	plot_fit(xr,fit,betaJ,sd_betaJ,chi2J,ax=ax[0,0],annotate=write_fit_info)
	plot_fit(xr,fit,betaCJ,sd_betaCJ,chi2CJ,ax=ax[0,1],annotate=write_fit_info)

	i=0
#ax0
	for rl,m in zip(ridgeLine,markers):
		ax[0,0].scatter(DistJ[i],fwhmJ[i],s=4,marker=m,label='{}'.format(label[i]))
		ax[0,0].errorbar(DistJ[i],fwhmJ[i],yerr=fwhmEJ[i],fmt='none',elinewidth=0.4,errorevery=1,label='_nolegend_',alpha=0.3)
		ax[0,1].scatter(DistCJ[i],fwhmCJ[i],s=4,label='{}'.format(label[i]),marker=m)
		ax[0,1].errorbar(DistCJ[i],fwhmCJ[i],yerr=fwhmECJ[i],fmt='none',elinewidth=0.4,errorevery=1,label='_nolegend_',alpha=0.3)	
		if plot_flux:
			ax[nf,0].scatter(DistJ[i],fluxJ[i],s=4,marker=m,label='{}'.format(label[i]))
			ax[nf,0].errorbar(DistJ[i],fluxJ[i],yerr=fluxJ[i]*flux_uncertainty,fmt='none',elinewidth=0.4,errorevery=1,label='_nolegend_',alpha=0.3)
			ax[nf,1].scatter(DistCJ[i],fluxCJ[i],s=4,label='{}'.format(label[i]),marker=m)
			ax[nf,1].errorbar(DistCJ[i],fluxCJ[i],yerr=fluxCJ[i]*flux_uncertainty,fmt='none',elinewidth=0.4,errorevery=1,label='_nolegend_',alpha=0.3)
		if plot_ang:
			angJ,angEJ = app_ang(fwhmJ[i],DistJ[i],[fwhmEJ[i]])
			angCJ,angECJ = app_ang(fwhmCJ[i],DistCJ[i],[fwhmECJ[i]])
			if write_tex:
				ascii.write([angJ,angEJ],saveFile+'Jet_app_angle.'+label[i]+'tex',names=['theta_E','Dtheta_E'],formats={'theta_E':'%.2f','Dtheta_E':'%.2f'}, overwrite=True,format='latex')
				ascii.write([angCJ,angECJ],saveFile+'CJet_app_angle'+label[i]+'.tex',names=['theta_W','Dtheta_W'],formats={'theta_W':'%.2f','Dtheta_W':'%.2f'}, overwrite=True,format='latex')

			ax[na,0].scatter(DistJ[i],angJ,s=4,label='{}'.format(label[i]),marker=m)
			ax[na,1].scatter(DistCJ[i],angCJ,s=4,label='{}'.format(label[i]),marker=m)
			ax[na,0].errorbar(DistJ[i],angJ,yerr=angEJ,fmt='none',elinewidth=0.4,errorevery=1,label='_nolegend_',alpha=0.3)	
			ax[na,1].errorbar(DistCJ[i],angCJ,yerr=angECJ,fmt='none',elinewidth=0.4,errorevery=1,label='_nolegend_',alpha=0.3)	


		if plot_res:
			if fit == 'Powerlaw':
				yymJ=ff.powerlaw(betaJ,DistJ[i])
				yymCJ=ff.powerlaw(betaCJ,DistCJ[i])
			elif fit=='brokenPowerlaw':
				yymJ=ff.broken_powerlaw(betaJ,DistJ[i])
				yymCJ=ff.broken_powerlaw(betaCJ,DistCJ[i])

			resJ = ff.resid(fwhmJ[i],yymJ)
			resCJ = ff.resid(fwhmCJ[i],yymCJ)

		#resJ = ff.fractdiff(fwhmJ[i],yymJ)
			#resCJ = ff.fractdiff(fwhmCJ[i],yymCJ)
			minr = min([min(resJ),min(resCJ)])
			maxr = max([max(resJ),max(resCJ)])
			if yrmin > minr: yrmin = minr
			if yrmax < maxr: yrmax = maxr
			ax[nr,0].scatter(DistJ[i],resJ,s=4,marker=m)
			ax[nr,1].scatter(DistCJ[i],resCJ,s=4,marker=m)
			ax[nr,1].plot(xr,np.linspace(0,0,len(xr)),'k:',lw=1)
			ax[nr,0].plot(xr,np.linspace(0,0,len(xr)),'k:',lw=1)
		i+=1

	if plot_hist:
		reswJ= ff.resid(np.concatenate(fwhmJ),ff.broken_powerlaw(betaJ,np.concatenate(DistJ)),np.concatenate(fwhmEJ))
		plotHist(reswJ,ax[nn,0],plot_norm=True,xlabel='Weighted Residuals',asize=asize)
		reswCJ= ff.resid(np.concatenate(fwhmCJ),ff.broken_powerlaw(betaCJ,np.concatenate(DistCJ)),np.concatenate(fwhmECJ))
		plotHist(reswCJ,ax[nn,1],plot_norm=True,xlabel='Weighted Residuals',asize=asize)



### finalize axis settings
	handles, labels = ax[0,0].get_legend_handles_labels()

	ax[0,0].legend(handles,labels,loc=2,markerscale=1.5,labelspacing=0.1)
	ax[0,1].legend(handles,labels,loc=0,markerscale=1.5,labelspacing=0.1)
	for axs in ax.flat:
		axs.set_xlim(xmin,xmax)
		axs.set_ylim(ymin,ymax)
	if plot_res:		
		ax[nr,0].set_ylim(yrmin,yrmax)
		ax[nr,1].set_ylim(yrmin,yrmax)

	if plot_flux:
		ymina = 1e-1
		ymaxa = 1e3*1.1
		ax[nf,0].set_ylim(ymina,ymaxa)
		ax[nf,1].set_ylim(ymina,ymaxa)
	if plot_ang:
		ax[na,0].set_ylim(0,175)
		ax[na,1].set_ylim(0,175)

	if plot_hist:
		ax[nr,0].autoscale()
		ax[nr,1].autoscale()
	saveFile+='width+flux_fit'
	saveFile = saveFile+'.'+fig_extension
	plt.savefig(saveFile,bbox_inches='tight')
	sys.stdout.write('Saved plot to file: {}\n'.format(saveFile))
####################################
#### plot both jets into one figure
	if fit_both_jets_together:
		yrmin = -1
		yrmax = 1
		Jcolor	= [0.1,0.4,0.6,1]
		CJcolor = [0.45,0.67,0.43,0.7]
		Fcolor	= 'k'
		if add_rl:
			xJ	= np.abs(np.concatenate(DistJ[:-1]))
			xCJ = np.abs(np.concatenate(DistCJ[:-1]))
			yJ	= np.concatenate(fwhmJ[:-1])
			yCJ	= np.concatenate(fwhmCJ[:-1])
			yJE	= np.concatenate(fwhmEJ[:-1])
			yCJE	= np.concatenate(fwhmECJ[:-1])
		else:
			xJ	= np.abs(np.concatenate(DistJ))
			xCJ = np.abs(np.concatenate(DistCJ))
			yJ	= np.concatenate(fwhmJ)
			yCJ	= np.concatenate(fwhmCJ)
			yJE	= np.concatenate(fwhmEJ)
			yCJE	= np.concatenate(fwhmECJ)
	
		if fig_size == 'aanda*':
			fig_size = 'aanda'
		f,ax=plt.subplots(nsub,1,figsize=(set_size(fig_size,subplots=(nsub,1))),sharex='col',sharey='row',gridspec_kw={'hspace': 0, 'wspace': 0})

		if plot_hist:
			add_subplot_unshare(ax[nh])

	
		axesWidthPlot(ax[0],secxax='Distance along jet [$R_\mathrm{S}$]',secyax='De-convolved Jet width [$R_\mathrm{S}$]')
		if plot_res:
			axesWidthPlot(ax[nr],ylabel='Residuals',yscale='linear')
		if plot_ang:
			axesWidthPlot(ax[na],ylabel=r'$\phi_{app}$ [deg]',yscale='linear')

		if fit=='Powerlaw':
			plot_fit(xr,'Powerlaw',beta,sd_beta,chi2,ax=ax[0],annotate=write_fit_info)
			yymJ=ff.powerlaw(beta,xJ)
			yymCJ=ff.powerlaw(beta,xCJ)
			ax[0].plot(xr,ff.powerlaw(beta,xr),color=Fcolor,lw=1)
		elif fit=='brokenPowerlaw':
			plot_fit(xr,'brokenPowerlaw',beta,sd_beta,chi2,ax=ax[0],annotate=write_fit_info)
			yymJ=ff.broken_powerlaw(beta,xJ)
			yymCJ=ff.broken_powerlaw(beta,xCJ)
			ax[0].plot(xr,ff.broken_powerlaw(beta,xr),color=Fcolor,lw=1)
	
	
		ax[0].scatter(xJ,yJ,color=Jcolor,s=4,marker='o',label='{}'.format('eastern Jet'))
		ax[0].errorbar(xJ,yJ,yerr=yJE,fmt='none',color=Jcolor,elinewidth=0.2,errorevery=1,label='_nolegend_',alpha=0.2)	
		ax[0].scatter(xCJ,yCJ,color=CJcolor,s=4,marker='x',label='{}'.format('western Jet'))
		ax[0].errorbar(xCJ,yCJ,yerr=yCJE,fmt='none',color=CJcolor,elinewidth=0.2,errorevery=1,label='_nolegend_',alpha=0.2)	
		if plot_ang:
			angJ,angJE = app_ang(yJ,xJ,[yJE])
			angCJ,angCJE = app_ang(yCJ,xCJ,[yCJE])

			ax[na].scatter(xJ,angJ,s=4,label='_nolegend',marker='o',color=Jcolor)
			ax[na].scatter(xCJ,angCJ,s=4,label='_nolegend',marker='x',color=CJcolor)
			ax[na].errorbar(xJ,angJ,yerr=angJE,fmt='none',color=Jcolor,elinewidth=0.4,label='_nolegend_',alpha=0.3)	
			ax[na].errorbar(xCJ,angCJ,yerr=angCJE,fmt='none',color=CJcolor,elinewidth=0.4,label='_nolegend_',alpha=0.3)	


		if plot_res:
			resJ = ff.resid(yJ,yymJ)
			resCJ = ff.resid(yCJ,yymCJ)
			resw= ff.resid(np.concatenate([yJ,yCJ]),np.concatenate([yymJ,yymCJ]),np.concatenate([yJE,yCJE]))
			#reswCJ= ff.resid(yCJ,yymCJ,yCJE)
		
		
			ax[nr].scatter(xJ,resJ,marker='o',color=Jcolor,s=4)
			ax[nr].scatter(xCJ,resCJ,marker='x',color=CJcolor,s=4)
			ax[nr].plot(xr,np.linspace(0,0,len(xr)),'k:',lw=1)
			if yrmin > np.min([min(resJ),min(resCJ)]): 
				yrmin = np.min([min(resJ),min(resCJ)])
			if yrmax < np.max([max(resJ),max(resCJ)]): 
				yrmax = np.max([max(resJ),max(resCJ)])
			ax[nr].set_ylim(yrmin,yrmax)
			ax[nr].set_xlim(xmin,xmax)
	
		if plot_hist:
			plotHist(resw,ax[ph],plot_norm=True,color='red',xlabel='Weighted Residuals',asize=asize)
		
		handles, labels = ax[0].get_legend_handles_labels()
		ax[0].legend(handles,labels,loc=0,markerscale=1.5)
		ax[0].set_ylim(ymin,ymax)
		ax[0].set_xlim(xmin,xmax)
		saveFile = saveFile.split('.')[0]+'_both_jets'
		saveFile = saveFile+'.'+fig_extension
		plt.savefig(saveFile,bbox_inches='tight')
		sys.stdout.write('Saved plot to file: {}\n'.format(saveFile))

#########################################
##################################
#################

if __name__=="__main__":
	''' think about maybe printing info about the modules itself '''
