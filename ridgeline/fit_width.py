#!/usr/bin/env python
# -*- coding: utf-8 -*-

from astropy.io import ascii
import os
from glob import glob
from VLBIana.ridgeline.ridgeLine import *
from VLBIana.modules.jet_calculus import *
import VLBIana.modules.fit_functions as ff
from VLBIana.modules.plot_functions import *


plt.ioff()

def rl_fit(mapFile,ridgeLine,saveFile,label,theta,logFile,shiftFile=False,fit='Powerlaw',asize=5,plot_hist=False,write_fit_info=False,plot_res=True,plot_ang=True,plot_flux=False,fig_size='aanda*',add_rl=False,fig_extension='pdf',incl=90,binRidgeLine=False,fit_both_jets_together=False,errorFile=False,write_tex=False,zcut=False,add_points=False,plot_lines=False):
    ''' To fit the width of a jet
    '''
    if type(ridgeLine)==str:
        mapFile     = [mapFile]
        ridgeLine   = [ridgeLine]
        label           = [label]
        logFile     = [logFile]
        if errorFile:
            errorFile       = [errorFile]

    else:
        pass

    if add_rl:
        mapFile.append(add_rl[0])
        ridgeLine.append(add_rl[1])
        logFile.append(add_rl[2])
        label.append(add_rl[3])
        if errorFile:
            errorFile.append(add_rl[4])
        n = len(ridgeLine)
        nn = 8
        colors[n-1]=[0.9,0.4,0.4,1]
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
    theta = theta*np.pi/180

    if shiftFile:
        shiftRA,shiftDEC = readShift(shiftFile)
        DistStr = 'Distshift'
        RAStr = 'RAshift'
        DecStr = 'Decshift'

        Ridge = [RidgeLine(RL,log,cmap,shift=(shiftra,shiftdec),incl=incl) for RL,log,cmap,shiftra,shiftdec in zip(ridgeLine,logFile,mapFile,shiftRA,shiftDEC)]
    else:
        DistStr = 'Dist'
        RAStr = 'RA'
        DecStr = 'Dec'
        Ridge = [RidgeLine(RL,log,cmap,incl=incl) for RL,log,cmap in zip(ridgeLine,logFile,mapFile)]

    if errorFile:
        sys.stdout.write('Use errorfile')
        ridgelines = [RL.readRidgeline(theta,widthErr=eF) for RL,eF in zip(Ridge,errorFile)]
    else:
        ridgelines = [RL.readRidgeline(theta) for RL in Ridge]
    if add_rl:
        ridgelines[-1]= ridgelines[-1][np.abs(ridgelines[-1][DistStr])>=0.2]
    if binRidgeLine:
        print('use binned ridge line')
        ridgelines = [RL.binRidgeLine() for RL in Ridge]
    i=0
    Jet,CJet = [],[]
    if zcut:
        for i,rl in enumerate(ridgelines):
            ridgelines[i] = rl[np.abs(rl[DistStr])<=2]

    for rl in ridgelines:
        Jet.append(rl[np.sign(rl[DistStr])==1])
        CJet.append(rl[np.sign(rl[DistStr])==-1])

        i+=1

    fwhmE= 'FWHMDeconvolvedErr'
    fwhm    = 'FWHMDeconvolved'
#
    DistJ,DistCJ,fwhmJ,fwhmCJ,fwhmEJ,fwhmECJ,fluxJ,fluxCJ=[],[],[],[],[],[],[],[]
    i=0
    for J,CJ in zip(Jet,CJet):
        DistJ.append(J[DistStr].copy())
        fwhmJ.append(J[fwhm].copy())
        fwhmEJ.append(J[fwhmE].copy())
        fluxJ.append(J['Peak'].copy())
        DistCJ.append(np.abs(CJ[DistStr].copy()))
        fwhmCJ.append(CJ[fwhm].copy())
        fwhmECJ.append(CJ[fwhmE].copy())
        fluxCJ.append(CJ['Peak'].copy())

        i+=1
#
########################################
#### make fits ###
    X   = [rl[DistStr] for rl in ridgelines]
    Y   = [rl[fwhm] for rl in ridgelines]
    DY= [rl[fwhmE] for rl in ridgelines]

    if fit=='Powerlaw':
        saveFile+='_powerlaw'
        if add_rl:
            sys.stdout.write('Not using additional ridgeline for fitting\n')
            betaJ,sd_betaJ,chi2J,out = ff.fit_pl(DistJ[:-1],fwhmJ[:-1],fwhmEJ[:-1])
            betaCJ,sd_betaCJ,chi2CJ,out = ff.fit_pl(DistCJ[:-1],fwhmCJ[:-1],fwhmECJ[:-1])
            beta,sd_beta,chi2,out = ff.fit_pl(X[:-1],Y[:-1],DY[:-1])
        else:
            betaJ,sd_betaJ,chi2J,out = ff.fit_pl(DistJ,fwhmJ,fwhmEJ)
            betaCJ,sd_betaCJ,chi2CJ,out = ff.fit_pl(DistCJ,fwhmCJ,fwhmECJ)
            beta,sd_beta,chi2,out = ff.fit_pl(X,Y,DY)
        ascii.write([['Eastern Jet','Western Jet','Both Jets'],[betaJ[1],betaCJ[1],beta[1]],[sd_betaJ[1],sd_betaCJ[1],sd_beta[1]],[chi2J,chi2CJ,chi2]],saveFile+'fit_Jet.tex',names=['Fit to','a','da','chi2'],formats={'a':'%.2f','da':'%.2f','chi2':'%.2f'}, overwrite=True,format='latex')
    elif fit=='brokenPowerlaw':
        saveFile+='_broken_powerlaw'
        if add_rl:
            sys.stdout.write('Not using additional ridgeline for fitting\n')
            betaJ,sd_betaJ,chi2J,out = ff.fit_bpl(DistJ[:-1],fwhmJ[:-1],fwhmEJ[:-1],x0=[0.3,0,1,3])
            betaCJ,sd_betaCJ,chi2CJ,out = ff.fit_bpl(DistCJ[:-1],fwhmCJ[:-1],fwhmECJ[:-1],x0=[0.1,0,1,3])
        else:
            betaJ,sd_betaJ,chi2J,out = ff.fit_bpl(DistJ,fwhmJ,fwhmEJ,x0=[0.3,0,1,2])
            betaCJ,sd_betaCJ,chi2CJ,out = ff.fit_bpl(DistCJ,fwhmCJ,fwhmECJ,x0=[0.3,0,1,2])
        beta,sd_beta,chi2,out = ff.fit_bpl(X,Y,DY,x0=[0.1,0,0.5,2])
        ascii.write([['Eastern Jet','Western Jet','Both Jets'],[betaJ[0],betaCJ[0],beta[0]],[sd_betaJ[0],sd_betaCJ[0],sd_beta[0]],[betaJ[1],betaCJ[1],beta[1]],[sd_betaJ[1],sd_betaCJ[1],sd_beta[1]],[betaJ[2],betaCJ[2],beta[2]],[sd_betaJ[2],sd_betaCJ[2],sd_beta[2]],[betaJ[3],betaCJ[3],beta[3]],[sd_betaJ[3],sd_betaCJ[3],sd_beta[3]],[chi2J,chi2CJ,chi2]],saveFile+'fit_Jet.tex',names=['Fit to','w0','dw0','ku','dku','kd','dkd','bp','dbp','chi2'],formats={'w0':'%.2f','dw0':'%.2f','ku':'%.2f','dku':'%.2f','kd':'%.2f','dkd':'%.2f','bp':'%.2f','dbp':'%.2f','chi2':'%.2f'}, overwrite=True,format='latex')


###########################################
### for plotting ###
    ymin = 6e-3
    ymax = max([max(f) for f in Y])*3
    xmin    = min (np.abs(np.concatenate(X)))
    xmax = max (np.abs(np.concatenate(X)))*1.1
    if add_points:
        if type(add_points)==dict:
            add_points = [add_points]
        for addP in add_points:
            if type(addP['dist'])==list:
                addMin = min(np.abs(addP['dist']))
            else:
                addMin = addP['dist']
            if addMin<xmin:
                xmin= addMin

    xmin = Rs2mas(100)

    xr=np.arange(xmin,xmax,0.01)
    yrmin = -1
    yrmax = 1
    figsize=(set_size(fig_size,subplots=(nsub,2)))

    f,ax=plt.subplots(nsub,2,sharex='col',sharey='row',gridspec_kw={'hspace': 0, 'wspace': 0},figsize=figsize)
    #####
    axes = ax.flatten()
    axesWidthPlot(axes[0],secxax=r'Distance from core [$R_\mathrm{S}$]')
    axesWidthPlot(axes[1],secxax=r'Distance from core [$R_\mathrm{S}$]',secyax=r'De-convolved width \small{[$R_\mathrm{S}$]}',ylabel=False)
    axes[0].annotate('Eastern Jet', xy=(0.2,0.85),xycoords='axes fraction',size=14)
    axes[1].annotate('Western Jet', xy=(0.2,0.85),xycoords='axes fraction',size=14)

    for axs in ax.flatten():
        axs.set_xscale('log')
        axs.set_yscale('log')

    if plot_flux:
        axesWidthPlot(ax[nf,0],ylabel=r'Peak Flux density \small{$\left[\mathrm{\frac{mJy}{beam}}\right]$}')
        axesWidthPlot(ax[nf,1],xlabel=r'Distance from the core [mas]')

    if plot_ang:
        axesWidthPlot(ax[na,0],ylabel=r'$\phi_\mathrm{app}$ [deg]',yscale='linear')
        axesWidthPlot(ax[na,1],yscale='linear',ylabel=False)

    if plot_res:
        axesWidthPlot(ax[nr,1],ylabel=False,yscale='linear')
        axesWidthPlot(ax[nr,0],ylabel='Residuals',yscale='linear')

    if plot_hist:
        plot_res = True
        add_subplot_unshare(ax[nh,0])
        add_subplot_unshare(ax[nh,1])


########### Now do the plotting ######
    i=0
#ax0
    if add_points:
        sys.stdout.write('plotting additional points on top')

    for rl,m in zip(ridgeLine,markers):
        print (rl)
        print(markers)
        axes[0].scatter(DistJ[i],fwhmJ[i],s=4,marker=m,label='{}'.format(label[i]))
        axes[0].errorbar(DistJ[i],fwhmJ[i],yerr=fwhmEJ[i],fmt=m,ms=0,linewidth=0,elinewidth=0.4,errorevery=1,label='_nolegend_',alpha=0.3)
        axes[1].scatter(DistCJ[i],fwhmCJ[i],s=4,marker=m)
        axes[1].errorbar(DistCJ[i],fwhmCJ[i],yerr=fwhmECJ[i],fmt=m,ms=0,linewidth=0,elinewidth=0.4,errorevery=1,label='_nolegend_',alpha=0.3)

        if plot_flux:
            ax[nf,0].scatter(DistJ[i],fluxJ[i],s=4,marker=m)
            ax[nf,0].errorbar(DistJ[i],fluxJ[i],yerr=fluxJ[i]*flux_uncertainty,fmt=m,ms=0,linewidth=0,elinewidth=0.4,errorevery=1,alpha=0.3)
            ax[nf,1].scatter(DistCJ[i],fluxCJ[i],s=4,marker=m)
            ax[nf,1].errorbar(DistCJ[i],fluxCJ[i],yerr=fluxCJ[i]*flux_uncertainty,fmt=m,ms=0,linewidth=0,elinewidth=0.4,errorevery=1,alpha=0.3)
        if plot_ang:
            angJ,angEJ = app_ang(fwhmJ[i],DistJ[i],[fwhmEJ[i]])
            angCJ,angECJ = app_ang(fwhmCJ[i],DistCJ[i],[fwhmECJ[i]])
            if write_tex:
                ascii.write([angJ,angEJ],saveFile+'Jet_app_angle.'+label[i]+'tex',names=['theta_E','Dtheta_E'],formats={'theta_E':'%.2f','Dtheta_E':'%.2f'}, overwrite=True,format='latex')
                ascii.write([angCJ,angECJ],saveFile+'CJet_app_angle'+label[i]+'.tex',names=['theta_W','Dtheta_W'],formats={'theta_W':'%.2f','Dtheta_W':'%.2f'}, overwrite=True,format='latex')

            ax[na,0].scatter(DistJ[i],angJ,s=4,marker=m)
            ax[na,1].scatter(DistCJ[i],angCJ,s=4,marker=m)
            ax[na,0].errorbar(DistJ[i],angJ,yerr=angEJ,fmt=m,ms=0,linewidth=0,elinewidth=0.4,errorevery=1,alpha=0.3)
            ax[na,1].errorbar(DistCJ[i],angCJ,yerr=angECJ,fmt=m,ms=0,linewidth=0,elinewidth=0.4,errorevery=1,alpha=0.3)

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
            if yrmax < maxr: yrmax = maxr+ maxr*0.1
            ax[nr,0].scatter(DistJ[i],resJ,s=4,marker=m)
            ax[nr,1].scatter(DistCJ[i],resCJ,s=4,marker=m)
            ax[nr,1].plot(xr,np.linspace(0,0,len(xr)),'k:',lw=1)
            ax[nr,0].plot(xr,np.linspace(0,0,len(xr)),'k:',lw=1)

        i+=1

    if add_points:
        addmarker=['<','2','d']
        addcolor = ['red','magenta','orange']
        for i,addP in enumerate(add_points):
            if addP['jet']=='Jet':
                axis=0
            elif addP['jet']=='CJet':
                axis=1
            elif addP['jet']=='ALL':
                axis='ALL'
            else:
                axis='BOTH'

            if axis=='BOTH':
                axes[0].scatter(addP['dist'],addP['fwhm'],marker=addmarker[i],color=addcolor[i],label='{}'.format(addP['label']),zorder=2)
                axes[1].scatter(addP['dist'],addP['fwhm'],marker=addmarker[i],color=addcolor[i],zorder=2)

            elif axis=='ALL':
                add_DistJ = np.array(addP['dist'])[np.sign(addP['dist'])==1]
                add_DistCJ = np.abs(np.array(addP['dist'])[np.sign(addP['dist'])==-1])
                add_fwhmJ = np.array(addP['fwhm'])[np.sign(addP['dist'])==1]
                add_fwhmCJ = np.array(addP['fwhm'])[np.sign(addP['dist'])==-1]
                ax[0,0].scatter(add_DistJ,add_fwhmJ,marker=addmarker[i],s=100,c=addcolor[i],label='{}'.format(addP['label']))
                ax[0,1].scatter(add_DistCJ,add_fwhmCJ,marker=addmarker[i],s=100,c=addcolor[i])
            else:
                ax[0,axis].scatter(addP['dist'],addP['fwhm'],marker=addmarker[i],color=addcolor[i])
            if plot_ang:
                if 'ang' in addP.keys():
                    ax[na,axis].scatter(addP['dist'],addP['ang'],marker=addmarker[i],color=addcolor[i])
            if plot_flux:
                if 'flux' in addP.keys():
                    if axis=='BOTH':
                        ax[0,0].scatter(addP['dist'],addP['flux'],marker=addmarker[i],color=addcolor[i])
                        ax[0,1].scatter(addP['dist'],addP['flux'],marker=addmarker[i],color=addcolor[i])
                    else:
                        ax[0,axis].scatter(addP['dist'],addP['flux'],marker=addmarker[i],color=addcolor[i])

    if plot_hist:
        reswJ= ff.resid(np.concatenate(fwhmJ),ff.broken_powerlaw(betaJ,np.concatenate(DistJ)),np.concatenate(fwhmEJ))
        plotHist(reswJ,ax[nn,0],plot_norm=True,xlabel='Weighted Residuals',asize=asize)
        reswCJ= ff.resid(np.concatenate(fwhmCJ),ff.broken_powerlaw(betaCJ,np.concatenate(DistCJ)),np.concatenate(fwhmECJ))
        plotHist(reswCJ,ax[nn,1],plot_norm=True,xlabel='Weighted Residuals',asize=asize)


    plot_fit(xr,fit,betaJ,sd_betaJ,chi2J,ax=axes[0],annotate=write_fit_info,asize=asize,ls='--',lw='0.8',label='broken power-law')
    plot_fit(xr,fit,betaCJ,sd_betaCJ,chi2CJ,ax=axes[1],annotate=write_fit_info,asize=asize,ls='--',lw='0.8')
    if plot_lines:
        plot_fit(np.arange(0.001,1.5,0.1),'Powerlaw',[2.5e-1,0.5],[0,0],2,ax=axes[1],color='#dede00',lw=1.5)#annotate=r'$k=0.5$',asize=asize,color='r')
        plot_fit(np.arange(1.3,18),'Powerlaw',[1e-1,1],[0,0],2,ax=axes[1],color='#e41a1c',lw=1.5)#annotate=r'$k=1$',asize=asize,color='b')
        plot_fit(np.arange(0.2,8),'Powerlaw',[0.24,0],[0,0],2,ax=axes[1],color='#ff7f00',lw=1.5)#annotate=r'$k=1$',asize=asize,color='b')
        plot_fit(np.arange(0.001,1.5,0.1),'Powerlaw',[2.5e-1,0.5],[0,0],2,ax=axes[0],color='#dede00',lw=1.5,label='parabolic')#annotate=r'$k=0.5$',asize=asize,color='r')
        plot_fit(np.arange(1.3,18),'Powerlaw',[1e-1,1],[0,0],2,ax=axes[0],color='#e41a1c',lw=1.5,label='conical')#annotate=r'$k=1$',asize=asize,color='b')
        plot_fit(np.arange(0.2,8),'Powerlaw',[0.24,0],[0,0],2,ax=axes[0],color='#ff7f00',lw=1.5,label='cylindric')#annotate=r'$k=1$',asize=asize,color='b')


### finalize axis settings
    handles, labels = axes[0].get_legend_handles_labels()

    if plot_lines:
        axes[0].legend(handles,labels,loc='lower left',ncol=4,markerscale=1,labelspacing=0.1,bbox_to_anchor=(0.0, 1.3,2.,.4),mode='expand',borderaxespad=0.,handletextpad=0.1)
    else:
        axes[0].legend(handles,labels,loc='lower left',ncol=7,markerscale=1,labelspacing=0.1,bbox_to_anchor=(0.0, 1.3,2.,.4),mode='expand',borderaxespad=0.,handletextpad=0.1)

    for axs in ax.flat:
        axs.set_xlim(xmin,xmax)
        axs.set_ylim(ymin,ymax)

    if plot_res:
        ax[nr,0].set_ylim(yrmin,yrmax)
        ax[nr,1].set_ylim(yrmin,yrmax)

    if plot_flux:
        ymina = 3e-1
        ymaxa = 1e3*1.1
        ax[nf,0].set_ylim(ymina,ymaxa)
        ax[nf,1].set_ylim(ymina,ymaxa)
    if plot_ang:
        ax[na,0].set_ylim(-2,175)
        ax[na,1].set_ylim(-2,175)

    if plot_hist:
        ax[nr,0].autoscale()
        ax[nr,1].autoscale()

    for aa in ax.flat:
        aa.minorticks_on()
        aa.tick_params(which='both',direction='inout')
        aa.label_outer()

    axes[0].invert_xaxis()
    #set_corrected_size(f,figsize)
    saveFile = saveFile+'.'+fig_extension
    plt.savefig(saveFile,bbox_inches='tight',transparent=True)
    plt.close()
    sys.stdout.write('Saved plot to file: {}\n'.format(saveFile))
####################################
#### plot both jets into one figure
    if fit_both_jets_together:
        yrmin = -1
        yrmax = 1
        Jcolor  = [0.1,0.4,0.6,1]
        CJcolor = [0.45,0.67,0.43,0.7]
        Fcolor  = 'k'
        if add_rl:
            xJ  = np.abs(np.concatenate(DistJ[:-1]))
            xCJ = np.abs(np.concatenate(DistCJ[:-1]))
            yJ  = np.concatenate(fwhmJ[:-1])
            yCJ = np.concatenate(fwhmCJ[:-1])
            yJE = np.concatenate(fwhmEJ[:-1])
            yCJE    = np.concatenate(fwhmECJ[:-1])
        else:
            xJ  = np.abs(np.concatenate(DistJ))
            xCJ = np.abs(np.concatenate(DistCJ))
            yJ  = np.concatenate(fwhmJ)
            yCJ = np.concatenate(fwhmCJ)
            yJE = np.concatenate(fwhmEJ)
            yCJE    = np.concatenate(fwhmECJ)

        if fig_size == 'aanda*':
            fig_size = 'aanda'
        f,ax=plt.subplots(nsub,1,sharex='col',sharey='row',gridspec_kw={'hspace': 0, 'wspace': 0})
        figsize=(set_size(fig_size,subplots=(nsub,1)))

        if plot_hist:
            add_subplot_unshare(ax[nh])


        axesWidthPlot(ax[0],secxax='Distance from the core [$R_\mathrm{S}$]',secyax='De-convolved Jet width [$R_\mathrm{S}$]')
        if plot_res:
            axesWidthPlot(ax[nr],ylabel='Residuals',yscale='linear')
        if plot_ang:
            axesWidthPlot(ax[na],ylabel=r'$\phi_\mathrm{app}$ [deg]',yscale='linear')

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
        ax[0].legend(handles,labels,loc=0,markerscale=1.5,prop={'size':asize})
        ax[0].set_ylim(ymin,ymax)
        ax[0].set_xlim(xmin,xmax)
        for aa in ax.flat:
            aa.minorticks_on()
            aa.tick_params(which='both',direction='inout')
            aa.label_outer()


        saveFile = saveFile.split('.')[0]+'_both_jets'
        saveFile = saveFile+'.'+fig_extension
        set_corrected_size(f,figsize)
        plt.savefig(saveFile,bbox_inches='tight',transparent=True)
        plt.close()
        sys.stdout.write('Saved plot to file: {}\n'.format(saveFile))

#########################################
##################################
#################

if __name__=="__main__":
    ''' think about maybe printing info about the modules itself '''
