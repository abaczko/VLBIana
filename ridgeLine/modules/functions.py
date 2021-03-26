#!/usr/bin/env python
'''
Collections of functions to use with ParselTongue.
By Anne baczko@mpifr-bonn.mpg.de
'''
######################
# Python libraries ###
######################
from io import StringIO
#import time
import os,sys,re,datetime
import aips_tasks as AT
import helper_functions as HF
import numpy as np
import __main__ as main
from itertools import islice
from astropy.table import Table
from astropy.time import Time
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
#print('OP for main calibration are loaded')
from extract_info_sn import *
import logging

#####################################
local_dir = os.getcwd()
now = datetime.datetime.now()
date =str(now.day)+'_'+str(now.month)+'_'+str(now.year)
logger = logging.getLogger(__name__)
aips_out=local_dir+'/aips_out_'+date

#
######################################
#

def print_history(uvdata,writeout=False):
	'''
	To print the history lines to a variable
	'''
	history=uvdata.history
	if writeout:
		of = open(aips_out+'history.txt','w')
		for row in history:
			of.write(row+'\n')
		of.close()
	else:
		return history
#
def data_exists(uvdata):
	'''
	Test whether data loaded really exists.
	AIPSUV() will not give an error itself if the data tried to be loaded does nto exist.
	'''
	if uvdata.exists()==True:
		logger.info('Data (%s,%s,%d,%d) exists. Assume it is already TB sorted (INDXR has been run).\n',uvdata.name,uvdata.klass,uvdata.disk,uvdata.seq)
	else:
		logger.error("Something went wrong Data does not exist")
		sys.exit()
	logger.info("Data Loaded: (%s, %s, %d, %d)\n", uvdata.name,uvdata.klass,uvdata.disk,uvdata.seq)
#

def print_basic_infos(uvdata):
	'''
	Will run a collection of tasks to get back basic information about the dataset:
	PRTAN,LISTR,SNPLT of CL1,IMHEAD,POSSM of all data using CL1
	'''
	logger.info('Basic Infos are plotted/printed using the following tasks/verbs:\n PRTAN\n LISTR \n SNPLT \n IMHEAD \n POSSM \n')
	AT.prtan(uvdata)
	AT.listr(uvdata)
	AT.snplt(uvdata,inv=1)
	AT.imhead(uvdata)
	AT.possm(uvdata)
#
#
def max_ch(uvdata):
	'''
	returns the maximum channel number
	'''
	# to return the channel number of the experiment
	header=AT.get_header_info(uvdata)
	return header['FREQ']['naxis']
#logfile.close()

def get_antenna_number(uvdata,antenna):
	'''
	returns the antenna number, when given the antenna shortcut. 
	E.G. antenna='EF'
	'''
	antennas=uvdata.antennas
	for i in range(len(antennas)):
		if antennas[i]==antenna:
			print ('{0}:{1}'.format(antennas[i],i+1))
			return i+1
#
def get_several_listings(uvdata,calibrators):
	'''
	Prints out MATX with different parameter setups.
	'''
	AT.listr(uvdata,optype='MATX',sources=calibrators,bif=1,eif=1,dparm=[3,0])
	AT.listr(uvdata,optype='MATX',sources=calibrators,bif=1,eif=1,dparm=[4,0])
	AT.listr(uvdata,optype='MATX',sources=calibrators,bif=1,eif=1,dparm=[5,0])
	AT.listr(uvdata,optype='MATX',sources=calibrators,bif=1,eif=1,dparm=[14,0])

def get_matx(uvdata,oneif=False):
	'''
	Returns information about the array observing for each scan as printed out in MATX.

	keywords:
	oneif=True : MATX info will only be returned for the first IF
	oneif=False : MATX for all IFs will be returned

	The returned structure is: x[IFnumber]=dict('IF':ifnumber,'scan':array of all scans in that IF)
	x[IF]['scan'][scannumber][n]=dict{'matrix','start','stop','source','scan'}
	'matrix': returns the Antenna Matrix of the scan in form of a nested list
	the other keywords save main information about each scan.
	'''
	outfile = aips_out+uvdata.name+'_listr_MATX.txt' 
	IF=[]
	if oneif:
		nifs=1
	else:
		nifs=8
	if not os.path.isfile(outfile):
		AT.listr(uvdata,optype='MATX',sources=[],bif=1,eif=nifs)

	for ii in range(nifs):
		starttime=''
		with open(outfile,'r') as infile:
			for row in infile:
				if row.find('IF') != -1:
					IF.append({'IF':row[row.find('IF')+5],'scan':[]})
					ii+=1
					ss=0
				if row.startswith('Time'):
					skip=False
					if row.split()[2]==starttime:
						skip=True
					else:
						starttime= row.split()[2]
						IF[ii-1]['scan'].append({'scan':ss+1})
						IF[ii-1]['scan'][ss]['start']= row.split()[2]
						IF[ii-1]['scan'][ss]['stop']= row.split()[4]
						IF[ii-1]['scan'][ss]['source']= row.split()[7]
				if row.startswith('Ant') and skip==False:
					antennas = [x.strip() for x in row.split('--')]
					IF[ii-1]['scan'][ss]['matrix']=[antennas]
					head=[row][0]
				try:
					if row[3]=='|':
						index=[head.find(x) for x in antennas]
						datarow=[row.split()[0][:-1]]
						[datarow.append(row[i-1:i+3].strip()) for i in index[1:]]
						datarow=[('0' if not x else x) for x in datarow]
						datarow[1:]=[int(x) for x in datarow[1:]]
						IF[ii-1]['scan'][ss]['matrix'].append(datarow)
						if row[1:3]==antennas[-1]:
							ss+=1
							continue
				except:
					continue
	return IF
#

def get_scans(uvdata):
	'''
	Get back an astropy Table of all scans
	'''
	if not os.path.isfile(aips_out+uvdata.name+'_listr_SCAN.txt'):
		AT.listr(uvdata)
	else:
		pass
	lines=[]
	with open(aips_out+uvdata.name+'_listr_SCAN.txt','r') as f:
		for row in islice(f,0,5,1):
			pass
		header=next(f)
		for row in f:
			if row.startswith('\x0c'):
				break
			lines.append(row)
	
	header1=header.split()
	header1[-1]=' '.join([header1[-2],header1[-1]])
	header1.pop(-2)
	header1[-2]=' '.join([header1[-3],header1[-2]])
	header1.pop(-3)
	ind=[0,5,21,29,37,41,66,73,81,-1]
	l1=[]
	for l in lines:
		l1.append([l[ind[i]:ind[i+1]].strip() for i in range(len(ind)-1)])
	scans=Table(rows=l1,names=header1,masked=True)
	scans['Time_start']=[tt.split('-')[0].strip() for tt in scans['Timerange']]
	scans['Time_end']=[tt.split('-')[1].strip() for tt in scans['Timerange']]
	del scans['Timerange']
	return scans

def scantime(uvdata,scan):
	'''
	returns the timer of the specified scan as used for e.g. fring.timer
	'''
	scans=get_scans(uvdata)
	time1,time2=[0],[0]
	time1.extend([int(x) for x in scans[scan-1]['Time_start'].split('/')[1].split(':')])
	time2.extend([int(x) for x in scans[scan-1]['Time_end'].split('/')[1].split(':')])
	time1.extend(time2)
	return time1
#
def derive_solint(uvdata,timer,refant,gainu,solint=[0.05,4],snr_cut=5,dparm=[1,0],plotname=''):
	'''
	Try to find the best solution interval for global fring.
	The function returns a plot of solution interval versus SNR.
	Probably it can be used to let a global fring run automatically.
	The function uses the additional module extract_info_sn.py

	keywords:
	timer: timerange on a short scan
	refant,gainu,snr_cut,dparm: parameters to be used for fring run
	solint: range of solution intervals to be tested
	plotname: name of the final plot file
	'''
	aparm=[2,0,0,0,1,2,snr_cut,0]
	solints=[float(x) for x in np.arange(solint[0],solint[1],0.05)]
	snv=[]
	for solint in solints:
		AT.fring_global(uvdata,timer=timer,cals=[],refant=refant,gainu=gainu,solint=solint,solsub=0,solmin=0,aparm=aparm,dparm=dparm)
		snv.append(uvdata.table_highver('SN'))
#	ant=[2,8,11,15]
	#snv=range(10,88)
	fring_sol=[extract_info_sn(uvdata,i,minSNR=5.0) for i in snv]
	fring_sol= [fring_sol[i] for i in range(len(fring_sol)) if fring_sol[i]]
	solints=[solints[i] for i in range(len(fring_sol)) if fring_sol[i]]
	#
	ant1=[x['Pol1']['antennas'] for x in fring_sol]
	ant2=[x['Pol2']['antennas'] for x in fring_sol]
	indx1 = [np.where(np.array(x['Pol1']['goodrefs'])) for x in fring_sol]
	indx2 = [np.where(np.array(x['Pol2']['goodrefs'])) for x in fring_sol]
	
	snr_pol1=[[] for s in range(len(solints))]
	snr_pol2=[[] for s in range(len(solints))]
	for s in range(len(solints)):
		for a in range(len(ant1[s])):
			snr_pol1[s].append([np.mean(fring_sol[s]['Pol1']['weight'][indx1[s][0][i]][indx1[s][1][j]][3]) for i,j in zip(range(len(indx1[s][0])),range(len(indx1[s][1]))) if indx1[s][0][i]==a])
		for a in range(len(ant2[s])):
			snr_pol2[s].append([np.mean(fring_sol[s]['Pol2']['weight'][indx2[s][0][i]][indx2[s][1][j]][3]) for i,j in zip(range(len(indx2[s][0])),range(len(indx2[s][1]))) if indx2[s][0][i]==a])
	
	ants1,ants2=[],[]
	[ants1.extend(x) for x in ant1]
	[ants2.extend(x) for x in ant2]
	ants1=np.unique(ants1)
	ants2=np.unique(ants2)
	cid = np.linspace(0,1,np.max([len(ants1),len(ants2)]))
	
	snr_pol1_a=[[] for i in range(len(ants1))]
	snr_pol2_a=[[] for i in range(len(ants2))]
	solint_pol1_a=[[] for i in range(len(ants1))]
	solint_pol2_a=[[] for i in range(len(ants2))]
	for s in range(len(snr_pol1)):
		for a in range(len(ant1[s])):
			for aa in range(len(ants1)):
				if ant1[s][a]==ants1[aa]:
					snr_pol1_a[aa].append(snr_pol1[s][a][0])
					solint_pol1_a[aa].append(solints[s])
	for s in range(len(snr_pol2)):
		for a in range(len(ant2[s])):
			for aa in range(len(ants2)):
				if ant2[s][a]==ants2[aa]:
					snr_pol2_a[aa].append(snr_pol2[s][a][0])
					solint_pol2_a[aa].append(solints[s])

	fig=plt.figure()
	ax= plt.subplot(111)
	for a in range(len(ants1)):
		ax.plot(solint_pol1_a[a],snr_pol1_a[a],"o",color=plt.cm.jet(cid[a]),markersize=6,label='Pol1 Ant {}'.format(ants1[a]))

	for a in range(len(ants2)):
		ax.plot(solint_pol2_a[a],snr_pol2_a[a],"^",color=plt.cm.jet(cid[a]),markersize=6,label='Pol2 Ant {}'.format(ants2[a]))
	
	ax.minorticks_on()
	box=ax.get_position()
	ax.set_position([box.x0,box.y0,box.width * 0.8,box.height])
	ax.legend(loc='upper left', bbox_to_anchor=(1,1.0))
	plt.xlabel('Solint [min]')
	plt.ylabel('SNR')
#	plt.title()
	plt.savefig(aips_out+'solint_test'+plotname+'.ps')
	for sn in snv:
		uvdata.zap_table('SN',sn)
#
def do_pccor(uvdata,antennas,calibrator,timer,refant):
	if OP.runpclod:
		uvdata.zap_table('PC',-1)
		logger.info('loading pc table\n')
		AT.pclod(uvdata,OP.pcfile)
		PC_out=AT.prtab(uvdata,'PC',1,return_of)

	cl_hv = uvdata.table_highver('CL')
	sn_hv= uvdata.table_highver('SN')
	logger.info('Running PCCOR\n')
	refant=get_antenna_number(uvdata,refant)
	AT.pccor(uvdata,refant=refant,timer=timer)
	AT.clcal(uvdata,sources=[],antennas=antennas,cals=[calibrator],gainv=cl_hv,gainu=cl_hv+1,snv=sn_hv+1,refant=refant,interpol='2PT')
	AT.possm(uvdata,sources=[calibrator],timer=timer,gainu=cl_hv+1)

def do_manual_phasecal(uvdata,calibrator,antennas_fring,antennas_clcal,refant,timer,aparm,dparm,solint,fgv,suba=0):
	'''
	Do a manual phasecal to allign between IFs.
	As input parameters list have to be given as noramlly several runs of fring are needed.
	For each lest item, the following steps are made: fring,clcal,possm of the scan before clcal and after clcal.
	'''
	cl_hv = uvdata.table_highver('CL')
	sn_hv= uvdata.table_highver('SN')
	#
	if suba==0:
		suba=[0 for i in range(len(refant))]
	for i in range(len(calibrator)):
		logger.info('#'*20)
		logger.info('\nStarting Fring for manual-phase-cal #%d\n',i)
		AT.fring_instr(uvdata,cals=[calibrator[i]],antennas=antennas_fring[i],timer=timer[i],refant=refant[i],aparm=aparm[i],dparm=dparm[i],solint=solint[i],suba=suba[i],fgv=fgv,gainu=cl_hv,snv=sn_hv+1)
		AT.clcal(uvdata,sources=[],cals=[calibrator[i]],antennas=antennas_clcal[i],gainv=cl_hv,gainu=cl_hv+1,snv=sn_hv+1,refant=refant[i],interpol='2PT',doblank=1)
		logger.info('Manual-phase-cal #%d done\nCurrent tables: CL%d, SN%d',i+1,cl_hv+1,sn_hv+1)
		cl_hv = uvdata.table_highver('CL')
		sn_hv= uvdata.table_highver('SN')
#		AT.possm(uvdata,sources=[calibrator[i]],timer=timer[i],antennas=antennas_fring[i],solint=solint[i],gainu= cl_hv,plotname='af')
#		AT.possm(uvdata,sources=[calibrator[i]],timer=timer[i],antennas=antennas_fring[i],solint=solint[i],gainu= cl_hv-1,plotname='bf')
#	AT.possm(uvdata,bpv=0,doband=-1,sources=OP.calibrators,gainu=cl_hv,aparm=[0,1,0,0,-200,200,0,0,1,0],plotname='AllCal',antennas=OP.possm_antennas)

#
def	do_global_fring(uvdata,timer,aparm,dparm,refant,cals=[[]],sources=[[]],antennas=[0],interpol=['2PT'],search=[0],solint=[2],gainu=[0],suba=[0],dofit=[0],solmin=[0],solsub=[0],dosnsmo=False,bpv=0,doband=-1,image=False,doblank=[1],dobtween=[0]):
	'''
	Applies everything normally done when running a global fring. 
	There are many things, that are done differently for each observation, so this function probably has to be adjusted.
	
	keywords:
	dosnmo=True: Smooth the SN table resulting from fring.
	'''
	cl_hv = gainu[0]
	cl_af = []
#	if get2n:
#		imfile=AIPSImage(OP.cmap_name,'CMAP',uvdata.disk,1)
#		if imfile.exists():
#			image=[OP.cmap_name,'CMAP',uvdata.disk,1]
#			logger.info('using clean image')
#		else:
#			catnr=AT.imlod(OP.cmap_name,uvdata.disk,OP.cmap_file)
#			image=AT.getndata(uvdata.disk,catnr)
	snv1 = uvdata.table_highver('SN')+1
	if type(refant) is str:
		for i in range(len(refant)):
			refant[i]=get_antenna_number(uvdata,refant[i])
	print (len(timer))
	for i,t in enumerate(timer):
		if type(refant[i]) is str:
			refant[i]=get_antenna_number(uvdata,refant[i])
		for j in range(len(search[i])):
			if type(search[i][j]) is str:
				search[i][j]=get_antenna_number(uvdata,search[i][j])
		logger.info('run Fring for timer={}'.format(i))
		if image[i]:
			logger.info('A clean model will be used')
			AT.fring_global(uvdata,antennas=antennas[i],timer=timer[i],cals=cals[i],dofit=dofit[i],refant=refant[i],suba=0,search=search[i],gainu=gainu[i],solint=solint[i],aparm=aparm[i],dparm=dparm[i],bpv=bpv,doband=doband,get2n=image[i],flux=1e-4)
		else:
			AT.fring_global(uvdata,antennas=antennas[i],timer=timer[i],cals=cals[i],dofit=dofit[i],refant=refant[i],suba=0,search=search[i],gainu=gainu[i],solint=solint[i],aparm=aparm[i],dparm=dparm[i],bpv=bpv,doband=doband)
		if dosnsmo:
			AT.snsmo(uvdata,doblank=OP.smooth_gf_doblank,dobtween=OP.smooth_gf_dobtween,smotype='VLDE',refant=refant[i],npiece=1,inv=uvdata.table_highver('SN'),outv=uvdata.table_highver('SN')+1,samptype='MWF',bparm=OP.smooth_gf_bparm,cparm=OP.smooth_gf_cparm)
		AT.snplt(uvdata,ine='SN',inv=uvdata.table_highver('SN'),optype='DELA',stokes='HALF',opcode='ALSI',plotname='dela_gf'+str(i)+'no_smooth')
		AT.snplt(uvdata,ine='SN',inv=uvdata.table_highver('SN'),optype='RATE',stokes='HALF',opcode='ALSI',plotname='rate_gf'+str(i)+'no_smooth')
		if dosnsmo:
			AT.snplt(uvdata,ine='SN',inv=uvdata.table_highver('SN'),optype='DELA',stokes='HALF',opcode='ALSI',plotname='dela_gf'+str(i)+'smooth')
			AT.snplt(uvdata,ine='SN',inv=uvdata.table_highver('SN'),optype='RATE',stokes='HALF',opcode='ALSI',plotname='rate_gf'+str(i)+'smooth')

	AT.clcal(uvdata,cals=[cals[0][0],cals[1][0]],sources=sources[0],dobtween=dobtween[0],doblank=doblank[0],timer=timer[0],antennas=dofit[0],suba=0,gainv=gainu[0],gainu=uvdata.table_highver('CL')+1,snv=snv1,inv=uvdata.table_highver('SN'),refant=refant[0],interpol=interpol[0],smotype='')

	try:
		if OP.ad_mpc_cals:
			for i in range(len(OP.ad_mpc_cals)):
				logger.info('#'*20)
				cl_hv=uvdata.table_highver('CL')
				logger.info('\nStarting Fring for manual-phase-cal #%d\n',i)
				AT.fring_instr(uvdata,cals=[OP.ad_mpc_cals[i]],antennas=OP.ad_mpc_antennas_fring[i],timer=OP.ad_mpc_timer[i],refant=OP.ad_mpc_refant[i],aparm=OP.ad_mpc_aparm[i],dparm=OP.ad_mpc_dparm[i],solint=OP.ad_mpc_solint[i],suba=0,fgv=0,gainu=cl_hv,snv=0,bpv=bpv,doband=doband)
				AT.clcal(uvdata,sources=OP.ad_mpc_sources[i],cals=[OP.ad_mpc_cals[i]],antennas=OP.ad_mpc_antennas_clcal[i],gainv=cl_hv,gainu=cl_hv+1,snv=uvdata.table_highver('SN'),refant=OP.ad_mpc_refant[i],interpol='2PT',timer=OP.ad_mpc_timer_clcal[i])
				logger.info('Manual-phase-cal #%d done\nCurrent tables: CL%d, SN%d',i+1,uvdata.table_highver('CL'),uvdata.table_highver('SN'))
	except:
		logger.warning('No additional manual phasecal done.\n Process without.')
		#sys.exit()
	return cl_af

def get_apcal_fit(uvdata,apcal_log):
	lines,ant,pol,trec,opac=[],[],[],[],[]
	with open(apcal_log,'r') as f:
		for line in f:
			if line.startswith('APCAL'):
				lines.append(line)

	for line in lines:
		if line.find('opac')!=-1:
			line=re.sub(' +',' ',line)
			temp=line.split()
			if len(temp[3])<4:
				temp[3]=temp[3]+temp[4]
				temp.pop(4)
			ant.append(temp[1])
			pol.append(temp[2])
			trec.append(temp[7])
			opac.append(temp[10].strip())
	apcal_fit=Table([ant,pol,trec,opac],names=['Antenna','Pol','Trec','Opac'],dtype=('S2','S2','f4','f4'))
	apcal_fit['Antenna_no']=[get_antenna_number(uvdata,x) for x in apcal_fit['Antenna']]
	nn=len(uvdata.antennas)
	tau0=nn*[0.1]
	trecvr=2*nn*[100]
	for i in range(nn):
		ind=np.where(apcal_fit['Antenna_no']==i+1)
		if len(ind[0])!=0:
			tau0[i]=float(np.round(np.mean(apcal_fit[ind]['Opac']),2))
			if tau0[i]==0.0:
				tau0[i]=0.1
			trecvr[2*i]=float(np.round(apcal_fit[ind]['Trec'][0],2))
			trecvr[(2*i)+1]=float(np.round(apcal_fit[ind]['Trec'][1],2))
	return trecvr,tau0


def do_apcal(uvdata,inv=0,aparm=[0],tyv=1,dofit=[0],tau0=[0],trecvr=[0],opcode='',calin='',savelog=True,repeat=2):
	apcal_log=aips_out+uvdata.name+'_APCAL_fit.log'
	temp =apcal_log.split('.')
	HF.delete_file(temp[0]+'*'+temp[1])
	apcal_log=HF.filename_append(apcal_log)
	AT.apcal(uvdata, aparm=aparm,inv=inv,tyv=tyv,dofit=dofit,tau0=tau0,trecvr=trecvr,opcode=opcode,calin=calin,savelog=apcal_log)
	trecvr,tau0 = get_apcal_fit(uvdata,apcal_log)
	trecvr_str=[','.join([str(x) for x in trecvr])]
	tau0_str=[','.join([str(x) for x in tau0])]

	if repeat>1:
		for i in range(repeat-1):
			trecvr,tau0 = get_apcal_fit(uvdata,apcal_log)
			apcal_log		= HF.filename_append(apcal_log)
			uvdata.zap_table('SN',uvdata.table_highver('SN'))
			#dofit=15*[1]
			AT.apcal(uvdata, aparm=aparm,inv=inv,tyv=tyv,dofit=dofit,tau0=tau0,trecvr=trecvr,opcode=opcode,calin=calin,savelog=apcal_log)
			trecvr2_str=[','.join([str(x) for x in trecvr])]
			tau02_str=[','.join([str(x) for x in tau0])]

	logger.info('''
	APCAL did run %i times. First with assuming typical values of 
	trecvr=100 and tau0=0.1 for all telescopes.
	For each following run the fit values for from the run of
	APCAL before was used as input for tau0 and trecvr.
					
	First run fit:
	%s
	%s

	Last run fit:
	%s
	%s
						
	SN%i produced.
	''' %(repeat,trecvr_str,tau0_str,trecvr2_str,tau02_str,uvdata.table_highver('SN')))
