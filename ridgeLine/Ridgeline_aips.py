#!/usr/bin/env ParselTongue
'''
Created on Friday Mar 26 2021

Derive the ridge-line, including width and peak flux of the jet. It only works well for fairly straight jets.

This module uses ParselTongue, which requires a working AIPS.

@author: Anne-Kathrin Baczko
'''
import sys,os,json,logging,logging.config,datetime
import modules.functions as AF
import modules.hrk as hrk
import modules.aips_tasks as AT
from scipy.stats import sem
from astropy.table import Table
from glob import glob
from astropy.io import ascii,fits
import numpy as np
import statistics

from AIPS import AIPS#, AIPSDisk # I am not sure if I need AIPSDisk
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage,AIPSCat
from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData
from AIPSTV import AIPSTV
from modules.helper_functions import *
#from logToTable import *
#logger = logging.getLogger(__name__)
local_dir = os.getcwd()+'/'

def start_logging(filename):
	filename = filename.split('/')[-1].split('.')[0]
	now			= datetime.datetime.now()
	date		= str(now.day)+'_'+str(now.month)+'_'+str(now.year)
	tasklog = local_dir+'PT_inputs_'+filename+'_'+date+'.log'

	delete_file(tasklog)
	logfile = open(tasklog,'a')
	#
	log_r = open(tasklog,'r')
	if len(log_r.read())==0:
		log_r.close()
		logfile.write('########################################################\n')
		logfile.write('# This is a ParselTongue Log file to keep track of the # \n')
		logfile.write('# inputs to AIPS Tasks in ParselTongue. 								#\n')
		logfile.write('########################################################\n')
		logfile.write('#\n'*2)
	else:
		log_r.close()
	#
	logfile.close()
	with open('logging_config.json','r') as f:
	  conf = json.load(f)
	
	conf["handlers"]["info_file_handler"]["filename"] = tasklog
	logging.config.dictConfig(conf)
	logger = logging.getLogger('__name__')

	aipslog = local_dir+'AIPS_messages_'+filename+'_'+date+'.log'
	delete_file(aipslog)
	sys.stdout.write('Will write AIPS messages to file \n{0} \n and task inputs to file \n{1}\n'.format(aipslog,tasklog))
	# initialize aips loggger
	return [logger,aipslog]

def get_blc_trc(imfile,image_fraction=1.):
	'''
	returns a set of blc and trc
	'''
	with fits.open(imfile) as hdulist:
		header = hdulist[0].header
	naxis1= header['NAXIS1']
	naxis2= header['NAXIS2']
	dx = naxis1*image_fraction/2.
	dy = naxis2*image_fraction/2.
	x1 = int(0+(naxis1/2.-dx))
	x2 = int(naxis1-(naxis1/2.-dx))
	y1 = int(0+(naxis2/2.-dy))
	y2 = int(naxis2-(naxis2/2.-dy))
	print('blc=[{},{}] and trc=[{},{}]'.format(x1,y1,x2,y2))
	return [x1,y1],[x2,y2]

def makeslices(imfile,source=False,freq=False,outdisk=1,**kwargs):
	logger = logging.getLogger('__name__')
	'''
	imfile: location to the image
	source: source name to use for naming of catalog in aips, if not set the imfile name will be used
	freq: frequency string to use for naming of catalog in aips, if not set the imfile name will be used
	'''
	args = {'ang'		: 0,
					'blc'		: False,
					'trc'		: False,
					'nfit'	: 0,
					'gmax'	: 1,
					'gwidth': 9,
					'gpos'	: 0,
					'step1' : True, #rotate image 
					'step2' :	True, #make slices
					'step3' : True, #fit slices
					'step4' :	True, #plot slices with fit and write output data table
					'ang_step': False,
					'image_fraction': 1,
					'del_imcat': True, # delete image files after program is finished
					'plot_slice' : False}
	args.update(kwargs)
	if not args['blc']:
		blc,trc = get_blc_trc(imfile,args['image_fraction'])
	else:
		blc		= args['blc']
		trc		= args['trc']
	nfit	= args['nfit']
	gmax	= args['gmax']
	gwidth= args['gwidth']
	gpos	= args['gpos']

	fitdir = 'SLFIT_{}degree'.format(args['ang'])
	if source:
		if freq:
			imname = source+'_'+freq
		else:
			imname = source
	else:
		imname = imfile.split('/')[-1].split('.')[0]
	if len(imname)>12:
		imname = imname[:12]
		logger.info('Image name too long, will use \'{}\'\n'.format(imname))
	imdata = AIPSImage(*[imname,'FITS',outdisk,1])
	if imdata.exists() == True:
		imdata.zap()
	imcat = AT.imlod(imname,outdisk,imfile,'FITS')
	imdata = AIPSImage(*AT.getndata(outdisk,imcat))
	if imdata.exists() == True:
		logger.info('Successfully loaded new image file.\n')
		logger.info('Data Loaded: (%s, %s, %d, %d)\n',imdata.name,imdata.klass,imdata.disk,imdata.seq)
	else:
		logger.info('Image file not loeaded. Check input.\n')
		return
	if args['ang_step']:
		raw_input("Press Enter to continue...")
	rotname = imname[:10]+'_R'
####
	if args['step1']:
		rotdata = AIPSImage(*[rotname,'FITS',outdisk,1])
		if rotdata.exists()==True:
			rotdata.zap()
		rotcat = AT.lgeom(imdata,rotname,outdisk,outseq=1,aparm=[0,0,args['ang'],0,0,0,0,0,0])
		rotdata = AIPSImage(*AT.getndata(outdisk,rotcat))
		
		if rotdata.exists()==True:
			logger.info('Rotated image exists.\n')
		else:
			logger.error('Rotated image not loaded.\n')
			sys.exit()
	else:
		rotdata = AIPSImage(*[rotname,'FITS',outdisk,1])

	sys.stdout.write('data that will be used {}'.format(rotdata))

###########3####
	if args['step2']:
		AT.extd(rotdata,'SL',-1)
		j=1
		for i in range(blc[0],trc[0]):
			AT.slices(rotdata,blc=[i,blc[1],1],trc=[i,trc[1],1])
			if args['plot_slice']:
				AT.sl2pl(rotdata,inv=j)
			j+=1
	
		if args['plot_slice']:
			outfile='%s_slice.pdf'%rotdata.name
			AT.lwpla(rotdata, outfile)
#############
	if args['step3']:
		if os.path.isdir(fitdir):
			for f in glob(fitdir+'/*'):
				os.remove(f)
		else:
			os.makedirs(fitdir)
		if type(gmax) != list:
			if type(nfit) == int:
				nn = np.arange(nfit,rotdata.table_highver('SL')+1)
			else:
				nn = np.concatenate([np.arange(n[0],n[1]+1) for n in nfit])

			for ii in nn:
				gmax1		= gmax
				gwidth1 = gwidth
				saveLog = fitdir+'/SLFIT_{}_{}.log'.format(rotdata.name,ii)
				AT.slfit(rotdata,inv=int(ii),gmax=gmax,gwidth=gwidth,gpos=gpos,savelog=saveLog)
				lines = []
				with open(saveLog,'r') as f:
					for line in f:
						if line.find('LOGFILE') == -1:
							lines.append(line)

				# automatically get the fitting results from this slice to use as starting values for the next one
				jj = 0
				for line in lines:
					if line.find('Results in physical units') != -1:
						gmax		= float(lines[jj+2].split()[3])
						peakunit= lines[jj+2].split()[5]
						gwidth	= float(lines[jj+4].split()[3])
						gwidthEr= float(lines[jj+4].split()[4])
						if peakunit == 'MicroJy/BEAM':
							gmax*=1e-3 #convert everything to MilliJy/beam
					#	elif peakunit == 'MilliJy/BEAM':
					#		gmax *= 1
						break
					jj+=1
				if np.logical_or.reduce((gwidth<gwidthEr,gmax<gmax1/5.,gmax>5*gmax1)):
					gmax	= gmax1
					gwidth= gwidth1
		else:
			for gma,gwi,nf in zip(gmax,gwidth,nfit):
				for ii in range(nf[0],nf[1]+1):
					AT.slfit(rotdata,inv=ii,gmax=gma,gwidth=gwi,savelog=fitdir+'/SLFIT_{}_{}.log'.format(rotdata.name,ii))

###########		
	if args['step4']:
		AT.extd(rotdata,'PL',-1)
		for i in range(0,trc[0]-blc[0]):
			AT.sl2pl(rotdata,inv=i+1,domod=1,dores=1)
		outfile='{}_fit_{}degree'.format(rotdata.name,args['ang'])
		AT.lwpla(rotdata,outfile+'.pdf',docol=1,plcolors=[[0,0,0],[0,0,0],[0,0,1],[0,1,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[1,1,1]])
		
		AT.slcol(rotdata,outfile+'.txt',inv=1,nfiles=rotdata.table_highver('SL'))

#################
	if args['del_imcat']:
		sys.stdout.write('Deleting loaded map catalog and rotated catalog.\n')
		rotdata.zap()
		imdata.zap()
###################
##############################################
def calc_stddev(logFiles,cut=5):
	ridgelines = [Table.read(LOG,format='ascii',data_start=1,names=['isl','m','n','Xpos','Ypos','Dist','Peak','FWHM']) for LOG in logFiles]
	RLslices = [Table(dtype=ridgelines[0]) for i in range(len(ridgelines[0]['isl']))]
	i=0
	RL_std=[]
	sliceDir = 'RLslices'
	if not os.path.isdir(sliceDir):
		os.makedirs(sliceDir)

	sys.stdout.write('Will remove FWHM> {}*median(FWHM), fits with Peaks<0 from error calculation.\n'.format(cut))
	for RLs in RLslices:
		isl = ridgelines[0]['isl'][i]
		for RL in ridgelines:
			RLs.add_row(RL[RL['isl']==isl][0])
		median = np.median(RLs['FWHM'])
		medianF=np.median(RLs['Peak'])
		mask = np.logical_and.reduce((np.abs(RLs['FWHM'])<cut*median,np.abs(RLs['FWHM'])>0.1*median,RLs['Peak']>0,RLs['Peak']<5*medianF))
		RLs = RLs[mask]
		print(RLs)
		RL_std.append(np.std(np.abs(RLs['FWHM']),ddof=1)*1e3)
		ascii.write(RLs,sliceDir+'/RL_slice_{}.txt'.format(isl),overwrite=True)
		i+=1
	ascii.write([ridgelines[0]['isl'],RL_std], 'RL_std.txt',names=['isl','RL_stddev'],overwrite=True)
