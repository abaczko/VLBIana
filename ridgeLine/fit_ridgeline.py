#!/usr/bin/env ParselTongue
#
'''
Created on Friday Mar 26 2021

This script uses aips through ParselTongue to read in an image fits file, rotate if if required to align the jet along the x-axis, makes slices for each pixel, and fit a single Gaussian to the slices. 

By providing an angle range (+-ang_range) and a step size (ang_step) the slicing and fitting is repeated for the final set of angles [ang-ang_range,ang+ang-range] and based on the fitted FWHM the standard deviation for each slie is derived and written to the file 'RL_std.txt'.

Further outputs:
	aips_out: Slice-fit-plots for each rotation angle
	'source_name'_'freq'_R_fit_xxxdegree.txt : File with the ridgeline fitting values if no source_name and freq is provided, the first 10-12 letters of the input file name are used
	SLFIT_'source_name'_'freq'_xxdegree.dat: AIPS output for slice fitting for each rotation
	RL_std.txt: standard deviation of FWHM
	PT_inputs_xxx.log : Logging file containing all inputs that were given to aips in addition to som einformative messages
	AIPS_messages_xxx.log: The aips message file

To define the box for slicing you can either define blc and trc by hand or let a function evaluates the trc and blc automatically assuming the source does cover the whole image or a part (then provide variable image_fraction). This function is in testing mode. Please check results carefully.
E.g.	1.) image size is 1024: blc=[0,0],trc=[1024,1024]
			2.) image_fraction=0.5 for images size 1024x512: blc=[226,128], trc=[768,384]

'nfit': set to 0 if you want to fit for each slice. Set to an array of arrays if only a subset of slices should be fit. It may be that after a first fitting some slices cannot be fitted.
E.g. slices 3-10 are not fit well. But everything until the last slice 150 is good: nfit=[[0,2],[11,15]]

You can either leave the other parameters (gwidth,gmax,gpos) to one value(float or int) or give an array here with the same number of arrays as in nfit. But if you do so, each parameter has to be an array or a single float, please don't mix by giving e.g. gwidth as array and gpos as int:
E.g.	gwidth=[1,0.5]
			gpos = [0,-1]
			gmax = [1,2]
@author: Anne-Kathrin Baczko
'''
import sys,os,json,logging,logging.config,datetime
#
import numpy as np
from Ridgeline_aips import *
import aips_tasks as AT
from helper_functions import *
from glob import glob
from logToTable import *
import shutil
#
##################################################################
# Setup variables for observation You only have to put in something in this section according to you image file.#
#######################
local_dir = os.getcwd()+'/'
aipsid		= xx #Put in your aips id
outdisk		= 1 #set appropriate if other aips disk should be used for writing temporary files.
imfile = local_dir+'xxx.fits' #Fits image input file 
ang = 24	# rotation angle in mathematically positive sense, so this is from top towards counter-clock-wise
ang_range = 1 
ang_step = 1
angles = np.arange(ang-ang_range,ang+ang_range+1,ang_step)
source_name =False #source_name and frequency give you the option to customize the file names. If you let it blank file names will be derived by using the first 10-12 letters of the input file name. 
frequency = False #if not False give a string e.g. 'K'
blc	= False
trc	= False
image_fraction = 1
nfit	= 0
gwidth= 1
gpos	= 0
gmax	= 1

# Setup logging ##
##################
logger,aipslog = start_logging(imfile)
	
################################
###############################################
# Starting up AIPS ############################
###############################################
#
from AIPS import AIPS#, AIPSDisk # I am not sure if I need AIPSDisk
from AIPSTask import AIPSTask #, AIPSList
from AIPSData import AIPSUVData, AIPSImage,AIPSCat
from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData
from AIPSTV import AIPSTV

logger.info('Will use local AIPS installation')
aips_out= local_dir+'aips_out/' #All outputs from AIPS (Plots, text-files) will be saved there
if not os.path.isdir(aips_out):
  os.makedirs(aips_out)

AIPS.log  = open(aipslog,'a')
#
AIPS.userno	= aipsid
sys.stdout.write('Printing AIPSCat.\n')
print AIPSCat(outdisk) #list the catalogs on the desired disk
###############################################
# SLICE analysis start ###########################
###############################################

for ff in glob('*.txt'):
	os.remove(ff)
for ff in glob('*.dat'):
	os.remove(ff)

jj=0
for ii in angles:
	if jj ==0:
		angstep = True
	else:
		angstep = False
	sys.stdout.write('Make SLFITs for rotation angle: {} degree.\n'.format(ii))
	makeslices(imfile,source=source_name,freq=frequency,blc=blc,trc=trc,nfit=nfit,gpos=gpos,gmax=gmax,gwidth=gwidth,ang=int(ii),ang_step=angstep,image_fraction=image_fraction)
	jj+=1

AIPS.log.close()
#########################
# Now the error is calculated
############################

logger.info('Will calculate std of width now.\n')
logFiles= glob(local_dir+'*fit*.txt')
calc_stddev(logFiles)

rotdir = glob('SLFIT*degree')
for rd in rotdir:
	logFiles = glob(rd+'/*.log')
	SLfits = log_to_table(logFiles)
	shutil.rmtree(rd,ignore_errors=True)
