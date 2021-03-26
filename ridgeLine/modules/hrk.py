#!/usr/bin/env ParselTongue
from AIPS import AIPS
from AIPSTask import AIPSTask
import os,warnings,glob,sys
from AIPSData import AIPSUVData, AIPSImage,AIPSCat
try:
    from StringIO import StringIO ## for Python 2
except ImportError:
    from io import StringIO ## for Python 3
from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData
import logging
#open logging
logger = logging.getLogger(__name__)
#############################################################################
### functions from Hans Kloeckner ###########################################
#############################################################################
def loadata(infile,outname,outclass='None',outdisk=1,outseq=0,wt=0,clint=0.1,digicor=-1,doconcat=-1):
	logger.info('Loading uvdata %s to catalog [%s,%s,%d,%d]\n',infile,outname,outclass,outdisk,outseq)

	"""
	Load a FITS file into the AIPS disk
    
	wt = Flagging threshold based on weights [range: 0 to 1 (flag all data)]
	"""
	#
	fitld = AIPSTask('FITLD')
	fitld.datain=str(infile)
	fitld.douvcomp  = -1
	fitld.digicor   = digicor
	fitld.doconcat	= doconcat
	fitld.wtthresh  = wt
	fitld.outname = outname
	fitld.clint = clint 
	if outclass != 'None':
		fitld.outclass  = outclass
		fitld.outdisk   = int(outdisk)
		fitld.outseq    = int(outseq)
		#fitld.msgkill   = -2
	fitld.go()
	return[fitld.outname,fitld.outclass,int(fitld.outdisk),int(fitld.outseq)]

def getheader(data,gethead='FREQ'):
	"""
	returns the value, referenc, increment, dimention of the 
	image of the header mostly used to get
	e.g. 'RA', 'DEC','STOKES', 'FREQ','IF'
	
	output would be for FREQ [61400000,0.5,125000.0,128]
	"""
	#
	if  datatype(data) == 'UV':
		from Wizardry.AIPSData import AIPSUVData as CHECKHEADER
	else:
		from Wizardry.AIPSData import AIPSImage as CHECKHEADER
	#
	chdata   = CHECKHEADER(data[0],data[1],int(data[2]),int(data[3]))
	chaxis   = chdata.header['ctype']
	chval    = -1
	chref    = -1
	chincr   = -1
	chdim    = -1
	
	for i in range(len(chaxis)):
		if chaxis[i] == gethead:
			if chaxis[i].startswith(gethead):
				chval  = chdata.header['crval'][i]
				chref  = chdata.header['crpix'][i]
				chincr = chdata.header['cdelt'][i]
				chdim  = chdata.header['naxis'][i]
	
	return[chval,chref,chincr,chdim]

def calpixval(data,axis='FREQ',pixel=0):
	"""
	Calculates the correct value of the pixel
	according to email exchange with Eric Geisen 
	freq = freq0 + (ch - refch) * increment
	"""
	headinfo = getheader(data,axis)
	
	if headinfo[0] != -1:
		if pixel == 0:
			value = headinfo[0] + (headinfo[3] - headinfo[1]) * headinfo[2]
		else:	headinfo = getheader(data,axis)

		value = headinfo[0] + (pixel - headinfo[1]) * headinfo[2]
		return(value)
	else:
		return(-1)

def datatype(data):
	"""
	return the datatype of the file
	"""
	#
	datatype = 'None'
	catalog_in = AIPSCat(int(data[2]))
	catalog    = str(catalog_in).split()
	#
	matchfname = []
	for i in range(len(catalog)):
		if catalog[i].find(data[0]) == 0:
			matchfname.append(i)
	matchfname.append(len(catalog))
	
	classindex = -1
	isname  = -1
	isclass = -1
	isdisk  = -1
	isseq   = -1
	i = 0
	while (isname != 1 or isclass != 1 or isdisk !=1 or isseq !=1 or classindex != -1) and i < len(matchfname)-1:
		isname  = -1
		isclass = -1
		isdisk  = -1
		isseq   = -1
		#print matchfname[i],matchfname[i+1]
		#print catalog[matchfname[i]:matchfname[i+1]]
		for j in range(matchfname[i],matchfname[i+1]):
		
			if catalog[j].replace('.','') == data[0].replace('.',''):
				isname = 1
			if catalog[j].replace('.','') == data[1]:
				isclass = 1
			if catalog[j].replace('.','') == str(int(data[2])):
				isdisk = 1
			if catalog[j].replace('.','') == str(int(data[3])):
				isseq = 1
			
			if catalog[j] == 'UV' and isname == 1 and isclass == 1 and isdisk == 1 and isseq ==1 and classindex == -1:
				classindex = j
			if catalog[j] == 'MA' and isname == 1 and isclass == 1 and isdisk == 1 and isseq ==1 and classindex == -1:
				classindex = j
		i += 1
	
	if classindex != -1:
		datatype = catalog[classindex]
	
	return(datatype)
 

def getntable(data,inext='CL'):
  """
  returns the number of speficied AIPS table 
  """
  if  datatype(data) == 'UV':
    from AIPSData import AIPSUVData as GETTAB
  else:
    from AIPSData import AIPSImage as GETTAB
  #
  tadata   = GETTAB(data[0],data[1],int(data[2]),int(data[3]))
  tables = tadata.tables
  ntable = -1
  for i in range(len(tables)):
    if tables[i][1] == 'AIPS '+inext:
      ntable = tables[i][0]
  return(int(ntable))
