#!/usr/bin/env ParselTongue
#
'''
This module contains several functions to pass commands to AIPS via ParselTongue.
These are mainly running several aips tasks with a given set of input parameters.
'''
#
import AIPS
from AIPSTV import AIPSTV
from AIPSTask import AIPSTask
import os,warnings,glob,sys
from AIPSData import AIPSUVData, AIPSImage,AIPSCat
#from string import split,find
from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData
import time 
import __main__ as main
import helper_functions as HF
import logging
logger = logging.getLogger(__name__)

if os.path.isfile('observation_parameters.py'):
	import observation_parameters as OP
	#date = OP.date
	#local_dir = OP.local_dir
	aips_out = OP.aips_out+'/'
else:
	sys.stdout.write('''Not using a local observation_parameters.py file.
	I asume you are using aips_taks.py as module
	and are not planning to use it in scope of the quasi-pipeline.
	A local directory aips_out will be created to place outputs there.''')
	aips_out=os.getcwd()+'/aips_out/'
#
'''
Functinos which mimic procedures and verbs in aips
'''
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
	AIPSUV() will not give an error itself if the data tried to be loaded does not exist.
	'''
	if uvdata.exists()==True:
		logger.info('Data (%s,%s,%d,%d) exists. Assume it is already TB sorted (INDXR has been run).\n',uvdata.name,uvdata.klass,uvdata.disk,uvdata.seq)
	else:
		logger.error("Something went wrong Data does not exist")
		sys.exit()
	logger.info("Data Loaded: (%s, %s, %d, %d)\n", uvdata.name,uvdata.klass,uvdata.disk,uvdata.seq)
#

def getndata(diskno,catno):
	'''
	similar to getn in aips to get parameters for the specified
	disk and catalog.

	returns [outname,outclass,outdisk,outseq] to be used by AIPSUVData
	'''
	proxy   = AIPS.disks[diskno].proxy()
	disk    = proxy.AIPSCat.cat(AIPS.disks[diskno].disk, AIPS.userno)
	catalog = [disk[i] for i in range(len(disk)) if disk[i]['cno']==catno]
	print ([catalog[0]['name'],catalog[0]['klass'],diskno,catalog[0]['seq']])
	return [catalog[0]['name'],catalog[0]['klass'],diskno,catalog[0]['seq']]
#
def extd(uvdata,ine=False,inv=False):
	if ine == False:
		sys.stdout.write('Please specify a table type\n')
	if inv == False:
		sys.stdout.write('No table version specified, will delete highest vers table\n')
		logger.info('No table version specified, will delete highest vers table \n')
		table= uvdata.table_highver('AIPS CL')
	else:
		table = inv
	uvdata.zap_table(ine,table)
	logger.info('deleted Table %s%d  \n',ine,inv)
	sys.stdout.write('deleted Table {0} {1}\n'.format(ine,inv))
#
def get_header_info(uvdata):
	"""
	Returnes a nested dictonary that looks the same as table output from imheader in AIPS 
	with rows ['COMPLEX','STOKES','FREQ','IF','RA','DEC',''] 
	and colums ['naxis','crval','crpix','cdelt','crota'
	"""
	header = uvdata.header
	naxis = header['naxis'] #complex, stokes, freq, if, ra, dec
	ctype = header['ctype']
	crpix = header['crpix']
	crval	= header['crval']
	cdelt = header['cdelt']
	crota	= header['crota']
	d=dict()
	for n in range(len(ctype)):
		d[ctype[n]]={}
		d[ctype[n]]['naxis']=naxis[n]
		d[ctype[n]]['crval']=crval[n]
		d[ctype[n]]['crpix']=crpix[n]
		d[ctype[n]]['cdelt']=cdelt[n]
		d[ctype[n]]['crota']=crota[n]
	return d
#
def max_ch(uvdata):
	'''
	returns the maximum channel number
	'''
	# to return the channel number of the experiment
	header=get_header_info(uvdata)
	return header['FREQ']['naxis']
##
#

'''
Calling tasks in aips with one line
'''

def accor(uvdata,solint=-10):
	accor = AIPSTask('accor')
	accor.indata = uvdata
	accor.solint = solint #scan-average
	HF.print_inp_to_log(accor,'ACCOR')
	accor()
#  return inputs
#
def acscl(uvdata,gainu=0,bpv=1,doband=1):
	acscl = AIPSTask('acscl')
	acscl.indata = uvdata
	#acscl.timer[1:] = [0]
	acscl.solint = 0
	acscl.docalib = 1
	acscl.gainuse =gainu
	acscl.bpver=bpv
	acscl.doband=doband
	HF.print_inp_to_log(acscl,'ACSCL')
	acscl()
	#return inputs
#
def antab(uvdata,infile,tyv=1,gcv=1):
	antab = AIPSTask('antab')
	antab.indata = uvdata
	antab.calin = infile
	antab.tyver = tyv
	antab.gcver = gcv
	HF.print_inp_to_log(antab,'ANTAB')
	antab()
	#return  inputs
#
def apcal(uvdata,tyv=0,gcv=0,opcode='',aparm=[],dofit=[],inv=0,trecvr=[],calin='',tau0=[0],savelog=False,solint=0):
	apcal = AIPSTask('apcal')
	if type(savelog)==str:
		apcal.log = open(savelog,'w')
	apcal.indata = uvdata
	apcal.tyver = tyv
	apcal.gcver = gcv
	apcal.dotv = -1
	apcal.solint = solint
	#apcal.timer = []
	apcal.opcode = opcode
	apcal.tau0[1:] = tau0
	apcal.invers = inv
	apcal.calin = calin
	apcal.dofit[1:] = dofit
	apcal.trecvr[1:] = trecvr
	apcal.aparm[1:] = aparm
	HF.print_inp_to_log(apcal,'APCAL')
	apcal()
	if apcal.opcode != '':
		outfile = uvdata.name+'APCAL_opacity_correction.ps'
		lwpla(uvdata, outfile)
	#return inputs
#
def bpass (uvdata,bpv=-1,cals=[],antennas=[],gainu=0,solint=0,refant=0,bpassprm=[0],ichansel=[],timer=[],fgv=1):
	bpass = AIPSTask('bpass')
	bpass.indata = uvdata
	bpass.calsour[1:] = cals
	bpass.antennas[1:] = antennas
	bpass.gainuse = gainu
	bpass.bpver	= bpv
	bpass.flagver = fgv
	bpass.solint = solint
	bpass.refant = refant
	bpass.docalib = 1
	bpass.timer[1:] = timer
	bpass.cmethod = 'GRID'
	for i in range(len(ichansel)):
		n=i+1
		bpass.ichansel[n][1:] = ichansel[i]
# form of 'ichansel': [[0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]]
	bpass.bpassprm[1:] = bpassprm
	#for row in inputs:
	# print row
	#raw_input("Press Enter to continue...")
	HF.print_inp_to_log(bpass,'BPASS')
	bpass()
	#return inputs
#
def clcal(uvdata,timer=[0,0,0,0],sources=[],cals=[],antennas =[],opcode='CALP',interpol='AMBG',gainv=0,gainu=0,snv=0,refant=0,doblank=0,dobtween=1,smotype='',samptype='',suba = -32000,inv=0):
	clcal = AIPSTask('clcal')
	clcal.indata = uvdata
	clcal.sources[1:] = sources
	clcal.calsour[1:] = cals
	clcal.timer[1:] = timer
	clcal.opcode = opcode
	clcal.subarray = suba
	clcal.interpol = interpol
	clcal.antennas[1:] = antennas
	clcal.snver = snv
	clcal.inver = inv
	clcal.gainver = gainv
	clcal.gainuse = gainu
	clcal.doblank = doblank
	clcal.dobtween = dobtween
	clcal.smotype = smotype
	clcal.samptype = samptype
	clcal.refant = refant
	HF.print_inp_to_log(clcal,'CLCAL')
	clcal()
	#return inputs
#
def eops(uvdata):
	if(os.path.isfile('./usno_finals.erp')):
		eops_file='./usno_finals.erp'
	else:
		valid = {"yes": True, "y": True, "ye": True, "no": False, "n": False}
		# According to http://www.aips.nrao.edu/index.shtml the goddad webserver is down. On the webside a copy can be downoaded. The path to the copy is set below.
		#eops_file = HF.download_file('https://gemini.gsfc.nasa.gov/solve_save/usno_finals.erp','./usno_finals.erp')
		eops_file = HF.download_file('ftp://ftp.lbo.us/pub/staff/wbrisken/EOP/usno_finals.erp','./usno_finals.erp')
		if eops_file == False:
			msg = "Run EOP with presonal 'usno_finals.erp'?" 
			choice = None
			while True:
				sys.stdout.write('EOP file download failed. Either provide EOP file yourself or skip clcor(eop)\n')
				choice = raw_input('{0} [yes/no]\n'.format(msg)).lower()
				if choice in valid:
					if valid[choice] is True:
						isfile=False
						while isfile is False:
							msg2='Please provied path to file:'
							eops_file = raw_input('{0}\n'.format(msg2))
							if(os.path.isfile(filename)):
								isfile ==True
							else:
								sys.stdout.write('Path  not found. Try again.')
					elif valid[choice] is False:
						logger.info('Skip EOP and procede with calibration./n')
						return
				else:
					sys.stdout.write("Please say 'yes' or 'no'")
	clcor = AIPSTask('clcor')
	clcor.indata = uvdata
	clcor.gainv = uvdata.table_highver('CL')
	clcor.gainu = uvdata.table_highver('CL')+1
	clcor.opcode = 'EOPS'
	clcor.infile = './usno_finals.erp'
	HF.print_inp_to_log(clcor,'CLCOR EPOS')
	clcor()
	#return inputs
#
def fittp(imdata,outfile=''):
	uvdata=AIPSUVData(*imdata)
	if outfile=='':
		outfile=uvdata.header['object']
	while os.path.isfile(aips_out+outfile):
		outname=outfile.split('.')[0]
		filetype=outfile.split('.')[-1]
		outname=outname.split('_')
		try:
			outname[-1]= str(int(outname[-1])+1)
		except ValueError:
			outname.append('2')
		outfile='_'.join(outname)+'.'+filetype

	fittp = AIPSTask('fittp')
	fittp.indata = uvdata
	fittp.dataout = aips_out+outfile
	HF.print_inp_to_log(fittp,'FITTP')
	fittp()
#
def fring_global(uvdata,cals=[],fgv=0,suba=0,refant=0,antennas=[],timer=[0,0,0,1],solint=4,solsub=0,solmin=0,aparm=[1,0,0,0,1,2,4.5,0,1,0],dparm=[1,400,400,0],snv=0,gainu=0,docal=1,search=[],echan=0,bchan=1,bpv=-1,doband=-1,dofit=[0],get2n=False,flux=0):
	fring = AIPSTask('fring')
	fring.indata = uvdata
	fring.flagver = fgv
	fring.calsour[1:] = cals
	fring.timerang[1:] = timer
	fring.solint = solint
	fring.solsub = solsub
	fring.solmin = solmin
	fring.subarray = suba
	fring.antennas[1:] = antennas
	fring.dofit[1:] = dofit
	fring.refant = refant
	fring.gainuse = gainu
	fring.snver = snv
	fring.aparm[1:] = aparm
	fring.dparm[1:] = dparm
	fring.docalib = docal
	fring.doband = doband
	fring.bpver = bpv
	fring.search[1:] = search
	fring.bchan = bchan
	fring.echan = echan
	if get2n:
	#	in2data=getndata(get2n[0],get2n[1])
		fring.in2data=AIPSImage(*get2n)
		fring.flux=flux
	fring.inputs()
	HF.print_inp_to_log(fring,'FRING')
	fring()
	#return inputs
#
def fring_instr(uvdata,fgv=0,suba=0,cals=[],refant=0,antennas=[],timer =[0,0,0,0],solint=-1,aparm=[2,0,0,0,0,2,0,0,0],dparm=[1,400,400,0,0,0,0,1,1],snv=0,gainu=0,docal=1,search=[],echan=0,bchan=1,bpv=-1,doband=-1):
	fring = AIPSTask('fring')
	fring.indata = uvdata
	fring.calsour[1:] = cals
	fring.timerang[1:] = timer
	fring.flagver = fgv
	fring.solint = solint
	fring.subarray = suba
	fring.antennas[1:] = antennas
	fring.refant = refant
	fring.gainuse = gainu
	fring.snver = snv
	fring.aparm[1:] = aparm
	fring.dparm[1:] = dparm
	fring.docalib = docal
	fring.doband = doband
	fring.bpver = bpv
	fring.search[1:] = search
	fring.bchan = bchan
	fring.echan = echan
	fring.inputs()
	HF.print_inp_to_log(fring,'FRING')
	fring()
	#return inputs
#
def imhead(uvdata):
	header = uvdata.header
	tables = uvdata.tables
	outfile = aips_out+uvdata.name+'_header.txt'
	HF.delete_file(outfile)
	header_file = open(aips_out+uvdata.name+'_header.txt','a')
	header_file.write('Header for uvdata {0}\n'.format(uvdata.name))
	for key in header.keys():
		header_file.write('{0}: {1}\n'.format(key,header[key]))
	for row in tables:
		header_file.write('{0}\n'.format(row))
	header_file.close()
#
def imlod(outn,outd,datain,outc='CMAP'):
	imlod	= AIPSTask('imlod')
	imlod.log= open('imlod.log','w')
	imlod.outname	=outn
	imlod.outclass= outc
	imlod.outseq	= 0
	imlod.outdisk	= outd
	imlod.datain	= datain
	HF.print_inp_to_log(imlod,'IMLOD')
	imlod()
	
	lines=[]
	with open('imlod.log','r') as f:
		for line in f:
			if line.startswith('IMLOD'):
				lines.append(line)
		for l in lines:
			if l.find('CNO')!=-1:
				catnr = l.split()[-1]
	return int(catnr)
#
def indxr(uvdata,clint=0.1):
	indxr = AIPSTask('indxr')
	indxr.indata = uvdata
	indxr.cparm[3] = clint
	HF.print_inp_to_log(indxr,'INDXR')
	indxr()
	#return inputs
	#
def lgeom(imdata,outn,outd,outseq=0,aparm=[0,0,0,0,0,0,0,0,0],blc=[0],trc=[0]):
	lgeom = AIPSTask('lgeom')
	lgeom.log = open('lgeom.log','w')
	lgeom.indata	= imdata
	lgeom.outname = outn
	lgeom.outclass= lgeom.inclass
	lgeom.outseq	= outseq
	lgeom.outdis	= outd
	lgeom.blc[1:]	= blc
	lgeom.trc[1:] = trc
	lgeom.aparm[1:] = aparm
	HF.print_inp_to_log(lgeom, 'LGEOM')
	lgeom()

	lines=[]
	with open('lgeom.log','r') as f:
		for line in f:
			if line.startswith('LGEOM'):
				lines.append(line)
		for l in lines:
			if l.find('Create')!=-1:
				catnr = l.split()[-1]
	os.remove('lgeom.log')
	print(catnr)
	return int(catnr)
#
def listr(uvdata,bif=0,eif=0,optype='SCAN',ine='CL',inv=1,sources=[],dparm=[]):
	if len(dparm)>0:
		outprint = aips_out+uvdata.name+'_listr_'+optype+'dparm1='+str(dparm[0])+'.txt'
	else:
		outprint = aips_out+uvdata.name+'_listr_'+optype+'.txt'
	HF.delete_file(outprint)
	listr = AIPSTask('listr')
	listr.indata = uvdata
	listr.optype = optype
	listr.eif = eif
	listr.bif = bif
	listr.dparm[1:] = dparm
	listr.sources[1:]= sources
	listr.inext = ine
	listr.inver = inv
	listr.docrt = -1
	listr.outprint = outprint 
	listr()
#
def lwpla(uvdata, outfile, inv=None, plv=1,docol=0,plcolors=False):
	HF.delete_file(aips_out+outfile)
	#Must set indata, outfile
  #assert (indata != None, outfile != None)

	lwpla = AIPSTask('lwpla')
	lwpla.indata = uvdata
	lwpla.outfile = aips_out+outfile
	lwpla.lpen = 1.
	lwpla.dparm[6] = 4
	lwpla.dparm[8] = 9
	lwpla.docolor = docol
	if docol==1:
		i=1
		for plc in plcolors:
			lwpla.plcolors[i]=[None]+plc
			i+=1
	lwpla.plver = plv
	if (inv == None or inv == 0):
		inv = uvdata.table_highver('AIPS PL')
	lwpla.inver = inv
	HF.print_inp_to_log(lwpla,'LWPLA')
	try:
		lwpla()
	except:
		sys.stdout.write('LWPLA did encounter an error.\n PL files will remain.\n Change input parameters and try again.')
	else:
		sys.stdout.write('PL files written to plotfile {0}.\n PL files will be deleted now.\n'.format(outfile))
		logger.info('PL files written to plotfile %s.\n PL files will be deleted now.\n',outfile)
		extd(uvdata,ine='PL',inv=-1)
	
#  return inputs
#
def msort(uvdata):
	sys.stdout.write('Sorting data by runing MSORT.')
	msort = AIPSTask('msort')
	msort.outdata = uvdata
	msort.outclass = 'MSORT'
	msort.indata = uvdata
	HF.print_inp_to_log(msort,'MSORT')
	msort()
	return AIPSUVData(uvdata.name,'MSORT',uvdata.disk,uvdata.seq)
#
def pang(uvdata,suba=0):
	clcor = AIPSTask('clcor')
	clcor.indata = uvdata
	clcor.opcode = 'PANG'
	clcor.subarray = suba
	clcor.gainv = uvdata.table_highver('CL')
	clcor.gainu = uvdata.table_highver('CL')+1
	clcor.clcorprm[1] = 1
	HF.print_inp_to_log(clcor,'CLCOR PANG')
	clcor()
	#return inputs
#
def pccor(uvdata,refant=0,cals=[''],timer=[]):
	pccor = AIPSTask('pccor')
	pccor.indata = uvdata
	pccor.timer[1:] = timer
	pccor.inver = 1
	pccor.refant = refant
	pccor.calsour[1:] = cals
	HF.print_inp_to_log(pccor,'PCCOR')
	pccor()
#
def pclod(uvdata,calin):
	pclod = AIPSTask('pclod')
	pclod.indata = uvdata
	pclod.calin = calin
	pclod()
#
def possm(uvdata, aparm=[0,1,0,0,-200,200,0,0,1,0], solint=-1, sources =[], timer =[0,0,0,0],
      codetype='A&P', nplots=9, stokes='HALF', antennas=[], baseline =[],
      docal=1, gainu=0, doband=0, bpv=0, fgv=0, dotv=-1,bchan=1,echan=0,plotname=None,**kwds):
	
	uvdata.zap_table('PL',-1)	
	sys.stdout.write('Produce possm plots for antennas = {0}'.format(antennas))
	sys.stdout.write('and Sources = {0}'.format(sources))
	possm = AIPSTask('possm')
	possm.indata = uvdata
	possm.solint = solint
	possm.stokes = stokes
	possm.sources[1:] = sources
	possm.timerang[1:] = timer
	possm.antennas[1:] = antennas
	possm.baseline[1:] = baseline
	possm.docalib = docal
	possm.gainuse = gainu
	possm.flagver = fgv
	possm.doband = doband
	possm.bpver = bpv
	possm.bchan = bchan
	possm.echan = echan
	if(aparm==None):
		possm.aparm[9]=1
	else:
		possm.aparm[1:] = aparm[0:]
	possm.codetype = codetype
	possm.nplots = nplots
	possm.dotv = dotv
	HF.print_inp_to_log(possm,'POSSM')
	possm()
	if gainu==0:
		CLtable=uvdata.table_highver('CL')
	else:
		CLtable=gainu
	if sources==[]:
		sourcelist='All'
	elif len(sources) >=1:
		sourcelist='_'.join(set(sources))
		#sourcelist='AC'
	if fgv==-1:
		outfile='{0}_possm_CL{1}_{2}.ps'.format(uvdata.name,CLtable,sourcelist)
	elif fgv>-1:
		if fgv==0:
			fgv=uvdata.table_highver('FG')
		outfile='{0}_possm_CL{1}_FG{3}_{2}.ps'.format(uvdata.name,CLtable,sourcelist,fgv)
	if bpv>0:
		outfile=outfile.split('.')[0]+'BP'+str(bpv)+'.ps'
	if plotname:
		outfile='_'.join(outfile.split('_')[0:4])+'_'+plotname+'.ps'
	lwpla(uvdata, outfile)
#
def prtab(uvdata,ine,inv,return_of=False):
	outfile = aips_out+'{0}_{1}{2}.txt'.format(uvdata.name,ine,inv)
	prtab = AIPSTask('prtab')
	prtab.indata = uvdata
	prtab.inext = ine
	prtab.invers = inv
	prtab.docrt =-1
	prtab.outprint = outfile
#	prtab.box = [0]
	prtab()
	if return_of:
		return outfile
#
def prtan(uvdata):
	outprint = aips_out+uvdata.name+'_prtan.txt'
	HF.delete_file(outprint)
	prtan = AIPSTask('prtan')
	prtan.indata = uvdata
	prtan.docrt = -1
	prtan.outprint = outprint
	prtan()
#
def setjy(uvdata):
	setjy = AIPSTask('setjy')
	setjy.inputs()
	setjy.indata = uvdata
	setjy.zerosp[1:] = [1,0]
	HF.print_inp_to_log(setjy,'SETJY')
	setjy()
#
def sl2pl(imdata,inv=0,domod=0,dores=0):
	sl2pl = AIPSTask('sl2pl')
	sl2pl.indata	= imdata
	sl2pl.invers	= inv
	sl2pl.domodel	= domod
	sl2pl.doslice	= 1
	sl2pl.doresid	= dores
	sl2pl.dotv		= -1
	#HF.print_inp_to_log(sl2pl,'SL2PL')
	sl2pl()
#
def slcol(imdata,outfile,inv=1,nfiles=0,pixxy=[0,0],opc='MODL', aparm=[1,1,0]):
	slcol = AIPSTask('slcol')
	slcol.indata		= imdata
	slcol.outtext		= outfile
	slcol.nfiles		= nfiles
	slcol.invers		= inv
	slcol.zinc			= 0
	slcol.pixxy[1:]	= pixxy
	slcol.opcode		= opc
	slcol.aparm[1:]	= aparm
	HF.print_inp_to_log(slcol,'SLCOL')
	slcol()
#
def slfit(imdata,inv,gmax,gwidth,gpos=0,ngauss=1,savelog=False):
	slfit = AIPSTask('slfit')
	if type(savelog)==str:
		slfit.log = open(savelog,'w')
	slfit.indata	= imdata
	slfit.invers	= inv
	slfit.ngauss	= ngauss
	slfit.gmax[1:]= [gmax]
	slfit.gpos[1:]= [[None,gpos]]
	slfit.gwidth[1:]= [[None,gwidth]]
	slfit.dopos[1:]	= [[None,1]]
	slfit.dowidth[1:]= [[None,1]]
	slfit.domax[1:] = [1]
	HF.print_inp_to_log(slfit,'SLFIT')
	slfit()
#
def slices(imdata,blc,trc,outtext=''):
 sli = AIPSTask('slice')
 sli.indata		= imdata
 sli.blc[1:]	= blc
 sli.trc[1:]	= trc
 sli.outtext	= outtext
 #HF.print_inp_to_log(sli,'SLICE')
 sli()
#
def snplt(uvdata,ine='CL',antennas=[],nplots=16,inv=0,optype='PHAS',stokes='',sources=[],opcode='ALIF',plotname=False):
	outn = uvdata.name
	if plotname !=False:
		outfile=outn+'_'+plotname+'_snplt_'+ine+str(inv)+'.ps'
	else:
		outfile=outn+'_snplt_'+ine+str(inv)+'.ps'
	snplt = AIPSTask('snplt')
	snplt.indata = uvdata
	snplt.intable = uvdata.table(ine,inv)
	snplt.optype = optype
	snplt.antennas[1:] = antennas
	snplt.stokes = stokes
	snplt.sources = sources
	snplt.dotv = -1
	snplt.do3col = 2
	snplt.opcode = opcode
	snplt.nplots = nplots
	HF.print_inp_to_log(snplt,'SNPLT')
	snplt()
	lwpla(uvdata, outfile)
#
def snplt_check_ty(uvdata):
	tv=AIPSTV()
	tv.start()
	snplt=AIPSTask('snplt')
	snplt.indata=uvdata
	snplt.nplots = 16
	snplt.dotv =1
	for antenna in range(1,len(uvdata.antennas)+1):
		snplt.antennas[1:] = [antenna]
		snplt.inext ='CL'
		snplt.invers =1
		snplt.optype = 'PHAS'
		snplt.grchan = 1
		snplt.inputs()
		snplt()
		snplt.inext ='TY'
		snplt.invers =1
		snplt.optype = 'TSYS'
		snplt.grchan = 2
		snplt.inputs()
		snplt()
		raw_input("Press Enter to continue...")
	tv.kill()
#
def snsmo(uvdata,antennas=[],npiece=0,smotype='BOTH',samptype='MWF',bparm=[],cparm=[],doblank=-1,inv=0,outv=0,refant=0,dobtween=0):
	snsmo = AIPSTask('snsmo')
	snsmo.indata		= uvdata
	snsmo.samptype	= samptype
	snsmo.antennas[1:] = antennas
	snsmo.npiece		= npiece
	snsmo.smotype		= smotype
	snsmo.bparm[1:] = bparm
	snsmo.cparm[1:] = cparm
	snsmo.doblank		= doblank #check if 1 or -1 is needed
	snsmo.inver			= inv
	snsmo.outver		= outv
	snsmo.refant		= refant
	snsmo.dobtween	= dobtween
	HF.print_inp_to_log(snsmo,'SNSMO')
	snsmo()
#
def split(uvdata,sources=[],outd=1,aparm=[2,1,1,1,0],antennas=[],fgv=0,bpv=0,doband=-1,gainu=0,bchan=1,echan=0,bif=1,eif=0,stokes='HALF'):
	split = AIPSTask('split')
	split.indata = uvdata
	split.antennas[1:]=antennas
	split.sources[1:] = sources
	split.docalib = 2
	split.aparm[1:] = aparm
	split.outdisk = outd
	split.bpv   = bpv
	split.flagver = fgv
	split.doband= doband
	split.bchan = bchan
	split.echan = echan
	split.bif   = bif
	split.eif   = eif
	split.stokes= stokes
	split.gainuse = gainu
	HF.print_inp_to_log(split,'SPLIT')
	split()
	#return  inputs
#
def swpol(uvdata,antennas):
	swpol = AIPSTask('swpol')
	swpol.indata = uvdata
	swpol.antennas[1:] = antennas
	swpol.outname = uvdata.name
	swpol.outclass = 'SWPOL'
	swpol.outseq = uvdata.seq
	swpol.outdisk = uvdata.disk
	HF.print_inp_to_log(swpol,'SWPOL')
	swpol()
	return AIPSUVData(uvdata.name,'SWPOL',uvdata.disk,uvdata.seq)
#
def tabed(uvdata,ine='an',optype='repl',keyvalue=[5,0],aparm=[5,0,0,4,4,14,0],inv=1,outv=1):
	tabed = AIPSTask('tabed')
	tabed.indata = uvdata
	tabed.outdata = uvdata
	tabed.inext = ine
	tabed.invers = inv
	tabed.outvers = outv
	tabed.aparm[1:]=aparm
	tabed.keyvalue[1:]=keyvalue
	tabed.optype = optype
	HF.print_inp_to_log(tabed,'TABED')
	tabed()
#
def tacop(uvdata,ine='SN',inv=0,outv=0):
	tacop = AIPSTask('tacop')
	tacop.indata = uvdata
	tacop.inext = ine
	tacop.invers = inv
	tacop.outvers = outv
	tacop.outdata = uvdata
	HF.print_inp_to_log(tacop,'TACOP')
	tacop()
#
def usuba(uvdata,opcode='SCAN'):
	usuba = AIPSTask('usuba')
	usuba.indata = uvdata
	usuba.opcode = opcode
	HF.print_inp_to_log(usuba,'USUBA')
	usuba()
#
def uvflg_flagfile(uvdata,intext,outfgv =1,opcode= 'FLAG'):
	uvflg = AIPSTask('uvflg')
	uvflg.indata =uvdata
	uvflg.intext = intext
	uvflg.outfgver = outfgv
	uvflg.opcode= opcode
	HF.print_inp_to_log(uvflg,'UVFLG')
	uvflg()
	#return inputs
#
def uvcop(uvdata,timer=[0,0,0,0],sources=[],antennas=[],uvcopprm=[1,0,0,1,0,0,0]):
	uvcop = AIPSTask('uvcop')
	uvcop.indata = uvdata
	uvcop.outname = uvdata.name
	uvcop.outclass = 'UVCOP'
	uvcop.outseq = uvdata.seq
	uvcop.outdisk = uvdata.disk
	uvcop.timer[1:] = timer
	uvcop.antennas[1:] = antennas
	uvcop.sources[1:] = sources
	uvcop.uvcopprm[1:] = uvcopprm
	uvcop()
	return AIPSUVData(uvdata.name,'UVCOP',uvdata.disk,uvdata.seq)
#
def uvflg(uvdata,stokes='',timer=[0,0,0,0],antennas=[],bchan=1,outfgv=0,echan=0,bif=1,eif=0,reason='edged'):
	uvflg = AIPSTask('uvflg')
	uvflg.indata = uvdata
	uvflg.antennas[1:] = antennas
	uvflg.timer[1:] = timer
	uvflg.bchan =bchan
	uvflg.echan = echan
	uvflg.stokes = stokes
	uvflg.bif =bif
	uvflg.eif =eif
	uvflg.outfgver=outfgv
	uvflg.opcode = 'FLAG'
	uvflg.reason =reason
	HF.print_inp_to_log(uvflg,'UVFLG')
	uvflg()
	#return inputs
#
def vlog(uvdata,infile,outfile):
	files=['PCAL','TSYS','FLAG','WX']
	of=[outfile+'.'+x for x in files]
	[os.remove(x) for x in of if os.path.isfile(x)]
	vlog = AIPSTask('vlog')
	vlog.indata = uvdata
	vlog.calin = infile
	vlog.outfile = outfile
	vlog.fqtol = 1e3
	vlog.prtlev = 0
	vlog()
