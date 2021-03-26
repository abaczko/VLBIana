#!/usr/bin/env python
'''
Collections of functions to use with ParselTongue.
By Anne baczko@mpifr-bonn.mpg.de
'''
######################
# Python libraries ###
######################

try:
    from StringIO import StringIO ## for Python 2
except ImportError:
    from io import StringIO ## for Python 3

#from cStringIO import StringIO
try:
	from urllib2 import urlopen
except:
	from urllib.request import urlopen
import datetime as dt
import __main__ as main
import sys,os
import logging
#####################################
logger = logging.getLogger(__name__)



class Capturing(list):
	'''
	Class to get strings which had been written to the terminal
	'''
	def __enter__(self):
		self._stdout = sys.stdout
		sys.stdout = self._stringio = StringIO()
		return self
	def __exit__(self, *args):
		self.extend(self._stringio.getvalue().splitlines())
		del self._stringio    # free up some memory
		sys.stdout = self._stdout
#
def return_inputs(task):
	'''
	To get back the inputs from tasks 
	to be able to write them to log file
	'''
	with Capturing() as output:
		task.inputs()
	return output
#
def print_inp_to_log(task,taskname):
	'''
	To print down the task inputs into the log file
	'''
	logger.info('#\n')
	logger.info('Inputs for %s:\n',taskname)
	for row in return_inputs(task):
	  logger.info(row+'\n')
#
def download_file(url,file):
	'''
	For downloading a file, e.g. usno.finals for eop correction
	'''
	try:
		filedata = urlopen(str(url))
	except:
		logger.info('File %s could not be downloaded', url)
		return False
	datatowrite = filedata.read()
	with open(file, 'wb') as f:
		f.write(datatowrite)
#
def delete_file(filename):
	'''
	Test whether a file exists or not. If so, it will be deleted.
	This is for example used, to delete files before lwpla writes 
	to a file.
	'''
	if(os.path.isfile(filename)):
		logger.info('File %s already exists. Will be deleted\n',filename)
		os.remove(filename)
	else:
		logger.info('File %s will be created.\n',filename)
#
def qu_delete_file(filename):
	'''
	A second version of testing the existence of a file.
	This time the function asks if the file should be deleted or not.
	'''
	valid = {"yes": True, "y": True, "ye": True, "no": False, "n": False}
	if(os.path.isfile(filename)):
		msg = 'Delete it?'
		choice = None
		while True:
			sys.stdout.write('File {0} already exists.\n'.format(filename))
			choice=raw_input('{0} [y/n] \n'.format(msg)).lower()
			if choice in valid:
				if valid[choice] is True:
					os.remove(filename)
					sys.stdout.write('File {0} deleted.\n'.format(filename))
				elif valid[choice] is False:
					warnings.warn('File {0} already exists.  AIPS will probably append to it.'.format(filename))
					logger.error('File %s already exists.  AIPS will probably append to it.',filename)
				return valid[choice]
			else:
				sys.stdout.write("Please say 'yes' or 'no \n")
#
def day_of_date(date):
	date=dt.datetime.strptime(date,'%Y-%m-%d')
	new_year_day=dt.datetime(date.year,1,1)
	return (date-new_year_day).days + 1

def filename_append(outfile):
	while os.path.isfile(outfile):
		temp=outfile.split('.')
		outname=temp[0]
		filetype=temp[1]
		outname=outname.split('_')
		try:
			outname[-1]= str(int(outname[-1])+1)
		except ValueError:
			outname.append('2')
		outfile='_'.join(outname)+'.'+filetype
	return outfile


