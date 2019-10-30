from astropy.table import Table
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as optimization
from itertools import groupby
from scipy import stats
from scipy.odr import *
from astropy.io import fits

def linar_fit(p,x):
	a,b=p
	return a+b*x

def power_law(p,x):
	a,b=p
	return a*(x**b)

def odr_fit(func,data,x0,fit_type=0,verbose=False,maxit=1e4):
	model=Model(func)
	if len(data)==3:
		x,y,dy=data
	elif len(data)==4:
		x,y,dy,dx=data
	else:
		print "Syntax: odr_log(func,[x,y,dy[,dx]],x0...)"
	x0=x0
	if 'dx' in vars() or 'dx' in globacls():
		fitdata=RealData(x,y,sx=dx,sy=dy)
	else:
		fitdata=RealData(x,y,sy=dy)
	myodr=ODR(fitdata,model,beta0=x0,maxit=maxit)
	myodr.set_job(fit_type=fit_type)
	if verbose == 2: 
		myodr.set_iprint(final=2)
	out=myodr.run()
	out.pprint()
	if out.stopreason[0] == 'Iteration limit reached':
		print '(WWW) poly_lsq: Iteration limit reached, result not reliable!'
	return out.beta[0],out.beta[1],out.sd_beta[0],out.sd_beta[1]

def read_fits(file,col1,col2,dcol1,dcol2):
	hdulist=fits.open(file)
	data=hdulist[1].data
	x=data[col1]
	y=data[col2]
	dx=data[dcol1]
	dy=data[dcol2]
	hdulist.close()
	return x,y,dx,dy

def log_error(x,dx):
	return 1./np.log(10)*dx/x

def keyfunc(s):
	  return [int(''.join(g)) if k else ''.join(g) for k, g in groupby(s, str.isdigit)]


