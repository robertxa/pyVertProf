######!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
##########################################################
#                                                        #  
#            Age-Elevation profile script                #
#                                                        #  
#                 By Xavier Robert                       #
#                Grenoble, 16/03/02                      #
#                   Lima, 18/01/10                       #
#                                                        #  
##########################################################


 This code is to compute regressions for points with their error
      using different regression methods


"""
# Do divisions with Reals, not with integers
# Must be at the beginning of the file
from __future__ import division

# Import necessary Python modules
modulesNames = ['sys', 'math', 'warnings', 'copy']
for module in modulesNames:
	try:
		# because we want to import using a variable, do it this way
		module_obj = __import__(module)
		# create a global object containging our module
		globals()[module] = module_obj
	except ImportError:
		#sys.exit(u"ERROR : Module " + module + " not present. \n\n Please, install it \
		raise ImportError(u"ERROR : Module " + module + " not present. \n\n Please, install it \
			      \n\n Edit the source code for more information")
from shutil import copyfile
try:
	import numpy as np                               # need version 1.7 or higher
	from numpy.random import normal,randint
except ImportError:
	raise ImportError(u"ERROR : Module Numpy not present. \n\n Please, install it \
		      \n\n Edit the source code for more information")
try:
	import os
	from os import path, access, R_OK, mkdir         # W_OK for write permission.
except ImportError:
	raise ImportError("ERROR : Module os not present. \n\n Please, install it \
		      \n\n Edit the source code for more information")
try:
	import scipy as sp                               # need version 0.12 or higher
except ImportError:
	raise ImportError(u"ERROR : Module scipy not present. \n\n Please, install it \
		      \n\n Edit the source code for more information")
try:
	import matplotlib.pyplot as plt                  # module to plot figures
	from matplotlib import cm
	from matplotlib.patches import Polygon
except ImportError:
	raise ImportError(u"ERROR : Module matplotlib not present. \n\n Please, install it \
		      \n\n Edit the source code for more information")
try:
	from pylab import savefig                        # Module to manage figures
except ImportError:
	raise ImportError(u"ERROR : Module pylab not present. \n\n Please, install it \
		      \n\n Edit the source code for more information")
try:
	from kapteyn import kmpfit
except ImportError:
	raise ImportError(u"ERROR : Module kapteyn not present. \n\n Please, install it \
		      \n\n Edit the source code for more information")

from statsfuncs import *

############################################################################
############################################################################


def vertprofile(datafnme = u'', work_dir = u'',  header = 1, struct = [1,2,3,4],
	labelx = 'to be completed', labely = 'to be completed', rangex = None, rangey = None,
	statstypes = [0,1,2,3], confprob = 95.0,
	fontsz = 10, fontleg = 9,
	output = 'graph'):

	"""
	This code is to compute regressions for points with their error
      using different regression methods
	
	USAGE
		>>>from vertprof import vertprofile
		>>>vertprofile(datafnme = u'', work_dir = u'',  header = 1, struct = [1,2,3,4],
						labelx = 'to be completed', labely = 'to be completed', rangex = None, rangey = None,
						statstypes = [0,1,2,3], confprob = 95.0,
						fontsz = 10, fontleg = 9,
						output = 'graph')
		
	INPUTS
		To use options or inputs, you need to set them as
		vertprof(option_name = option_value, [...])
	
		1. work_dir: Working directory (local path)
				ex: work_dir = 'Purgatorio3/'
				Default None
		2. datafnme: name of text data file
	           should be in the form : Alt - Age - Age_Error
				ex: datafnme = 'Purgatorio3.txt'
		3. header: number of lines of the header in the data
				ex: header = 1 (Default)
		4. struct: Structure of the data
	         Struct = [Xdata, Xerr, Ydata, Yerr]
	         where Xdata, Xerr, Ydata, Yerr give their respective columns in the data file
	         If one of the data type do not exist, set the corresponding column to NONE
			ex : struct = [1,2,3,4] (Default)
		5. output: name of output
				ex: output = 'graph' (Default)
		6. labelx/labely: label x-axis/y-axis
				ex: labelx = 'Exposure age (ka)'
					labely ='Distance to the top of the scarp (m)'
		7. rangex/rangey: Set the x- and y-range
	               Depends on type of data
					ex: rangex = [0,8]
						rangey = [10,4] (in that case, the axis is inverted)
		8. statstypes: Type(s) of the stats to plot
					0 = kmpfit effective variance
					1 = kmpfit unweighted
					2 = Williamson
					3 = Cl relative weighting in X &/or Y
					ex: statstype = [0,1,2,3] (Default)
						statstype = [1,3]
		9. fontsz: labels fontsize
				ex: fontsz = 10 (Default)
		10. fontleg: legend fontsize
				ex: fontleg = 9 (Default)
		11. confprob: the confidence interval probabilty (in %)
				ex: confprob = 95.0 (Default)

	OUTPUTS
		The module outputs a pdf graph, and a text file with the stats outputs for each method asked
	
	CONTACT
		If needed, do not hesitate to contact the author. 
		Please, use https://isterre.fr/spip.php?page=contact&id_auteur=303

	LICENCE
		This package is licenced with `CCby-nc` (https://creativecommons.org/licenses/by-nc/2.0/)
	"""
	# If the minima input data are not given print the help file
	if not datafnme: 
		help(vertprofile)
		raise TypeError(color.RED + u'ERROR : ' + color.END + u'Name of Input file not given...')
	
	# check if the working directory exists, if not, create it
	if work_dir.strip():
		if work_dir[-1:] != u'/': work_dir = work_dir + u'/'
		if path.exists(work_dir) == False:
			os.mkdir(work_dir)
		else:
			warnings.warn(color.YELLOW + u"\nWARNING : " + color.END + u"Working directory already exists and will be erased\n")
		# Check if the file exists in the main folder
		#       if not, ckeck in the working folder
		#          if only in the working folder, copy it to the main folder
		if not path.isfile(datafnme):
			if not path.isfile(work_dir + datafnme):
				#sys.exit(color.RED + u"ERROR : " + color.END + u"Input file %s does not exist..." % datafnme)
				raise NameError(color.RED + u"ERROR : " + color.END + u"Input file %s does not exist..." % datafnme)
			else:
				copyfile(work_dir + datafnme, datafnme)	
		if not path.isfile(work_dir + datafnme):
			copyfile(datafnme, work_dir + datafnme)
		datapath = work_dir + datafnme
	else:
		if not path.isfile(work_dir + datafnme):
			#sys.exit(color.RED + u"ERROR : " + color.END + u"Input file %s does not exist..." % datafnme)
			raise NameError(color.RED + u"ERROR : " + color.END + u"Input file %s does not exist..." % datafnme)
		else:
			datapath = datafnme
			work_dir = u''
		
		
	# read data
	data = np.loadtxt(datapath, skiprows = header)
	y = data[:,struct[0]-1]
	if struct[1] != 0 and struct[1] != None:
		erry = data[:,struct[1]-1]
	else:
		erry = 0
	x = data[:,struct[2]-1]
	if struct[3] != 0 and struct[3] != None:
		errx = data[:,struct[3]-1]
	else:
		errx = 0
	N = len(x)
	
	# Open new file to print the results
	f1w = open(work_dir + u'results_' + datafnme, 'w')
	string = u"Results of weigthed least square from file " + datafnme + u"\n"
	print(color.RED + color.BOLD + string + color.END); f1w.write(string + u"\n")
	string = u"xmax/min = " + str(x.max()) + u", " + str(x.min())
	print(string); f1w.write(string + u"\n")
	string = u"ymax/min = " + str(y.max()) + u", " + str(y.min())
	print(string); f1w.write(string + u"\n")
	string = u"x/y mean = " + str(x.mean()) + u", " + str(y.mean())
	print(string); f1w.write(string + u"\n")
	
	string = (u"\nMethods used are: \n" +
	          u"   - kmpfit:  kapteyn method, from " +
	          u"https://www.astro.rug.nl/software/kapteyn/kmpfittutorial.html \n" +
	          u"              with error on X and Y or Y only or none \n" +
	          u"   - ODR: Orthogonal Distance Regression \n" +
	          u"   - Williamson: least square fitting with errors in X and Y " +
	          u"according to \n" +
	          u"                 Williamson (Canadian Journal of Physics, 46, 1845-1847, 1968) \n")
	print(color.CYAN + string + color.END); f1w.write(string)
	
	beta0 = [5.0, 1.0]         # Initial estimates
	if struct[3] != 0 and struct[3] != None:
		if struct[1] != 0 or struct[1] != None:
			# Prepare fit routine
			fitobj = kmpfit.Fitter(residuals=residuals, data=(x, y, errx, erry))
			fitobj.fit(params0=beta0)
			string = (u"\n\n======== Results kmpfit: weights for both coordinates =========")
			print(color.BLUE + color.BOLD + string + color.END)
			f1w.write(string + "\n")
			string = (u"Fitted parameters:      " + str(fitobj.params) + u"\n"
			         u"Covariance errors:      " + str(fitobj.xerror) + u"\n"
			         u"Standard errors         " + str(fitobj.stderr) + u"\n"
			         u"Chi^2 min:              " + str(fitobj.chi2_min) + u"\n"
			         u"Reduced Chi^2:          " + str(fitobj.rchi2_min) + u"\n"
			         u"Iterations:             " + str(fitobj.niter))
			print(string); f1w.write(string + u"\n")
		else:
			# Prepare fit routine
			fitobj2 = kmpfit.Fitter(residuals=residuals2, data=(x, y, erry))
			fitobj2.fit(params0=beta0)
			string = (u"\n\n======== Results kmpfit errors in Y only =========")
			print(color.BLUE + color.BOLD + string + color.END)
			f1w.write(string + u"\n")
			string = (u"Fitted parameters:      " + str(fitobj2.params) + u"\n"
			          u"Covariance errors:      " + str(fitobj2.xerror) + u"\n"
			          u"Standard errors         " + str(fitobj2.stderr) + u"\n"
			          u"Chi^2 min:              " + str(fitobj2.chi2_min) + u"\n"
			          u"Reduced Chi^2:          " + str(fitobj2.rchi2_min))
			print(string); f1w.write(string)
	
	# Unweighted (unit weighting)
	fitobj3 = kmpfit.Fitter(residuals=residuals2, data=(x, y, np.ones(N)))
	fitobj3.fit(params0=beta0)
	string = (u"\n\n======== Results kmpfit unit weighting =========")
	print(color.BLUE + color.BOLD + string + color.END)
	f1w.write(string + "\n")
	string = (u"Fitted parameters:      " + str(fitobj3.params) + u"\n"
	          u"Covariance errors:      " + str(fitobj3.xerror) + u"\n"
	          u"Standard errors         " + str(fitobj3.stderr) + u"\n"
	          u"Chi^2 min:              " + str(fitobj3.chi2_min) + u"\n"
	          u"Reduced Chi^2:          " + str(fitobj3.rchi2_min))
	print(string); f1w.write(string)

	from scipy.odr import Data, Model, ODR, RealData, odr_stop
	# Compare result with ODR
	linear = Model(model)
	if struct[3] != 0 and struct[3] != None:
		if struct[1] != 0 and struct[1] != None:
			mydata = RealData(x, y, sx=errx, sy=erry)
		else:
			mydata = RealData(x, y, sy=erry)
	else:
		mydata = RealData(x, y)
	
	myodr = ODR(mydata, linear, beta0=beta0, maxit=5000, sstol=1e-14)
	myoutput = myodr.run()
	string = (u"\n\n======== Results ODR =========")
	print(color.BLUE + color.BOLD + string + color.END)
	f1w.write(string + u"\n")
	string = (u"Fitted parameters:      " + str(myoutput.beta) + u"\n"
	          u"Covariance errors:      " + str(np.sqrt(myoutput.cov_beta.diagonal())) + u"\n"
	          u"Standard errors:        " + str(myoutput.sd_beta) + u"\n"
	          u"Minimum chi^2:          " + str(myoutput.sum_square) + u"\n"
	          u"Minimum (reduced)chi^2: " + str(myoutput.res_var))
	print(string); f1w.write(string)
	
	# Compute Williamson regression
	a_will, b_will, siga, sigb, beta, x_av, x__mean = williamson(x, y, errx, erry, struct, myoutput)
	
	if struct[3] != 0 and struct[3] != None:
		if struct[1] != 0 and struct[1] != None:
			chi2 = (residuals(fitobj.params,(x,y,errx,erry))**2).sum()
		else:
			chi2 = (residuals2(fitobj2.params,(x,y,erry))**2).sum()
	else:
		chi2 = (residuals2(fitobj3.params,(x,y,np.ones(N)))**2).sum()

	string = (u"\n\n======== Results Williamson =========")
	print(color.BLUE + color.BOLD + string + color.END)
	f1w.write(string + u"\n")
	string = (u"Average x weighted:     " + str(x_av) + u"\n"
	          u"Average x unweighted:   " + str(x__mean) + u"\n"
	          u"Fitted parameters:      " + u"[" + str(a_will) + u"," + str(b_will) + u"]\n"
	          u"Covariance errors:      " + u"[" + str(siga) + u"," + str(sigb) + u"]\n"
	          u"Minimum chi^2:          " + str(chi2))
	print(string); f1w.write(string)	
	
	string = (u"\n====================================="
	         u"\nPractical results:                     a       +       b     *   x")
	print(color.BLUE + color.BOLD + string + color.END)
	f1w.write(u"\n" + string + u"\n")
	string = (u"kmpfit unweighted:           %13.4f    %13.4f"
	         %(fitobj3.params[0], fitobj3.params[1]))
	print(string)
	f1w.write(string + u"\n")
	string = (u"   Exhumation Rate: %13.4f" %(1/fitobj3.params[1]))
	print(string)
	f1w.write(string + u"\n")
	string = (u"   Min. Exh. Rate: %13.4f" %(1/(fitobj3.params[1] - np.sqrt(fitobj3.stderr[1]))))
	print(string)
	f1w.write(string + u"\n")
	string = (u"   Max. Exh. Rate: %13.4f" %(1/(fitobj3.params[1] + np.sqrt(fitobj3.stderr[1]))))
	print(string)
	f1w.write(string + "\n")
	if struct[3] != 0 and struct[3] != None:
		if struct[1] != 0 and struct[1] != None:
			string = u"kmpfit effective variance:    %13.4f    %13.4f"%(fitobj.params[0], fitobj.params[1])
			print(string)
			f1w.write(string + u"\n")
			string = (u"   Exhumation Rate: %13.4f" %(1/fitobj.params[1]))
			print(string)
			f1w.write(string + u"\n")
			string = (u"   Min. Exh. Rate: %13.4f" %(1/(fitobj.params[1] - np.sqrt(fitobj.stderr[1]))))
			print(string)
			f1w.write(string + u"\n")
			string = (u"   Max. Exh. Rate: %13.4f" %(1/(fitobj.params[1] + np.sqrt(fitobj.stderr[1]))))
			print(string)
			f1w.write(string + u"\n")
		else:
			string = u"kmpfit weights in Y only:     %13.4f    %13.4f"%(fitobj2.params[0], fitobj2.params[1])
			print(string)
			f1w.write(string + u"\n")
			string = (u"   Exhumation Rate: %13.4f" %(1/fitobj2.params[1]))
			print(string)
			f1w.write(string + u"\n")
			string = (u"   Min. Exh. Rate: %13.4f" %(1/(fitobj2.params[1] - np.sqrt(fitobj2.stderr[1]))))
			print(string)
			f1w.write(string + u"\n")
			string = (u"   Max. Exh. Rate: %13.4f" %(1/(fitobj2.params[1] + np.sqrt(fitobj2.stderr[1]))))
			print(string)
			f1w.write(string + u"\n")
	string = u"ODR:                          %13.4f    %13.4f"%(beta[0], beta[1])
	print(string)
	f1w.write(string + u"\n")
	string = (u"   Exhumation Rate: %13.4f" %(1/beta[1]))
	print(string)
	f1w.write(string + u"\n")
	string = (u"   Min. Exh. Rate: %13.4f" %(1/(beta[1] - np.sqrt(myoutput.sd_beta[1]))))
	print(string)
	f1w.write(string + u"\n")
	string = (u"   Max. Exh. Rate: %13.4f" %(1/(beta[1] + np.sqrt(myoutput.sd_beta[1]))))
	print(string)
	f1w.write(string + u"\n")
	string = u"Williamson:                   %13.4f    %13.4f"%(a_will, b_will)
	print(string)
	f1w.write(string + u"\n")
	string = (u"   Exhumation Rate: %13.4f" %(1/b_will))
	print(string)
	f1w.write(string + u"\n")
	# For the next lines, not sure that we have to use the chi2 as error...
	string = (u"   Min. Exh. Rate: %13.4f" %(1/(b_will - chi2)))
	print(string)
	f1w.write(string + u"\n")
	string = (u"   Max. Exh. Rate: %13.4f" %(1/(b_will + chi2)))
	print(string)
	f1w.write(string + u"\n")
	print(u"\n")
	
	# calcul of confidence intervals
	dfdp = [1, x, x**2]
	#confprob = 95.0
	if struct[3] != 0 and struct[3] != None:
		if struct[1] != 0 and struct[1] != None:
			ydummy, upperband, lowerband = confidence_band(x, dfdp, confprob, fitobj, model)
			labelconf = u"CI (" + str(confprob) + u"%) relative weighting in X & Y"
			string = (u"y = a + b * x with a = %6.4f +/- %6.4f and b = %6.4f +/- %6.4f " 
			         %(fitobj.params[0], np.sqrt(fitobj.stderr[0]), 
			           fitobj.params[1], np.sqrt(fitobj.stderr[1])))
		else:
			ydummy, upperband, lowerband = confidence_band(x, dfdp, confprob, fitobj2, model)
			labelconf = u"CI (" + str(confprob) + u"%) relative weighting in Y"
			string = (u"y = a + b * x with a = %6.4f +/- %6.4f and b = %6.4f +/- %6.4f " 
			         %(fitobj2.params[0], np.sqrt(fitobj2.stderr[0]), 
			           fitobj2.params[1], np.sqrt(fitobj2.stderr[1])))
	else:
		ydummy, upperband, lowerband = confidence_band(x, dfdp, confprob, fitobj3, model)
		labelconf = u"CI (" + str(confprob) + u"%) with no weighting"
		string = (u"y = a + b * x with a = %6.4f +/- %6.4f and b = %6.4f +/- %6.4f " 
		         %(fitobj3.params[0], np.sqrt(fitobj3.stderr[0]), 
		           fitobj3.params[1], np.sqrt(fitobj3.stderr[1])))
	#verts = zip(x, lowerband) + zip(x[::-1], upperband[::-1])
	# Because we need to invert X and Y for the graph
	verts = zip(lowerband, x) + zip(upperband[::-1], x[::-1])
	
	# Close output file
	f1w.close()

	# Some plotting
	plt.rc(u'font', size = fontsz)
	plt.rc(u'legend', fontsize = fontleg)
	# Open figure
	fig = plt.figure(1)
	fig.clear()
	#d = (x.max() - x.min())/10
	d = (y.max() - y.min())/10
	X = np.linspace(x.min()-d, x.max()+d, 50)
	Y = np.linspace(y.min()-d, y.max()+d, 50)
	frame = fig.add_subplot(1,1,1, aspect=1, adjustable=u'datalim')
	# Plot the data and the least squared weigthed model 
	#      depending if we have errors on X/Y or only Y

	if struct[1] != 0 and struct[1] != None:
		if struct[3] != 0 and struct[3] != None:
			frame.errorbar(y, x, xerr=erry, yerr=errx,  fmt='o', label=u"Data")
			if 0 in statstypes:
				frame.plot(model(fitobj.params,X), X, 'g', ls='--', lw=2, label=u"kmpfit effective variance")
		else:
			frame.errorbar(y, x, xerr=erry,  fmt='o', label=u"Data")
			if 0 in statstypes:
				frame.plot(model(fitobj2.params,X), X, 'g', label=u"kmpfit errors in y only")
	else:
		frame.plot(y, x, fmt = 'o', label=u"Data")
	# Plot the classic regression
	if 1 in statstypes:
		frame.plot(model(fitobj3.params,X), X, 'r', label=u"kmpfit unweighted")
	# Plot the Williamson regression
	if 2 in statstypes:
		frame.plot(model((a_will,b_will),x), x, 'y', label=u"Williamson")
	if struct[3] != 0 and struct[3] != None and 3 in statstypes:
		# plot the confidence intervals
		poly = Polygon(verts, closed=True, fc='c', ec='c', alpha=0.3, #ls = 'dashed',
		           label="CI (" + str(confprob) + u"%) relative weighting in X & Y")
		frame.add_patch(poly)
	# Set graph characteristics/attributes
	frame.set_xlabel(labelx)
	frame.set_ylabel(labely)
	frame.set_title(string)
	frame.grid(True)
	#leg = frame.legend(loc=1)
	leg = frame.legend(loc="best")
	if rangex != None: plt.xlim(rangex[0],rangex[1])
	if rangey != None: plt.ylim(rangey[0],rangey[1])
	fig.savefig(work_dir + output + u'.pdf')
	plt.close()

   
############################################################################
############################################################################

# Main code
if __name__ == u'__main__':
	
	print 'zoup'
	# Run the code
	vertprofile()
	
# END of code
