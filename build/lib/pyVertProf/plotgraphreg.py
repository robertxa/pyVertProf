######!/usr/bin/env python
# -*- coding: utf-8 -*-

def plotgraphreg(frame, struct, x, y, errx, erry, X, Y, labeldata, statstypes, a_will, b_will, verts, confprob):
	
	"""
	Function to...
	
	INPUT
		frame: 
		struct: 
		x, y: 
		errx, erry: 
		X, Y: 
		statstypes: 
		a_will, b_will: 
		verts: 
		confprob:
		
	OUTPUT
		frame
	
	"""
	# Plot the data and the least squared weigthed model 
	#      depending if we have errors on X/Y or only Y
	if struct[1] != 0 and struct[1] != None:
		if struct[3] != 0 and struct[3] != None:
			frame.errorbar(y, x, xerr=erry, yerr=errx,  fmt='o', label=labeldata)
			if 0 in statstypes:
				frame.plot(model(fitobj.params,X), X, 'g', ls='--', lw=2, label=u"kmpfit effective variance")
		else:
			frame.errorbar(y, x, xerr=erry,  fmt='o', label=labeldata)
			if 0 in statstypes:
				frame.plot(model(fitobj2.params,X), X, 'g', label=u"kmpfit errors in y only")
	else:
		frame.plot(y, x, fmt = 'o', label=labeldata)
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
	
	return frame
	