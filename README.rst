pyVertProf
========

This code is to compute regressions for points with their error using different regression methods

It uses a DEM that is projected in UTM (best input), but will also work with geographic coordinates (lat-long).

This module may contain numerous bugs, and could probably be ameliorate and/or optimized. If you have any comments, do not hesitate to add a new branch or to contact the author.
To know the history, please check the file History.txt

Install
-------

To install it :

.. code-block:: bash

	pip install pyVertProf

To update it :

.. code-block:: bash

	pip install -U pyVertProf

If during the update you get a problem with the update of a dependency, you may try :

.. code-block:: bash

	pip install -U --no-deps pyVertProf


The module has been written and tested with Python 2.7, but not with Python 3.

Dependencies
------------
This code needs the following python modules and their dependencies, you may install them before the installation of the **pyswath** module:
	- kapteyn

Usage
-----

Inside a (i)python environnement:

To import the module:

.. code-block:: python

>>> from pyVertProf import vertprofile
	
To plot a swath profile [A,B] through the raster 'DEM/dem.tif':

.. code-block:: python

    >>> vertprofile(datafnme = u'', work_dir = u'',  header = 1, struct = [1,2,3,4], labelx = 'to be completed', labely = 'to be completed', rangex = None, rangey = None, statstypes = [0,1,2,3], confprob = 95.0, fontsz = 10, fontleg = 9, output = 'graph')

Options/inputs
--------------

To use options or inputs, you need to set them as	

.. code-block:: python

    >>> vertprof(option_name = option_value, [...])
	
Options/inputs are (option_names):
	1. work_dir: Working directory (local path)
	
				ex: ``work_dir = 'Purgatorio3/'``
	
				Default ``None``
	2. datafnme: name of text data file
	           should be in the form : Alt - Age - Age_Error
	
				ex: ``datafnme = 'Purgatorio3.txt'`
	
	3. header: number of lines of the header in the data
				
				ex: ``header = 1`` (Default)
				
	4. struct: Structure of the data
	         Struct = [Xdata, Xerr, Ydata, Yerr]
	         where Xdata, Xerr, Ydata, Yerr give their respective columns in the data file
	         If one of the data type do not exist, set the corresponding column to NONE
			ex : ``struct = [1,2,3,4] (Default)``
	
	5. output: name of output
			
			ex: ``output = 'graph'`` (Default)
	
	6. labelx/labely: label x-axis/y-axis
				
				ex: ``labelx = 'Exposure age (ka)'``
					``labely ='Distance to the top of the scarp (m)'``
	
	7. rangex/rangey: Set the x- and y-range
	               Depends on type of data
					
					ex: ``rangex = [0,8]``
						``rangey = [10,4]`` (in that case, the axis is inverted)
	
	8. statstypes: Type(s) of the stats to plot
					0 = kmpfit effective variance : `kapteyn method <https://www.astro.rug.nl/software/kapteyn/kmpfittutorial.html>'_ with error on X and Y or Y only or none
					1 = kmpfit unweighted : Orthogonal Distance Regression
					2 = Williamson : least square fitting with errors in X and Y according to Williamson (Canadian Journal of Physics, 46, 1845-1847, 1968)
					3 = Cl relative weighting in X &/or Y
					
					ex: ``statstype = [0,1,2,3]`` (Default)
						``statstype = [1,3]``
	
	9. fontsz: labels fontsize
				
				ex: ``fontsz = 10`` (Default)
	
	10. fontleg: legend fontsize
				
				ex: ``fontleg = 9`` (Default)
	
	11. confprob: the confidence interval probabilty (in %)
				
				ex: ``confprob = 95.0 (Default)``

Help files
----------

To get help in your (i)python environnement:

.. code-block:: python

	>>> help(vertprofile)

or simply:

.. code-block:: python

	>>> vertprofile()

Examples
--------

.. code-block:: python

    >>> vertprofile(datafnme = u'test.txt', work_dir = u'test',  header = 1, struct = [1,2,3,4], labelx = u'Ages (Ka)', labely = u'Depth (m)', rangex = [0,8], rangey = [10,4], statstypes = [0,1,2,3], confprob = 95.0)

The previous line permits to build the graph :

.. image:: https://github.com/robertxa/pyVertProf/tree/master/pyVertProf/graph.png?raw=true
   :scale: 100 %
   :align: center
			
Outputs
-------

The code build two files :
	- Graph.pdf: Graphical results of the computation
	- results_datafnme.txt: Output of the fitting methods

Contact
-------

If needed, do not hesitate to add a new branch or to contact the author. 
Please, use `https://isterre.fr/spip.php?page=contact&id_auteur=303 <https://isterre.fr/spip.php?page=contact&id_auteur=303>`_

Licence
-------

This package is licenced with `CCby-nc <https://creativecommons.org/licenses/by-nc-sa/3.0/>`_
