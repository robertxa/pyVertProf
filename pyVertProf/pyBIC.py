#!/usr/bin/env python3.9
# -*- coding: utf-8 -*-
"""
This Python script calculates the best fitting AER using least-square fitting taking into account 
measurement errors. The number of slopes is evaluated using the Bayesian Information Criteria (BIC). 
Original Matlab Code written by C. Glotzbach in 2010. Please cite the EPSL-Paper: Glotzbach et al. (2011).
Please, if you use this module, cite also: Robert X., pyAER, a python BIC Implementation for AER, 2022; Zenodo : 

# Copyright (c) 2022 Xavier Robert <xavier.robert@ird.fr>
# SPDX-License-Identifier: GPL-3.0-or-later


INPUTS:
    

OUTPUT:
    

USAGE:
    

INSTALL:
    In the top folder, run in your terminal:
    ~$ python setup.py install
    
    or use pip:
    pip install pyBIC
"""

import time, os, warnings
from copy import copy
import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations
from alive_progress import alive_bar              # https://github.com/rsalmei/alive-progress

from .statsfuncs import *

####################################################
####################################################

def readData(filename_data, systems = ['AHE', 'AFT', 'ZHE', 'ZFT'],
             names = ['LAT','LONG','HEIGHT', 'AHE', 'DAHE', 'AFT', 'DAFT',
                      'ZHE', 'DZHE', 'ZFT', 'DZFT', 'KAR', 'DKAR']):
    """_summary_

    Args:
        filename_data (string)  : Filename where data are stored.
        systems (list, optional): Systems to analyze and plot. 
                                  Defaults to ['AHE', 'AFT', 'ZHE', 'ZFT'].
        names (list, optional)  : Header of the data file.
                                  Defaults to ['LAT','LONG','HEIGHT', 'AHE', 'DAHE', 'AFT', 'DAFT', 
                                               'ZHE', 'DZHE', 'ZFT', 'DZFT', 'KAR', 'DKAR'].

    Raises:
        NameError: Input file does not exists.

    Returns:
        fdata (Pandas dataframe): Dataframe of data.
        nobs (dictionnary)      : Dictionnary of the number of data for each system.
    """

    # Check if the input file exists
    if not os.path.isfile(filename_data):
        raise NameError(color.RED + u"ERROR : " + color.END + u"Input file %s does not exist..." %(filename_data))
    # Read data file
    if filename_data[-4:] == '.txt':
        # In the case of the old Pecube (v2.x) structure 
        #   (where the first line correspond to the number of data)
        # Add the corresponding headers - agescols; Change it to pd dataframe?
        fdata = pd.read_csv(filename_data, 
                            #sep = "\t|;|:|,| |",
                            sep = None,
                            engine='python',
                            skiprows = 1, 
                            names = names)
        # Add a column with the name of samples
        junk = np.arange(len(fdata))
        samp = list(np.zeros(len(fdata)))
        for i in range (len(fdata)):
            samp[i] = 'sample_' + str(junk[i])
        fdata.insert(loc=0, column='SAMPLE', value=samp, allow_duplicates=False)

    elif filename_data[-4:] == '.csv':
        # If input data are csv format with a header,
        # or as Pecube v4.x input file (.csv)
        # with a ',' as separator
        fdata = pd.read_csv(filename_data, sep = ',')
    
    elif filename_data[-4:] == '.xls' or filename_data[-4:] == 'xlsx':
        # If data are in an Excel sheet (.xls or .xlsx)
        fdata = pd.read_excel(filename_data)
    else:
        # if not any case from above, the format is not taken in account,
        # Raise an error !
        raise NameError(color.RED + u"ERROR : " + color.END + u'input file %s not txt, csv or xlsx format...' %(filename_data))
    # Find sample names
    sample = list(fdata['SAMPLE'])
    # Find the obsservation type present in the data file
    obstype = fdata.columns.tolist()
    # Sort data with increasing elevations
    fdata = fdata.sort_values("HEIGHT", axis = 0, ascending = True)
    # Add an Index to the DataFrame, using the sample names
    fdata = fdata.set_index("SAMPLE", drop = False)
    # change -9999 and 0 values to np.nan
    #fdata = fdata.replace([-9999, 0], np.nan)
    fdata = fdata.replace(-9999, np.nan)

    # Find the number of lines of the data
    #nd = len(fdata)
    # Find the number of AHe, AFT, ZHe, ZFT
    # Check the header to find wich system are in the data
    #if 'AHE' in obstype:
    #    a = len(fdata[~np.isnan(fdata["AHE"])]) # number of AHe ages
    #else:
    #    a = None
    #if 'AFT' in obstype:
    #    b = len(fdata[~np.isnan(fdata["AFT"])]) # number of AFT ages
    #else:
    #    b = None
    #if 'ZFT' in obstype:
    #    c = len(fdata[~np.isnan(fdata["ZFT"])]) # number of ZFT ages
    #else:
    #    c = None
    #if 'ZHE' in obstype:
    #    d = len(fdata[~np.isnan(fdata["ZHE"])]) # number of ZHe ages
    #else:
    #    d = None
    
    nobs = {}
    for item in obstype:
        nobs[item] = len(fdata[~np.isnan(fdata[item])])

    # fdata is the table,
    # nobs is a dictionnary of the number of data
    #nobs = {'AHE' : a, 'AFT' : b, 'ZHE' : d, 'ZFT' : c}
    return fdata, nobs


#########################################################################################
def calculateBIC(filename_data, fdata, nobs, systems = ['AHE', 'AFT', 'ZHE', 'ZFT'], nslopes = None):
    """
    Function to compute the BIC criterion.    

    Args:
        filename_data (string)    : Input data file.
        fdata (Pandas dataframe)  : Dataframe of data.
        nobs (dictionnary)        : Dictionnary of the number of data for each system.
        systems (list, optional)  : System(s) to analyse/plot.
                                    Defaults: ['AHE', 'AFT', 'ZHE', 'ZFT'].
        nslope (integer, optional): number of slopes to test; 
                                    If not given, the default value ('None') imposes a search of the
                                    number of slopes
    """

    if not os.path.isdir('Graphs/'):
        # Do a directory where to store the graphs    
        os.mkdir('Graphs/')
    else:
        print(color.YELLOW + '\tWARNING: ' + color.END + '"Graphs/" folder already exists, it will be eraised')

    if not nslopes:
        # Find the maximum number of slopes (=nslopes)
        nslopes = 0
        maxN = 0
        for system in systems:
            if maxN < len(fdata[~np.isnan(fdata[system])][system]):
                maxN = len(fdata[~np.isnan(fdata[system])][system])
        nslopes = len(np.arange(0, maxN, 3)) - 1
    print("\tTotal number of slopes to process: %s for %s systems" %(str(nslopes), str(len(systems))))

    with alive_bar(len(systems) * nslopes, title = "\x1b[32;1m- Processing BIC\x1b[0m", length = 30) as bar:
        # Do a loop on each system to plot and on each number of slope
        # Open a text file to write the results
        f1w = open('Graphs/resultsBIC.txt', 'w')
        f1w.write('BIC computation results Computed with the Glotzbach et al., 2011 method upgraded by Xavier Robert (2022)\n')
        f1w.write('\t\tUsing file: %s \n\n' %(filename_data))
        f1w.write('I am analyzing %s systems: %s\n\n' %(len(systems), systems))

        # Open a new figure with suplots
        fig, ax = plt.subplots(nrows = len(systems), ncols = nslopes, 
                               sharex = False, sharey = True, 
                               figsize = (14/3 * nslopes, len(systems) * 6))
                               #figsize = (14, len(systems) * 6))

        # Do a loop on the system
        for system in systems:
            f1w.write('%s system:\n' %(system))
            if nobs[system] > 3:
                # take the rows that correspond to the system only, remove the rows with nan in the system
                data = fdata[~np.isnan(fdata[system])]
                # Do a test to know if the number of required slopes is lower 
                #   than the maximum number of slopes possible
                nslopes_max = len(np.arange(0, len(data[system]), 3)) - 1
                if nslopes < nslopes_max:
                    rangeslopes = nslopes + 1
                else:
                    if nslopes != nslopes_max:
                        print(color.YELLOW + '\tWARNING: ' + color.END + 'Total number of slopes possible (%s) for the system %s is lower than the number of slopes (%s) asked...' %(nslopes_max, system, nslopes))
                        print('\t\tI am computing only %s slopes' %(nslopes_max))
                    rangeslopes = len(np.arange(0, len(data[system]), 3))
                
                for nslope in range (1, rangeslopes):
                    # Do a loop on the number of slopes, to compute every possibility    
                    # Initiate parameters
                    BIC_best = 1E10
                    rw_best = np.zeros(nslope)
                    exh_mean_best = np.zeros(nslope)
                    intercept_mean_best = np.zeros(nslope)
                    exh_max_best = np.zeros(nslope)
                    intercept_max_best = np.zeros(nslope)
                    exh_min_best = np.zeros(nslope)
                    intercept_min_best = np.zeros(nslope)
                    knick_point_best = np.zeros(nslope + 1, int)
                    knick_point = np.zeros(nslope + 1, int)
                    
                    if nslope > 1:
                        # Find knickpoint combinaisons
                        combi_tot = find_combinaisons(data, system, nslope)
                    else:
                        combi_tot = [0]
                    exh_mean = np.zeros(nslope)
                    intercept_mean = np.zeros(nslope)
                    exh_max = np.zeros(nslope)
                    exh_min = np.zeros(nslope)
                    intercept_max = np.zeros(nslope)
                    exh_min_best = np.zeros(nslope)
                    intercept_min = np.zeros(nslope)

                    for k in range(len(combi_tot)):
                        # Compute the different slopes age-limits
                        # intitiate parameters
                        log_L = 0
                        rw = np.zeros(nslope)
                        mod = np.zeros(len(data[system]))
                        # Build the knickpoint vector
                        knick_point[0] = int(0)
                        knick_point[nslope] = int(nobs[system] - 1)
                        if nslope > 1 :
                            if nslope == 2 :
                                knick_point[1] = combi_tot[k][0]
                            else:
                                for i in range (1, nslope):
                                    knick_point[i] = combi_tot[k][i-1]
                        
                        for i in range (0, nslope):
                            # Compute the stats
                            mw, bw, rw[i], smw, sbw, xw, yw = lsqfityw(data['HEIGHT'][knick_point[i]:knick_point[i+1]+1],
                                                                        data[system][knick_point[i]:knick_point[i+1]+1], 
                                                                        data['D' + system][knick_point[i]:knick_point[i+1]+1])
                            exh_mean[i] = 1 / mw
                            intercept_mean[i] = -bw / mw
                            mod[knick_point[i]:knick_point[i+1]+1] = (data['HEIGHT'][knick_point[i]:knick_point[i+1]+1] - intercept_mean[i]) / exh_mean[i]
                
                            log_L = compute_log_L([knick_point[i], knick_point[i+1]+1], data, system, mod, log_L)
            
                            exh_max[i] = 1 / (mw - smw)
                            exh_min[i] = 1 / (mw + smw)
                            intercept_max[i] = -(bw + sbw)/(mw - smw)
                            intercept_min[i] = -(bw - sbw)/(mw + smw)
                        
                        # Compute the BIC
                        BIC = 2 * nslope * np.log(nobs[system] - 1) - 2 * log_L

                        if BIC < BIC_best:
                            BIC_best = BIC.copy()
                            log_L1 = log_L.copy()
                            exh_mean_best = exh_mean.copy()
                            intercept_mean_best = intercept_mean.copy()
                            exh_max_best = exh_max.copy()
                            intercept_max_best = intercept_max.copy()
                            exh_min_best = exh_min.copy()
                            intercept_min_best = intercept_min.copy()
                            knick_point_best = knick_point.copy()
                            rw_best = rw.copy()

                    # Plot the sub-graph
                    fig, ax = plotBIC(systems, system, data, nobs,
                                      BIC_best, log_L1, exh_mean_best, exh_min_best, exh_max_best, 
                                      intercept_mean_best, intercept_min_best, intercept_max_best,
                                      knick_point_best, rw_best,
                                      subplot = nslope, fig = fig, ax = ax, f1w = f1w)
                    # Update the progress-bar
                    bar()
            else:
                print(color.YELLOW + '\tWARNING: ' + color.END + 'Not enougth data to do BIC computation with the %s system...' %(system))
                print('\t\tI am skipping it...')
                f1w.write('\tNot enougth data to do BIC computation...\n')  
                # Plot the sub-graph
                for n in range(1, nslopes + 1):
                    fig, ax = plotBIC(systems, system, data, nobs,
                                      BIC = None, log_L1 = None, 
                                      exh_mean = None, exh_min = None, exh_max = None, 
                                      intercept_mean = None, intercept_min = None, intercept_max = None,
                                      knick_point = None, rw = None,
                                      subplot = n, fig = fig, ax = ax, f1w = f1w)
                    # Update the progress-bar
                    bar()  
            f1w.write('----------------------------------------------------------------------\n')
            
    # Save the fig in a pdf file in a Graph directory
    print("Saving plot...\n")
    plt.savefig('Graphs/BIC.pdf')
    # Close the fig
    plt.close()
    # Close the text file
    f1w.close()
    print('**********************************************************************')
        
    return


#########################################################################################
def plotBIC(systems, system, data, nobs,
            BIC, log_L1, exh_mean, exh_min, exh_max, 
            intercept_mean, intercept_min, intercept_max,
            knick_point, rw,
            subplot, fig, ax, f1w = None):
    """_summary_

    Args:
        system (str)           : System used (AFT, AHE, ZFT, ZHE, ...)
        data (pandas dataframe): Data to plot
        nobs (dictionnary)     : number of samples
        BIC (float)            : BIC value
        log_L1 (float)         : Log(L)
        exh_mean (array)       : Mean exhumation rate
        exh_min (array)        : Min exhumation rate
        exh_max (array)        : Max exhumation rate
        knick_point (array)    : 4 values array that gives the kinks when several slopes
                                    if subplot = 1: knick_point = [0, nobs[system]-1, 0, 0]
                                    if subplot = 2: knick_point = [0, knick1, nobs[system]-1, 0]
                                    if subplot = 2: knick_point = [0, knick1, knick2, nobs[system]-1]
        interceps (array)      : Intercepts
        rw (array)             : r2 value for the regression
        subplot (int, optional): _description_. Defaults to 1.
        fig (fig, optional)    : Current Matplotlib figure. Defaults to fig.
        ax (axe, optional)     : Current axes for the Matplotlib figure fig. Defaults to ax.
        f1w (file)             : File where to write the BIC results

    Returns:
        fig (fig)    : Current Matplotlib figure.
        ax (axe)     : Current axes for the Matplotlib figure fig.
    """
    if BIC == None:
        pddata = pd.DataFrame({'Y' : [min(data['HEIGHT']), max(data['HEIGHT'])],
                               'X' : [0, 10]})
        axi = sns.lineplot(data = pddata, x = 'X', y = 'Y', 
                     color = (1,0,0), linewidth = 4, 
                     ax = ax[systems.index(system), subplot-1])
        pddata = pd.DataFrame({'Y' : [min(data['HEIGHT']), max(data['HEIGHT'])],
                               'X' : [10, 0]})
        axi = sns.lineplot(data = pddata, x = 'X', y = 'Y', 
                     color = (1,0,0), linewidth = 4, 
                     ax = ax[systems.index(system), subplot-1])
    else:
        if subplot == 1:
            f1w.write('\t' + str(subplot) + ' slope computation -- BIC = ' + f"{BIC:.2f}" + '; Log(L) = ' + f"{log_L1:.2f}" + "\n")
        else:
            f1w.write('\t' + str(subplot) + ' slope(s) computation -- BIC = ' + f"{BIC:.2f}" + '; Log(L) = ' + f"{log_L1:.2f}" + "\n")
    
        # Do a loop on the number of segments/slopes (that correspond to the number of the subplot)
        for i in range (0, subplot):
            # Make an elevation vector to compute the ages from the linear relationship
            Y = np.arange(min(data['HEIGHT'][knick_point[i]:knick_point[i+1]]), 
                          #max(data['HEIGHT'][knick_point[i]:knick_point[i+1]]+ 1))
                          max(data['HEIGHT'][knick_point[i]:knick_point[i+1]+1]+ 1))
            pddata = pd.DataFrame({'Y' : Y,
                                   'X' : (Y - intercept_mean[i]) / exh_mean[i],
                                   'Xerrmax' : (Y - intercept_max[i]) / exh_max[i],
                                   'Xerrmin' : (Y - intercept_min[i]) / exh_min[i]})
            # Plot the regression line
            sns.lineplot(data = pddata, x = 'X', y = 'Y', 
                         color = (0,0,0), dashes = False, linewidth = 2, 
                         ax = ax[systems.index(system), subplot-1])
        
            # Plot the error bars
            sns.lineplot(data = pddata, x = 'Xerrmax', y = 'Y', 
                         color = (0.5, 0.5, 0.5), dashes = True, linewidth = 2, 
                         ax = ax[systems.index(system), subplot-1])
            sns.lineplot(data = pddata, x = 'Xerrmin', y = 'Y', 
                         color = (0.5, 0.5, 0.5), dashes = True, linewidth = 2, 
                         ax = ax[systems.index(system), subplot-1])

            # Write the results in an associated text file
            AER_text3 = '\t\t' + str(i+1) + '. Slope : ' + f"{exh_mean[i]*1e-3:.2f}" + ' km/Ma\n'
            f1w.write(AER_text3)
            AER_text4 = '\t\t  min = ' + f"{exh_min[i]*1e-3:.2f}" + ' km/Ma -- max = ' f"{exh_max[i]*1e-3:.2f}" + ' km/Ma\n'
            f1w.write(AER_text4)
            if subplot > 1:
                AER_text5 = '\t\t  between ' + f"{data[system][knick_point[i]]:.2f}" + ' Ma and ' + f"{data[system][knick_point[i+1]]:.2f}" + ' Ma\n'
                f1w.write(AER_text5)
            AER_text6 = '\t\t  R^2 = ' + f"{rw[i]:.2f}" + '\n\n'
            f1w.write(AER_text6)

            # Plot the real observations
            axi = sns.scatterplot(x = system, y = 'HEIGHT', data = data, 
                              markers = 'ok', ci = None,
                              ax = ax[systems.index(system), subplot-1])
        axi.errorbar(data[system], data['HEIGHT'], xerr = data['D' + system], fmt = 'ok')

        # Plot a vertical red dashed segments to mark the knickpoints
        for i in range (1, len(knick_point) - 1):
            axi.axvline(x = data[system][knick_point[i]], 
                        ymax = (data['HEIGHT'][knick_point[i]] - axi.get_ylim()[0])/(axi.get_ylim()[1] - axi.get_ylim()[0]), 
                        ls = '--', color = 'red', lw = 0.5)

    # Set the legends
    if max(data['HEIGHT']) < 10 :
        axi.set(xlabel = '%s age [Ma]' %(system), ylabel = "Elevation [km]")
        factor = 1  
    else:
        axi.set(xlabel = '%s age [Ma]' %(system), ylabel = "Elevation [m]")
        factor = 1000
    axi.set_xlim(0)

    if BIC != None:
        # Add the BIC value on the graph
        AER_text1 = 'BIC: ' + f"{BIC:.2f}" + '\n' + 'logL: ' + f"{log_L1:.2f}"
        props = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.5)
        axi.text(0.5 * max(data[system])/10, max(data['HEIGHT']) - factor * 0.2, AER_text1, bbox = props)  

    return fig, ax


#############################################################################
def BIC(filename_data = None, systems = None, nslopes = None):

    #warnings.filterwarnings("ignore")
    print('***********************************************************************')
    print('This is a BIC computation for thermochronologic vertical profiles')
    print('\tWritten by Xavier Robert (ISTerre - IRD, France)')
    print('\tport from the Matlab file used and described in Glotzbach et al., 2011')
    print('**********************************************************************\n')
    # Read the data
    print('Reading data...')
    fdata, nobs = readData(filename_data, systems)
    # Compute BIC and plot the graphs
    #print('Computing BIC...')
    calculateBIC(filename_data, fdata, nobs, systems, nslopes)
    print('\n')

    return
