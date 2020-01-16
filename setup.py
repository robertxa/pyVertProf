######!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2020 Xavier Robert <xavier.robert@ird.fr>
# SPDX-License-Identifier: GPL-3.0-or-later


from setuptools import setup, find_packages

# Import of the lib pyVertProf
import pyVertProf

def readme():
	with open('README.rst') as f:
		return f.read()

setup(name='pyVertProf',
	#version='0.1.1',
	version=pyVertProf.__version__,
	description='package that provide tools to perform stats on vertical profils',
	long_descritpion=open('README.rst').read(),
	url='https://github.com/robertxa/pyVertProf',
	download_url='https://github.com/robertxa/pyVertProf/archive/master.zip',
	author='Xavier Robert',
	author_email='xavier.robert@univ-grenoble-alpes.fr',
	license='GPL-V3.0',
	packages=find_packages(),
	install_requires=[
	      'kapteyn'
	],
	classifiers=[
		#"Programming language :: Python",
		"Operating System :: OS Independent",
		#"Programming Language :: Python :: 3",
		"Topic :: Scientific/Engineering :: Visualization"
	#	"Topic :: Scientific/Engineering :: GIS"
	],
	include_package_data=True,
	zip_safe=False)
      