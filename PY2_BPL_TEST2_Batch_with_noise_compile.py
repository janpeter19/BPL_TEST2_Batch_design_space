# Figure - Simulation of batch reactor 
#          with functions added to facilitate explorative simulation work
#
# MIT License
# Copyright (c) 2022, Jan Peter Axelsson, All rights reserved.
#
# Author: Jan Peter Axelsson
# 2012-10-16 - Created
# 2012-10-25 - Modified 
# 2013-08-31 - Modified 
# 2016-09-15 - Modified
# 2018-01-29 - Tested with JModelica version: 2.1 and OK
# 2018-02-08 - Adjusted graph to be similar to parameter sweep diagrams
# 2018-02-28 - Include functions that facilitate interactive usage
# 2018-03-04 - Now we have: newplot(), par(), init(),  simu() - all Simnon style!
# 2018-03-24 - Change to the structure of BPL and include disp()
# 2018-05-03 - Harmonize mode 'Continued' in simu() with the other explore-scripts
# 2018-05-04 - Added chocie of phase plane plot and made plotted variables global
# 2018-05-11 - Changed default parVal to nan to allow also to set parameter zero
# 2018-05-17 - Added function describe() that prints description of parameters
# 2018-05-19 - Tested for JModelica version 2.2 and OK
# 2018-05-23 - Added help function xRange to simplify parameter sweep
# 2018-08-24 - Added automatic change of linetype for multiple plots
# 2018-10-04 - Polish for br5.mo
# 2018-10-05 - Further polish for br5.mo
# 2018-10-12 - Modified newplot() so that the title can be an argument
# 2018-10-13 - Modified newplot() and handling of update of linetype
# 2018-10-15 - Create stateDict automatically
# 2018-10-24 - Updated describe() with reading unit
# 2019-01-18 - Start work to comply with subset of MSL media handling
# 2019-01-21 - Added information in package medium
# 2019-01-22 - Changed reactor to vector formulation instead - stateDict handled
# 2019-01-23 - Addedinstance of Medium in them main model just for Python access
# 2019-01-24 - Started work to make reactor general for different cultures - S,X
# 2019-01-25 - Added function describeMedium()
# 2019-02-01 - Adapted for br5e-version where only ReactorType needs Medium explicitely
# 2019-02-02 - Try to divide up the code in library and application
# 2019-02-17 - Resume work on dividing code
# 2019-02-18 - Changed to X=1 from X=2 as default
# 2019-03-11 - Adpated to br5g and use of LiquidCon and EnvironmentCon as parameters
# 2019-03-12 - Adpated to br5h and - inner and outer for culture
# 2019-05-21 - Adaptation to BR5j - changes following DEMO22
# 2019-05-24 - Adaptation to BR5k
#--------------------------------------------------------------------------------------
# 2020-02-27 - Python2 script for compilation and only for Windows
#            - Added system_info() that prints system information
#            - Change newplot and simu using objectoriented diagrams 
#            - Simplified handling of simulation results
#            - Tested with JModelica 2.14 and seems ok
#            - Adjusted describeBroth and changed its name
# 2020-03-16 - Adapted for BR5m with stream etc
# 2020-03-16 - Indluced in system_info() information if FMU is ME or CS
#--------------------------------------------------------------------------------------
# 2020-07-13 - Start BP6a_test2 from BR5m
# 2020-07-14 - Changed name model_file to application_file
# 2020-10-10 - Simplified Yxs to Y
# 2020-11-21 - Adapted to ReactorType with n_inlets, n_outlets and n_ports
#--------------------------------------------------------------------------------------
# 2021-02-10 - Adapt for BPL_v2
# 2021-03-19 - Change of inner/outer connection - Stéphande Veluts recomendation
# 2021-09-08 - Run it with BPL version 2.0.7 beta
# 2022-09-03 - Run it with BPL verion 2.1.0 beta and extended for TEST2B
#--------------------------------------------------------------------------------------

# Setup framework
import platform
import locale
import numpy as np 
import matplotlib.pyplot as plt
from itertools import cycle
from pymodelica import compile_fmu 
from pyfmi import load_fmu
from pyfmi.fmi import FMUException
import pyfmi as pyfmi                    # just to get version number actually
from importlib_metadata import version   # included in future Python 3.8

# Describe framework
def system_info():
    """Print system information"""
    FMU_type = model.__class__.__name__
    print
    print 'System information'
    print ' -OS:', platform.system()
    print ' -Python:', platform.python_version()
    print ' -PyFMI:', pyfmi.__version__
    print ' -FMU by:', model.get_generation_tool()
    print ' -FMI:', model.get_version()
    print ' -Type:', FMU_type
    print ' -Name:', model.get_name()
    print ' -Description:', model.get_description()  
    

# Define model file name and class name 
model_name = 'BPL_TEST2.BatchWithNoise' 
application_file = 'BPL_TEST2.mo'
library_file = 'Z:BPL/package.mo'

# Compile model
fmu_model = compile_fmu(model_name, [application_file, library_file], target='cs')

# Load model
global model; model = load_fmu(fmu_model)

