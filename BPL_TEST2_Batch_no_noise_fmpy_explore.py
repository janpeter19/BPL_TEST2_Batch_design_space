# Figure - Simulation of batch reactor 
#          with a set of functions and global variables added to facilitate explorative simulation work.
#          The general part of this code is called FMU-explore and is planned to be avaialbe as a separate package.
#
# MIT License
# Copyright (c) 2022, Jan Peter Axelsson, All rights reserved.
#------------------------------------------------------------------------------------------------------------------
# 2020-01-30 - Adpated to Python 3 and tested using Docker-Ubuntu FMU
#            - Introduced locale and setting of it - important for OpenModelica
#            - Introduced check of platform to adapt code for Windows/Linux in dialog
#            - print()
#            - np.nan
# 2020-01-31 - Object oriented diagrams
#            - Simplified handling of simulation results 
# 2020-02-01 - Adjusted describeMeium() to get the look better
# 2020-02-02 - Now only for Python 3, i.e. special Python 2 script for compilation
# 2020-02-04 - Update describe() to include time
# 2020-02-20 - Modified the script to be like for the BPL - student version
# 2020-02-20 - Modified the script for Jupyter: simu() no plt.show() and use opts
#------------------------------------------------------------------------------------------------------------------
# 2020-07-13 - Start for BP6a_test2 inherited from br5
# 2020-07-27 - Introduce choice of Linux FMU - JModelica or OpenModelica
# 2020-07-27 - Change of simu('cont') and handling of stateDict and model.get..
# 2020-10-01 - Upddated with new BP6a from BP6c
# 2020-10-10 - Simplified Yxs to Y
# 2020-11-21 - Adapted to ReactorType with n_inlets, n_outlets and n_ports
# 2021-02-04 - Adjust describe() for change to liquidphase
#------------------------------------------------------------------------------------------------------------------
# 2021-02-10 - Adapted for BPL_v2
# 2021-03-20 - Adapted for BPL ver 2.0.3
# 2021-04-25 - Adapted for BPL ver 2.0.4
#------------------------------------------------------------------------------------------------------------------
# 2021-05-06 - Adapted for BPL ver 2.0.5 and updated BPL interface
# 2021-05-10 - Tidy up indentations
# 2021-05-17 - Enable also to run on Linux
# 2021-05-19 - Corrected indentation in describe() and default value give for disp()
# 2021-05-22 - Improved decribe() to include short names in par[]
# 2021-05-24 - Shortened the command in simu() to update parameters and initial values
# 2021-05-25 - Change to parDict and parLocation and use par() for parDict.update()
# 2021-05-26 - Introduced init() again that is useful for clarity but par() can still be used 
# 2021-06-05 - Updated init() to handle also strings
# 2021-06-24 - Moved to BPL_dev and test on integration of newplot() and simu()
# 2021-06-25 - For newplot() added on alternative "PhasePlane" again and took away safetey around eval in simu()
# 2021-08-05 - Introduced describe_parts() and corrected disp() to handle number of displayed decimals 
# 2021-09-13 - Tested with BPL ver 2.0.7
# 2021-10-02 - Updated system_info() with FMU-explore version
# 2022-01-25 - Updated to FMU-explore 0.8.8
# 2022-02-01 - Updated to FMU-explore 0.8.9
# 2022-03-26 - Updated to FMU-explore 0.9.0 - model.reset(), and par(), init()
# 2022-04-29 - Updated to FMU-explore 0.9.1 
# 2022-05-28 - Introduce variable mu in parLocation for use in describe() but also disp()
# 2022-09-03 - Adapt for design space studies
# 2022-09-04 - FMU-explore extended exception list in describe_parts but could perhaps improve MSL.usage further
# 2022-09-05 - Introduced describe_MSL() to handle MSL information, both version and usage
# 2022-09-06 - Introduced new plotTypes for newplot()
# 2022-09-07 - Adjust for no noise system
# 2022-09-22 - Updated for FMU-explore 0.9.4 -include scipy version when used in the notebook
# 2022-10-05 - Updated for FMU-explore 0.9.5 with disp() that do not include extra parameters with parLocation
# 2023-02-09 - Updated to FMU-explore 0.9.6e
# 2023-02-13 - Consolidate FMU-explore to 0.9.6 and means parCheck and par() udpate and simu() with opts as arg
# 2023-02-24 - Corrected MSL-usage information for OpenModelica Linux
# 2023-03-29 - Update FMU-explore 0.9.7
# 2023-03-30 - Update FMU-explore 0.9.7 for FMPy
# 2023-05-31 - Adjusted to from importlib.meetadata import version
# 2024-05-27 - Update of FMU-explore for FMPy to 1.0.0
#------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------
#  Framework
#------------------------------------------------------------------------------------------------------------------

# Setup framework
import sys
import platform
import locale
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.image as img
import zipfile 

from fmpy import simulate_fmu
from fmpy import read_model_description
import fmpy as fmpy

from itertools import cycle
from importlib.metadata import version

# Set the environment - for Linux a JSON-file in the FMU is read
if platform.system() == 'Linux': locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')

#------------------------------------------------------------------------------------------------------------------
#  Setup application FMU
#------------------------------------------------------------------------------------------------------------------

# Provde the right FMU and load for different platforms in user dialogue:
global fmu_model, model
if platform.system() == 'Windows':
   print('Windows - run FMU pre-compiled JModelica 2.14')
   flag_vendor = 'JM'
   flag_type = 'CS'
   fmu_model ='BPL_TEST2_BatchNoNoise_windows_jm_cs.fmu'        
   model_description = read_model_description(fmu_model)   
elif platform.system() == 'Linux':
   flag_vendor = 'OM'
   flag_type = 'ME'
   if flag_vendor in ['OM','om']:
      print('Linux - run FMU pre-compiled OpenModelica 1.23.0-dev') 
      if flag_type in ['CS','cs']:         
         fmu_model ='BPL_TEST2_BatchNoNoise_linux_om_cs.fmu'    
         model_description = read_model_description(fmu_model)  
      if flag_type in ['ME','me']:         
         fmu_model ='BPL_TEST2_BatchNoNoise_linux_om_me.fmu'    
         model_description = read_model_description(fmu_model)  
   else:    
      print('There is no FMU for this platform')

# Provide various opts-profiles
if flag_type in ['CS', 'cs']:
   opts_std = {'ncp': 500}
elif flag_type in ['ME', 'me']:
   opts_std = {'ncp': 500}
else:    
   print('There is no FMU for this platform')
  
# Provide various MSL and BPL versions
if flag_vendor in ['JM', 'jm']:
#   MSL_usage = model.get('MSL.usage')[0]
   constants = [v for v in model_description.modelVariables if v.causality == 'local'] 
   MSL_usage = [x[1] for 
      x in [(constants[k].name, constants[k].start) for k in range(len(constants))] if 'MSL.usage' in x[0]][0]   
   MSL_version = [x[1] for 
      x in [(constants[k].name, constants[k].start) for k in range(len(constants))] if 'MSL.version' in x[0]][0]
   BPL_version = [x[1] for 
      x in [(constants[k].name, constants[k].start) for k in range(len(constants))] if 'BPL.version' in x[0]][0] 
elif flag_vendor in ['OM', 'om']:
   MSL_usage = '3.2.3 - used components: none' 
   MSL_version = '3.2.3'
   BPL_version = 'Bioprocess Library version 2.2.0' 
else:    
   print('There is no FMU for this platform')

# Simulation time
global simulationTime; simulationTime = 5.0
global prevFinalTime; prevFinalTime = 0

# Dictionary of time discrete states
timeDiscreteStates = {} 

# Define a minimal compoent list of the model as a starting point for describe('parts')
component_list_minimum = ['bioreactor', 'bioreactor.culture']

# Provide process diagram on disk
fmu_process_diagram ='BPL_TEST2_Perfusion_process_diagram_om.png'

#------------------------------------------------------------------------------------------------------------------
#  Specific application constructs: stateDict, parDict, diagrams, newplot(), describe()
#------------------------------------------------------------------------------------------------------------------

# Create stateDict that later will be used to store final state and used for initialization in 'cont':
global stateDict; stateDict =  {}
stateDict = {variable.derivative.name:None for variable in model_description.modelVariables \
                                            if variable.derivative is not None}
stateDict.update(timeDiscreteStates) 

global stateDictInitial; stateDictInitial = {}
for key in stateDict.keys():
    if not key[-1] == ']':
         if key[-3:] == 'I.y':
            stateDictInitial[key] = key[:-10]+'I_start'
         elif key[-3:] == 'D.x':
            stateDictInitial[key] = key[:-10]+'D_start'
         else:
            stateDictInitial[key] = key+'_start'
    elif key[-3] == '[':
        stateDictInitial[key] = key[:-3]+'_start'+key[-3:]
    elif key[-4] == '[':
        stateDictInitial[key] = key[:-4]+'_start'+key[-4:]
    elif key[-5] == '[':
        stateDictInitial[key] = key[:-5]+'_start'+key[-5:] 
    else:
        print('The state vector has more than 1000 states')
        break

global stateDictInitialLoc; stateDictInitialLoc = {}
for value in stateDictInitial.values():
    stateDictInitialLoc[value] = value

# Create dictionaries parDict[] and parLocation[]
global parDict; parDict = {}
parDict['V_start'] = 1.0
parDict['VX_start'] = 1.0
parDict['VS_start'] = 10.0

parDict['Y'] = 0.5
parDict['qSmax'] = 1.0
parDict['Ks'] = 0.1

parDict['S_min'] = 1.0
parDict['time_final_max'] = 6.0
parDict['X_final_min'] = 5.0

global parLocation; parLocation = {}
parLocation['V_start'] = 'bioreactor.V_start'
parLocation['VX_start'] = 'bioreactor.m_start[1]' 
parLocation['VS_start'] = 'bioreactor.m_start[2]' 

parLocation['Y'] = 'bioreactor.culture.Y'
parLocation['qSmax'] = 'bioreactor.culture.qSmax'
parLocation['Ks'] = 'bioreactor.culture.Ks'

parLocation['S_min'] = 'monitor.S_min'
parLocation['time_final_max'] = 'monitor.time_final_max'
parLocation['X_final_min'] = 'monitor.X_final_min'

# Parameter value check - especially for hysteresis to avoid runtime error
global parCheck; parCheck = []
parCheck.append("parDict['Y'] > 0")
parCheck.append("parDict['qSmax'] > 0")
parCheck.append("parDict['Ks'] > 0")
parCheck.append("parDict['V_start'] > 0")
parCheck.append("parDict['VX_start'] >= 0")
parCheck.append("parDict['VS_start'] >= 0")

# Extended list of parameters and variables only for describe and not change
global key_variables; key_variables = []
parLocation['mu'] = 'bioreactor.culture.mu'; key_variables.append(parLocation['mu'])
parLocation['bioreactor.c[1]'] = 'bioreactor.c[1]'; key_variables.append(parLocation['bioreactor.c[1]'])
parLocation['bioreactor.c[2]'] = 'bioreactor.c[2]'; key_variables.append(parLocation['bioreactor.c[2]'])
parLocation['bioreactor.culture.q[1]'] = 'bioreactor.culture.q[1]'; key_variables.append(parLocation['bioreactor.culture.q[1]'])
parLocation['monitor.time_final'] = 'monitor.time_final'; key_variables.append(parLocation['monitor.time_final'])
parLocation['monitor.X_final'] = 'monitor.X_final'; key_variables.append(parLocation['monitor.X_final'])
parLocation['monitor.batch_evaluation'] = 'monitor.batch_evaluation'; key_variables.append(parLocation['monitor.batch_evaluation'])
parLocation['monitor.S_min'] = 'monitor.S_min'; key_variables.append(parLocation['monitor.S_min'])
parLocation['sensor.out.c[2]'] = 'sensor.out.c[2]'; key_variables.append(parLocation['sensor.out.c[2]'])

# Create list of diagrams to be plotted by simu()
global diagrams
diagrams = []

# Define standard diagrams
def newplot(title='Batch cultivation', plotType='TimeSeries'):
   """ Standard plot window
        title = ''
       two possible diagrams
        diagram = 'TimeSeries' default
        diagram = 'PhasePlane' """
    
   # Reset pens
   setLines()

   # Transfer of global axes to simu()
   global ax, ax1, ax2, ax3, ax4, ax5

   # Plot diagram 
   if plotType == 'TimeSeries':

      plt.figure()
      ax1 = plt.subplot(4,1,1)
      ax2 = plt.subplot(4,1,2)
      ax3 = plt.subplot(4,1,3)
      ax4 = plt.subplot(4,1,4)

      ax1.set_title(title)  
      ax1.set_ylabel('X and S [g/L]'); ax1.grid()     
      ax2.set_ylabel('S measured [g/L]'); ax2.grid()      
      ax3.set_ylabel('mu [1/h]'); ax3.grid()
      ax4.set_ylabel('Batch evaluation'); ax4.set_xlabel('Time [h]'); ax4.grid()      
          
      # List of commands to be executed by simu() after a simulation  
      diagrams.clear()
      diagrams.append("ax1.plot(sim_res['time'],sim_res['bioreactor.c[1]'],color='r',linestyle=linetype)")
      diagrams.append("ax1.plot(sim_res['time'],sim_res['bioreactor.c[2]'],color='b',linestyle=linetype)")   
      diagrams.append("ax1.legend(['X','S'])")   
      diagrams.append("ax2.plot(sim_res['time'],sim_res['sensor.out.c[2]'],color='b',linestyle=linetype)")   
      diagrams.append("ax2.plot([0, simulationTime], [model.get('monitor.S_min'), model.get('monitor.S_min')],color='g',linestyle='--')")     
      diagrams.append("ax3.plot(sim_res['time'],sim_res['bioreactor.culture.q[1]'],color='r',linestyle=linetype)")   
      diagrams.append("ax4.step(sim_res['time'],sim_res['monitor.batch_evaluation'],where='post',color='b',linestyle=linetype)") 

   elif plotType == 'TimeSeries_2':

      plt.figure()
      ax1 = plt.subplot(4,1,1)
      ax2 = plt.subplot(4,1,2)
      ax3 = plt.subplot(4,1,3)
      ax4 = plt.subplot(4,1,4)

      ax1.set_title(title)  
      ax1.set_ylabel('X [g/L]'); ax1.grid()     
      ax2.set_ylabel('S [g/L]'); ax2.grid()    
      ax3.set_ylabel('mu [1/h]'); ax3.grid()   
      ax4.set_ylabel('Batch evaluation'); ax4.set_xlabel('Time [h]'); ax4.grid()      
          
      # List of commands to be executed by simu() after a simulation  
      diagrams.clear()
      diagrams.append("ax1.plot(sim_res['time'],sim_res['bioreactor.c[1]'],color='b',linestyle=linetype)")
      diagrams.append("ax2.plot(sim_res['time'],sim_res['bioreactor.c[2]'],color='b',linestyle=linetype)")  
      diagrams.append("ax2.plot([0, simulationTime], [model_get('monitor.S_min'), model_get('monitor.S_min')],color='g',linestyle='--')") 
      diagrams.append("ax3.plot(sim_res['time'],sim_res['bioreactor.culture.q[1]'],color='b',linestyle=linetype)")      
      diagrams.append("ax4.step(sim_res['time'],sim_res['monitor.batch_evaluation'],where='post',color='b',linestyle=linetype)") 

   elif plotType == 'TimeSeries_2_diagrams':

      plt.figure()
      ax1 = plt.subplot(4,1,1)
      ax2 = plt.subplot(4,1,2)
      ax3 = plt.subplot(4,1,3)
      ax4 = plt.subplot(4,1,4)

      ax1.set_title(title)  
      ax1.set_ylabel('X [g/L]'); ax1.grid()     
      ax2.set_ylabel('S [g/L]'); ax2.grid()    
      ax3.set_ylabel('mu [1/h]'); ax3.grid() 
      ax4.set_ylabel('Batch evaluation'); ax4.set_xlabel('Time [h]'); ax4.grid()      
          
      # List of commands to be executed by simu() after a simulation  
      diagrams.clear()

   elif plotType == 'PhasePlane':
    
      plt.figure()
      ax = plt.subplot(1,1,1)
    
      ax.set_title(title)
      ax.grid()
      ax.set_ylabel('S [g/L]')
      ax.set_xlabel('X [g/L]')

      # List of commands to be executed by simu() after a simulation         
      diagrams.clear()
      diagrams.append("ax.plot(sim_res['bioreactor.c[1]'],sim_res['bioreactor.c[2]'],color='b',linestyle=linetype)")
             
   else:
      print("Plot window type not correct")

# Define describtions partly coded here and partly taken from the FMU
def describe(name, decimals=3):
   """Look up description of culture, media, as well as parameters and variables in the model code"""
        
   if name == 'culture':
      print('Simplified text book model - only substrate S and cell concentration X')      
 
   elif name in ['broth', 'liquidphase', 'media']: 
      """Describe medium used"""
      X = model.get('liquidphase.X')[0] 
      X_description = model.get_variable_description('liquidphase.X') 
      X_mw = model.get('liquidphase.mw[1]')[0]
         
      S = model.get('liquidphase.S')[0] 
      S_description = model.get_variable_description('liquidphase.S')
      S_mw = model.get('liquidphase.mw[2]')[0]
         
      print()
      print('Reactor broth substances included in the model')
      print()
      print(X_description, '    index = ', X, 'molecular weight = ', X_mw, 'Da')
      print(S_description, 'index = ', S, 'molecular weight = ', S_mw, 'Da')
  
   elif name in ['parts']:
      describe_parts(component_list_minimum)

   elif name in ['MSL']:
      describe_MSL()

   else:
      describe_general(name, decimals)
      
#------------------------------------------------------------------------------------------------------------------
#  General code 
FMU_explore = 'FMU-explore for FMPy version 1.0.0'
#------------------------------------------------------------------------------------------------------------------

# Define function par() for parameter update
def par(parDict=parDict, parCheck=parCheck, parLocation=parLocation, *x, **x_kwarg):
   """ Set parameter values if available in the predefined dictionaryt parDict. """
   x_kwarg.update(*x)
   x_temp = {}
   for key in x_kwarg.keys():
      if key in parDict.keys():
         x_temp.update({key: x_kwarg[key]})
      else:
         print('Error:', key, '- seems not an accessible parameter - check the spelling')
   parDict.update(x_temp)
   
   parErrors = [requirement for requirement in parCheck if not(eval(requirement))]
   if not parErrors == []:
      print('Error - the following requirements do not hold:')
      for index, item in enumerate(parErrors): print(item)

# Define function init() for initial values update
def init(parDict=parDict, *x, **x_kwarg):
   """ Set initial values and the name should contain string '_start' to be accepted.
       The function can handle general parameter string location names if entered as a dictionary. """
   x_kwarg.update(*x)
   x_init={}
   for key in x_kwarg.keys():
      if '_start' in key: 
         x_init.update({key: x_kwarg[key]})
      else:
         print('Error:', key, '- seems not an initial value, use par() instead - check the spelling')
   parDict.update(x_init)

# Define fuctions similar to pyfmi model.get(), model.get_variable_descirption(), model.get_variable_unit()
def model_get(parLoc, model_description=model_description):
   """ Function corresponds to pyfmi model.get() but returns just a value and not a list"""
   par_var = model_description.modelVariables
   for k in range(len(par_var)):
      if par_var[k].name == parLoc:
         try:
            if par_var[k].name in start_values.keys():
                  value = start_values[par_var[k].name]
            elif par_var[k].variability in ['constant', 'fixed']:        
                  value = float(par_var[k].start)     
            elif par_var[k].variability == 'continuous':
               try:
                  timeSeries = sim_res[par_var[k].name]
                  value = timeSeries[-1]
               except (AttributeError, ValueError):
                  value = None
                  print('Variable not logged')
            else:
               value = None
         except NameError:
            print('Error: Information available after first simution')
            value = None
   return value

def model_get_variable_description(parLoc, model_description=model_description):
   """ Function corresponds to pyfmi model.get_variable_description() but returns just a value and not a list"""
   par_var = model_description.modelVariables
#   value = [x[1] for x in [(par_var[k].name, par_var[k].description) for k in range(len(par_var))] if parLoc in x[0]]
   value = [x.description for x in par_var if parLoc in x.name]   
   return value[0]
   
def model_get_variable_unit(parLoc, model_description=model_description):
   """ Function corresponds to pyfmi model.get_variable_unit() but returns just a value and not a list"""
   par_var = model_description.modelVariables
#   value = [x[1] for x in [(par_var[k].name, par_var[k].unit) for k in range(len(par_var))] if parLoc in x[0]]
   value = [x.unit for x in par_var if parLoc in x.name]
   return value[0]
      
# Define function disp() for display of initial values and parameters
def disp(name='', decimals=3, mode='short'):
   """ Display intial values and parameters in the model that include "name" and is in parLocation list.
       Note, it does not take the value from the dictionary par but from the model. """
   
   def dict_reverser(d):
      seen = set()
      return {v: k for k, v in d.items() if v not in seen or seen.add(v)}
   
   if mode in ['short']:
      k = 0
      for Location in [parLocation[k] for k in parDict.keys()]:
         if name in Location:
            if type(model_get(Location)) != np.bool_:
               print(dict_reverser(parLocation)[Location] , ':', np.round(model_get(Location),decimals))
            else:
               print(dict_reverser(parLocation)[Location] , ':', model_get(Location))               
         else:
            k = k+1
      if k == len(parLocation):
         for parName in parDict.keys():
            if name in parName:
               if type(model_get(Location)) != np.bool_:
                  print(parName,':', np.round(model_get(parLocation[parName]),decimals))
               else: 
                  print(parName,':', model_get(parLocation[parName])[0])

   if mode in ['long','location']:
      k = 0
      for Location in [parLocation[k] for k in parDict.keys()]:
         if name in Location:
            if type(model_get(Location)) != np.bool_:       
               print(Location,':', dict_reverser(parLocation)[Location] , ':', np.round(model_get(Location),decimals))
         else:
            k = k+1
      if k == len(parLocation):
         for parName in parDict.keys():
            if name in parName:
               if type(model_get(Location)) != np.bool_:
                  print(parLocation[parName], ':', dict_reverser(parLocation)[Location], ':', parName,':', 
                     np.round(model_get(parLocation[parName]),decimals))

# Line types
def setLines(lines=['-','--',':','-.']):
   """Set list of linetypes used in plots"""
   global linecycler
   linecycler = cycle(lines)

# Show plots from sim_res, just that
def show(diagrams=diagrams):
   """Show diagrams chosen by newplot()"""
   # Plot pen
   linetype = next(linecycler)    
   # Plot diagrams 
   for command in diagrams: eval(command)

# Define simulation
def simu(simulationTime=simulationTime, mode='Initial', options=opts_std, diagrams=diagrams):
   """Model loaded and given intial values and parameter before, and plot window also setup before."""   
   
   # Global variables
   global sim_res, prevFinalTime, stateDict, stateDictInitial, stateDictInitialLoc, start_values
   
   # Simulation flag
   simulationDone = False
   
   # Internal help function to extract variables to be stored
   def extract_variables(diagrams):
       output = []
       variables = [v for v in model_description.modelVariables if v.causality == 'local']
       for j in range(len(diagrams)):
           for k in range(len(variables)):
               if variables[k].name in diagrams[j]:
                   output.append(variables[k].name)
       return output

   # Run simulation
   if mode in ['Initial', 'initial', 'init']: 
      
      start_values = {parLocation[k]:parDict[k] for k in parDict.keys()}
      
      # Simulate
      sim_res = simulate_fmu(
         filename = fmu_model,
         validate = False,
         start_time = 0,
         stop_time = simulationTime,
         output_interval = simulationTime/options['ncp'],
         record_events = True,
         start_values = start_values,
         fmi_call_logger = None,
         output = list(set(extract_variables(diagrams) + list(stateDict.keys()) + key_variables))
      )
      
      simulationDone = True
      
   elif mode in ['Continued', 'continued', 'cont']:
      
      if prevFinalTime == 0: 
         print("Error: Simulation is first done with default mode = init'")
         
      else:         
         # Update parDictMod and create parLocationMod
         parDictRed = parDict.copy()
         parLocationRed = parLocation.copy()
         for key in parDict.keys():
            if parLocation[key] in stateDictInitial.values(): 
               del parDictRed[key]  
               del parLocationRed[key]
         parLocationMod = dict(list(parLocationRed.items()) + list(stateDictInitialLoc.items()))
   
         # Create parDictMod and parLocationMod
         parDictMod = dict(list(parDictRed.items()) + 
            [(stateDictInitial[key], stateDict[key]) for key in stateDict.keys()])      

         start_values = {parLocationMod[k]:parDictMod[k] for k in parDictMod.keys()}
  
         # Simulate
         sim_res = simulate_fmu(
            filename = fmu_model,
            validate = False,
            start_time = prevFinalTime,
            stop_time = prevFinalTime + simulationTime,
            output_interval = simulationTime/options['ncp'],
            record_events = True,
            start_values = start_values,
            fmi_call_logger = None,
            output = list(set(extract_variables(diagrams) + list(stateDict.keys()) + key_variables))
         )
      
         simulationDone = True
   else:
      
      print("Error: Simulation mode not correct")

   if simulationDone:
      
      # Plot diagrams from simulation
      linetype = next(linecycler)    
      for command in diagrams: eval(command)
   
      # Store final state values in stateDict:        
      for key in stateDict.keys(): stateDict[key] = model_get(key)  
         
      # Store time from where simulation will start next time
      prevFinalTime = sim_res['time'][-1]
      
   else:
      print('Error: No simulation done')
            
# Describe model parts of the combined system
def describe_parts(component_list=[]):
   """List all parts of the model""" 
       
   def model_component(variable_name):
      i = 0
      name = ''
      finished = False
      if not variable_name[0] == '_':
         while not finished:
            name = name + variable_name[i]
            if i == len(variable_name)-1:
                finished = True 
            elif variable_name[i+1] in ['.', '(']: 
                finished = True
            else: 
                i=i+1
      if name in ['der', 'temp_1', 'temp_2', 'temp_3', 'temp_4', 'temp_5', 'temp_6', 'temp_7']: name = ''
      return name
    
#   variables = list(model.get_model_variables().keys())
   variables = [v.name for v in model_description.modelVariables]
        
   for i in range(len(variables)):
      component = model_component(variables[i])
      if (component not in component_list) \
      & (component not in ['','BPL', 'Customer', 'today[1]', 'today[2]', 'today[3]', 'temp_2', 'temp_3']):
         component_list.append(component)
      
   print(sorted(component_list, key=str.casefold))

# Describe MSL   
def describe_MSL(flag_vendor=flag_vendor):
   """List MSL version and components used"""
   print('MSL:', MSL_usage)
 
# Describe parameters and variables in the Modelica code
def describe_general(name, decimals):
  
   if name == 'time':
      description = 'Time'
      unit = 'h'
      print(description,'[',unit,']')
      
   elif name in parLocation.keys():
      description = model_get_variable_description(parLocation[name])
      value = model_get(parLocation[name])
      try:
         unit = model_get_variable_unit(parLocation[name])
      except FMUException:
         unit =''
      if unit =='':
         if type(value) != np.bool_:
            print(description, ':', np.round(value, decimals))
         else:
            print(description, ':', value)            
      else:
        print(description, ':', np.round(value, decimals), '[',unit,']')
                  
   else:
      description = model_get_variable_description(name)
      value = model_get(name)
      try:
         unit = model_get_variable_unit(name)
      except FMUException:
         unit =''
      if unit =='':
         if type(value) != np.bool_:
            print(description, ':', np.round(value, decimals))
         else:
            print(description, ':', value)     
      else:
         print(description, ':', np.round(value, decimals), '[',unit,']')

# Plot process diagram
def process_diagram(fmu_model=fmu_model, fmu_process_diagram=fmu_process_diagram):   
   try:
       process_diagram = zipfile.ZipFile(fmu_model, 'r').open('documentation/processDiagram.png')
   except KeyError:
       print('No processDiagram.png file in the FMU, but try the file on disk.')
       process_diagram = fmu_process_diagram
   try:
       plt.imshow(img.imread(process_diagram))
       plt.axis('off')
       plt.show()
   except FileNotFoundError:
       print('And no such file on disk either')
         
# Describe framework
def BPL_info():
   print()
   print('Model for bioreactor has been setup. Key commands:')
   print(' - par()       - change of parameters and initial values')
   print(' - init()      - change initial values only')
   print(' - simu()      - simulate and plot')
   print(' - newplot()   - make a new plot')
   print(' - show()      - show plot from previous simulation')
   print(' - disp()      - display parameters and initial values from the last simulation')
   print(' - describe()  - describe culture, broth, parameters, variables with values/units')
   print()
   print('Note that both disp() and describe() takes values from the last simulation')
   print('and the command process_diagram() brings up the main configuration')
   print()
   print('Brief information about a command by help(), eg help(simu)') 
   print('Key system information is listed with the command system_info()')

def system_info():
   """Print system information"""
#   FMU_type = model.__class__.__name__
   constants = [v for v in model_description.modelVariables if v.causality == 'local']
   
   print()
   print('System information')
   print(' -OS:', platform.system())
   print(' -Python:', platform.python_version())
   try:
       scipy_ver = scipy.__version__
       print(' -Scipy:',scipy_ver)
   except NameError:
       print(' -Scipy: not installed in the notebook')
   print(' -FMPy:', version('fmpy'))
   print(' -FMU by:', read_model_description(fmu_model).generationTool)
   print(' -FMI:', read_model_description(fmu_model).fmiVersion)
   if model_description.modelExchange is None:
      print(' -Type: CS')
   else:
      print(' -Type: ME')
   print(' -Name:', read_model_description(fmu_model).modelName)
   print(' -Generated:', read_model_description(fmu_model).generationDateAndTime)
   print(' -MSL:', MSL_version)    
   print(' -Description:', BPL_version)   
   print(' -Interaction:', FMU_explore)
   
#------------------------------------------------------------------------------------------------------------------
#  Startup
#------------------------------------------------------------------------------------------------------------------

BPL_info()
