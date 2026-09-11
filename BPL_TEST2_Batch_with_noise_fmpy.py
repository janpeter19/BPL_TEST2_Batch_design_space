# setup data TEST2_Batch_with_noise_fmpy 
# Author: Jan Peter Axelsson
#------------------------------------------------------------------------------------------------------------------
# 2026-09-11 - Created from earlier script that originates from the fall 2022. Today FMU uses BPL 2.3.2.
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
from fmpy import simulate_fmu
from fmpy import read_model_description

# Set the environment - for Linux a JSON-file in the FMU is read
if platform.system() == 'Linux': locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')

#------------------------------------------------------------------------------------------------------------------
#  Setup application FMU
#------------------------------------------------------------------------------------------------------------------

# Provde the right FMU and load for different platforms in user dialogue:
if platform.system() == 'Windows':
   print('Windows - run FMU pre-compiled JModelica 2.14')
   fmu_model ='BPL_TEST2_BatchWithNoise_windows_jm_cs.fmu'        
   model_description = read_model_description(fmu_model)  
   flag_vendor = 'JM'
   flag_type = 'CS'
elif platform.system() == 'Linux':
   flag_vendor = 'OM'
   flag_type = 'ME'
   if flag_vendor in ['OM','om']:
      print('Linux - run FMU pre-compiled OpenModelica') 
      if flag_type in ['CS','cs']:      
         fmu_model ='BPL_TEST2_BatchWithNoise_linux_om_cs.fmu'       
         model_description = read_model_description(fmu_model)  
      if flag_type in ['ME','me']:         
         fmu_model ='BPL_TEST2_BatchWithNoise_linux_om_me.fmu'  
#        fmu_model ='BPL_TEST2_BatchWithNoise_linux_2404_om_me.fmu'   
         model_description = read_model_description(fmu_model)  
   else:    
      print('There is no FMU for this platform')   
   
# Provide various opts-profiles
if flag_type in ['CS', 'cs']:
   opts_std = {'NCP': 500}
elif flag_type in ['ME', 'me']:
   opts_std = {'NCP': 500}
else:    
   print('There is no FMU for this platform')
  
# Provide various MSL and BPL versions
if flag_vendor in ['JM', 'jm']:
   constants = [v for v in model_description.modelVariables if v.causality == 'local'] 
   MSL_usage = [x[1] for x in [(constants[k].name, constants[k].start) \
                     for k in range(len(constants))] if 'MSL.usage' in x[0]][0]   
   MSL_version = [x[1] for x in [(constants[k].name, constants[k].start) \
                       for k in range(len(constants))] if 'MSL.version' in x[0]][0]
   BPL_version = [x[1] for x in [(constants[k].name, constants[k].start) \
                       for k in range(len(constants))] if 'BPL.version' in x[0]][0] 
elif flag_vendor in ['OM', 'om']:
   MSL_usage = '4.1.0 - used components: Noise.NormalNoise' 
   MSL_version = '4.1.0'
   BPL_version = 'Bioprocess Library version 2.3.2' 
else:    
   print('There is no FMU for this platform')

#------------------------------------------------------------------------------------------------------------------

# Simulation time
simulationTime = 5.0

# Dictionary of time discrete states
timeDiscreteStates = {} 

# Create stateValue that later will be used to store final state and used for initialization in 'cont':
stateValue =  {}
stateValue = {variable.derivative.name:None for variable in model_description.modelVariables \
                                            if variable.derivative is not None}
stateValue.update(timeDiscreteStates) 

global stateValueInitial; stateValueInitial = {}
for key in stateValue.keys():
    if not key[-1] == ']':
         if key[-3:] == 'I.y':
            stateValueInitial[key] = key[:-10]+'I_start'
         elif key[-3:] == 'D.x':
            stateValueInitial[key] = key[:-10]+'D_start'
         else:
            stateValueInitial[key] = key+'_start'
    elif key[-3] == '[':
        stateValueInitial[key] = key[:-3]+'_start'+key[-3:]
    elif key[-4] == '[':
        stateValueInitial[key] = key[:-4]+'_start'+key[-4:]
    elif key[-5] == '[':
        stateValueInitial[key] = key[:-5]+'_start'+key[-5:] 
    else:
        print('The state vector has more than 1000 states')
        break

stateValueInitialLoc = {}
for value in stateValueInitial.values():
    stateValueInitialLoc[value] = value

# Define a minimal compoent list of the model as a starting point for describe('parts')
component_list_minimum = ['bioreactor', 'bioreactor.culture']

# Provide process diagram on disk
fmu_process_diagram = 'BPL_TEST2_Batch_with_noise_om.png'

#------------------------------------------------------------------------------------------------------------------
#  Specific application constructs: stateValue, parValue, parLocation, parCheck, diagrams, ax, lines
#------------------------------------------------------------------------------------------------------------------

# Create dictionaries parValue[] and parLocation[]
parValue = {}
parValue['V_start'] = 1.0
parValue['VX_start'] = 1.0
parValue['VS_start'] = 10.0

parValue['Y'] = 0.5
parValue['qSmax'] = 1.0
parValue['Ks'] = 0.1

parValue['S_min'] = 1.0
parValue['time_final_max'] = 6.0
parValue['X_final_min'] = 5.0
parValue['sigma'] = 0.48

parValue['samplePeriod'] = 0.1
parValue['seed'] = 1
parValue['useGlobalSeed'] = False
parValue['useAutomaticLocalSeed'] = False

parLocation = {}
parLocation['V_start'] = 'bioreactor.V_start'
parLocation['VX_start'] = 'bioreactor.m_start[1]' 
parLocation['VS_start'] = 'bioreactor.m_start[2]' 

parLocation['Y'] = 'bioreactor.culture.Y'
parLocation['qSmax'] = 'bioreactor.culture.qSmax'
parLocation['Ks'] = 'bioreactor.culture.Ks'

parLocation['S_min'] = 'monitor.S_min'
parLocation['time_final_max'] = 'monitor.time_final_max'
parLocation['X_final_min'] = 'monitor.X_final_min'
parLocation['sigma'] = 'sensor.sigma'

parLocation['samplePeriod'] = 'sensor.samplePeriod'
parLocation['seed'] = 'sensor.noise.fixedLocalSeed'
parLocation['useGlobalSeed'] = 'sensor.noise.useGlobalSeed'
parLocation['useAutomaticLocalSeed'] = 'sensor.noise.useAutomaticLocalSeed'

# Parameter value check - especially for hysteresis to avoid runtime error
parCheck = []
parCheck.append("parValue['Y'] > 0")
parCheck.append("parValue['qSmax'] > 0")
parCheck.append("parValue['Ks'] > 0")
parCheck.append("parValue['V_start'] > 0")
parCheck.append("parValue['VX_start'] >= 0")
parCheck.append("parValue['VS_start'] >= 0")

# Extended list of parameters and variables only for describe and not change
keyVariables = []
parLocation['mu'] = 'bioreactor.culture.mu'; keyVariables.append(parLocation['mu'])
parLocation['bioreactor.c[1]'] = 'bioreactor.c[1]'; keyVariables.append(parLocation['bioreactor.c[1]'])
parLocation['bioreactor.c[2]'] = 'bioreactor.c[2]'; keyVariables.append(parLocation['bioreactor.c[2]'])
parLocation['bioreactor.culture.q[1]'] = 'bioreactor.culture.q[1]'; keyVariables.append(parLocation['bioreactor.culture.q[1]'])
parLocation['monitor.time_final'] = 'monitor.time_final'; keyVariables.append(parLocation['monitor.time_final'])
parLocation['monitor.X_final'] = 'monitor.X_final'; keyVariables.append(parLocation['monitor.X_final'])
parLocation['monitor.batch_evaluation'] = 'monitor.batch_evaluation'; keyVariables.append(parLocation['monitor.batch_evaluation'])
parLocation['monitor.S_min'] = 'monitor.S_min'; keyVariables.append(parLocation['monitor.S_min'])
parLocation['sensor.out.c[2]'] = 'sensor.out.c[2]'; keyVariables.append(parLocation['sensor.out.c[2]'])

# Create list of diagrams to be plotted by simu()
diagrams = []

# Create an empty list axes to be defined in newplot() and plotted by simu() or show()
ax = []

# Create list of pens for the diagrams
lines = ['-','--',':','-.']
