# setup application functions BPL_TEST2_Batch_with_noise, dependent on previous import of functions from fmu_explore 
# Author: Jan Peter Axelsson
#------------------------------------------------------------------------------------------------------------------
# 2026-08-28 - Created
#------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------
#  Specific application functions: newplot(), describe()
#------------------------------------------------------------------------------------------------------------------

# Define standard diagrams
def newplot(title='Batch cultivation - noise on S only', plotType='TimeSeries'):
   """ Standard plot window
        title = ''
       two possible diagrams
        diagram = 'TimeSeries' default
        diagram = 'PhasePlane' """
    
   # Reset pens
   resetPen()

   # Plot diagram 
   if plotType == 'TimeSeries':

      ax1 = plt.subplot(4,1,1)
      ax2 = plt.subplot(4,1,2)
      ax3 = plt.subplot(4,1,3)
      ax4 = plt.subplot(4,1,4)

      ax.clear()
      ax.append(ax1)
      ax.append(ax2)
      ax.append(ax3)
      ax.append(ax4)

      ax[0].set_title(title)  
      ax[0].set_ylabel('X and S [g/L]'); ax[0].grid()     
      ax[1].set_ylabel('S measured [g/L]'); ax[1].grid()      
      ax[2].set_ylabel('mu [1/h]'); ax[2].grid()
      ax[3].set_ylabel('Batch evaluation'); ax[3].set_xlabel('Time [h]'); ax[3].grid()      
          
      # List of commands to be executed by simu() after a simulation  
      diagrams.clear()
      diagrams.append("ax[0].plot(sim_res['time'],sim_res['bioreactor.c[1]'],color='r',linestyle=linetype)")
      diagrams.append("ax[0].plot(sim_res['time'],sim_res['bioreactor.c[2]'],color='b',linestyle=linetype)")   
      diagrams.append("ax[0].legend(['X','S'])")   
      diagrams.append("ax[1].plot(sim_res['time'],sim_res['sensor.out.c[2]'],color='b',linestyle=linetype)")   
      diagrams.append("ax[1].plot([0, simulationTime], [model.get('monitor.S_min'), model.get('monitor.S_min')],color='g',linestyle='--')")     
      diagrams.append("ax[2].plot(sim_res['time'],sim_res['bioreactor.culture.q[1]'],color='r',linestyle=linetype)")   
      diagrams.append("ax[3].step(sim_res['time'],sim_res['monitor.batch_evaluation'],where='post',color='b',linestyle=linetype)") 

   elif plotType == 'TimeSeries_2':

      ax1 = plt.subplot(5,1,1)
      ax2 = plt.subplot(5,1,2)
      ax3 = plt.subplot(5,1,3)
      ax4 = plt.subplot(5,1,4)
      ax5 = plt.subplot(5,1,5)

      ax.clear()
      ax.append(ax1)
      ax.append(ax2)
      ax.append(ax3)
      ax.append(ax4)
      ax.append(ax5)

      ax[0].set_title(title)  
      ax[0].set_ylabel('X [g/L]'); ax[0].grid()     
      ax[1].set_ylabel('S [g/L]'); ax[1].grid()    
      ax[2].set_ylabel('mu [1/h]'); ax[2].grid()   
      ax[3].set_ylabel('S measured [g/L]'); ax[3].grid()      
      ax[4].set_ylabel('Batch evaluation'); ax[4].set_xlabel('Time [h]'); ax[4].grid()      
          
      # List of commands to be executed by simu() after a simulation  
      diagrams.clear()
      diagrams.append("ax[0].plot(sim_res['time'],sim_res['bioreactor.c[1]'],color='b',linestyle=linetype)")
      diagrams.append("ax[1].plot(sim_res['time'],sim_res['bioreactor.c[2]'],color='b',linestyle=linetype)")  
      diagrams.append("ax[2].plot(sim_res['time'],sim_res['bioreactor.culture.q[1]'],color='b',linestyle=linetype)")     
      diagrams.append("ax[3].plot(sim_res['time'],sim_res['sensor.out.c[2]'],color='b',linestyle=linetype)")   
      diagrams.append("ax[3].plot([0, simulationTime], [model.get('monitor.S_min'), model.get('monitor.S_min')],color='g',linestyle='--')")     
      diagrams.append("ax[4].step(sim_res['time'],sim_res['monitor.batch_evaluation'],where='post',color='b',linestyle=linetype)") 
     
      
   elif plotType == 'TimeSeries_2_diagrams':

      ax1 = plt.subplot(4,1,1)
      ax2 = plt.subplot(4,1,2)
      ax3 = plt.subplot(4,1,3)
      ax4 = plt.subplot(4,1,4)
      
      ax.clear()
      ax.append(ax1)
      ax.append(ax2)
      ax.append(ax3)
      ax.append(ax4)

      ax[0].set_title(title)  
      ax[0].set_ylabel('X [g/L]'); ax[0].grid()     
      ax[1].set_ylabel('S [g/L]'); ax[1].grid()    
      ax[2].set_ylabel('mu [1/h]'); ax[2].grid() 
      ax[3].set_ylabel('Batch evaluation'); ax[3].set_xlabel('Time [h]'); ax[3].grid()      
          
      # List of commands to be executed by simu() after a simulation  
      diagrams.clear()

   elif plotType == 'PhasePlane':
    
      ax1 = plt.subplot(1,1,1)
      
      ax.clear()
      ax.append(ax1)
    
      ax[0].set_title(title)
      ax[0].grid()
      ax[0].set_ylabel('S [g/L]')
      ax[0].set_xlabel('X [g/L]')

      # List of commands to be executed by simu() after a simulation         
      diagrams.clear()
      diagrams.append("ax[0].plot(sim_res['bioreactor.c[1]'],sim_res['bioreactor.c[2]'],color='b',linestyle=linetype)")
                          
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
#  Startup
#------------------------------------------------------------------------------------------------------------------

FMU_explore_info()
