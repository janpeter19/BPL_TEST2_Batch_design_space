# Setup application functions BPL_TEST2_Batch_no_noise, dependent on previous import of functions from fmu_explore
# Author: Jan Peter Axelsson
# ------------------------------------------------------------------------------------------------------------------
# 2026-08-28 - Created from script that originates back in 2020. Today FMU uses BPL 2.3.2.
# 2026-09-28 - Brought back numpy since quite often used in the notebooks
# 2026-09-28 - Change indentaiton from 3 spaces to 4 using black
# 2026-10-19 - Modify for oms and add a clean-up section at the end
# ------------------------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------
#  Framework
# --------------------------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------------------------------
#  Specific application functions: newplot(), describe()
# --------------------------------------------------------------------------------------------------


def newplot(title="Batch cultivation", plotType="TimeSeries"):
    """Standard plot window
     title = ''
    two possible diagrams
     diagram = 'TimeSeries' default
     diagram = 'PhasePlane'"""

    # Reset pens
    resetPen()

    # Plot diagram
    if plotType == "TimeSeries":

        ax1 = plt.subplot(4, 1, 1)
        ax2 = plt.subplot(4, 1, 2)
        ax3 = plt.subplot(4, 1, 3)
        ax4 = plt.subplot(4, 1, 4)

        ax.clear()
        ax.append(ax1)
        ax.append(ax2)
        ax.append(ax3)
        ax.append(ax4)

        ax[0].set_title(title)
        ax[0].set_ylabel("X and S [g/L]")
        ax[0].grid()
        ax[1].set_ylabel("S measured [g/L]")
        ax[1].grid()
        ax[2].set_ylabel("mu [1/h]")
        ax[2].grid()
        ax[3].set_ylabel("Batch evaluation")
        ax[3].set_xlabel("Time [h]")
        ax[3].grid()

        # List of commands to be executed by simu() after a simulation
        diagrams.clear()
        diagrams.append(
            "ax[0].plot(sim_res['time'],sim_res['Batch.bioreactor.c[1]'],color='r',linestyle=linetype)"
        )
        diagrams.append(
            "ax[0].plot(sim_res['time'],sim_res['Batch.bioreactor.c[2]'],color='b',linestyle=linetype)"
        )
        diagrams.append("ax[0].legend(['X','S'])")
        diagrams.append(
            "ax[1].plot(sim_res['time'],sim_res['sensor.out.c[2]'],color='b',linestyle=linetype)"
        )
#        diagrams.append(
#            "ax[1].plot([0, simulationTime], [model.get('monitor.S_min'), model.get('monitor.S_min')],color='g',linestyle='--')"
#        )
        diagrams.append(
            "ax[2].plot(sim_res['time'],sim_res['Batch.bioreactor.culture.q[1]'],color='r',linestyle=linetype)"
        )
        diagrams.append(
            "ax[3].step(sim_res['time'],sim_res['monitor.batch_evaluation'],where='post',color='b',linestyle=linetype)"
        )

    elif plotType == "TimeSeries_2":

        ax1 = plt.subplot(4, 1, 1)
        ax2 = plt.subplot(4, 1, 2)
        ax3 = plt.subplot(4, 1, 3)
        ax4 = plt.subplot(4, 1, 4)

        ax.clear()
        ax.append(ax1)
        ax.append(ax2)
        ax.append(ax3)
        ax.append(ax4)

        ax[0].set_title(title)
        ax[0].set_ylabel("X [g/L]")
        ax[0].grid()
        ax[1].set_ylabel("S [g/L]")
        ax[1].grid()
        ax[2].set_ylabel("mu [1/h]")
        ax[2].grid()
        ax[3].set_ylabel("Batch evaluation")
        ax[3].set_xlabel("Time [h]")
        ax[3].grid()

        # List of commands to be executed by simu() after a simulation
        diagrams.clear()
        diagrams.append(
            "ax[0].plot(sim_res['time'],sim_res['Batch.bioreactor.c[1]'],color='b',linestyle=linetype)"
        )
        diagrams.append(
            "ax[1].plot(sim_res['time'],sim_res['Batch.bioreactor.c[2]'],color='b',linestyle=linetype)"
        )
#        diagrams.append(
#            "ax[1].plot([0, simulationTime], [model.get('monitor.S_min'), model.get('monitor.S_min')],color='g',linestyle='--')"
#        )
        diagrams.append(
            "ax[2].plot(sim_res['time'],sim_res['Batch.bioreactor.culture.q[1]'],color='b',linestyle=linetype)"
        )
        diagrams.append(
            "ax[3].step(sim_res['time'],sim_res['monitor.batch_evaluation'],where='post',color='b',linestyle=linetype)"
        )

    elif plotType == "TimeSeries_2_diagrams":

        ax1 = plt.subplot(4, 1, 1)
        ax2 = plt.subplot(4, 1, 2)
        ax3 = plt.subplot(4, 1, 3)
        ax4 = plt.subplot(4, 1, 4)

        ax.clear()
        ax.append(ax1)
        ax.append(ax2)
        ax.append(ax3)
        ax.append(ax4)

        ax[0].set_title(title)
        ax[0].set_ylabel("X [g/L]")
        ax[0].grid()
        ax[1].set_ylabel("S [g/L]")
        ax[1].grid()
        ax[2].set_ylabel("mu [1/h]")
        ax[2].grid()
        ax[3].set_ylabel("Batch evaluation")
        ax[3].set_xlabel("Time [h]")
        ax[3].grid()

        # List of commands to be executed by simu() after a simulation
        diagrams.clear()

    elif plotType == "PhasePlane":

        ax1 = plt.subplot(1, 1, 1)

        ax.clear()
        ax.append(ax1)

        ax[0].set_title(title)
        ax[0].grid()
        ax[0].set_ylabel("S [g/L]")
        ax[0].set_xlabel("X [g/L]")

        # List of commands to be executed by simu() after a simulation
        diagrams.clear()
        diagrams.append(
            "ax[0].plot(sim_res['Batch.bioreactor.c[1]'],sim_res['Batch.bioreactor.c[2]'],color='b',linestyle=linetype)"
        )

    else:
        print("Plot window type not correct")


def describe(name, decimals=3):
    """Look up description of culture, media, as well as parameters and variables in the model code"""

    if name == "culture":
        print("Simplified text book model - only substrate S and cell concentration X")

    elif name in ["MSL"]:
        describe_MSL()

    else:
        Application.describe_general(name, decimals)


# --------------------------------------------------------------------------------------------------
#  Startup
# --------------------------------------------------------------------------------------------------

FMU_explore_info()


#------------------------------------------------------------------------------------------------------------------
#  Clean-up
#------------------------------------------------------------------------------------------------------------------

del parValue, parLocation, parCheck
del fmu_model, fmu_process_diagram
del MSL_usage, MSL_version, BPL_version
del options, simulationTime, timeDiscreteStates
del flag_type, flag_vendor
del lines

