Tutorial
=========

The building and simulation of a model structure is performed in the pyFSE
package as follows:

Import the appropriate modules and packages for the analysis::

	import os
	import numpy as np
	import matplotlib.pyplot as plt

Load the pyFUSE package (make sure it is installed or added to the Python path)::
	
	import pyFUSE as pf
	
Load the example data of the Nete catchment::

	datapath='.\exampledata'
	rain=np.loadtxt(os.path.join(datapath,'Rain_Cal_warm_1jan02_31dec05'))
	evapo=np.loadtxt(os.path.join(datapath,'ET_Cal_warm_1jan02_31dec05'))
	Flow_obs=np.loadtxt(os.path.join(datapath,'Flow_Cal_Meas_13aug02_31dec05')) #data in m3/s
	ID=rain.size-Flow_obs.size #to exclude warming up period

To build a model, 3 setup files are needed:
	* Model structure file
	* Parameter file
	* constant values file

General functions can be used to load the files::	
	
	topt=pf.Set_options(os.path.join(datapath,'pyFUSE_strinput.txt'))
	tcnt=pf.Set_cnts(os.path.join(datapath,'pyFUSE_ctninput.txt'))
	tpar=pf.Set_pars(os.path.join(datapath,'pyFUSE_parinput.txt'))

Create the model instance::
	
	tmod = pf.pyFUSE_Model('tutorialmodel',topt,tpar,rain,evapo,tcnt)

Run the model and create output::

	output, info = tmod.run(run_id='Testrun')

Make a plot of the model run, with observations and simulation in m3/s::

	plt.plot(tmod.total_outflow(run_id='Testrun')[ID:],label='Simulation')
	plt.plot(Flow_obs,'--',label='Observation')
	plt.xlabel('Time (hour)')
	plt.ylabel(r'Flow ($m^3/s$)')
	plt.legend()
	plt.grid()
	plt.show()

.. image:: exampleoutput.png
	:width: 753px
	:height: 308px
	
Run 100 Monte Carlo Run of the model::

	nMC=100
	tmod.run_MC(nMC)
   
End analysis and close the hdf5 file::

	tmod.close_h5()