# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 13:30:00 2012

@author: VHOEYS
"""

import sys
import os

import numpy as np
import matplotlib.pyplot as plt

sys.path.append('D:\Modellen\Version2012\Main\Model\Lumped')

import pyFUSE as pf

###########################################################
# DATA NETE Rain_Cal_warm_1jan02_31dec05
###########################################################
datapath=".\exampledata"
rain=np.loadtxt(os.path.join(datapath,'Rain_Cal_warm_1jan02_31dec05'))
evapo=np.loadtxt(os.path.join(datapath,'ET_Cal_warm_1jan02_31dec05'))
Flow_obs=np.loadtxt(os.path.join(datapath,'Flow_Cal_Meas_13aug02_31dec05'))
ID=rain.size-Flow_obs.size

rain = rain[:150]
evapo = evapo[:150]

###########################################################
# Construct model
###########################################################

topt=pf.Set_options(os.path.join(datapath,'pyFUSE_strinput.txt'))
tcnt=pf.Set_cnts(os.path.join(datapath,'pyFUSE_ctninput.txt'))
tpar=pf.Set_pars(os.path.join(datapath,'pyFUSE_parinput.txt'))

tmod = pf.pyFUSE_Model('mdtt.hdf5',topt,tpar,rain,evapo,tcnt)

###########################################################
# Run the model
###########################################################

tmod.run(run_id='optimized')

###########################################################
# Plot output
###########################################################
#plt.plot(tmod.total_outflow(run_id='optimized')[ID:],label='Simulation')
#plt.plot(Flow_obs,'--',label='Observation')	
#plt.xlabel('Time (hour)')
#plt.ylabel(r'Flow ($m^3/s$)')
#plt.legend()
#plt.grid()
#plt.show()

#100 Monte Carlo Run of the model
#tmod.run_MC(5)
##plot the different outputs
#for i in range(1,3):
#    print str(i)
#    plt.plot(tmod.ofile['mdtt/MCrun'+str(i)+'/FLUX/FLUX_q_all'])
   
#for name, value in tmod.ofile['mdtt'].attrs.iteritems():
#    print name+":", value    
#tmod.close_h5()

######################
##Sobol variance base sensitivity
######################

#compare two pars
testpars={}
testpars['S1max'] = ModPar('S1max',51.,400.,150.,'randomUniform')
testpars['ku'] = ModPar('ku',0.004,2.,0.5,'randomUniform')
testpars['timeo'] = ModPar('timeo',4.,24.,20.,'randomUniform')
testpars['timei'] = ModPar('timei',4.,240.,200.,'randomUniform')
testpars['timeb'] = ModPar('timeb',400.,2400.,2000.,'randomUniform')

nbaseruns = 5
sens1 = SobolVariance(testpars, ModelType = 'pyFUSE')
sens1.SobolVariancePre(nbaseruns)
sens1.run_pyFUSE(tmod)

#average flow is not best comparison
output = np.zeros((nbaseruns*(2+sens1.ndim),1))
for i in range(nbaseruns*(2+sens1.ndim)):
    output[i,0] = tmod.total_outflow(run_id='SobolRun'+str(i))[:].max()
    plt.plot(tmod.total_outflow(run_id='SobolRun'+str(i))[:])
    
sens1.SobolVariancePost(output)
sens1.plotSi()
sens1.plotSTi()



#quick delete
#for i in range(nbaseruns*(2+sens1.ndim)):
#    del tmod.ofile['mdtt/SobolRun'+str(i)]






