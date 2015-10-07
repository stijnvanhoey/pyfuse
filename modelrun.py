#-------------------------------------------------------------------------------
# Name:        FUSE style FLEXIBLE HYDROLOGICAL MODEL
# Purpose:     extend the functionalities of the FORTRAN implementation by Clark
#
# Author:      vhoeys
#
# Created:
# Copyright:   (c) vhoeys 2011
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

#import os
#import sys
#
#import time
#import random
#import math

import matplotlib.pyplot as plt

#import fileinput
import numpy as np
import h5py
import time

from scipy.interpolate import interp1d
from scipy.integrate import odeint

try:
    import odespy
    odespy_import = True
except:
    odespy_import = False
    print("Odespy was not found, ODE integration is limited to 'odeint'!")

from pyFUSE.distributions import *
from pyFUSE.fluxes import *
from pyFUSE.modelsetup import *
from pyFUSE.parameter import *

###############################################################################
###   MODEL ENVIRONMENT
###############################################################################

class pyFUSE_Model():
    '''
    Create pyFUSE hydrological model, a python version and extention of the FORTRAN
    FUSE model environment by Clarke, 2008 [1]
    This is not a wrapper of the Fortran implementation of Clark, but a
    complete rewrite of the original code in order to make further extensions
    easier.

    Parameters
    ------------
    name: str
        A given name for the constructed model structure
    OPTIONS: dict
        Dictionary with the model options for construction
    PARS: dict
        Dictionary of parameters used for the model evaluation, coming from
        input parameter file
    RAIN: array
        numpy array containting the rainfall data
    EVAPO: array
        numpy array containting the potential evapotranspiration data,
        with same time frame and resolution as the rain
    CONST: dict
        dictionary of neede constant values for the model configuration and
        calculation
    INITFRAC: float
        fraction of maximum storage used as initial condition

    References
    -------------

    [1] Clark, Martyn P., A. G. Slater, D. E. Rupp, R. A. Woods, Jasper A.
    Vrugt, H. V. Gupta, Thorsten Wagener, and L. E. Hay. Framework for
    Understanding Structural Errors (FUSE): A modular framework to diagnose
    differences between hydrological models. Water Resources Research 44 (2008): 14.

    '''

    def __init__(self,name,OPTIONS,PARS,RAIN,EVAPO,CONST,INITFRAC=0.1,oldfile=False):
        '''
        initiate the model structure:
            get pars extraction
            get constants
            get rain/evapo data

        for initial conditions: automated fraction of stores, with frac-par

        RAIN en EVAPO: timeseries from same length: defining the simulation times

        See literature:
        Clark, Martyn P., A. G. Slater, D. E. Rupp, R. A. Woods, Jasper A.
        Vrugt, H. V. Gupta, Thorsten Wagener, and L. E. Hay. Framework for
        Understanding Structural Errors (FUSE): A modular framework to diagnose
        differences between hydrological models. Water Resources Research 44
        (2008): 14.

        Original fortran code from Clark, Martyn P.
        '''

        #cONTROL INPUT INSTANCES
        ##################################

        if isinstance(name,str):
            self.name = name

        if isinstance(OPTIONS,dict):
            self.OPTIONS=OPTIONS

        #LOAD PARAMETERS
        self.mcpossible=True
        if isinstance(PARS,str):
            print 'Parameterfile used'
            self.parameters = Set_pars(PARS) #parameters: dictionary with all information (dict of Modpars)
            print self.parameters
            self.PARS = Set_pars_for_run(self.parameters)
            print self.PARS
        elif isinstance(PARS.values()[0],ModPar):
            print 'Parameter library with distribution info'
            self.parameters=PARS
            self.PARS = Set_pars_for_run(self.parameters)
        elif isinstance(PARS.values()[0],float):
            print 'Parameter only values used'
            self.PARS = PARS
            self.mcpossible=False
        else:
            raise Exception('Bad parameter input type')

        if isinstance(CONST,dict):
            self.CONST=CONST

        if RAIN.size <> EVAPO.size:
            raise Exception('Timeserie length of RAIN and EVAPO need to be the same')
        else:
            self.RAIN=RAIN
            self.EVAPO=EVAPO
            print 'Simulation of the model will run for %d timesteps' %RAIN.size


        #Define the constants
        ##################################
        self.area = np.float64(self.CONST['area'])     #catchment area
        self.dt = np.float64(self.CONST['timestep'])  #relative to the frequency of the rain input
        self.totn = self.RAIN.size

        #Define time-steps to run output for
        # output to save is without last value
        ##################################
        self.totaltime=np.float(self.totn)
        self.totaltimesteps=self.totaltime/self.dt
        self.time= np.linspace(0,self.totaltime-1,self.totaltimesteps) #eg. hours if datainput hours and timestep=1.0; geeft mogelijkheid om met daginput ook uur-output te verkrijgen. NIET met uur ook dag (teveel gezever ivm hoe uurlijke info omzetten in dag, kan beter preprocessing zijn)
#        self.time=np.arange(0,int(self.totaltime),int(self.dt))

        #Prepare input data timeseries for linear interpolation in substepping
        ##################################
        self.time_int = np.linspace(0.0,self.totaltime-1,self.totaltime)

        #add dummy value to rain and ET to make sure calculation runs enough timesteps
#        self.RAIN=np.hstack((self.RAIN,np.zeros(2)))
#        self.EVAPO=np.hstack((self.EVAPO,np.zeros(2)))

        self.rain_int = interp1d(self.time_int,self.RAIN,kind='nearest',bounds_error=False,fill_value=0.0)
        self.evapo_int = interp1d(self.time_int,self.EVAPO,kind='nearest',bounds_error=False,fill_value=0.0)

        #LOAD PARLIB and CREATE RELEVANT ONE BASED ON OPTIONS -> get .optguess
        ##################################
        #self.map_PARS()

        #LOAD FLUXLIB and CREATE RELEVANT ONE BASED ON OPTIONS
        ##################################
        self.map_ROUT() #volgorde belangrijk!
        self.map_FLUXES()

        #LOAD STATESLIB and CREATE RELEVANT ONE BASED ON OPTIONS
        ##################################
        self.map_STATES()
        print self.STATES

        # get fraction of total volume to start with
        self.INITFRAC = INITFRAC

        #Create hdf5 instance to put all outputs in
        #the idea: a model creates a file; different conditions (length calc, parlength ...) different group -> given in attr
        #each flux and state is different dataset in the group
        #################################
        if oldfile == False:
            if self.name[-5:]=='.hdf5':
                self.ofile = h5py.File(self.name,'w')
                self.name=self.name[:-5]
                self.hdf5main=self.ofile.create_group(self.name)
            else:
                self.ofile = h5py.File(self.name+'.hdf5','w')
                self.hdf5main=self.ofile.create_group(self.name)

            for opt in self.OPTIONS:
                self.hdf5main.attrs[opt] = self.OPTIONS[opt]
        else:
             self.ofile = h5py.File(self.name+'.hdf5', 'r+')
             self.hdf5main = self.ofile[self.name]


        #MINI number to make sure some variables are never zero or below
        self.xmini=1e-6


    def check_impossible_options(self):
        '''
        Control method to check if the combined options are reasonable and
        physically meaningful

        TODO: make automatic and propose alternatives
        '''

#                ! don't allow a lower tension tank when there are two upper ones
#        IF (LIST_ARCH1(ISW_ARCH1)%MCOMPONENT(1:10).EQ.'tension2_1'.AND. &
#            LIST_ARCH2(ISW_ARCH2)%MCOMPONENT(1:10).EQ.'tens2pll_2') CYCLE
#        ! don't allow percolation below field capacity if there are multiple upper tanks
#        IF (LIST_ARCH1(ISW_ARCH1)%MCOMPONENT(1:10).NE.'onestate_1'.AND. &
#            LIST_QPERC(ISW_QPERC)%MCOMPONENT(1:10).EQ.'perc_w2sat') CYCLE
#            idem voor tension1_1 en perc_w2sat !!
#        freefirst evaporation not with tension2_1 => gelijke verdeling?!? dan niet
#        onestate_1 en sequential is zelfde als onestate_1 en freefirst
#
#        linear ET enkel bij onelayer
#                perc-nodrain only possible with overland-hymod en single linear baseflow?!?

        pass

    def check_impossible_pars(self):
        '''
        Control method to adapt parameters when values are not in correspondence
        to eachother, e.g. both samples from distributions for Monte Carlo analysis

        TODO: resample until good value
        '''

        pass

    def create_datasets(self,run_id,custom_period=None):
        '''
        Based on the relevant fluxes and variables, the datasets for saving
        in the hdf5 binary format are created.

        '''
        if custom_period:  #adjust length of the datasets
            raise Exception('not yet implemented')
        else:
            #-1 timesteps, as timestep is not saved for flux, only calculated as dummy
            data_length_state=self.totaltimesteps.copy()
            data_length_flux=self.totaltimesteps.copy()

        #CREATE GROUP FOR REQUESTED OUTPUT TIMESTEPS
        self.rungroup=self.hdf5main.create_group(run_id)
        self.stategroup=self.rungroup.create_group('STATE')
        for state in self.state_str:
            self.stategroup.create_dataset('STATE_'+state, (int(data_length_state),), '=f8')

        self.fluxgroup=self.rungroup.create_group('FLUX')
        for flux in self.flux_str:
            self.fluxgroup.create_dataset('FLUX_'+flux, (int(data_length_flux),), '=f8')

        self.pargroup=self.rungroup.create_group('PAR')
        for par in self.PARS:
            if par == 'frac_future':
                self.pargroup.create_dataset('PAR_'+par , (int(self.PARS['frac_future'].size),), '=f8')
            else:
                self.pargroup.create_dataset('PAR_'+par , (int(1),), '=f8')


        #CREATE SEPERATE GROUP FOR ALL OUTPUT CALCULATION STEPS
        #        self.rungroup=self.ofile.create_group(run_id+'_ALLSTEPS')

    def statetoh5(self,run_id,y_out):
        '''
        Save outputs in hdf5 arrays

        Order of working is the same in init, odeint and state_str, making sure values are correct
        states to save is [INIT t1... tn], so timesteps+1 length
        '''
        cnt=0
        for st in self.state_str:
            self.ofile[self.name+'/'+run_id+'/STATE/STATE_'+st][:] = y_out[:,cnt]
            cnt+=1

    def partoh5(self,run_id):
        '''
        Save parameters of run in hdf5 arrays
        '''
        for par in self.PARS:
            if par =='frac_future':
                self.ofile[self.name+'/'+run_id+'/PAR/PAR_'+par][:] = self.PARS[par][:]
            else:
                self.ofile[self.name+'/'+run_id+'/PAR/PAR_'+par][:] = self.PARS[par]

    def fluxtoh5(self,run_id,y_out):
        '''
        Save fluxes to h5

        (implementation is not ideal and sensitive to mistakes)
        fluxes to save is [t1.. tn], so timesteps length (1 shorter than length states)
        '''
        data_length_flux=self.totaltimesteps.copy()

        for tstep in range(int(data_length_flux)):
            #update of the states
            self.state_lib_update(y_out[tstep,:])
            #update of the rain,evapo
            self.FLUXES['rain']=self.RAIN[tstep]
            self.FLUXES['pet']=self.EVAPO[tstep]

            self.calc_fluxes()
            for fl in self.flux_str:
                if 'rout' in fl:
                    self.ofile[self.name+'/'+run_id+'/FLUX/FLUX_'+fl][tstep] = self.FLUXES[fl][-1]
                elif 'qgamma' in fl:
                    pass
                elif 'qfuture' in fl:
                    pass
                else:
                    self.ofile[self.name+'/'+run_id+'/FLUX/FLUX_'+fl][tstep] = self.FLUXES[fl]

    def map_PARS(self):
        '''
        Only select and work with the relevant parameters for the model under consideration
        '''
        pass

    def map_ROUT(self):
        """
        Map only the relevant states to work with in the model
        If routing is linear: the 3-integer string option need to be interpretedand translated in different structural options
        If zero routing reservoirs, no extra fluxes are used in the model environment
        """
        self.ROUTLIB={}
        if self.OPTIONS['routing'] == 'rout_ind':
            #read the integer and translate in 3 options
#            if type(self.OPTIONS['reservoirs']) is 'str' and len(self.OPTIONS['reservoirs']) == 3:
            if isinstance(self.OPTIONS['reservoirs'], str) and len(self.OPTIONS['reservoirs']) == 3:
                try:
                    dummy = int(self.OPTIONS['reservoirs'])
                except:
                    raise Exception('the reservoir option is not representing a combination of three integer values')
            else:
                raise Exception('Make sure the reservoir option is a 3-integer combination representing the number of reservoirs in the subflows')


            if self.OPTIONS['reservoirs'][0] == '0':
                print 'overland flow is not routed by extra linear reservoirs'
            else:
                self.ROUTLIB['routover'] = int(self.OPTIONS['reservoirs'][0])


            if self.OPTIONS['reservoirs'][2] == '0':
                print 'baseflow is not routed by extra linear reservoirs'
            else:
                self.ROUTLIB['routbase'] = int(self.OPTIONS['reservoirs'][2])


            if self.OPTIONS['interflow'] == 'intflwnone':
                print 'no interflow routing, since interflow is not added in the model; second integer in reservoirs option is f no relevance'
            else:
                if self.OPTIONS['reservoirs'][1] == '0':
                    print 'interflow is not routed by extra linear reservoirs'
                else:
                    self.ROUTLIB['routinter'] = int(self.OPTIONS['reservoirs'][1])


    def map_STATES(self):
        '''
        Translates the STATES LIB towards the array of variables used by the solver
        and create the lib of specific variables used in the calculation

       Returns
       --------
       self.STATES : dict
           dict with model specific states
       self.state_str : list
           list of strings with variable names (specific order to help mapping in odeint function)
       self.state_size : int
           number of states in the model
        '''

#        if self.OPTIONS['soilstorage'] == 'onelayer':
#            self.state_size=0
#            self.state_str=[]
#
#            self.STATES={}
#            #WE FOCUSSEN ONS HIEROP OM HET TE TESTEN!!
#            if self.OPTIONS['uplayer'] == 'easytest':
#                self.state_str.append('S')
#                self.STATES['S']=0.
#                self.state_size+=1
#                print self.state_size
#                print self.state_str

        if self.OPTIONS['soilstorage'] == 'twolayer':
            self.state_size=0
            self.state_str=[]

            self.STATES={}
            if self.OPTIONS['uplayer']== 'onestate_1':
                self.state_str.append('S1')
                self.STATES['S1']=0.
                self.state_size+=1

            elif self.OPTIONS['uplayer']== 'tension1_1' or \
            self.OPTIONS['uplayer']== 'surface1_1':
                self.state_str.append('S1T')
                self.STATES['S1T']=0.
                self.state_str.append('S1F')
                self.STATES['S1F']=0.
                self.state_size+=2

            elif self.OPTIONS['uplayer']== 'tension2_1':
                self.state_str.append('S1TA')
                self.STATES['S1TA']=0.
                self.state_str.append('S1TB')
                self.STATES['S1TB']=0.
                self.state_str.append('S1F')
                self.STATES['S1F']=0.
                self.state_size+=3

        #op basis van identifier t -> all fluxen seven in een dataformaat!!
        ####################################################################

            if self.OPTIONS['lowlayer_baseflow'] == 'tens2pll_2':
                self.state_str.append('S2T')
                self.STATES['S2T']=0.
                self.state_str.append('S2FA')
                self.STATES['S2FA']=0.
                self.state_str.append('S2FB')
                self.STATES['S2FB']=0.
                self.state_size+=3

            elif self.OPTIONS['lowlayer_baseflow'] == 'unlimfrc_2' or \
            self.OPTIONS['lowlayer_baseflow'] == 'unlimpow_2' or \
            self.OPTIONS['lowlayer_baseflow'] == 'topmdexp_2' or \
            self.OPTIONS['lowlayer_baseflow'] == 'fixedsiz_2':
                self.state_str.append('S2')
                self.STATES['S2']=0.
                self.state_size+=1

    def map_FLUXES(self):
        '''
        Translates the STATES LIB towards the array of variables used by odeint
        and create the lib of specific variables used in the calculation

        Returns
        ---------

        self.FLUXES: dict
            dict with model specific states
        self.flux_str : str
            list of strings with variable names (specific order to help mapping in odeint function)
        self.flux_size : int
            number of states in the model

        !also pet (potential evaporation )is used as flux to get a total output of all info
        '''
        #rain and pet are present in every model
        self.flux_str=[]
        self.flux_str.append('rain')
        self.flux_str.append('pet')

        self.flux_size=0
        self.flux_size+=2

        self.FLUXES={}
        self.FLUXES['rain']=0.
        self.FLUXES['pet']=0.


        #ONELAYER
#        if self.OPTIONS['soilstorage'] == 'onelayer':
#
#            #WE FOCUSSEN ONS HIEROP OM HET SYSTEEM TE TESTEN!!
#            if self.OPTIONS['uplayer'] == 'easytest':
#                self.flux_str.append('e')
#                self.FLUXES['e']=0.
#                self.flux_str.append('qsx')
#                self.FLUXES['qsx']=0.
#                self.flux_size+=2

        #TWOLAYER
        if self.OPTIONS['soilstorage'] == 'twolayer':

            #uplayer
            if self.OPTIONS['uplayer']== 'onestate_1':
                self.flux_str.append('e1')
                self.FLUXES['e1']=0.
                self.flux_str.append('qufof')
                self.FLUXES['qufof']=0.
                self.flux_str.append('q12')
                self.FLUXES['q12']=0.
                self.flux_str.append('qif')
                self.FLUXES['qif']=0.
                self.flux_str.append('qsx')
                self.FLUXES['qsx']=0.
                self.flux_size+=5

            elif self.OPTIONS['uplayer']== 'surface1_1':
                self.flux_str.append('e1A')
                self.FLUXES['e1A']=0.
                self.flux_str.append('e1B')
                self.FLUXES['e1B']=0.
                self.flux_str.append('e1')
                self.FLUXES['e1']=0.
                self.flux_str.append('qstof')
                self.FLUXES['qstof']=0.
                self.flux_str.append('q12')
                self.FLUXES['q12']=0.
                self.flux_str.append('qif')
                self.FLUXES['qif']=0.
                self.flux_str.append('qsx')
                self.FLUXES['qsx']=0.
                self.flux_size+=7

            elif self.OPTIONS['uplayer']== 'tension1_1':
                self.flux_str.append('e1')
                self.FLUXES['e1']=0.
                self.flux_str.append('qutof')
                self.FLUXES['qutof']=0.
                self.flux_str.append('qufof')
                self.FLUXES['qufof']=0.
                self.flux_str.append('q12')
                self.FLUXES['q12']=0.
                self.flux_str.append('qif')
                self.FLUXES['qif']=0.
                self.flux_str.append('qsx')
                self.FLUXES['qsx']=0.
                self.flux_size+=6

            elif self.OPTIONS['uplayer']== 'tension2_1':
                self.flux_str.append('e1A')
                self.FLUXES['e1A']=0.
                self.flux_str.append('e1B')
                self.FLUXES['e1B']=0.
                self.flux_str.append('e1')
                self.FLUXES['e1']=0.
                self.flux_str.append('qurof')
                self.FLUXES['qurof']=0.
                self.flux_str.append('qutof')
                self.FLUXES['qutof']=0.
                self.flux_str.append('qufof')
                self.FLUXES['qufof']=0.
                self.flux_str.append('q12')
                self.FLUXES['q12']=0.
                self.flux_str.append('qif')
                self.FLUXES['qif']=0.
                self.flux_str.append('qsx')
                self.FLUXES['qsx']=0.
                self.flux_size+=9

            #LOWLAYER
            if self.OPTIONS['lowlayer_baseflow'] == 'tens2pll_2':
                self.flux_str.append('e2')
                self.FLUXES['e2']=0.
                self.flux_str.append('qbA')
                self.FLUXES['qbA']=0.
                self.flux_str.append('qbB')
                self.FLUXES['qbB']=0.
                self.flux_str.append('qb') #sum of qba and qbb
                self.FLUXES['qb']=0.
                self.flux_str.append('qstof')
                self.FLUXES['qstof']=0.
                self.flux_str.append('qsfofA')
                self.FLUXES['qsfofA']=0.
                self.flux_str.append('qsfofB')
                self.FLUXES['qsfofB']=0.
                self.flux_str.append('qsfof') #sum of qsfofa and qsfofB
                self.FLUXES['qsfof']=0.
                self.flux_size+=8


            elif self.OPTIONS['lowlayer_baseflow'] == 'unlimfrc_2' or \
            self.OPTIONS['lowlayer_baseflow'] == 'unlimpow_2' or \
            self.OPTIONS['lowlayer_baseflow'] == 'topmdexp_2' or \
            self.OPTIONS['lowlayer_baseflow'] == 'fixedsiz_2':
                self.flux_str.append('e2')  #zero in certain options
                self.FLUXES['e2']=0.
                self.flux_str.append('qsfof')
                self.FLUXES['qsfof']=0.
                self.flux_str.append('qb')
                self.FLUXES['qb']=0.
                self.flux_size+=3

        #the extra routing components
        if self.OPTIONS['routing'] == 'rout_ind':
            for key, value in self.ROUTLIB.items():
#                print key,value,'routinginfo'
                self.flux_str.append(key)
                self.FLUXES[key]=np.float64(0.0)*np.ones(value) #array with number of fluxes the number of basins
                self.flux_size+=1
#            print self.FLUXES
        elif self.OPTIONS['routing'] == 'rout_all1':
            self.flux_str.append('qgamma')
            self.FLUXES['qgamma']=0.
            self.flux_str.append('qfuture')
            self.FLUXES['qfuture']=np.zeros(self.PARS['frac_future'].size)   #ARRAY!!
            self.flux_size+=2


        self.flux_str.append('q_all')
        self.FLUXES['q_all']=0.  #sum of qsx,qb and qif coming out of storage model
        self.flux_str.append('e_all')
        self.FLUXES['e_all']=0.  #sum of e1 and e2 coming out of storage model
        self.flux_size+=2


    def state_lib_update(self,cal_states):
        '''
        Update the states in the dictionary

        Var-names in string (fixed places comparable with vars in odeint) used to map towards dictionary
        + calculate derived states
        '''

        ##  GET VALUES FROM ODEINT AND PUT IN LIBRARY
        cnt=0
        for st in self.state_str:
            self.STATES[st]=cal_states[cnt]
#            print cal_states[cnt]
            cnt+=1


        ## Calculate derived STATES, needed to correspond to the flux calculations
        # this makes the STATE-lib larger, but the state_str remaines and is the most important for the administration
        if self.OPTIONS['soilstorage'] == 'twolayer':
            #UPLAYER
            if self.OPTIONS['uplayer']== 'onestate_1':
                '''
                FSTATE%TENS_1  = MIN(FSTATE%WATR_1, DPARAM%MAXTENS_1)      ! tension storage
                FSTATE%FREE_1  = MAX(XMIN, FSTATE%WATR_1 - DPARAM%MAXTENS_1) ! free storage
                '''
                self.STATES['S1T'] = min(self.STATES['S1'],self.PARS['S1Tmax'])
                self.STATES['S1F'] = max(self.xmini,self.STATES['S1'] - self.PARS['S1Tmax'])

            elif self.OPTIONS['uplayer']== 'surface1_1':
                '''
                NAM implementation of the upper layer

                '''
                self.STATES['S1'] = self.STATES['S1T'] + self.STATES['S1F']

            elif self.OPTIONS['uplayer']== 'tension1_1':
                '''
                FSTATE%WATR_1  = FSTATE%TENS_1 + FSTATE%FREE_1             ! total storage
                '''
                self.STATES['S1'] = self.STATES['S1T'] + self.STATES['S1F']

            elif self.OPTIONS['uplayer']== 'tension2_1':
                '''
                 FSTATE%TENS_1  = FSTATE%TENS_1A + FSTATE%TENS_1B           ! tension storage
                 FSTATE%WATR_1  = FSTATE%TENS_1  + FSTATE%FREE_1            ! total storage
                '''
                self.STATES['S1T'] = self.STATES['S1TA'] + self.STATES['S1TB']
                self.STATES['S1'] = self.STATES['S1T'] + self.STATES['S1F']

            else:
                raise Exception('wrong option for upper layer')

            #LOWLAYER
            if self.OPTIONS['lowlayer_baseflow'] == 'tens2pll_2':
                '''
                FSTATE%FREE_2  = FSTATE%FREE_2A + FSTATE%FREE_2B              ! free storage
                FSTATE%WATR_2  = FSTATE%TENS_2  + FSTATE%FREE_2               ! total storage
                '''
                self.STATES['S2F'] = self.STATES['S2FA'] + self.STATES['S2FB']
                self.STATES['S2'] = self.STATES['S2T'] + self.STATES['S2F']

            elif self.OPTIONS['lowlayer_baseflow'] == 'unlimfrc_2' or \
            self.OPTIONS['lowlayer_baseflow'] == 'unlimpow_2' or \
            self.OPTIONS['lowlayer_baseflow'] == 'topmdexp_2' or \
            self.OPTIONS['lowlayer_baseflow'] == 'fixedsiz_2':
                '''
                FSTATE%TENS_2  = MIN(FSTATE%WATR_2, DPARAM%MAXTENS_2)         ! tension storage
                FSTATE%FREE_2  = MAX(0._sp, FSTATE%WATR_2 - DPARAM%MAXTENS_2) ! free storage
                '''
                self.STATES['S2T'] = min(self.STATES['S2'], self.PARS['S2Tmax'])
                self.STATES['S2F'] = max(0.,self.STATES['S2'] - self.PARS['S2Tmax'])


    def calc_fluxes(self):
        '''
        Calculate fluxes, based on the current state values,
        follow up of the different flux calculations is important
        '''

        if self.OPTIONS['soilstorage']=='twolayer':
            self.FLUXES = self.Evaporation.calc(self.STATES,self.FLUXES)
            self.FLUXES = self.Surface.calc(self.STATES,self.FLUXES) #MOET VOOR BASEFLOW EN PERCOLATION GEBEUREN!!!!
            self.FLUXES = self.Percolation.calc(self.STATES,self.FLUXES)
            self.FLUXES = self.Interflow.calc(self.STATES,self.FLUXES)
            self.FLUXES = self.Baseflow.calc(self.STATES,self.FLUXES)
            self.FLUXES = self.Misscell.calc(self.STATES,self.FLUXES)
            self.FLUXES = self.Routing.calc(self.STATES,self.FLUXES,self.ROUTLIB)  #moet laatst



    def deriv_mod(self,cal_states,t,raint,evapo,PARS,OPTIONS):  #odeint
        '''
        Calculate derivatives
        gets libraries: STATES,PARS and OPTIONS
        gets input-timeseries taint,evapo

        TODO: control if fluxes to large or too small -> repair

        '''
        ##--------------------------------------------------------------------
        #NEW IDEA!: also an update fluxes - combined with options and based on the mapped fluxes
        # -> give all STATES and calculate the fluxes (dict), based on OPTIONS and the fluxes in it
        # run evaporation
        # run ...
        #dan ook routing analytisch in 1 zwoeng te doen! en geen trucskes met updaten en zo
        ##--------------------------------------------------------------------

        ##--------------------------------------------------------------------
        #update rain and et
        rains=raint(t)  #needs to be interp1d object
        self.FLUXES['rain']=rains
        pet=evapo(t)    #needs to be interp1d object
        self.FLUXES['pet']=pet #!!! ZP vooraleer de ander fluxen te gebeuren!!
        ##--------------------------------------------------------------------


        savetimefreq=500.
        if t > savetimefreq*self.timecounter:
            self.timecounter+=1
            print 'Time step ',str(savetimefreq*self.timecounter),' calculated...'
            print '    '

        ##--------------------------------------------------------------------
        #update library of STATES with current values to calculate the fluxes
        #the derived STATES are already calculated
        self.state_lib_update(cal_states)
        ##--------------------------------------------------------------------

        ##--------------------------------------------------------------------
        #update library of FLUXES with current values to calculate the NEW STATES
        self.calc_fluxes()
        ##--------------------------------------------------------------------

        ##--------------------------------------------------------------------
        ## Seperate library to map later on with real library
        state_deriv={}
        ##
        dStates = np.zeros(self.state_size)
        ##--------------------------------------------------------------------

        if OPTIONS['soilstorage'] == 'onelayer':
#            if OPTIONS['uplayer'] == 'VHMstyle':
#                S = (self.rain_int - qsx) - self.Evaporation.calc(self.STATES,pet,lay='up') -qif -qufof #hier qufof nodig? TE VERANDEREN
#            if OPTIONS['uplayer'] == 'logSPM':
#                S = self.rain_int - q_quick - q_rge - q_et


            #WE FOCUSSEN ONS HIEROP OM HET TE TESTEN!!
            if OPTIONS['uplayer'] == 'easytest':
#                S = rains - self.Evaporation.calc(self.STATES,pet) - self.Surface.testeasy(self.STATES)
                dS = self.FLUXES['rain'] - self.FLUXES['e'] - self.FLUXES['qsx']
                state_deriv['dS']=dS

            else:
                raise Exception('Not yet implemented')


        elif OPTIONS['soilstorage'] == 'twolayer':  #! voor ET: gebruik lay='up' and lay='low'
            #UPLAYER
            if OPTIONS['uplayer']== 'onestate_1':
                state_deriv['dS1'] = (self.FLUXES['rain'] - self.FLUXES['qsx']) \
                - self.FLUXES['e1'] - self.FLUXES['q12'] - self.FLUXES['qif'] \
                - self.FLUXES['qufof']

            elif OPTIONS['uplayer']== 'surface1_1':
                state_deriv['dS1F'] = (self.FLUXES['rain'] - self.FLUXES['qsx']) \
                - self.FLUXES['e1A'] - self.FLUXES['q12'] - self.FLUXES['qif'] \
                - self.FLUXES['qstof']

                state_deriv['dS1T'] = self.FLUXES['qstof'] - self.FLUXES['e1B']

            elif OPTIONS['uplayer']== 'tension1_1':
                state_deriv['dS1T'] = (self.FLUXES['rain'] - self.FLUXES['qsx']) \
                - self.FLUXES['e1'] - self.FLUXES['qutof']

                state_deriv['dS1F'] = self.FLUXES['qutof'] - self.FLUXES['q12'] \
                - self.FLUXES['qif'] - self.FLUXES['qufof']

            elif OPTIONS['uplayer']== 'tension2_1':
                state_deriv['dS1TA'] = (self.FLUXES['rain'] - self.FLUXES['qsx']) \
                - self.FLUXES['e1A'] - self.FLUXES['qurof']

                state_deriv['dS1TB'] = self.FLUXES['qurof'] - self.FLUXES['e1B']  \
                - self.FLUXES['qutof']

                state_deriv['dS1F'] = self.FLUXES['qutof'] - self.FLUXES['q12'] \
                - self.FLUXES['qif'] - self.FLUXES['qufof']
            else:
                raise Exception('wrong option for upper layer')

            #LOWLAYER
            if OPTIONS['lowlayer_baseflow'] == 'tens2pll_2':
                '''
                DY_DT%TENS_2  = M_FLUX%QPERC_12*(1._SP-MPARAM%PERCFRAC) - M_FLUX%EVAP_2 - M_FLUX%TENS2FREE_2
                DY_DT%FREE_2A = M_FLUX%QPERC_12*(MPARAM%PERCFRAC/2._SP) + M_FLUX%TENS2FREE_2/2._SP - M_FLUX%QBASE_2A &
                  - M_FLUX%OFLOW_2A
                DY_DT%FREE_2B = M_FLUX%QPERC_12*(MPARAM%PERCFRAC/2._SP) + M_FLUX%TENS2FREE_2/2._SP - M_FLUX%QBASE_2B &
                  - M_FLUX%OFLOW_2B
                '''
                state_deriv['dS2T'] = (1. - self.PARS['kappa']) * self.FLUXES['q12'] \
                - self.FLUXES['e2'] - self.FLUXES['qstof']
                state_deriv['dS2FA'] = self.FLUXES['q12'] * (self.PARS['kappa']/2.) \
                + self.FLUXES['qstof']/2. - self.FLUXES['qbA'] - self.FLUXES['qsfofA']
                state_deriv['dS2FB'] = self.FLUXES['q12'] * (self.PARS['kappa']/2.) \
                + self.FLUXES['qstof']/2. - self.FLUXES['qbB'] - self.FLUXES['qsfofB']

            elif OPTIONS['lowlayer_baseflow'] == 'unlimfrc_2' or \
            OPTIONS['lowlayer_baseflow'] == 'unlimpow_2' or \
            OPTIONS['lowlayer_baseflow'] == 'topmdexp_2' or \
            OPTIONS['lowlayer_baseflow'] == 'fixedsiz_2':
#                if OPTIONS['uplayer']== 'surface1_1':
#                    state_deriv['dS2']  = self.FLUXES['qstof'] - self.FLUXES['e2']
#                else:
                state_deriv['dS2']  = self.FLUXES['q12'] - self.FLUXES['e2'] - self.FLUXES['qb'] - self.FLUXES['qsfof']

        ##--------------------------------------------------------------------
        ## check if state_deriv and state_str are of same size
        if len(state_deriv) <> len(self.state_str):
            print 'state deriv is'
            print state_deriv
            print 'state str is'
            print state_str
            raise Exception('Mismatch between calculated derivatives and number of state variables')

        ## Map states with self.state_str and put in right order
        cnt=0
        for st in self.state_str:
            dst='d'+st
            if dst in state_deriv:
                dStates[cnt] = state_deriv[dst]
            else:
                raise Exception('Mismatch between calculated derivatives and number of state variables')
            cnt+=1

        ##--------------------------------------------------------------------
        return dStates


    def get_all_runids(self):
        '''
        Get list of all used run_ids saved in the hdf5 file
        '''
        return self.ofile.keys()

    def clean_outputs(self,todelete):
        '''
        Delete the none interesting output groups and the related datasets

        Parameters
        --------------
        todelete: list
            list of strings with the groups to delete
        '''
        for ids in todelete:
            print ids
            try:
                del self.ofile[ids]
            except:
                print '%s is not an existing group, so not used as run_id for this model' %ids

    def update_pars(self,dic2update):
        '''
        Takes existing library and only updates the chosen parameters from
        the given dictionary
        '''
        for key, value in dic2update.items():
            if key in self.PARS:
                print 'Par value of parameter ', key,' changed: old value was ', self.PARS[key],', the new value is ',value
                self.PARS[key] = value
            else:
                print 'Parameter ',key,' not updated, since not in current dictionary'

    def load_new_pars(self,inff=None):
        '''
        Load a parameter set from a input parameter file or dict and put in parameter
        dictionary or load them from the pars given in dict

        Currently only from file is supported, directly passing a dictionary should also be possible and will be implemented

        Parameters
        -------------
        inff: str
            name of the input parameter file
        '''
        if isinstance(inff,str):
            pp=Set_pars(inff)
            ppr=Set_pars_for_run(pp)
            if isinstance(ppr,dict):
                self.PARS=pprn
        elif isinstance(inff,dict):
            #assumes ready for run
            #pars updated
#            self.PARS = inff  #old implementation: entire pardict
            self.update_pars(inff) #new, only relevant
            self.PARS = Derived_pars(self.PARS)

    def get_list_ode_integrators(self):
        return odespy.solvers.list_available_solvers()


    def run(self,custom_period=None,new_pars=False,run_id='testrun'):
        '''
        Run the model!

        Parameters
        ------------
        custom_period:
            calculate the model for a specific subperiod of the total data lenght
        new_pars: str
            textfile of the new parameters used to (re)run the model or
            a dictionary
        run_id: str
            used to identify the outputs of the specific modelrun,
            if nothing given, a testrun-group is added
        '''

        if new_pars:
            self.load_new_pars(new_pars)
            #check if run_id is is renewed
            if run_id in self.ofile.keys():
                raise Exception('Give new run id when parameter values are changed')
            else:
                print 'New parameterization saved under the run_id',run_id
        else:
            print 'Model simulation with given set of parameters'


        #GET METHODS FOR FLUXES -> create instance with the model specific inputs
        ##################################
        self.Evaporation = Evaporation(self.PARS,self.OPTIONS)
        self.Surface = Surface(self.PARS,self.OPTIONS)
        self.Percolation = Percolation(self.PARS,self.OPTIONS)
        self.Interflow = Interflow(self.PARS,self.OPTIONS)
        self.Baseflow = Baseflow(self.PARS,self.OPTIONS)
        self.Routing = Routing(self.PARS,self.OPTIONS)
        self.Misscell = Misscell(self.PARS,self.OPTIONS)
        print ' '
        print 'Model fluxes are loaded'
        print ' '

        #----------------------------------------------------------------------
        # CHECK CALCULATION TIME
        #----------------------------------------------------------------------
        e0 = time.time() # elapsed time since the epoch
        c0 = time.clock() # total CPU time spent in the script so far
        self.timecounter=1.
        #----------------------------------------------------------------------

        #----------------------------------------------------------------------
        #Create the state and flux arrays to save
        #----------------------------------------------------------------------
#        self.ofile.create_group(run_id)
        if not custom_period:
            self.create_datasets(run_id)
        else:
            raise Exception('this functionality is not yet implemented in the model')
        #----------------------------------------------------------------------

        #----------------------------------------------------------------------
        #initiate timestep-counter
        #self.current_timestep=1
        #----------------------------------------------------------------------

        #----------------------------------------------------------------------
        #Prepare initial conditions AND MAKE RELEVANT FOR  OPTIONS ==> MAPPERFUNCIONS
        #TODO: make it INITFRAC from max values
        ##################################
        #prepare initial values of different variables based on number of variables
        self.INIT=np.ones(self.state_size)*1.
        #specially for NAM test:
        #tmself.INIT=[100., 0.5, 0.0]
        #specially for PDM test:
        #self.INIT=[20., 0.5, 0.0]
        print 'Initial values are loaded'
        print ' '

        #----------------------------------------------------------------------

        #----------------------------------------------------------------------
        #Solve with LSODA scheme -> odeint
        # more freedom in the solver method needed, but for the first version, odeint is ok
        #----------------------------------------------------------------------
        print 'Model simulation started...'
        print ' '

        #odespy solvers
        use_odespy = True
        self.ode_integrator = 'RK4'

        if use_odespy and odespy_import:
            solver = eval("odespy." + self.ode_integrator + "(self.deriv_mod)")
            solver.set_initial_condition(self.INIT)
            solver.set(f_args = (self.rain_int, self.evapo_int,
                                 self.PARS, self.OPTIONS,))
            state_out, t = solver.solve(self.time)
            infodict = 'odespy does not support convergence infodict'

            print 'Using' + solver.quick_description + 'solver'

        else:
            state_out, infodict = odeint(self.deriv_mod, self.INIT, self.time,
                                        full_output=1, printmessg=True,
                                        args=(self.rain_int, self.evapo_int,
                                              self.PARS,self.OPTIONS), hmax=1.0)

     # ODESPY stijl: solver = eval("odespy." + self.ode_integrator + "(self._fun_ODE)")
     # https://github.ugent.be/biomath/biointense/blob/master/biointense/ode_generator.py  lijn 1057
        # todo -> user solver.get() + change options!

        print '...Model simulation ends'
        print ' '

        #print infodict
        #----------------------------------------------------------------------

        #----------------------------------------------------------------------
        #SAVE STATES TO hdf5 folder
        #----------------------------------------------------------------------
        self.statetoh5(run_id,state_out)
        #----------------------------------------------------------------------
        print 'State results saved...'
        print ' '
        #----------------------------------------------------------------------
        #CALCULATE FLUXES for specific timesteps and save to hdf5 folder
        #----------------------------------------------------------------------
        self.fluxtoh5(run_id,state_out)
        print 'Fluxes calculated and saved...'
        print ' '
        #alternatief: alles saven en op basis geg van de infodict de fluxen per tijdstap berekenen
        #----------------------------------------------------------------------

        #----------------------------------------------------------------------
        #PUT PARS IN hdf5 file
        #----------------------------------------------------------------------
        self.partoh5(run_id)
        #Add constant values as attribute
        for cnt in self.CONST:
            self.rungroup.attrs[cnt] = self.CONST[cnt]
        #----------------------------------------------------------------------

        #TODO:  opleten voor de trunc ations en andere files in de originerle FUSE: opnemen
        #TODO: CHECK filesto what end this is needed?: fix_states, limit_xtry
        #TODO: raise Exception for impossible options
        #TODO:  Add attributes to the h5py files

        #----------------------------------------------------------------------
        ##CHECK TIME TO CALCULATE
        elapsed_time = time.time() - e0
        cpu_time = time.clock() - c0
        print 'The elapsed time for the model run was %f seconds' %elapsed_time
        print 'The elapsed cpu-time for the model run was %f seconds' %cpu_time
        #----------------------------------------------------------------------

        print 'Do not forget to use te close_h5 when finished working on the hdf5 file'
        return state_out, infodict

    def run_MC(self,nruns):
        '''
        Do a Monte Carlo run of the specified model, using the given distributions

        Parameters
        ------------
        nruns: int
            number of monte carlo runs to perform of the model

        Returns
        ---------
        output: hdf5 datasets
            An update of the self.ofile is made, with the different runs saved in
            the hdf5 format
        '''
        if self.mcpossible == False:
            raise Exception('No monte carlo possible, since only one value for each parameter')
        for i in range(nruns):
            #sample value
            print 'Simulation ',str(i+1),' starts...'
            parin={}
            for par in self.parameters:
                parin[par]=self.parameters[par].aValue()
            parin = Derived_pars(parin)

            #run model
            self.run(new_pars=parin,run_id='MCrun'+str(i+1))

    def run_EE(self):
        '''
        Explicit Euler for solving the equation and comparison

        !only valid for simple testmodel for comparison
        '''

        state=np.zeros(self.totaltimesteps)
        state[0]=self.INIT

        for tstep in range(1,int(self.totaltimesteps+1)):
            state[tstep]=state[tstep-1]+self.deriv_mod([state[tstep-1]],tstep,self.rain_int,self.evapo_int,self.PARS,self.OPTIONS)
        return state

    def control_balance(self,runid='testrun'):
        dS=testModel.ofile[self.name+'/'+run_id+'/STATE_S'][1:]-testModel.ofile[self.name+'/'+run_id+'/STATE_S'][:][:-1]

        errbal= self.ofile[self.name+'/'+run_id+'/FLUX_rain'][:][:-1] - self.ofile[self.name+'/'+run_id+'/FLUX_q_all'][:][:-1] \
        - self.ofile[self.name+'/'+run_id+'/FLUX_e_all'][:][:-1] - dS

        plt.figure()
        plt.plot(errbal)
        print np.sum(errbal)
        return errbal

    def total_outflow(self,run_id='testrun'):
        '''
        Outflow of the catchment in m3/s when hourly timestep

        Parameters
        ------------
        run_id: str
            name of the run_id to get the flow from

        Returns
        ----------
        totalout: array
            array of the flow output
        '''
        v = np.float64(self.area * 1000.0 / (60.0 * 60.0))
        self.totalout=self.ofile[self.name+'/'+run_id+'/FLUX/FLUX_q_all'][:]*v
        return self.totalout

    def array_output(self,outname,outtype='FLUX',run_id='testrun'):
        '''
        Get output array of the selected outname; fluxes in mm/hour, states in mm and pars in the given units

        Parameters
        ------------
        outtype: str
            one of these values: 'FLUX', 'STATE' or 'PAR'
        outname: str
            the specific flux/stae of par needed
        run_id: str
            name of the run_id to get the flow from

        Returns
        ----------
        output: array
            output of the specified model output
        '''
        try:
            return self.ofile[self.name+runid+'/'+outtype+'/'+outtype+'_'+outname][:]
        except:
            raise Exception('Choose output name from model')

    def close_h5(self):
        '''
        close the hdf5 file to stop the model analysis
        '''
        self.ofile.close()

    def get_const_info(self,run_id='testrun1'):
        '''
        Get overview of the constant values used in the run
        '''
        for name, value in self.ofile[self.name+'/'+run_id].attrs.iteritems():
            print name+":", value

    def get_model_info(self):
        '''
        Get overview of the structure options working with
        '''
        for name, value in self.ofile[self.name].attrs.iteritems():
            print name+":", value


def Load_model(hdffile, parfile):
    '''
    Load an old model run and set up the model

    Parameters
    -----------
    hdffile: HDF5 file
        previous simulation outputs file to re-use

    parfile: parameter textfile
        parameterfile of the model under consideration

    Returns
    --------
    pyFUSE_Model: model instance
        loaded model to do calculations

    See Also
    ----------
    pyFUSE.Set_pars:
        Loading in parameter sets example of input parameter file

    '''
    #GET model name
    if hdffile[-5:]=='.hdf5':
        name=hdffile[:-5]
    else:
        name=hdffile

    #LOAD hdf5 file
    try:
        ff = h5py.File(hdffile, 'r')
    except:
        raise Exception('Wrong hdf5 file')

    modelname = ff.keys()[0]

    #GET model OPTIONS
    options={}
    for name, value in ff[modelname].attrs.iteritems():
        options[name] = value

    #GET model constant values
    cnt={}
    arun=ff[modelname].keys()[0] #assumes all runs with same constant values
    for name, value in ff[modelname+'/'+ arun].attrs.iteritems():
        cnt[name] = value

    #GET model aprameters => use parfile
    #from file

    #Get rain
    rain = ff[modelname+'/'+ arun+'/FLUX/FLUX_rain'][:]

    #Get evapotranspiration
    evapo = ff[modelname+'/'+ arun+'/FLUX/FLUX_pet'][:]
    ff.close()
    Model=pyFUSE_Model(modelname,options,parfile,rain,evapo,cnt,oldfile=True)

    return Model



