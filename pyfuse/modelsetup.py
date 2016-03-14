# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 16:51:28 2012

@author: VHOEYS
"""

import re
import pprint

import numpy as np

from infodicts import *
from parameter import *
from extrafunctions import *

###############################################################################
###   START-UP FUNCTIONS TO SELECT MODELS/PARS FOR CREATING MODEL STRUCTURE
###############################################################################

def Set_options(filename=False,default=False):
    '''
    Read options from file or if no file take standard implemented here

    Parameters
    -----------
    filename: str
        name of the Options containing file, following the given template

    Returns
    ----------
    Dictionary, used to set-up the model structure
    '''

    if not filename==False:
        Options={}
        f=open(filename,'r')
        for line in f.readlines():
            if line[0]=='#':
                continue
            else:
                line.split()
                optioninputs=re.split(r'\s+',line.strip())

                if not optioninputs[0] in all_options():
                    raise Exception('Option name is not matching one of the FUSE model decisions,\
                                    information in all_options()')
                if not optioninputs[2] in all_options()[optioninputs[0]]:
                    raise Exception('Model option value is not matching one of the FUSE option values,\
                                    information in all_options()')
                Options[optioninputs[0]]=optioninputs[2]
        f.close()
        Options['soilstorage']= 'twolayer'
        #to delete, but still old version
        Options['reservoirs']='220' #TODO veranderen!!!

    elif not default==False:
        if default == 'NAM':
            print 'NAM model options are selected'
            Options={}

            #NAM-STYLE
            Options['soilstorage']= 'twolayer'#FIXED OPTION: ALWAYS ON!!
            Options['uplayer']= 'surface1_1'#'#'onestate_1'#'surface1_1'
            Options['lowlayer_baseflow']= 'unlimfrc_2'
            Options['surface']= 'oflwtresh'
            Options['percolation']='perc_tresh'
            Options['evaporation']= 'sequential'
            Options['interflow']='intflwtresh'
            Options['routing']= 'rout_ind' #linear
            Options['reservoirs']='220' #only when routing is selected

        elif default == 'PDM':
            print 'PDM model options are selected'

            Options={}
            Options['soilstorage']= 'twolayer'
            Options['uplayer']= 'onestate_1'#'#'onestate_1'#'surface1_1'
            Options['lowlayer_baseflow']= 'unlimfrc_2'
            Options['surface']= 'arno_x_vic'
            Options['percolation']='perc_w2sat'
            Options['evaporation']= 'sequential'
            Options['interflow']='intflwnone'
            Options['routing']= 'rout_ind' #linear
            Options['reservoirs']='200' #HERWERKING: 0111 (alle deelstromen) of 1000 (volledig op totaal) => wel par geworden!!

        elif default == 'Hymod':
            Options={}
            Options['soilstorage']= 'twolayer'
            Options['uplayer']= 'onestate_1'#'#'onestate_1'#'surface1_1'
            Options['lowlayer_baseflow']= 'unlimfrc_2'
            Options['surface']= 'arno_x_vic'
            Options['percolation']='perc_nodrain'
            Options['evaporation']= 'sequential'
            Options['interflow']='intflwnone'
            Options['routing']= 'rout_ind' #linear
            Options['reservoirs']='300'

        else:
            raise Exception('NAM, PDM and hymod are the default model structures')

    else:
        raise Exception('use default structure name or give options file')

    print 'The selected options are'
    pprint.pprint(Options)
    return Options

def Set_pars(filename):
    '''
    Read info from file and put in parameter dictionary
    and calculate derived parameter values

    Parameters
    ------------
    filename: str
        name of the parameters containing file, following the given template

    Returns
    ----------
    Dictionary, used to parameterize the model structure

    '''
    set_par={}

    #READ FILE and match the pars in the general oveview
    f=open(filename,'r')
    for line in f.readlines():
#        print line
#    for line in fileinput.input([filename]):
        if line[0]=='#':
            continue
        else:
            line.split()
            parinputs=re.split(r'\s+',line.strip())
            if not parinputs[0] in all_pars_info():
                raise Exception('Parameter name is not matching one of the FUSE parameter names,\
                                information in all_pars_info()')
            try:
                parmin=np.float(parinputs[1])
                parmax=np.float(parinputs[2])
                paroptguess=np.float(parinputs[3])

            except:
                raise Exception('Parameter %s in input file not well defined'%parinputs[0])

            if not parinputs[4] in ['randomUniform', 'randomTriangular', 'randomTrapezoidal', 'randomNormal', 'randomLogNormal']:
                raise Exception('selected distribution not available')

            if parinputs[4] in ['randomUniform']:
                set_par[parinputs[0]]=ModPar(parinputs[0],parmin,parmax,paroptguess,parinputs[4])
            elif parinputs[4] in ['randomTriangular']:
                set_par[parinputs[0]]=ModPar(parinputs[0],parmin,parmax,paroptguess,parinputs[4],np.float(parinputs[5]))
            else:
                set_par[parinputs[0]]=ModPar(parinputs[0],parmin,parmax,paroptguess,parinputs[4],np.float(parinputs[5]),np.float(parinputs[6]))

            #            set_par[parinputs[0]]=ModPar(parinputs[0],parmin,parmax,paroptguess,parinputs[4])
    f.close()

    return set_par

def Set_pars_for_run(Parlib):
    '''
    Retrieve optguess values of parameters to calculate a simulation

    Parameters
    ------------
    Parlib: dict
        Dictionary containing the parameter characteristics, with distribution
        and extra information

    Returns
    --------
    pars: dict
        Dictionary with parnames and their repsective values to simulate
    '''
    pars={}
    for par in Parlib.iterkeys():
        pars[par]=Parlib[par].optguess
    Derived_pars(pars)
    return pars

def Derived_pars(set_par_lib):
    '''
    Help function to calculate parameters derived by the given par values
    calculated and added to the library

    When performing Monte Carlo, needs to be recalculated each time
    '''
    set_par_lib['S1Tmax']=set_par_lib['S1max'] * set_par_lib['fitens']  #bucketsize file
    set_par_lib['S2Tmax']=set_par_lib['S2max'] * set_par_lib['fitens']
    set_par_lib['S1Fmax']=set_par_lib['S1max'] * (1.-set_par_lib['fitens'])
    set_par_lib['S2Fmax']=set_par_lib['S2max'] * (1.-set_par_lib['fitens'])
    set_par_lib['S1TAmax']=set_par_lib['S1Tmax'] * set_par_lib['firchr']
    set_par_lib['S1TBmax']=set_par_lib['S1Tmax'] * (1.-set_par_lib['firchr'])
    set_par_lib['S2FAmax']=set_par_lib['S2Fmax'] * set_par_lib['fibase']
    set_par_lib['S2FBmax']=set_par_lib['S2Fmax'] * (1.-set_par_lib['fibase'])
    set_par_lib['r2']=1. - set_par_lib['r1']

#    set_par_lib['Smax']= set_par_lib['S1max']  #Smax gets max from value S1max
    #calculate mean of power transformed topographic index
    set_par_lib = calc_meantipow(set_par_lib)  #adds 2 parameters!
    set_par_lib = qtimedelay(set_par_lib,deltim=1.)
    #TODO: fracfuture for gamma-runoff
    return set_par_lib

def Set_cnts(filename=False):
    '''
    Read info from file and put in constants dictionary

    Parameters
    ------------
    filename: str
        name of the constants containing file, following the given template

    Returns
    ----------
    Dictionary, used to simulate the model

    '''
    if filename:
        const={}
        f=open(filename,'r')
        for line in f.readlines():
            if line[0]=='#':
                continue
            else:
                line.split()
                ctninputs=re.split(r'\s+',line.strip())
                print ctninputs[0]

                if not ctninputs[0] in all_cnt_info():
                    raise Exception('constant name is not matching one of the FUSE model constants,\
                                    information in all_cnt_info()')
                const[ctninputs[0]]=np.float(ctninputs[2])
        f.close()
    else:
        const={}
        const['area']=np.float(362.)
        const['timestep']=np.float(1.0)    #relative to the frequency of the rain input
    return const