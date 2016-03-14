# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 16:48:33 2012

@author: VHOEYS
"""

import pprint
import numpy as np

##############################################################################
###   DICTIONARY FOR INFORMATION ABOUT MODEL FRAMEWORK OF ALL VARS,PARS,...
###############################################################################

def all_options():
    """
    Returns a print out of all possible fluxes of all model combinations possible to select
    
    (1) upper-layer architecture
    
        tension1_1:
            upper layer broken up into tension and free storage
        tension2_1:
            tension storage sub-divided into recharge and excess
        onestate_1:
            upper layer defined by a single state variable
        surface1_1:
            upper layer defined by a surface storage representing and a tension reservoir
        

    (2) lower-layer architecture and baseflow
    
        tens2pll_2:
            tension reservoir plus two parallel tanks
        unlimfrc_2:
            baseflow resvr of unlimited size (0-HUGE), frac rate
        unlimpow_2:
            baseflow resvr of unlimited size (0-HUGE), power recession
        fixedsiz_2:
            baseflow reservoir of fixed size

    (3) surface runoff
    
        arno_x_vic:
            ARNO/Xzang/VIC parameterization (upper zone control)
        prms_varnt:
            PRMS variant (fraction of upper tension storage)
        tmdl_param:
            TOPMODEL parameterization (only valid for TOPMODEL qb)
        oflwtresh:
            threshold based overland flow generation

    (4) percolation
        
        perc_f2sat:
            water from (field cap to sat) avail for percolation
        perc_w2sat:
            water from (wilt pt to sat) avail for percolation
        perc_lower:
            perc defined by moisture content in lower layer (SAC)
        perc_tresh:
            threshold based percolation
        perc_nodrain:
            percolation represents the baseflow routing

    (5) evaporation
    
        sequential:
            sequential evaporation model
        rootweight:
            root weighting

    (6) interflow
    
        intflwnone:
            no interflow
        intflwsome:
            linear interflow
        intflwtresh:
            threshold based interflow generation
    
    (7) time delay in runoff
    
        rout_all1:
            use a Gamma distribution with shape parameter = 2.5 to rout combined flow
        no_rout:
            no routing
        rout_ind:
            rout subflows individual    
    """
    All_options={}
#    All_options['soilstorage'] = 'onelayer', 'twolayer' #AFSCHAFFEN DIE HANDEL => altijd 2 layers
    #idea one layer is an upper layer:
    All_options['uplayer']= 'onestate_1', 'tension1_1','tension2_1','easytest','surface1_1'
    All_options['lowlayer_baseflow']= 'tens2pll_2','unlimfrc_2','unlimpow_2','fixedsiz_2'
    All_options['surface']= 'arno_x_vic','prms_varnt','tmdl_param','testeasy','oflwtresh' 
    All_options['percolation']='perc_f2sat','perc_w2sat','perc_lower','perc_nodrain','perc_tresh'
    All_options['evaporation']= 'sequential','rootweight'#, 'freefirst'  
    All_options['interflow']='intflwnone','intflwsome','intflwtresh'    
    All_options['routing']= 'rout_all1','no_rout','rout_ind'
#    All_options['reservoirs']='Any combination of tree integers (<10) with respectively overland, interflow and base flow: eg 122'
         
#    pprint.pprint(All_options)
    return All_options

def all_states():
    """
    Returns a print of all possible states of all model combinations
    """
    
    all_states_info={}
    all_states_info['S1']='Total water content in the upper soil layer'
    all_states_info['S1T']='Tension water content in the upper soil layer'
    all_states_info['S1TA']='Primary tension water content in the upper soil layer'
    all_states_info['S1TB']='Secondary tension water content in the upper soil layer'
    all_states_info['S1F']='Free water content in the upper soil layer'
    
    all_states_info['S2']='Total water content in the lower soil layer'
    all_states_info['S2T']='Tension water content in the lower soil layer'
    all_states_info['S2FA']='Free water content in the primary base flow reservoir'
    all_states_info['S2FB']='Free water content in the secondary base flow reservoir'
    
    all_states_info['S']='Total water content in the single soil layer'        
         
    return all_states_info

def all_fluxes():
    """
    Returns a print of  all possible fluxes of all model combinations possible to select
    """
    all_fluxes_info={}
    all_fluxes_info['rain']='Precipitation'
    all_fluxes_info['pet']='Potential evaporation'
    all_fluxes_info['e']='otal evaporation from both soil layers'        
    all_fluxes_info['e1']='Evaporation from the upper soil'
    all_fluxes_info['e2']='Evaporation from the lower soil'   
    all_fluxes_info['e1A']='Evaporation from the primary tension store'
    all_fluxes_info['e1B']='Evaporation from the secondary tension store'
    all_fluxes_info['qsx']='Surface runoff'  
    all_fluxes_info['q12']='Percolation of water from the upper to the lower layer'       #water from upper to lower storage
    all_fluxes_info['qif']=  'Interflow'
    all_fluxes_info['qb']= 'Base flow'
    all_fluxes_info['qbA']= 'Base flow from the primary reservoir'
    all_fluxes_info['qbB']= 'Base flow from the secondary reservoir'
    all_fluxes_info['qurof']='Overflow of water from the primary tension in the upper soil layer'
    all_fluxes_info['qutof']='Overflow of water from tension storage in the upper soil layer'
    all_fluxes_info['qufof']= 'Overflow of water from free storage in the upper soil layer'
    all_fluxes_info['qstof']= 'Overflow of water from tension storage in the lower soil layer'
    all_fluxes_info['qsfof']='Overflow of water from free storage in the lower soil layer'
    all_fluxes_info['qsfofA']='Overflow of water from primary base flow storage in the lower soil layer'
    all_fluxes_info['qsfofB']='Overflow of water from secondary base flow storage in the lower soil layer'
    all_fluxes_info['routov']='Overland flow after routing'
    all_fluxes_info['routif']='Interflow flow after routing'
    all_fluxes_info['routb']='Baseflow flow after routing'    

    return all_fluxes_info  

def all_pars_info():
    """
    Returns a print of all possible parameters of all model combinations 
    """
    all_pars_info={}
    all_pars_info['S1max']='Maximum storage in the upper layer'
    all_pars_info['S2max']='Maximum storage in the lower layer '
    all_pars_info['fitens']='Fraction total storage as tension storage'
    all_pars_info['firchr']='Fraction of tension storage in primary zone (upper layer) '
    all_pars_info['fibase']='Fraction of free storage in primary reservoir (lower layer) '
    all_pars_info['r1']='Fraction of roots in the upper layer'
    all_pars_info['ku']='Percolation rate'
    all_pars_info['c']='Percolation exponent in SAC'
    all_pars_info['alfa']='Percolation multiplier for the lower layer cfr SAC'
    all_pars_info['psi']='Percolation exponent for the lower layer'
    all_pars_info['kappa']='Fraction of percolation to tension storage in the lower layer '
    all_pars_info['ki']='Interflow rate'
    all_pars_info['ks']='Base flow rate'
    all_pars_info['n']='Base flow exponent'
    all_pars_info['v']='Base flow depletion rate for single reservoir'
    all_pars_info['vA']='Base flow depletion rate for primary reservoir'
    all_pars_info['vB']='Base flow depletion rate for secondary reservoir'    
    all_pars_info['Acmax']='Maximum saturated area (fraction)'    
    all_pars_info['b']= 'ARNO/VIC b exponent'
    all_pars_info['loglambda']='Mean of the LOG-transformed topographic index distribution '    
    all_pars_info['chi']='Shape parameter defining the topographic index distribution'    
    all_pars_info['mut']='Time delay in runoff'
    all_pars_info['alfah']='Splitting parameter when modelling hymod model option'
    all_pars_info['tg']='Threshold value for groundwater recharge'
    all_pars_info['tif']= 'interflow treshold function'
    all_pars_info['tof']= 'overland flow treshold function'
    all_pars_info['ko']= 'overland flow runoff coefficient'   
    all_pars_info['timeo']= 'overland flow time constant'    
    all_pars_info['timei']= 'interflow flow time constant'
    all_pars_info['timeb']= 'base flow time constant'
    all_pars_info['nreso']= 'overland flow number of reservoirs'    
    all_pars_info['nresi']= 'interflow flow number of reservoirs'
    all_pars_info['nresb']= 'base flow number of reservoirs'
        
        
    all_pars_info['S1Tmax']='Maximum tension storage in the upper layer (derived)'
    all_pars_info['S2Tmax']='Maximum tension storage in the lower layer (derived)'
    all_pars_info['S1Fmax']='Maximum free storage in the upper layer (derived)'
    all_pars_info['S2Fmax']='Maximum free storage in the lower layer (derived)'
    all_pars_info['S1TAmax']='Maximum storage in the primary tension reservoir (derived)'
    all_pars_info['S1TBmax']='Maximum storage in the secondary tension reservoir (derived)'
    all_pars_info['S2FAmax']='Maximum storage in the primary base flow reservoir (derived)'     
    all_pars_info['S2FBmax']='Maximum storage in the secondary base flow reservoir (derived)'         
    all_pars_info['r2']= 'Root fraction in the lower soil layer (derived)'
    all_pars_info['powlambda']='Mean of the POWER-transformed topographic index (derived)' 
    all_pars_info['maxpow']='max value of power-transformed topographic index, m**(1/n)'
    all_pars_info['frac_future']='Fraction of runoff for future timesteps (array!)'

    #! when sampling    
    all_pars_info['be']='Evaporation exponent to differentiate between linear and non-linear dependence'
    return all_pars_info

def all_cnt_info():
    """
    Returns a print of all constants needed
    """
    all_cnt_info={}
    all_cnt_info['area']='Catchment size in square km'
    all_cnt_info['timestep']='Time step for the model relative to an hourly timestep' 
    return all_cnt_info
