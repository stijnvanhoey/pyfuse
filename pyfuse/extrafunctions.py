# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 16:45:07 2012

@author: VHOEYS
"""

import numpy as np
from scipy import special


###############################################################################
###   USEFUL FUNCTIONS  -> need to go to sperate folder later on
###############################################################################

def Logistic(State,Statemax,Psmooth=0.01):
    ''' 
    Logistic function to smooth the threshold at the overflow of a bucket
    
    Parameters
    ------------
    State: float
        Current state value
    Statemax: float
        Maximum allowed value of the state value
    Psmooth: float 
        smoothing parameter to choose amount of smoothing default value 0.01
    
    Returns
    --------  
    logismooth: float
        smoothing multiplier
        
    References
    ------------
    [1] Clark, Martyn P., A. G. Slater, D. E. Rupp, R. A. Woods, Jasper A. 
    Vrugt, H. V. Gupta, Thorsten Wagener, and L. E. Hay. Framework for
    Understanding Structural Errors (FUSE): A modular framework to diagnose 
    differences between hydrological models. Water Resources Research 44 (2008): 14.
    Original code from Clark, Martyn P.
    
    [2] Kavetski, Dmitri, and George Kuczera. “Model Smoothing Strategies to 
    Remove Microscale Discontinuities and Spurious Secondary Optima in 
    Objective Functions in Hydrological Calibration.” Water Resources Research 43, no. 3 (March 8, 2007): 1–9. http://www.agu.org/pubs/crossref/2007/2006WR005195.shtml.
    '''
    epsilon = 5.0       #Multiplier to ensures storagde is always less than capacity
    
    Asmooth = Psmooth * Statemax        #actual smoothing
    logismooth = 1. / ( 1. + np.exp(-(State - (Statemax - Asmooth * epsilon) ) / Asmooth) )
    return logismooth


def calc_meantipow(set_par):
    '''
    Calculate mean of power transformed topographic index
    
    Implementation here is the 3 parameter gamma distribution
    version following [1], with the values for chi, psi and loglambda. 
    In the [2] version, the chi and psi 
    are interchanged and the 1 parameter version is applied
       

    more info: utilities -> gammadist
    
    needs par-library as input!
    
    References
    ------------
    [1] Sivapalan, M., K. Beven, and F Eric Wood. “On Hydrologic Similarity 
    2. A Scaled Model of Storm Runoff Production.” Water Resources Research 
    23, no. 12 (1987): 2266–2278. 
    
    [2] Clark, Martyn P., A. G. Slater, D. E. Rupp, R. A. Woods, Jasper A. 
    Vrugt, H. V. Gupta, Thorsten Wagener, and L. E. Hay. Framework for
    Understanding Structural Errors (FUSE): A modular framework to diagnose 
    differences between hydrological models. Water Resources Research 44 (2008): 14.
    Original code from Clark, Martyn P.
    '''
    Ti_off=3.
    phi= (set_par['loglambda']-Ti_off)/set_par['chi']    #Chi -- loglamb is the first parameter (mean)

    # loop through the frequency distribution
    LOWERV = 0.
    LOWERP = 0.
    AVELOG = 0.  #testing
    AVEPOW = 0.
    
    Nbins=2000
    Ti_max=50.
    
    for ibin in range (1,Nbins):
        # get probability for the current bin
        UPPERV = (float(ibin)/Nbins) * Ti_max               # upper value in frequency bin
        GMARG2 = max(0., UPPERV - Ti_off) / set_par['chi']          # 2nd argument to the Gamma function  Ti_arg = max(0., Ti_log - Ti_off) / chi
        UPPERP = special.gammainc(phi, GMARG2)           # GAMMP is the incomplete Gamma function GAMMP(Ti_shp, GMARG2)
        PROBIN = UPPERP-LOWERP                              # probability of the current bin
        # get the scaled topographic index value
        LOGVAL = 0.5*(LOWERV+UPPERV)                        # log-transformed index for the current bin
        POWVAL = (np.exp(LOGVAL))**(1./set_par['n'])        # power-transformed index for the current bin
        AVELOG = AVELOG + LOGVAL*PROBIN                    # ! average log-transformed index (testing)
        AVEPOW += POWVAL*PROBIN                             # average power-transformed index
#        print LOWERV, UPPERV, LOGVAL, POWVAL, AVEPOW
#        !write(*,'(7(f9.3,1x))') lowerv, upperv, logval, powval, avelog, avepow
        # save the lower value and probability
        LOWERV = UPPERV                                     # lower value for the next bin
        LOWERP = UPPERP                                     # cumulative probability for the next bin

    set_par['maxpow']=POWVAL
    set_par['powlambda']=AVEPOW
    
    return set_par

def qtimedelay(set_par,deltim=1.):
    '''
    Gamma-function based weight function to control the runoff delay, this
    function calculates the fractions to derive the fluxes, defined by the
    frac_future parameter
       
    Parameters
    ------------
    set_par: dict
        The dictionary with the model parameters
    
    Returns
    ---------
    set_par: dict
        updated dictionary with the frac_future parameter added        
    '''
    alpha = 2.5
    alamb = alpha/set_par['mut']
    psave = 0.0
    set_par['frac_future'] = np.zeros(100.)  #Parameter added
    ntdh = set_par['frac_future'].size
    
    deltim = deltim 
#    print 'qtimedelay is calculated with a unit of',deltim,'hours to have parameter values comparable to Clarke, 2008'
    
    for jtim in range(ntdh):
#        print jtim
        tfuture = jtim * deltim
#        print alamb*tfuture
        cumprob = special.gammainc(2.5, alamb*tfuture)# hoeft niet want verschil wordt genomen: /special.gamma(alpha)
#        print cumprob
        set_par['frac_future'][jtim]=max(0.,cumprob-psave)
        psave = cumprob
    
    if cumprob < 0.99:
        print 'not enough bins in the frac_future'

    #make sure sum to one
    set_par['frac_future'][:]=set_par['frac_future'][:]/set_par['frac_future'][:].sum()
    
    return set_par
   

###############################################################################
###   LINEAR RESERVOIRS
###############################################################################

def linres(n_res,q_init,co,k):
    '''
    Recursive calculation of a cascade of linear reservoirs, with fluxes
    defined in mm
    
    Parameters
    -----------
    n_res: int
        number of reservoirs in the cascade
    q_init: narray
        initial fluxes of the different reservoirs in array
    co: float
        incoming flow
    k: float
        residence time constant for the reservoir
    
    Returns
    ---------
    q_init: narray
        narray with the new fluxes ofr each reservoir

    TODO: improve by incorporating previous timestep incoming flow        
    '''
    if n_res==1:
        q_init[0]=q_init[0]*np.exp(-1/k) + co*(1 - np.exp(-1/k))
        return q_init
    else:
        q_init[n_res-1]=q_init[n_res-1]* np.exp(-1/k)+linres(n_res-1,q_init,co,k)[n_res-2]*(1 - np.exp(-1/k))
        return q_init

def linresv(n_res,q_init,co,v,k):
    '''
    Recursive calculation of a cascade of linear reservoirs, with incoming
    fluxes defined in mm and outgoing flow in m3/s
    
    See Also
    ---------
    pyFUSE.linres
    '''    
    if n_res==1:
        q_init[0]=q_init[0]*np.exp(-1/k) + co*(1 - np.exp(-1/k))*v
        return q_init
    else:
        q_init[n_res-1]=q_init[n_res-1]* np.exp(-1/k)+linres(n_res-1,q_init,co,k)[n_res-2]*(1 - np.exp(-1/k))*v
        return q_init 