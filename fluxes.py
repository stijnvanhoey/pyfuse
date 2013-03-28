# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 16:49:56 2012

@author: VHOEYS
"""

from pyFUSE.extrafunctions import *

###############################################################################
###   FLUXES OF THE MODEL
###############################################################################

class Flux_base():
    """
    Base class for the fluxes to be used by all fluxes used in the model representation
    
    Basically a flux is a calculated process giving an output based on the states on a particular moment
    This base-class is used add specific methods for extraction of information (todo)
    
    Since for one model run, parameters and options are constant, 
    when setting-up the model, all fluxes are loaded and pars dictionary is directed
    
    Information when adding a new flux option: 
        - need to add elif in the calc of the flux
        - create new method (+ adding in general list!)
    
    Parameters
    ----------
    PARS : dictionary
        Dictionary containing all the model parameters
    OPTIONS: dictionary 
        Dictionary containting all the model structure options
    """
    
    def __init__(self,PARS,OPTIONS):
        if isinstance(PARS,dict):
            self.PARS=PARS
            #            print 'parameters to calculate model fluxes loaded'
            #todo: control here only 1 float value for each parameter!
        else:
            raise Exception('Pars to use in fluxes must be given in dictionary format')

        if isinstance(OPTIONS,dict):
            self.OPTIONS=OPTIONS
            #            print 'parameters to calculate model fluxes loaded'
        else:
            raise Exception('OPTIONS to use in fluxes must be given in dictionary format')            

class Evaporation(Flux_base):
    """
    Class for Evapotranspiration fluxes calculation, with 2 implemented options:
        1. sequential: upper layer is first evaporating and the remaining is extracted from the lower layer
        2. root-weight: the amount of roots (parameter r1) defines the relative evaporation

    Information when adding a new flux option: 
        * need to add elif in the calc 
        * create new method (+ adding in general list!)
        
    Parameters
    ----------
    PARS : dictionary
        Dictionary containing all the model parameters
    OPTIONS: dictionary 
        Dictionary containting all the model structure options  

    Methods
    ----------
    calc(STATES,FLUXES)
        The calc method picks the correct definition based on the OPTIONS dictionary defined by the mode structure        
    sequential(STATES,FLUXES) 
        Sequential calculation of the evapotranspiration. First the upper layer is evaporating,
        wherafter the lower layer is evaporating. In the case of two tension reservoirs in the 
        upper layer, e1 and e2 are extracted from those two storages. 
        For the surface option, e1 is evaporation from the surface storage, e2, the upper tension storage    
        
        Model States: dependent of the upper and lower configuration, one or more of:
            S1TA, S1TB, S1T, S2T, S1F, S2T        
        Model Parameters:            
            be, S1Tmax, S2Tmax
        Fluxes returned:
            e1, e2, e_all
            
    root-weight(STATES,FLUXES)
        root-weight calculation of the evapotranspiration. First the upper layer is evaporating
        based on the r1 rootweight and the lower layer is evaporating based on (1-r1).
        In the case of two tension reservoirs in the upper layer, 
        e1 and e2 are extracted from those two storages. 
        For the surface option, r1 parameter is related to the surface storage
        and (1-r1) to the upper tension storage               
        
        Model States: dependent of the upper and lower configuration, one or more of:
            S1TA, S1TB, S1T, S2T, S1F, S2T        
        Model Parameters:            
            be, S1Tmax, S2Tmax, r1 (r2 = 1-r1)
        Fluxes returned:
            e1, e2, e_all            

    Attributes
    ------------
    STATES: dict
            Dictionary containing all the model states of the current calculation time step          
    FLUXES: dict
            Dictionary containting all the fluxes structure options
        
    """

    def __init__(self,PARS,OPTIONS):
        Flux_base.__init__(self,PARS,OPTIONS)
#        print 'ET flux loaded'
        print self.OPTIONS
            
    def calc(self,STATES,FLUXES):
        if self.OPTIONS['evaporation'] == 'sequential':
            return self.sequential(STATES,FLUXES)

        elif self.OPTIONS['evaporation'] == 'rootweight':
            return self.rootweight(STATES,FLUXES)
            
        else:
            raise Exception('No valid option chosen for evaporation')            

    def sequential(self,STATES,FLUXES):
        if self.OPTIONS['soilstorage'] == 'twolayer':
            if self.OPTIONS['uplayer'] == 'tension2_1':
                '''
                M_FLUX%EVAP_1A = MFORCE%PET * TSTATE%TENS_1A/DPARAM%MAXTENS_1A
                M_FLUX%EVAP_1B = (MFORCE%PET - M_FLUX%EVAP_1A) * TSTATE%TENS_1B/DPARAM%MAXTENS_1B
                M_FLUX%EVAP_1  = M_FLUX%EVAP_1A + M_FLUX%EVAP_1B
                '''  
#                e1A = FLUXES['pet'] * (STATES['S1TA']/self.PARS['S1TAmax'])**self.PARS['be']
                e1A = FLUXES['pet'] * (1.-((self.PARS['S1TAmax']-STATES['S1TA'])/self.PARS['S1TAmax'])**self.PARS['be'])
#                e1B = (FLUXES['pet'] - e1A)  * (STATES['S1TB']/self.PARS['S1TBmax'])**self.PARS['be']
                e1B = (FLUXES['pet'] - e1A) * (1.-((self.PARS['S1TBmax']-STATES['S1TB'])/self.PARS['S1TBmax'])**self.PARS['be'])
                e1 = e1A + e1B
                FLUXES['e1A']=e1A
                FLUXES['e1B']=e1B
                FLUXES['e1']=e1
                FLUXES['e_all']=e1
            
            elif  self.OPTIONS['uplayer'] == 'tension1_1' or \
            self.OPTIONS['uplayer'] == 'onestate_1':
                '''
                M_FLUX%EVAP_1A = 0._sp
                M_FLUX%EVAP_1B = 0._sp
                M_FLUX%EVAP_1  = MFORCE%PET * TSTATE%TENS_1/DPARAM%MAXTENS_1
                '''
#                e1 = FLUXES['pet'] * (STATES['S1T']/self.PARS['S1Tmax'])**self.PARS['be']
                e1 = FLUXES['pet'] * (1.-((self.PARS['S1Tmax']-STATES['S1T'])/self.PARS['S1Tmax'])**self.PARS['be'])                
                FLUXES['e1']=e1
                FLUXES['e_all']=e1
            
            elif self.OPTIONS['uplayer'] == 'surface1_1':
                '''
                surface style
                '''
                Psmooth = 0.01
                e1A = FLUXES['pet']*(1.-np.exp(-(STATES['S1F'])/Psmooth)) #\cite{Kavetski2007}
#                e1B = (FLUXES['pet']-e1A) * (STATES['S1T']/self.PARS['S1Tmax'])**self.PARS['be']
                e1B = (FLUXES['pet']-e1A) * (1.-((self.PARS['S1Tmax']-STATES['S1T'])/self.PARS['S1Tmax'])**self.PARS['be'])                
                FLUXES['e1A']=e1A
                FLUXES['e1B']=e1B 
                FLUXES['e1']=e1A+e1B                
                FLUXES['e_all']=e1A+e1B

            else: 
                raise Exception('Two layer system must have firts layer option of tension2_1, tension1_1 or onestate_1')
                
        #TWO LAYERS ->lower layer                
        if self.OPTIONS['soilstorage'] == 'twolayer':
            #CHECK LOWER LAYER CHOICE                
            if self.OPTIONS['lowlayer_baseflow'] == 'tens2pll_2' or self.OPTIONS['lowlayer_baseflow'] =='fixedsiz_2':
                #CHECK THE UPPER LAYER CHOICE
                if self.OPTIONS['uplayer'] == 'tension1_1'or self.OPTIONS['uplayer'] == 'onestate_1':
                    '''
                    M_FLUX%EVAP_2 = (MFORCE%PET-M_FLUX%EVAP_1) * (TSTATE%TENS_2/DPARAM%MAXTENS_2)
                    '''
#                    e2 = (FLUXES['pet'] - FLUXES['e1']) * (STATES['S2T']/self.PARS['S2Tmax'])**self.PARS['be']
                    e2 = (FLUXES['pet'] - FLUXES['e1']) * (1.-((self.PARS['S2Tmax']-STATES['S2T'])/self.PARS['S2Tmax'])**self.PARS['be'])                
                    FLUXES['e2']=e2
                    FLUXES['e_all']=FLUXES['e_all']+e2
                    
                elif self.OPTIONS['uplayer'] == 'tension2_1' or self.OPTIONS['uplayer'] == 'surface1_1':
                    e2 = 0.0
                    FLUXES['e2']=e2
                else:
                    raise Exception('must be tension1_1, onestate_1 or tension2_1')
                    
            elif self.OPTIONS['lowlayer_baseflow'] == 'unlimfrc_2' or self.OPTIONS['lowlayer_baseflow'] == 'unlimpow_2':
                e2 = 0.0
                FLUXES['e2']=e2
            else: 
                raise Exception('Lower layer must be tens2pll_2,fixedsiz_2, unlimfrc_2,unlimpow_2,topmdexp_2')   
        return FLUXES                
                
    def rootweight(self,STATES,FLUXES):
        #TWO LAYER ->upper layer
        if self.OPTIONS['soilstorage'] == 'twolayer':
            if self.OPTIONS['uplayer'] == 'tension2_1':
                '''
                M_FLUX%EVAP_1A = MFORCE%PET * MPARAM%RTFRAC1 * TSTATE%TENS_1A/DPARAM%MAXTENS_1A
                M_FLUX%EVAP_1B = MFORCE%PET * DPARAM%RTFRAC2 * TSTATE%TENS_1B/DPARAM%MAXTENS_1B
                M_FLUX%EVAP_1  = M_FLUX%EVAP_1A + M_FLUX%EVAP_1B
                '''
                e1A = FLUXES['pet'] * self.PARS['r1'] * (STATES['S1TA']/self.PARS['S1TAmax'])**self.PARS['be']
                e1B = FLUXES['pet'] * self.PARS['r2'] * (STATES['S1TB']/self.PARS['S1TBmax'])**self.PARS['be']
                e1 = e1A + e1B
                FLUXES['e1']=e1
                FLUXES['e_all']=e1
                
            elif  self.OPTIONS['uplayer'] == 'tension1_1' or \
            self.OPTIONS['uplayer'] == 'onestate_1':
                '''
                M_FLUX%EVAP_1A = 0._sp
                M_FLUX%EVAP_1B = 0._sp
                M_FLUX%EVAP_1  = MFORCE%PET * MPARAM%RTFRAC1 * TSTATE%TENS_1/DPARAM%MAXTENS_1
                '''
                e1 = FLUXES['pet'] * self.PARS['r1'] * (STATES['S1T']/self.PARS['S1Tmax'])**self.PARS['be']
                FLUXES['e1']=e1
                FLUXES['e_all']=e1
                
            elif self.OPTIONS['uplayer'] == 'surface1_1':
                '''
                surface style
                '''
                Psmooth = 0.01
                e1A = FLUXES['pet']* self.PARS['r1'] *(1.-np.exp(-(STATES['S1F'])/Psmooth)) #\cite{Kavetski2007}
                e1B =  FLUXES['pet']* self.PARS['r2'] * (STATES['S1T']/self.PARS['S1Tmax'])**self.PARS['be']
                FLUXES['e1A']=e1A
                FLUXES['e1B']=e1B 
                FLUXES['e1']=e1A+e1B                
                FLUXES['e_all']=e1A+e1B               
                
            else: 
                raise Exception('Two layer system must have firts layer option of tension2_1, tension1_1 or onestate_1')                
        
        #TWO LAYER ->lower layer            
        if self.OPTIONS['soilstorage'] == 'twolayer':
            #CHECK LOWER LAYER CHOICE                
            if self.OPTIONS['lowlayer_baseflow'] == 'tens2pll_2' or self.OPTIONS['lowlayer_baseflow'] =='fixedsiz_2':
                #CHECK THE UPPER LAYER CHOICE
                if self.OPTIONS['uplayer'] == 'tension1_1'or self.OPTIONS['uplayer'] == 'onestate_1':
                    '''
                    M_FLUX%EVAP_2 = MFORCE%PET * DPARAM%RTFRAC2 * (TSTATE%TENS_2/DPARAM%MAXTENS_2)
                    '''
                    e2 = FLUXES['pet'] * self.PARS['r2'] * (STATES['S2T']/self.PARS['S2Tmax'])**self.PARS['be']
                    FLUXES['e2']=e2
                    FLUXES['e_all']=FLUXES['e_all']+e2                    
#                    return FLUXES
                    
                elif self.OPTIONS['uplayer'] == 'tension2_1' or self.OPTIONS['uplayer'] =='surface1_1':
                    e2 = 0.0
                    FLUXES['e2']=e2
#                    return FLUXES
                else:
                    raise Exception('must be tension1_1, onestate_1 or tension2_1')
                    
            elif self.OPTIONS['lowlayer_baseflow'] == 'unlimfrc_2' or self.OPTIONS['lowlayer_baseflow'] == 'unlimpow_2' or self.OPTIONS['lowlayer_baseflow'] == 'topmdexp_2':
                e2 = 0.0
                FLUXES['e2']=e2
                
            else: 
                raise Exception('Lower layer must be tens2pll_2,fixedsiz_2, unlimfrc_2,unlimpow_2,topmdexp_2')   
        return FLUXES
            
#    def freefirst(self,STATES,FLUXES): 
#        '''
#        NAM style linear evapotranspiration with first evaporating free storage
#        water evaporates only when upper is one big storage or got at least a free storage
#        similar for one and twolayer configuration        
#        
#        FLUXES is dictionary wth all relevant fluxes of the model
#        
#        VHM here depreciated: uevap-concept is translated into tension-system of others!
#        '''    
#        #water extraction from the upper layer (smoothed)
#
#        #TWO LAYER ->upper layer
#        if self.OPTIONS['soilstorage'] == 'twolayer':
#            if self.OPTIONS['uplayer'] == 'tension2_1' or \
#            self.OPTIONS['uplayer'] == 'tension1_1':
#                '''
#                M_FLUX%EVAP_1A = 0._sp
#                M_FLUX%EVAP_1B = 0._sp
#                M_FLUX%EVAP_1  = MFORCE%PET * MPARAM%RTFRAC1 * TSTATE%TENS_1/DPARAM%MAXTENS_1
#                '''
#                Psmooth = 0.01
#                e1 = FLUXES['pet']*(1.-np.exp(-(STATES['S1F'])/Psmooth))                
#                FLUXES['e1']=e1
#                FLUXES['e_all']=e1
#
#            elif self.OPTIONS['uplayer'] == 'onestate_1':
#                '''
#                M_FLUX%EVAP_1A = 0._sp
#                M_FLUX%EVAP_1B = 0._sp
#                M_FLUX%EVAP_1  = MFORCE%PET * MPARAM%RTFRAC1 * TSTATE%TENS_1/DPARAM%MAXTENS_1
#                '''
#                Psmooth = 0.01
#                e1 = FLUXES['pet']*(1.-np.exp(-(STATES['S1'])/Psmooth))                
#                FLUXES['e1']=e1
#                FLUXES['e_all']=e1                
#                
#            else: 
#                raise Exception('Two layer system must have firts layer option of tension2_1, tension1_1 or onestate_1')                
#        
#        #TWO LAYER ->lower layer            
#        if self.OPTIONS['soilstorage'] == 'twolayer':
#            #CHECK LOWER LAYER CHOICE                
#            if self.OPTIONS['lowlayer_baseflow'] == 'tens2pll_2' or self.OPTIONS['lowlayer_baseflow'] =='fixedsiz_2':
#                #CHECK THE UPPER LAYER CHOICE
#                if self.OPTIONS['uplayer'] == 'tension1_1'or self.OPTIONS['uplayer'] == 'onestate_1':
#                    '''
#                    M_FLUX%EVAP_2 = MFORCE%PET * DPARAM%RTFRAC2 * (TSTATE%TENS_2/DPARAM%MAXTENS_2)
#                    '''
#                    e2 = (FLUXES['pet']-e1) * (STATES['S2T']/self.PARS['S2Tmax'])**self.PARS['be']
#                    FLUXES['e2']=e2
#                    FLUXES['e_all']=FLUXES['e_all']+e2                    
##                    return FLUXES
#                    
#                elif self.OPTIONS['uplayer'] == 'tension2_1':
#                    e2 = 0.0
#                    FLUXES['e2']=e2
##                    return FLUXES
#                else:
#                    raise Exception('must be tension1_1, onestate_1 or tension2_1')
#                    
#            elif self.OPTIONS['lowlayer_baseflow'] == 'unlimfrc_2' or self.OPTIONS['lowlayer_baseflow'] == 'unlimpow_2' or self.OPTIONS['lowlayer_baseflow'] == 'topmdexp_2':
#                e2 = 0.0
#                FLUXES['e2']=e2
#                
#            else: 
#                raise Exception('Lower layer must be tens2pll_2,fixedsiz_2, unlimfrc_2,unlimpow_2,topmdexp_2')   
#
#        return FLUXES         


class Percolation(Flux_base):
    """
    Class for Percolation fluxes calculation, with 5 implemented options:
        1. perc_w2sat: percolation from the entire upper layer storage
        2. perc_f2sat: percolation from the free upper layer storage
        3. perc_lower: percolation based on the lower layer storage
        4. perc_nodrain: no percolation, but infiltration (direct splitting overland; Hymod)
        5. perc_tresh:  percolation based on both layers
        
    Parameters
    -----------
    PARS : dictionary
        Dictionary containing all the model parameters
    OPTIONS: dictionary 
        Dictionary containting all the model structure options  

    Methods
    ---------- 
    calc(STATES,FLUXES)
        The calc method picks the correct definition based on the OPTIONS dictionary defined by the mode structure            
    
    perc_tresh(STATES,FLUXES)
        Treshold calculation of the percolation. When using surface storage,
        percolation is linear related to the surface storage and the upper layer
        tension storage, in all other cases the tension storage of the lower layer
        influences the amount of percolation
        
        Model States:
            S1F, S2T (or S1T when surface1_1 option in upper layer) 

        Model Parameters:
            tg, S1Fmax, S2Tmax (or S1Tmax when surface1_1 option in upper layer) 
        
        Fluxes updated:  
            q12 

    perc_nodrain(STATES,FLUXES)
        Hymod approach, with separating of the excess runoff in baseflow and 
        surface component. The baseflow component is here identified by the 
        percolation and the lower layer is conceptualized as baseflow routing.

        Model States:
            none, since purely dependent from qsx-FLUX

        Model Parameters:
            alfah
        
        Fluxes updated:
            q12
        
    perc_w2sat(STATES,FLUXES)
        Percolation in function of the upper layer soil storage, with in general large 
        values for parameter c to limit drainage below field capacity.

        Model States:
            S1

        Model Parameters:
            S1max, c
        
        Fluxes updated:
            q12

    perc_f2sat(STATES,FLUXES)
        Percolation in function of the free upper layer soil storage, with in general  
        values for parameter c close to unity. 

        Model States:
            S1F

        Model Parameters:
            ku,S1Fmax,c
        
        Fluxes updated:
            q12
           
    perc_lower(STATES,FLUXES)
        Percolation in function of the free upper layer soil storage, with in general  
        values for parameter c close to unity.        

        Model States:
            S1F, S2

        Model Parameters:
            S1Fmax, alfa, S2max, psi
            qbsat is calculated in function of the selected lower layer option
        
        Fluxes updated:
            q12

    Attributes
    ------------
    STATES: dict
            Dictionary containing all the model states of the current calculation time step          
    FLUXES: dict
            Dictionary containting all the fluxes structure options            
    """    
    def __init__(self,PARS,OPTIONS):
        Flux_base.__init__(self,PARS,OPTIONS)
#        print 'Percolation flux loaded'

    def calc(self,STATES,FLUXES):
        if self.OPTIONS['percolation'] == 'perc_w2sat':
            return self.perc_w2sat(STATES,FLUXES)

        elif self.OPTIONS['percolation'] == 'perc_f2sat':
            return self.perc_f2sat(STATES,FLUXES)

#        elif self.OPTIONS['percolation'] == 'perc_pdm':
#            return self.perc_pdm(STATES,FLUXES)

        elif self.OPTIONS['percolation'] == 'perc_lower':
            return self.perc_lower(STATES,FLUXES)  
            
        elif self.OPTIONS['percolation'] == 'perc_nodrain':
            return self.perc_nodrain(STATES,FLUXES)
        
        elif self.OPTIONS['percolation'] == 'perc_tresh':
            return self.perc_tresh(STATES,FLUXES)            
        
        else:
            raise Exception('No valid option chosen for percolation')

    def perc_tresh(self,STATES,FLUXES):           
       
        if self.OPTIONS['uplayer'] == 'surface1_1':     
            wf=Logistic(STATES['S1F'],self.PARS['S1Fmax'],Psmooth=0.01)
            wfg=Logistic(STATES['S1T']/self.PARS['S1Tmax'],self.PARS['tg'],Psmooth=0.01) 
            q12 = ((STATES['S1T']/self.PARS['S1Tmax'] - self.PARS['tg'])/(1.-self.PARS['tg']))*STATES['S1F']*wf*wfg
            FLUXES['q12']=q12
        else:
            wf=Logistic(STATES['S1F'],self.PARS['S1Fmax'],Psmooth=0.01)
            wfg=Logistic(STATES['S2T']/self.PARS['S2Tmax'],self.PARS['tg'],Psmooth=0.01) 
            q12 = ((STATES['S2T']/self.PARS['S2Tmax'] - self.PARS['tg'])/(1.-self.PARS['tg']))*STATES['S1F']*wf*wfg
            FLUXES['q12']=q12
            
        return FLUXES        
   
    def perc_nodrain(self,STATES,FLUXES):                  
        q12 = (1.-self.PARS['alfah'])*FLUXES['qsx']
        FLUXES['q12']=q12
        return FLUXES  
    
    def perc_w2sat(self,STATES,FLUXES): #cfr. gravity drainage term in Richard's equation- c large to prevent drainage below field capacity ! water from (wilt pt to sat) avail for percolation
        q12 = self.PARS['ku'] *(STATES['S1']/self.PARS['S1max'])**self.PARS['c']
        FLUXES['q12']=q12
        return FLUXES

    def perc_f2sat(self,STATES,FLUXES):  #c close to unity - ! water from (field cap to sat) avail for percolation
        q12 = self.PARS['ku'] *(STATES['S1F']/self.PARS['S1Fmax'])**self.PARS['c']
        FLUXES['q12']=q12
        return FLUXES

#    def perc_pdm(self,STATES,FLUXES):  #c close to unity - ! water from (field cap to sat) avail for percolation
#        '''
#        '''
#        Stau=0.0
#        wf = Logistic(STATES['S1'],Stau,Psmooth=0.01)
#        q12 = self.PARS['ku']* wf *(STATES['S1']-Stau)**self.PARS['c']
#        FLUXES['q12']=q12
#        return FLUXES        

    def perc_lower(self,STATES,FLUXES): #! perc defined by moisture content in lower layer (SAC)      
        d_lz = 1.0 + self.PARS['alfa']*(STATES['S2']/self.PARS['S2max'])**self.PARS['psi']
#        print self.calc_Qbsat()
#        print d_lz,'dlz'
        q12 = self.calc_Qbsat() * d_lz*(STATES['S1F']/self.PARS['S1Fmax'])
        FLUXES['q12']=q12
        return FLUXES
    
    def calc_Qbsat(self):
        '''
        Help function of the percolation class,
        calculates the maximum baseflow in function of the lower layer options

        Parameters
        ----------
        
        Returns
        ----------
        qbsat (q0 in model description)        
        '''
        
        if self.OPTIONS['lowlayer_baseflow'] == 'tens2pll_2':
            qbsat = self.PARS['vA']* self.PARS['S2FAmax'] + self.PARS['vB'] * self.PARS['S2FBmax']
    
        elif self.OPTIONS['lowlayer_baseflow'] == 'unlimfrc_2':
            qbsat = self.PARS['v']* self.PARS['S2max']
    
        elif self.OPTIONS['lowlayer_baseflow'] == 'unlimpow_2':
            #      ! This is a bit tricky.  The capacity of the aquifer is m*n, where m is a scaling
            #      ! parameter.  We have the capacity, i.e., MPARAM%MAXWATR_2/1000., and need the
            #      ! TOPMODEL "m" parameter
            TOPmdm = (self.PARS['S2max']/1000.)/self.PARS['n'] # NOTE: mm --> m
            #      ! ...and, compute baseflow
            qbsat = self.PARS['ks'] * (TOPmdm/(self.PARS['powlambda']**self.PARS['n']))
            
#        elif self.OPTIONS['lowlayer_baseflow'] == 'topmdexp_2':
#            #          ! for simplicity we use the CAPACITY as the TOPMODEL scaling parameter
#            TOPmdm = self.PARS['S2max']/1000.                   # NOTE: mm --> m
#            #          ! ..., and compute baseflow
#            qbsat = self.PARS['ks'] * TOPmdm * np.exp(self.PARS['loglambda'])
    
        elif self.OPTIONS['lowlayer_baseflow'] == 'fixedsiz_2':
            qbsat = self.PARS['ks']
        
        else:
            raise Exception('Must be of tens2pll_2,unlimfrc_2,unlimpow_2,topmdexp_2 or fixedsiz_2')
        
        return qbsat

class Interflow(Flux_base):
    """
    Class for Percolation fluxes calculation, with 3 implemented options:
        1. intflwnone: no interflow (hypodermic flow)
        2. intflwsome: linear extraction of free storage
        3. intflwtresh: threshold based interflow of upper layer
 
    Parameters
    ----------
    PARS : dictionary
        Dictionary containing all the model parameters
    OPTIONS: dictionary 
        Dictionary containting all the model structure options  

    Methods
    ---------- 
    calc(STATES,FLUXES)
        The calc method picks the correct definition based on the OPTIONS dictionary defined by the mode structure            
    
    intflwnone(STATES,FLUXES)
        No interflow, like in TOPMODEL and ARNO/VIC
        
        Model States:
            none

        Model Parameters:
            none
        
        Fluxes updated:  
            qif = 0.0
    
    intflwsome(STATES,FLUXES)
        Linear interflow conceptualization from free storage
        
        Model States:
            S1F

        Model Parameters:
            ki, S1Fmax
        
        Fluxes updated:  
            qif
            
    intflwtresh(STATES,FLUXES)
        Linear interflow conceptualization from free storage if above threshold
        
        Model States:
            S1F, S2T (or S1T if surface1_1 option for upper layer)

        Model Parameters:
            ki, tif, S2Tmax (or S1Tmax if surface1_1 option for upper layer)
        
        Fluxes updated:  
            qif

    Attributes
    ------------
    STATES: dict
            Dictionary containing all the model states of the current calculation time step          
    FLUXES: dict
            Dictionary containting all the fluxes structure options

    """

    def __init__(self,PARS,OPTIONS):
        Flux_base.__init__(self,PARS,OPTIONS)
#        print 'Interflow flux loaded'

    def calc(self,STATES,FLUXES):
        if self.OPTIONS['interflow'] == 'intflwnone':
            return self.intflwnone(STATES,FLUXES)

        elif self.OPTIONS['interflow'] == 'intflwsome':
            return self.intflwsome(STATES,FLUXES)
            
        elif self.OPTIONS['interflow'] == 'intflwtresh':
            return self.intflwsome(STATES,FLUXES) 
            
        else:
            raise Exception('No valid option chosen for interflow')

    def intflwnone(self,STATES,FLUXES):        
        qif = 0.0
        FLUXES['qif']=qif
        return FLUXES

    def intflwsome(self,STATES,FLUXES):       
        qif = self.PARS['ki'] * (STATES['S1F']/self.PARS['S1Fmax'])
        FLUXES['qif']=qif
        return FLUXES

    def intflwtresh(self,STATES,FLUXES):     
        if self.OPTIONS['uplayer'] == 'surface1_1':     
            wfif=Logistic(STATES['S1']/self.PARS['S1Tmax'],self.PARS['tif'],Psmooth=0.01)
            qif = self.PARS['ki'] * wfif* STATES['S1F'] *((STATES['S1T']/self.PARS['S1Tmax']-self.PARS['tif'])/(1.-self.PARS['tif']))
            FLUXES['qif']=qif
        else:
            wfif=Logistic(STATES['S2']/self.PARS['S2Tmax'],self.PARS['tif'],Psmooth=0.01)
            qif = self.PARS['ki'] * wfif* STATES['S1F'] *((STATES['S2T']/self.PARS['S2Tmax']-self.PARS['tif'])/(1.-self.PARS['tif']))
            FLUXES['qif']=qif
            
        return FLUXES
        

class Surface(Flux_base):
    """
    Class for Surface fluxes calculation, with 4 implemented options:
        1. arno_x_vic: Probability distribution style (cfr.Variable Infiltration Concept)
        2. prms_varnt: PRMS model concept
        3. tmdl_param: power law transmissivity profile TOPMODEL
        4. oflwtresh: surface flow when above threshold
        (5. testeasy:  only for testing purposes)        
    
    Parameters
    ----------
    PARS : dictionary
        Dictionary containing all the model parameters
    OPTIONS: dictionary 
        Dictionary containting all the model structure options  

    Methods
    ---------- 
    calc(STATES,FLUXES)
        The calc method picks the correct definition based on the OPTIONS dictionary defined by the mode structure            
    
    arno_x_vic(STATES,FLUXES)
        VIC conceptualization for surface storage (or Pareto distribution in Moore-PDM concept)
        
        Model States:
            S1

        Model Parameters:
            S1max, b
            alfah if percolation based on hymod concep
        
        Fluxes updated:  
            qsx 
    
    prms_varnt(STATES,FLUXES)       
        PRMS conceptualization for surface storage, linea dependence,
        catchment conceptualized as one reservoir
   
        Model States:
            S1T

        Model Parameters:
            S1Tmax, Acmax
            alfah if percolation based on hymod concept
        
        Fluxes updated:  
            qsx    
   
    tmdl_param(STATES,FLUXES) 
        TOPMODEL conceptualization for surface storage, based on the topographic
        distribution function

        Model States:
            S2

        Model Parameters:
            S2max, maxpow (derived), loglambda, chi, n
            alfah if percolation based on hymod concept
        
        Fluxes updated:  
            qsx        

        Implementation for FUSE [1] is not really 3-par version, but interpretation of the 1 parameter
        as described in [2], but are essentially the same, see utilities -> gammadistr.py file  
          
    oflwtresh(STATES,FLUXES)    
        Threshold conceptualization for surface storage (cfr. NAM model), 
        catchment conceptualized as one reservoir

        Model States:
            S2T (S1T when surface1_1 upper layer option)

        Model Parameters:
            ko, tof, S2Tmax (S1Tmax when surface1_1 upper layer option)
            alfah if percolation based on hymod concept    
        
        Fluxes updated:  
            qsx   

    Attributes
    ------------
    STATES: dict
            Dictionary containing all the model states of the current calculation time step          
    FLUXES: dict
            Dictionary containting all the fluxes structure options

    """

    def __init__(self,PARS,OPTIONS):
        Flux_base.__init__(self,PARS,OPTIONS)
#        print 'Surface flux loaded'

    def calc(self,STATES,FLUXES):
        if self.OPTIONS['surface'] == 'arno_x_vic':
            return self.arno_x_vic(STATES,FLUXES)

        elif self.OPTIONS['surface'] == 'prms_varnt':
            return self.prms_varnt(STATES,FLUXES) 
            
        elif self.OPTIONS['surface'] == 'tmdl_param':
            return self.tmdl_param(STATES,FLUXES) 
            
        elif self.OPTIONS['surface'] == 'testeasy':
            return self.testeasy(STATES,FLUXES)   

        elif self.OPTIONS['surface'] == 'oflwtresh':
            return self.oflwtresh(STATES,FLUXES)                
        
            
        else:
            raise Exception('No valid option chosen for surface flow')        

        
    def arno_x_vic(self,STATES,FLUXES):       
        Sat_area = 1. - (1. - min(STATES['S1']/self.PARS['S1max'], 1.))**self.PARS['b']
        qsx = Sat_area * FLUXES['rain']
        if self.OPTIONS['percolation'] == 'perc_nodrain':
            FLUXES['qsx']=qsx*self.PARS['alfah']
        else:
            FLUXES['qsx']=qsx            
        return FLUXES
    
    def prms_varnt(self,STATES,FLUXES):
        Sat_area = min(STATES['S1T']/self.PARS['S1Tmax'],1.) * self.PARS['Acmax']
        qsx = Sat_area * FLUXES['rain']
        if self.OPTIONS['percolation'] == 'perc_nodrain':
            FLUXES['qsx']=qsx*self.PARS['alfah']
        else:
            FLUXES['qsx']=qsx  
        return FLUXES        
    
    def tmdl_param(self,STATES,FLUXES):
        nozerodivide=1.e-8 #prevent zero dividing 
        Ti_sat = self.PARS['powlambda']/(STATES['S2']/(self.PARS['S2max']+nozerodivide))
        
        if Ti_sat > self.PARS['maxpow']:
            Sat_area = 0.0
        else:
            Ti_log = np.log(Ti_sat**self.PARS['n'])
            Ti_off=3.0
    #        Ti_chi = (loglambda-Ti_off)/chi
            Ti_arg = max(0., Ti_log - Ti_off) / self.PARS['chi']
            phi=(self.PARS['loglambda']-Ti_off)/self.PARS['chi']
            Sat_area = 1.0 - special.gammainc(phi, Ti_arg)

        qsx = Sat_area * FLUXES['rain']
        
        if self.OPTIONS['percolation'] == 'perc_nodrain':
            FLUXES['qsx']=qsx*self.PARS['alfah']
        else:
            FLUXES['qsx']=qsx  
        return FLUXES 

    def oflwtresh(self,STATES,FLUXES):
        if self.OPTIONS['uplayer'] == 'surface1_1': 
            wfof=Logistic(STATES['S1T']/self.PARS['S1Tmax'],self.PARS['tof'],Psmooth=0.01)
            Sat_area = self.PARS['ko'] * wfof* ((STATES['S1T']/self.PARS['S1Tmax']-self.PARS['tof'])/(1.-self.PARS['tof']))
            qsx = Sat_area * FLUXES['rain']
        
        else:
            wfof=Logistic(STATES['S2T']/self.PARS['S2Tmax'],self.PARS['tof'],Psmooth=0.01)
            Sat_area = self.PARS['ko'] * wfof* ((STATES['S2T']/self.PARS['S2Tmax']-self.PARS['tof'])/(1.-self.PARS['tof']))
            qsx = Sat_area * FLUXES['rain']            
        
        if self.OPTIONS['percolation'] == 'perc_nodrain':
            FLUXES['qsx']=qsx*self.PARS['alfah']
        else:
            FLUXES['qsx']=qsx  
        return FLUXES
    
    def testeasy(self,STATES,FLUXES):
        '''
        Only for testing purposes!
        '''
        qsx = STATES['S']*0.2
        if self.OPTIONS['percolation'] == 'perc_nodrain':
            FLUXES['qsx']=qsx*self.PARS['alfah']
        else:
            FLUXES['qsx']=qsx  
        return FLUXES


class Baseflow(Flux_base):
    """
    Class for Base Flow and lower layer calculation, with 4 implemented options:
        1. tens2pll_2: Two parallel reservoirs 
        2. unlimfrc_2: Storage of unlimited size
        3. unlimpow_2: TOPMODEL adaptive version
        4. fixedsiz_2: Storage of fixed size
        
    Parameters
    ----------
    PARS : dictionary
        Dictionary containing all the model parameters
    OPTIONS: dictionary 
        Dictionary containting all the model structure options  

    Methods
    ---------- 
    calc(STATES,FLUXES)
        The calc method picks the correct definition based on the OPTIONS dictionary defined by the mode structure            
    
    tens2pll_2(STATES,FLUXES)      
        Two parallel linear reservoirs used in conjunction with 2 parallel 
        reservoirs in state equations
        
        Model States:
            S2FA, S2FB

        Model Parameters:
            vA, vB
        
        Fluxes updated:  
            qbA, qbB, qb

    unlimfrc_2(STATES,FLUXES)
        Linear or non-linear reservoir used in combination with a 
        single reservoir of infinite size 
        
        Model States:
            S2

        Model Parameters:
            n, v
        
        Fluxes updated:  
            qb    
    
    unlimpow_2(STATES,FLUXES)
        Non-linear reservoir used in combination with a 
        single reservoir of infinite size; used to conceptualize the 
        TOPMODEL parameterization of the power law transmissivity profile
        
        Model States:
            S2

        Model Parameters:
            n, S2max
            qbsat is calculated in function of the selected lower layer option
        
        Fluxes updated:  
            qb               

    fixedsiz_2(STATES,FLUXES)
        Linear or Non-linear reservoir used in combination with a 
        single reservoir of fixed size
        
        Model States:
            S2

        Model Parameters:
            ks, n, S2max
        
        Fluxes updated:  
            qb 

    Attributes
    ------------
    STATES: dict
            Dictionary containing all the model states of the current calculation time step          
    FLUXES: dict
            Dictionary containting all the fluxes structure options
            
    """

    def __init__(self,PARS,OPTIONS):
        Flux_base.__init__(self,PARS,OPTIONS)
#        print 'Baseflow flux loaded'

    def calc(self,STATES,FLUXES):
        #'lowlayer_baseflow']= 'tens2pll_2','unlimfrc_2','unlimpow_2','fixedsiz_2','topmdexp_2'

        if self.OPTIONS['lowlayer_baseflow'] == 'tens2pll_2':
            return self.tens2pll_2(STATES,FLUXES)

        elif self.OPTIONS['lowlayer_baseflow'] == 'unlimfrc_2':
            return self.unlimfrc_2(STATES,FLUXES) 
            
        elif self.OPTIONS['lowlayer_baseflow'] == 'unlimpow_2':
            return self.unlimpow_2(STATES,FLUXES) 
            
        elif self.OPTIONS['lowlayer_baseflow'] == 'fixedsiz_2': #nonlinear
            return self.fixedsiz_2(STATES,FLUXES)            

#        elif self.OPTIONS['lowlayer_baseflow'] == 'fromupper': #baseflow from upper storage (infiltration)
#            return self.fromupper(STATES,FLUXES)  
        
#        elif self.OPTIONS['lowlayer_baseflow'] == 'topmdexp_2':
#            return self.topmdexp_2(STATES,FLUXES) 
            
        else:
            raise Exception('No valid option chosen for baseflow')        
        
    def tens2pll_2(self,STATES,FLUXES):       
        qbA = self.PARS['vA'] * STATES['S2FA']
        qbB = self.PARS['vB'] * STATES['S2FB']
        qb = qbA + qbB
        FLUXES['qb']=qb
        FLUXES['qbA']=qbA
        FLUXES['qbB']=qbB
        return FLUXES
    
    def unlimfrc_2(self,STATES,FLUXES):       
        qb = self.PARS['v'] * (STATES['S2'])**self.PARS['n']
        FLUXES['qb']=qb
        return FLUXES
    

    def unlimpow_2(self,STATES,FLUXES):        
        qb = self.calc_Qbsat() * ( STATES['S2']/self.PARS['S2max'])**self.PARS['n']
        FLUXES['qb']=qb
        return FLUXES        


    def fixedsiz_2(self,STATES,FLUXES):       
        qb = self.PARS['ks'] * (STATES['S2']/self.PARS['S2max'])**self.PARS['n']
        FLUXES['qb']=qb
        return FLUXES        
        

    def calc_Qbsat(self):        
        if self.OPTIONS['lowlayer_baseflow'] == 'tens2pll_2':
            qbsat = self.PARS['vA']* self.PARS['S2FAmax'] + self.PARS['vB'] * self.PARS['S2FBmax']
    
        elif self.OPTIONS['lowlayer_baseflow'] == 'unlimfrc_2':
            qbsat = self.PARS['v']* self.PARS['S2max']
    
        elif self.OPTIONS['lowlayer_baseflow'] == 'unlimpow_2':
            #      ! This is a bit tricky.  The capacity of the aquifer is m*n, where m is a scaling
            #      ! parameter.  We have the capacity, i.e., MPARAM%MAXWATR_2/1000., and need the
            #      ! TOPMODEL "m" parameter
            TOPmdm = (self.PARS['S2max']/1000.)/self.PARS['n'] # NOTE: mm --> m
            #      ! ...and, compute baseflow
            qbsat = self.PARS['ks'] * (TOPmdm/(self.PARS['powlambda']**self.PARS['n']))
               
        elif self.OPTIONS['lowlayer_baseflow'] == 'fixedsiz_2': #this is redundant
            qbsat = self.PARS['ks']
        
        else:
            raise Exception('Must be of tens2pll_2,unlimfrc_2,unlimpow_2,topmdexp_2 or fixedsiz_2')
        
        return qbsat
   
class Routing(Flux_base):
    """
    Class for Routingof the relevant subflows with 3 options:
        1. rout_all1: routing the subflows combined
        2. no_rout: no routing
        3. rout_ind: routing the subflow indeividual
        
    Parameters
    ----------
    PARS : dictionary
        Dictionary containing all the model parameters
    OPTIONS: dictionary 
        Dictionary containting all the model structure options      

    Methods
    -----------
    calc(STATES,FLUXES,ROUTLIB)
        The calc method picks the correct definition based on the OPTIONS dictionary defined by the mode structure            
        
    rout_all1(STATES,FLUXES)      
        Routing of the sum of all subflows
        
        Model States:
            none, analytical solution used based on fluxes only

        Model Parameters:
            frac_future (derived from number of reservoirs and residence parameter)
        
        Fluxes updated:  
            qgamma, qfuture, q_all

    no_rout(STATES,FLUXES)  
        No routing, subflows are summed for each time stepq_all

        Model States:
            none

        Model Parameters:
            none
        
        Fluxes updated:  
            q_all

    rout_ind(STATES,FLUXES,ROUTLIB)
        Individual routing of the subflows             
    
        Model States:
            none, analytical solution used based on fluxes only

        Model Parameters:
            timeo,timei,timeb
        
        Fluxes updated:  
             q_all, routover, routinter, routbase

    Attributes
    ------------
    STATES: dict
            Dictionary containing all the model states of the current calculation time step          
    FLUXES: dict
            Dictionary containting all the fluxes structure options
    ROUTLIB: dict
        Dictionary containting all the routing characteristics when using linear reservoirs convolution 

    Notes
    -------
    TODO: make also gamma distribution based to enable also non-integer reservoir value
    """
    def __init__(self,PARS,OPTIONS):
        Flux_base.__init__(self,PARS,OPTIONS)
#        print 'Routing flux loaded'    


    def calc(self,STATES,FLUXES,ROUTLIB):
        if self.OPTIONS['routing'] == 'rout_all1':
            return self.rout_all1(STATES,FLUXES)
        elif self.OPTIONS['routing'] == 'no_rout':
            return self.no_rout(STATES,FLUXES)
        elif self.OPTIONS['routing'] == 'rout_ind':
            return self.rout_ind(STATES,FLUXES,ROUTLIB)
        else:
            raise Exception('No valid option chosen for routing')         
        
    def rout_all1(self,STATES,FLUXES):
        #sum fluxes:
        if 'qufof' in FLUXES:
            #additional runof
#            FLUXES['q_all'] =  FLUXES['q_all'] + FLUXES['qufof']
            FLUXES['qsx']  = FLUXES['qsx'] + FLUXES['qufof']                            
        if 'qsfof' in FLUXES:
            #additional baseflow!!
#            FLUXES['q_all'] =  FLUXES['q_all'] + FLUXES['qsfof'] 
            FLUXES['qb'] = FLUXES['qb'] + FLUXES['qsfof'] 
            
        if 'qsx' in FLUXES:
            FLUXES['q_all'] =  FLUXES['qsx']            
        if 'qif' in FLUXES:
            FLUXES['q_all'] =  FLUXES['q_all'] + FLUXES['qif']            
        if 'qb' in FLUXES:
            FLUXES['q_all'] =  FLUXES['q_all'] + FLUXES['qb'] 
        
        #Calculate routing
        ntdh =  self.PARS['frac_future'].size
        for jtim in range(ntdh):
            FLUXES['qfuture'][jtim] = FLUXES['qfuture'][jtim] + FLUXES['q_all']* self.PARS['frac_future'][jtim]
        
        #outflow of the moment
        FLUXES['qgamma'] = FLUXES['qfuture'][0]
        #move array back
        for jtim in range(1,ntdh):
            FLUXES['qfuture'][jtim-1]=FLUXES['qfuture'][jtim]
        FLUXES['qfuture'][ntdh-1] = 0.0 
        
        return FLUXES

    def rout_ind(self,STATES,FLUXES,ROUTLIB):
        #Prepare the output subflows first
        if 'qufof' in FLUXES:
            #additional runof
            FLUXES['qsx']  = FLUXES['qsx'] + FLUXES['qufof']                            
        if 'qsfof' in FLUXES:
            #additional baseflow!!
            FLUXES['qb'] = FLUXES['qb'] + FLUXES['qsfof']  

        if 'routover' in FLUXES:
            FLUXES['routover'] = linres(ROUTLIB['routover'],FLUXES['routover'],FLUXES['qsx'],self.PARS['timeo'])
            FLUXES['q_all'] = FLUXES['routover'][-1]
        else:
            FLUXES['q_all'] = FLUXES['qsx']            

        if 'routinter' in FLUXES:
            FLUXES['routinter'] = linres(ROUTLIB['routinter'],FLUXES['routinter'],FLUXES['qif'],self.PARS['timei'])
            FLUXES['q_all'] =  FLUXES['q_all'] + FLUXES['routointer'][-1]
        else:
            if 'qif' in FLUXES:
                FLUXES['q_all'] = FLUXES['q_all'] + FLUXES['qif']            
            
        if 'routbase' in FLUXES:
            FLUXES['routbase'] = linres(ROUTLIB['routbase'],FLUXES['routbase'],FLUXES['qb'],self.PARS['timeb'])
            FLUXES['q_all'] =  FLUXES['q_all'] + FLUXES['routbase'][-1]
        else:
            FLUXES['q_all'] = FLUXES['q_all'] + FLUXES['qb']  

        return FLUXES            
     
#        #Even trishen en S2 bypassen door q12 te routen en e2 op 0 te zetten
#        if 'routbase' in FLUXES:
#            FLUXES['routbase'] = linres(ROUTLIB['routbase'],FLUXES['routbase'],FLUXES['q12'],self.PARS['timeb'])
#            FLUXES['q_all'] =  FLUXES['q_all'] + FLUXES['routbase'][-1]
#        else:
#            FLUXES['q_all'] = FLUXES['q_all'] + FLUXES['qb']  
    
    def no_rout(self,STATES,FLUXES):
#        FLUXES['q_all'] = FLUXES['qsx'] + FLUXES['qif'] + FLUXES['qb']
#        print FLUXES

        if 'qufof' in FLUXES:
            #additional runof
#            FLUXES['q_all'] =  FLUXES['q_all'] + FLUXES['qufof']
            FLUXES['qsx']  = FLUXES['qsx'] + FLUXES['qufof']                            
        if 'qsfof' in FLUXES:
            #additional baseflow!!
#            FLUXES['q_all'] =  FLUXES['q_all'] + FLUXES['qsfof'] 
            FLUXES['qb'] = FLUXES['qb'] + FLUXES['qsfof'] 
            
        if 'qsx' in FLUXES:
            FLUXES['q_all'] =  FLUXES['qsx']            
        if 'qif' in FLUXES:
            FLUXES['q_all'] =  FLUXES['q_all'] + FLUXES['qif']            
        if 'qb' in FLUXES:
            FLUXES['q_all'] =  FLUXES['q_all'] + FLUXES['qb']     

        return FLUXES


class Misscell(Flux_base):
    """
    Class for miscellaneous bucket overflow fluxes as an automatic consequence 
    of the selected options in upper and lower layer, using logistic smoothing
    functions to prevent from dicontinuities.
        
    Parameters
    -----------
    PARS : dictionary
        Dictionary containing all the model parameters
    OPTIONS: dictionary 
        Dictionary containting all the model structure options  

    Attributes
    ------------
    STATES: dict
            Dictionary containing all the model states of the current calculation time step          
    FLUXES: dict
            Dictionary containting all the fluxes structure options
            
    See Also
    ----------
    pyFUSE.Logistic: For the logistic smoother used
    """
    
    def __init__(self,PARS,OPTIONS,solver=True):
        Flux_base.__init__(self,PARS,OPTIONS)
#        print 'Miscellaneous fluxes flux loaded'


    def calc(self,STATES,FLUXES):
        '''
        The calc method picks the correct definition based on the OPTIONS dictionary defined by the mode structure
                   
        Parameters
        ----------
        STATES : dictionary
            Dictionary containing all the model states of the current calculation time step
        FLUXES: dictionary 
            Dictionary containting all the fluxes structure options, updated by the flux calculation            

        Returns
        -------
        FLUXES: dictionary
            updated values
        '''
        
        #UPPER LAYER
        if self.OPTIONS['uplayer']== 'onestate_1':
            return self.onestate_1(STATES,FLUXES)    
        elif self.OPTIONS['uplayer']== 'tension1_1':
            return self.tension1_1(STATES,FLUXES)    
        elif self.OPTIONS['uplayer']== 'tension2_1':
            return self.tension2_1(STATES,FLUXES) 
        elif self.OPTIONS['uplayer']== 'surface1_1':
            return self.surface1_1(STATES,FLUXES)             
        else:
            raise Exception('No valid option chosen for baseflow upper layer')

        #LOWER LAYER
        if self.OPTIONS['lowlayer_baseflow'] == 'tens2pll_2':
            return self.tens2pll_2(STATES,FLUXES)

        elif self.OPTIONS['lowlayer_baseflow'] == 'unlimfrc_2':
            return self.no_limit(STATES,FLUXES) 
            
        elif self.OPTIONS['lowlayer_baseflow'] == 'unlimpow_2':
            return self.no_limit(STATES,FLUXES) 

        elif self.OPTIONS['lowlayer_baseflow'] == 'topmdexp_2':
            return self.no_limit(STATES,FLUXES) 
            
        elif self.OPTIONS['lowlayer_baseflow'] == 'fixedsiz_2':
            return self.fixedsiz_2(STATES,FLUXES)            
        
        else:
            raise Exception('No valid option chosen for baseflow lower layer')
        
    
    def tension2_1(self,STATES,FLUXES):
        '''
        States
            S1TA,S1TB,S1F
        
        Parameters
            S1TAmax, S1TBmax, S1Fmax
        
        FLUXES: dictionary
            updated by the percolation calculation: qurof, qutof, qufof
        '''
        
        w_func = Logistic(STATES['S1TA'],self.PARS['S1TAmax'])
        qurof = w_func * (FLUXES['rain'] - FLUXES['qsx'])
        FLUXES['qurof']=qurof
        w_func = Logistic(STATES['S1TB'],self.PARS['S1TBmax'])
        qutof = w_func * FLUXES['qurof']
        FLUXES['qutof']=qutof
        w_func = Logistic(STATES['S1F'],self.PARS['S1Fmax'])
        qufof = w_func * FLUXES['qutof']
        FLUXES['qufof']=qufof  
        
        return FLUXES
    
    def tension1_1(self,STATES,FLUXES):
        '''
        States
            S1T,S1F
        
        Parameters
            S1Tmax, S1Fmax
    
        FLUXES: dictionary
            updated by the percolation calculation: qutof, qufof
        '''
        w_func = Logistic(STATES['S1T'],self.PARS['S1Tmax'])
        qutof = w_func * (FLUXES['rain'] - FLUXES['qsx'])
        FLUXES['qutof']=qutof        
        w_func = Logistic(STATES['S1F'],self.PARS['S1Fmax'])
        qufof = w_func * FLUXES['qutof']
        FLUXES['qufof']=qufof         
        
        return FLUXES

    def onestate_1(self,STATES,FLUXES):
        '''
        States
            S1
            
        Parameters
            S1max
        
        FLUXES: dictionary
            updated by the percolation calculation: qufof
        '''
        w_func = Logistic(STATES['S1'],self.PARS['S1max'])
        qufof = w_func * (FLUXES['rain'] - FLUXES['qsx'])
        FLUXES['qufof']=qufof                
        return FLUXES
    
    def surface1_1(self,STATES,FLUXES):
        '''
        States
            S1F
        
        Parameters
            S1Fmax
        
        FLUXES: dictionary
            updated by the percolation calculation: qstof
        '''
        wfr=Logistic(STATES['S1F'],self.PARS['S1Fmax'],Psmooth=0.01)
        qstof = FLUXES['rain']*wfr
        FLUXES['qstof']=qstof                
        return FLUXES       
        
    def tens2pll_2(self,STATES,FLUXES):
        '''
        States
            S2T, S2FA, S2FB
        
        Parameters
            S2Tmax, S2FAmax, S2FBmax, kappa
        
        FLUXES: dictionary
            updated by the percolation calculation: qstof, qsfofa, qsfofb, qsfof
        '''
        # compute flow from tension storage to free storage (mm s-1)
        w_func = Logistic(STATES['S2T'],self.PARS['S2Tmax'])
        qstof = w_func * FLUXES['q12'] * (1. - self.PARS['kappa']) #! in paper niet 1- gezet
        FLUXES['qstof']=qstof
        # compute over-flow of free water in the primary reservoir
        w_func = Logistic(STATES['S2FA'],self.PARS['S2FAmax'])
        qsfofa = w_func * (FLUXES['q12'] * (self.PARS['kappa']/2.) + FLUXES['qstof']/2.)
        FLUXES['qsfofa']=qsfofa
        # compute over-flow of free water in the secondary reservoir
        w_func = Logistic(STATES['S2FB'],self.PARS['S2FBmax'])
        qsfofb = w_func * (FLUXES['q12'] * (self.PARS['kappa']/2.) + FLUXES['qstof']/2.)
        FLUXES['qsfofb']=qsfofb
        # compute total overflow
        FLUXES['qsfof']=qsfofa+qsfofb
        return FLUXES
    
    def no_limit(self,STATES,FLUXES):
        '''
        States
            none
                
        Parameters
            none

        FLUXES: dictionary
            updated by the percolation calculation: qsfof=0.0
        '''
        
        FLUXES['qsfof']=0.0
        return FLUXES
        
    def fixedsiz_2(self,STATES,FLUXES):
        '''
        States
            S2
        
        Parameters
            S2
        
        FLUXES: dictionary
            updated by the percolation calculation: qsfof
        '''
        w_func = Logistic(STATES['S2'],self.PARS['S2max'])        
        FLUXES['qsfof'] = w_func * FLUXES['q12']
        return FLUXES