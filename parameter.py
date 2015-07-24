# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 16:37:52 2012

@author: VHOEYS
"""

import numpy as np

#from pyFUSE.distributions import *
from distributions import *

class ModPar(object):    
    """
    Model Parameter class
    
    The parameter is defined by his distribution, boundaries and optimal guess
    
    Parameters
    ----------
    name : string
        Name of the parameter
    minval : float
        Minimum value of the parameter distribution
    maxval :  float
        Maximum value of the parameter distribution
    optguess : float
        Optimal guess of the parameter, must be between min and max value
    pardistribution : string
        choose a distributionfrom: randomUniform, randomTriangular, randomTrapezoidal, randomNormal, randomLogNormal
    *kargs  : 
        Extra arguments necessaty for the chosen distribution
    
    Examples
    ----------
    >>> import pyFUSE as pf
    >>> par1=pf.ModPar('parameter1',0.0,3.0,2.5,'randomUniform')
    >>> par1.aValue()
    >>> 2.20263437200365 #random
    >>> par1.optguess
    >>> 2.5
    >>> par1.pardistribution
    >>> 'randomUniform
    """
    def __init__(self,name,minval,maxval,optguess,pardistribution,*kargs):
        if not isinstance(pardistribution,str):
            raise ValueError('Parameter name must be a string')
        else:
            self.name=name
        if not isinstance(minval,float):
            raise ValueError('Minvalue must be Float')
        else:
            self.min = minval
        if not isinstance(minval,float):
            raise ValueError('Maxvalue must be Float')
        else:
            self.max = maxval
        self.bound = maxval-minval
        if self.bound < 0.0:
            raise Exception('Minimal and maximum value of bounds interchanged')
        if not isinstance(optguess,float):
            raise ValueError('Maxvalue must be Float')
        elif not minval<optguess<maxval:
            raise Exception('optimal guess must be between min and max value')
        else:
            self.optguess = optguess
        if not isinstance(pardistribution,str): #todo: aangeven welke types mogelijk
            raise ValueError('Pardistribution must be a string describing the Priori estimated distribution; choose from: randomUniform, randomTriangular, randomTrapezoidal, randomNormal, randomLogNormal')
        else:
            self.pardistribution = pardistribution
            
        if self.pardistribution =='randomUniform':
            if not len(kargs)==0:
                raise Exception('No extra arguments added when using randomuniform')
            self.mode='NaN'
            self.mode1='NaN'
            self.mode2='NaN'
            self.mu='NaN'
            self.sigma='NaN'

#        elif self.pardistribution =='discreteUniform':
#            if not len(kargs)==0:
#                raise Exception('No extra arguments added when using discreteuniform')
#            self.mode='NaN'
#            self.mode1='NaN'
#            self.mode2='NaN'
#            self.mu='NaN'
#            self.sigma='NaN'
            
        elif self.pardistribution =='randomTriangular':
            if not len(kargs)==1:
                raise Exception('randomTraingular needs one extra argument: mode (center of triangle)')
            if not minval<kargs[0]<maxval:
                raise Exception('mode must be between min and max value')
            if not isinstance(kargs[0],float):
                raise ValueError('mode must be Float')
            else:
                self.mode=kargs[0]
            self.mode1='NaN'
            self.mode2='NaN'
            self.mu='NaN'
            self.sigma='NaN'
        elif self.pardistribution =='randomTrapezoidal':
            if not len(kargs)==2:
                raise Exception('randomTrapezoidal needs two extra argument: mode1 and mode2')
            if not isinstance(kargs[0],float):
                raise ValueError('mode 1 must be Float')
            if not isinstance(kargs[1],float):
                raise ValueError('mode 2 must be Float')
            if not kargs[0]<kargs[1]:
                raise Exception('mode 1 must have smaller value then mode 2')
            if not minval<kargs[0]<maxval:
                raise Exception('mode 1 must be between min and max value')
            if not minval<kargs[1]<maxval:
                raise Exception('mode 2 must be between min and max value')
            else:
                self.mode1=kargs[0]
                self.mode2=kargs[1]
            self.mu='NaN'
            self.sigma='NaN'
        elif self.pardistribution =='randomNormal':
            if not len(kargs)==2:
                raise Exception('randomNormal needs two extra argument: mean and var')
            if not 2*kargs[1]+kargs[0]<maxval:
                raise Exception('maxval smaller then average+2*std')
            if not kargs[0]-2*kargs[1]<maxval:
                raise Exception('minval higher then average-2*std')
            else:
                self.mu=kargs[0]
                self.sigma=kargs[1]
        elif self.pardistribution =='randomLogNormal':
            if not len(kargs)==2:
                raise Exception('randomNormal needs two extra argument: mode1 and mode2')
            if not 2*kargs[1]+kargs[0]<maxval:
                raise Exception('maxval smaller then average+2*std')
            if not kargs[0]-2*kargs[1]<maxval:
                raise Exception('minval higher then average-2*std')
            else:
                self.mu=kargs[0]
                self.sigma=kargs[1]             
               
                
        #to add: LHSuniform, LHSnormal, pseudorandomUniform, pseudirandomNormal
        else:
            raise Exception('Wrong ditribution name! choose from: randomUniform, randomTriangular, randomTrapezoidal, randomNormal, randomLogNormal')

    def MCSample(self,nruns):
        '''
        Give a sample of nMC samples from the par distribution 
        
        Parameters
        ----------
        nruns: int
            number of Monte Carlo samples to take
        
        Returns
        --------
        mcsample: array
            numpy array with the specified number of Monte Carlo samples
        
        '''
        if self.pardistribution =='randomUniform':
            return randomUniform(left=self.min,right=self.max,rnsize=nruns)
#        if self.pardistribution =='discreteUniform':
#            return randomUniform(left=self.min,right=self.max,rnsize=nruns)            
        elif self.pardistribution =='randomTriangular':
            return randomTriangular(left=self.min, mode=self.mode, right=self.max, rnsize=nruns)
        elif self.pardistribution =='randomTrapezoidal':
            return randomTrapezoidal(left=self.min,mode1=self.mode1,mode2=self.mode2,right=self.max,rnsize=nruns)
        elif self.pardistribution =='randomNormal':
            return randomNormal(mu=self.mu, sigma=self.sigma, rnsize=nruns)
        elif self.pardistribution =='randomLogNormal':
            return randomLogNormal(mu=self.mu, sigma=self.sigma, rnsize=nruns)

    def aValue(self):
        '''
        Sample 1 value of the pardistribution
        '''
        if self.pardistribution =='randomUniform':
            return randomUniform(left=self.min,right=self.max,rnsize=1)[0]
        elif self.pardistribution =='randomTriangular':
            return randomTriangular(left=self.min, mode=self.mode, right=self.max, rnsize=1)[0]
        elif self.pardistribution =='randomTrapezoidal':
            return randomTrapezoidal(left=self.min,mode1=self.mode1,mode2=self.mode2,right=self.max,rnsize=1)[0]
        elif self.pardistribution =='randomNormal':
            return randomNormal(mu=self.mu, sigma=self.sigma, rnsize=1)[0]
        elif self.pardistribution =='randomLogNormal':
            return randomLogNormal(mu=self.mu, sigma=self.sigma, rnsize=1)[0]
    
    def set_optguess(self,value):
        '''
        Change the optimal guess value of the parameter
        
        Parameters
        ----------
        value: float
            New value for optimal guess
        '''
        self.optguess=value
    
    def optguess_from_optimization(self):
        raise Exception('Not provided yet since integration with total model is not yet provided')

    def ahist(self,nruns,nbins=30,saveit='show',*args,**kwargs):
        '''
        Returns a histogram from the current parameter distribution
        
        Parameters
        ------------
        nruns: int
            number of Monte Carlo samples
        nbins: int
            number of bins to use in the histogram
        saveit; 'show' or True
            if True, figure is saved, otherwise the picture is shown
        *args:
            matplotlib histogram keyword arguments
        **kwargs:
            kwargs are used to update the properties of the
            class:`~matplotlib.patches.Patch` instances returned by *hist*
                
        Returns
        --------
        fig: 
            figure with the histogram
        '''
        plt.hist(self.MCSample(nruns),bins=nbins,*args,**kwargs)
        plt.xlabel(self.name)
        if saveit==True:
            plt.savefig('par_'+self.name+'_hist_'+str(nruns)+'runs.png', dpi=300, facecolor='w', edgecolor='w',orientation='portrait', papertype=None,transparent=False)
        elif saveit=='show':
            plt.show()
        else:
            return fig

    def LatinH(self,nruns):
        '''
        Return a sample of nMC samples from the par distribution 
        with Latin Hypercube sampling (always randomUniform)
        
        Parameters
        -----------
        nruns: int
            number of Latin HYpercube samples to take
               
        '''
        if not self.pardistribution =='randomUniform':
            raise Exception('Latin hypercube only supported for uniform distribution')
            
        pranges=[self.name]
        low=self.min
        high=self.max
        delta=(high-low)/float(nruns)
        for j in range(0,nruns):
            pranges.append(random.uniform(low+j*delta,low+(j+1)*delta))

        s=range(0,nruns)
        result=[]
        for i in range(0,nruns):
            a = random.sample(s,1)[0]
            result.append(a)
            s.remove(a)

        sample=[]
        for j in range(0,len(result)):
            sample.append(pranges[result[j]+1])

        return np.array(sample)


def reScale(arr,vmin,vmax):
    '''
    Rescale the sampled values between 0 and 1 towards the real boundaries of the pars
    
    Parameters
    -----------
    arr: array
        array of the sampled values
    vmin: float
        minimal value to rescale to
    vmax: float
        maximum value to rescale to
    '''
    arrout=(vmax-vmin)*arr+vmin
    return arrout
   
def Sobol(ParsIn,nruns,seed=1):
    '''
    Return a sobol sampling of the parameter space; Sobol is always performed on the entire
    set of parameters used in the analysis.
    
    Parameters
    ------------    
    ParsIn: list of ModPar instances
        List with all the parameters to sample from
    nruns: int
        number of samples
    seed: int
        seed to start from, change this when performing multiple samples or to make sure
        the values are continued by using the last seed
    
    Returns
    --------
    Pars: narray 
        2D array with the rows the different runs and the pars in the columns
    '''
    ndim=len(ParsIn)
    Pars=np.zeros((nruns,ndim))
    
    for i in xrange(1,nruns+1):        
        [r, seed_out] = i4_sobol(ndim, seed)
        Pars[i-1,:]=r        
        seed = seed_out
    for i in range(ndim):
        Pars[:,i]=reScale(Pars[:,i],ParsIn[i].min,ParsIn[i].max)
    print 'The seed to continue this sampling is',seed
        
    return Pars