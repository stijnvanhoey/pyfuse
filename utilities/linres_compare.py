# -*- coding: utf-8 -*-
"""
Created on Sun Sep 23 16:14:13 2012

@author: VHOEYS
"""

import os
import numpy as np

from scipy.interpolate import interp1d
from scipy.integrate import odeint

from scipy.interpolate import interp1d
from scipy import arange, array, exp
from scipy import special

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator,MaxNLocator,LinearLocator,FixedLocator

mpl.rcParams['font.size'] = 12
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['lines.color'] = 'k'
mpl.rcParams['xtick.labelsize'] = 20


def extrap1d(interpolator):
    xs = interpolator.x
    ys = interpolator.y

    def pointwise(x):
        if x < xs[0]:
            return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
        elif x > xs[-1]:
            return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
        else:
            return interpolator(x)

    def ufunclike(xs):
        return array(map(pointwise, array(xs)))

    return ufunclike


def linres(n_res,q_init,cov,co,k):  #nog VHM gewijs met cov te doen van vorige tijdstap
    if n_res==1:
#        print q_init[0],'qinit1'
#        print np.exp(-1/k),'ewp'
        q_init[0]=q_init[0]*np.exp(-1./k) + (co+cov)*(1 - np.exp(-1/k))/2
#        print q_init[0],'qinit2'
        return q_init
    else:
        q_init[n_res-1]=q_init[n_res-1]* np.exp(-1/k)+linres(n_res-1,q_init,cov,co,k)[n_res-2]*(1 - np.exp(-1/k))
        return q_init



def deriv(u,t,rain,const):
#    print t
    k=0.2
    RAIN=rain(t)
#    print RAIN
    du1 = RAIN-k*u[0]
    du2 = k*u[0] -k*u[1]
    du3 = k*u[1] -k*u[2]
    
    return np.array([du1,du2,du3])        
    

def NAM_LSODA(const,InitCond,rain):
    
    totn = rain.size
#    print totn

    #Define time-steps to give output for
    ##################################
    totaltime=np.float(totn)
    time=np.arange(0,totaltime,1.)
#    print time
#    time= np.linspace(0.0,totaltime,totaltime) #eg. hours if datainput hours and timestep=1.0; geeft mogelijkheid om met daginput ook uur-output te verkrijgen. NIET met uur ook dag (teveel gezever ivm hoe uurlijke info omzetten in dag, kan beter preprocessing zijn)
#    print time
    
    #Prepare timeseries for linear interpolation in substepping
    ##################################
    rain_int = interp1d(time,rain,bounds_error=False)
#    rain_int2 = extrap1d(rain_int)
       
#    f_i = interp1d(x, y)
#    f_x = extrap1d(f_i)

    

    #Solve with LSODA scheme
    ##################################
    y=odeint(deriv,InitCond,time,full_output=0, printmessg=True, args=(rain_int,const),hmax=1.)

    #Calculate fluxes and prepare outputs
    ##################################
#    v = np.float64(area * 1000.0 / (60.0 * 60.0))

    return y
    

#datapath="D:\Modellen\Version2012\HPC"
#Rain=np.loadtxt(os.path.join(datapath,'Rain_Cal_warm_1jan02_31dec05'))
#Rain=Rain[:200]
#
##ODE SOLVER
#InitCond=np.array([0.,0.,0.])
#const=0.
#trg=NAM_LSODA(const,InitCond,Rain) 
#k=0.2
#plt.plot(k*trg)   
#
#
##LINRES SOLUTION CHOW
#nn=3
#ffs=np.zeros((201,nn))
#k=0.2
#ff=np.ones(nn)*0.0
#ffs[0]=ff
#
#for t in range(200):
#    ff=linres(nn,ff,Rain[t-1],Rain[t],1./k)
#    ffs[t]=ff
#    
#plt.plot(ffs)


#GAMMA DISTRIBUTIE
def gammat(nn,k):
    '''
    gamma-function based weight function to control the runoff delay => Chow, 1988 
    '''
    
    frac_future=np.zeros(50.)  #Parameter added
    ntdh = frac_future.size
    
    deltim=1. 
#    print 'qtimedelay is calculated with a unit of',deltim,'hours to have parameter values comparable to Clarke, 2008'
    
    for jtim in range(ntdh):
        tfuture=jtim*deltim
        prob=(1./(k*special.gamma(nn)))*(tfuture/k)**(nn-1) * np.exp(-tfuture/k)
        
        frac_future[jtim]=max(0.,prob)
    
#    if cumprob < 0.99:
#        print 'not enough bins in the frac_future'

    #make sure sum to one
    frac_future[:]=frac_future[:]/frac_future[:].sum()
    
    return frac_future
    
def gamma_tdelay(set_par):
    '''
    gamma-function based weight function to control the runoff delay => Chow, 1988 
    '''
    nn=set_par['nres']
    
    frac_future=np.zeros(500.)  #Parameter added
    ntdh = frac_future.size
    
    deltim=1. 
#    print 'qtimedelay is calculated with a unit of',deltim,'hours to have parameter values comparable to Clarke, 2008'
    
    for jtim in range(ntdh):
        tfuture=jtim*deltim
        prob=(1./(k*special.gamma(nn)))*(tfuture/k)**(nn-1) * np.exp(-tfuture/k)
        
        frac_future[jtim]=max(0.,prob)
    
#    if cumprob < 0.99:
#        print 'not enough bins in the frac_future'

    #make sure sum to one
    frac_future[:]=frac_future[:]/frac_future[:].sum()
    
    return frac_future    


def qtimedelay(set_par,deltim=1.):
    '''
    gamma-function based weight function to control the runoff delay
    '''
    alpha=3.
    alamb = alpha/set_par['mut']
    psave=0.0
    set_par['frac_future']=np.zeros(50.)  #Parameter added
    ntdh = set_par['frac_future'].size
    
    deltim=deltim 
    print 'qtimedelay is calculated with a unit of',deltim,'hours to have parameter values comparable to Clarke, 2008'
    
    for jtim in range(ntdh):
#        print jtim
        tfuture=jtim*deltim
#        print alamb*tfuture
        cumprob= special.gammainc(alpha, alamb*tfuture)# hoeft niet want verschil wordt genomen: /special.gamma(alpha)
#        print cumprob
        set_par['frac_future'][jtim]=max(0.,cumprob-psave)
        psave = cumprob
    
    if cumprob < 0.99:
        print 'not enough bins in the frac_future'

    #make sure sum to one
    set_par['frac_future'][:]=set_par['frac_future'][:]/set_par['frac_future'][:].sum()
    
    return set_par    

#nn=2.
#k=2.0
#
#frr=gammat(nn,k)
#qfuture = np.zeros(frr.size)
#qoo=np.zeros(Rain.size)
#
#for t in range(200):
#    #Calculate routing
#    ntdh = frr.size
#    for jtim in range(ntdh):
#        qfuture[jtim] = qfuture[jtim] + Rain[t] * frr[jtim]
#    
#    #outflow of the moment
#    qoo[t]=qfuture[0]
#    if qfuture[0] > 0.0:
#        print qfuture[0]
#    
#    #move array back
#    for jtim in range(1,ntdh):
#        qfuture[jtim-1]=qfuture[jtim]
    qfuture[ntdh-1] = 0.0 

#plt.plot(qoo,'--')    
#
#set_par={}
#set_par['mut']=1/0.2
#pars=qtimedelay(set_par,deltim=1.)
#tt=pars['frac_future']
#
#plt.plot(tt)
#plt.plot(frr)
#
#
#
#####PHD PLOT
#plt.figure() 
#plt.subplots_adjust(wspace = 0.05)
#
#ax1=plt.subplot(121)   
#nn=2.
#k= [2.5,5.,10.]
#liness=['k-','k--','k-.']
#
#cnt=0
#for ll in k:
#    frr=gammat(nn,ll)
#    ax1.plot(frr,liness[cnt],label=r'$k$ = '+str(ll))
#    cnt+=1    
#      
##plt.xlabel(r'$\zeta$ ($ln(\alpha / \tan \beta $)')
#plt.xlabel(r'time')
#plt.ylabel(r'$h(t)$')
#plt.legend()
#
#majorLocator1= MaxNLocator(nbins=3,integer=True)
#ax1.xaxis.set_major_locator(majorLocator1)
#
#majorLocator2= MaxNLocator(nbins=3)
#ax1.yaxis.set_major_locator(majorLocator2)
#
#ax2=plt.subplot(122,sharey=ax1)   
#ax2.get_yaxis().set_visible(False)
#nn=[2.,3.5, 5.]
#k=4.
#liness=['k-','k--','k-.']
#
#cnt=0
#for ll in nn:
#    frr=gammat(ll,k)
#    ax2.plot(frr,liness[cnt],label=r'$n$ = '+str(ll))
#    cnt+=1
#
#ax2.xaxis.set_major_locator(majorLocator1)
#      
##plt.xlabel(r'$\zeta$ ($ln(\alpha / \tan \beta $)')
#plt.xlabel(r'time')
##plt.ylabel(r'$\frac{A_c}{A}$')
#plt.legend()

#plt.savefig('Rout_pars.png')
#plt.savefig('Rout_pars.pdf')    
#plt.savefig('Rout_pars.eps')    


def Logistic1(State,Statemax,Psmooth=0.01):
    ''' 
    Uses a logistic function to smooth the threshold at the top of a bucket
    
    See literature:    
    Clark, Martyn P., A. G. Slater, D. E. Rupp, R. A. Woods, Jasper A. Vrugt, H. V. Gupta, Thorsten Wagener, and L. E. Hay. Framework for Understanding Structural Errors (FUSE): A modular framework to diagnose differences between hydrological models. Water Resources Research 44 (2008): 14.
    Original code from Clark, Martyn P.
    '''
    epsilon = 5.0       #Multiplier to ensures storagde is always less than capacity
    
    Asmooth = Psmooth * Statemax        #actual smoothing
#    LOGISMOOTH = 1. / ( 1. + np.exp(-(State - (Statemax - Asmooth * epsilon) ) / Asmooth) )
    LOGISMOOTH = 1. / ( 1. + np.exp(-(State - (Statemax) ) / Asmooth) )
    return LOGISMOOTH
    
def Logistic2(State,Statemax,Psmooth=0.01):
    ''' 
    Uses a logistic function to smooth the threshold at the top of a bucket
    
    See literature:    
    Clark, Martyn P., A. G. Slater, D. E. Rupp, R. A. Woods, Jasper A. Vrugt, H. V. Gupta, Thorsten Wagener, and L. E. Hay. Framework for Understanding Structural Errors (FUSE): A modular framework to diagnose differences between hydrological models. Water Resources Research 44 (2008): 14.
    Original code from Clark, Martyn P.
    '''
    epsilon = 5.0       #Multiplier to ensures storagde is always less than capacity
    
    Asmooth = Psmooth * Statemax        #actual smoothing
    LOGISMOOTH = 1. / ( 1. + np.exp(-(State - (Statemax - Asmooth * epsilon) ) / Asmooth) )
#    LOGISMOOTH = 1. / ( 1. + np.exp(-(State - (Statemax) ) / Asmooth) )
    return LOGISMOOTH    
    
    
plt.figure()  
plt.subplots_adjust(wspace = 0.05)
ax1=plt.subplot(121)  
Smax=50.
S=np.arange(0.,100,0.1)
ll=np.zeros(S.size)
smooth=[0.0001,0.01,0.1]

tresh1=np.zeros(500)
tresh2=np.ones(500)
tresh=np.hstack((tresh1,tresh2))
ax1.plot(S,tresh,'k',label='step',linewidth=2.)
liness=['k:','k--','k-.']

cnt=0
for smo in smooth:
    smo=smo*Smax
    for i in range(S.size):
        ll[i]=Logistic1(S[i],Smax,Psmooth=smo)
    ax1.plot(S,ll,liness[cnt],label=r'$\omega$='+str(smo))
    cnt+=1

ax1.set_ylim([0.0,1.02])
ax1.set_ylabel(r'$\Phi(S,S_{max},\omega)$')

ax1.set_xticks([50.])
ax1.set_xticklabels([r'$S_{max}$'])
ax1.set_yticks([0.,0.5,1.])
#ax1.set_ylabel()
ax1.legend(loc=4,frameon=False,bbox_to_anchor=(1.05,0.001))
#ax1.legend(loc=4,frameon=True,bbox_to_anchor=(1.2,0.001))

#secodnd subplot
ax2=plt.subplot(122)  
Smax=50.
S=np.arange(0.,55.,0.1)
ll=np.zeros(S.size)
smooth=[0.0001,0.01,0.1]

tresh1=np.zeros(500.)
tresh2=np.ones(10*5)
tresh=np.hstack((tresh1,tresh2))
ax2.plot(S,tresh,'k',label='step',linewidth=2.)
liness=['k:','k--','k-.']

cnt=0
for smo in smooth:
    smo=smo*Smax
    for i in range(S.size):
        ll[i]=Logistic2(S[i],Smax,Psmooth=smo)
    ax2.plot(S,ll,liness[cnt],label=r'$\omega$='+str(smo))
    cnt+=1

ax2.get_yaxis().set_visible(False)
ax2.set_ylim([0.0,1.02])
ax2.set_xlim([0.0,55.])

ax2.set_xticks([50.])
ax2.set_xticklabels([r'$S_{max}$'])
#ax2.set_yticks([0.,0.5,1.])
#ax1.set_ylabel()
#ax2.legend(loc=8)
#leg=ax2.legend(loc='upper center', bbox_to_anchor=(-0.025, 1.12),mode="expand",frameon=False, shadow=False,ncol=4)





plt.savefig('Smooth_pars.png')
plt.savefig('Smooth_pars.pdf')    
plt.savefig('Smooth_pars.eps') 

    
    
    
    
    
    




