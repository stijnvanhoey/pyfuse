# -*- coding: utf-8 -*-
"""
Created on Wed Sep 19 14:40:44 2012

@author: VHOEYS
"""
import numpy as np
import scipy as sp
from scipy import special
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['font.size'] = 14
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['lines.color'] = 'k'

def calc_meantipow(off,loglam,chi,n):
    '''
    calculate mean of power transformed topographic index
    
    See literature:    
    Clark, Martyn P., A. G. Slater, D. E. Rupp, R. A. Woods, Jasper A. Vrugt, H. V. Gupta, Thorsten Wagener, and L. E. Hay. Framework for Understanding Structural Errors (FUSE): A modular framework to diagnose differences between hydrological models. Water Resources Research 44 (2008): 14.
    Original code from Clark, Martyn P.
    '''
    Ti_off=off
    Ti_shp=chi                          #shape of the Gamma distribution  chi eigenlijk
    Ti_chi= (loglam-Ti_off)/Ti_shp    #Chi -- loglamb is the first parameter (mean) phi eigenlijk
    print Ti_chi,'tichi'
    
    nn=n
    
    # values for testing (Sivapalan et al., WRR, December 1987)
#    TI_OFF = 3.82_SP  ! TI_OFF = 2.92_SP
#    TI_SHP = 2.48_SP  ! TI_SHP = 3.52_SP
#    TI_CHI = 1.00_SP  ! TI_CHI = 0.742_SP

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
        GMARG2 = max(0., UPPERV - Ti_off) / Ti_chi          # 2nd argument to the Gamma function
        UPPERP = special.gammainc(Ti_shp, GMARG2)           # GAMMP is the incomplete Gamma function GAMMP(Ti_shp, GMARG2)
        PROBIN = UPPERP-LOWERP                              # probability of the current bin
        # get the scaled topographic index value
        LOGVAL = 0.5*(LOWERV+UPPERV)                        # log-transformed index for the current bin
        POWVAL = (np.exp(LOGVAL))**(1./nn)        # power-transformed index for the current bin
        AVELOG = AVELOG + LOGVAL*PROBIN                     #! average log-transformed index (testing)
        AVEPOW += POWVAL*PROBIN                             # average power-transformed index
#        print LOWERV, UPPERV, LOGVAL, POWVAL, AVEPOW
#        !write(*,'(7(f9.3,1x))') lowerv, upperv, logval, powval, avelog, avepow
        # save the lower value and probability
        LOWERV = UPPERV                                     # lower value for the next bin
        LOWERP = UPPERP                                     # cumulative probability for the next bin

    return POWVAL,AVEPOW,AVELOG
   
def calc_meantipow_2(off,loglam,chi,nn):
    '''
    calculate mean of power transformed topographic index
    See literature:    
    Clark, Martyn P., A. G. Slater, D. E. Rupp, R. A. Woods, Jasper A. Vrugt, H. V. Gupta, Thorsten Wagener, and L. E. Hay. Framework for Understanding Structural Errors (FUSE): A modular framework to diagnose differences between hydrological models. Water Resources Research 44 (2008): 14.
    Original code from Clark, Martyn P.
    '''
    mu=off
    loglambda=loglam                    
    phi= (loglambda-mu)/chi    #Chi -- loglamb is the first parameter (mean) phi eigenlijk
    n=nn
    
    # values for testing (Sivapalan et al., WRR, December 1987)
#    TI_OFF = 3.82_SP  ! TI_OFF = 2.92_SP
#    TI_SHP = 2.48_SP  ! TI_SHP = 3.52_SP
#    TI_CHI = 1.00_SP  ! TI_CHI = 0.742_SP

    # loop through the frequency distribution
    avelog = 0.  #testing
    avepow = 0.
    
    Nbins=20000
    Ti_max=50.
    width=Ti_max/Nbins
    
    for ibin in range (1,Nbins):
        # get probability for the current bin
        zeta = (float(ibin)/Nbins) * Ti_max               # upper value in frequency bin
#        print zeta
        
        temp=max(0.,(zeta-mu)/chi)
        fzeta=width *(1./(chi*special.gamma(phi))) * temp**(phi-1) * np.exp(-temp)  
        powval = (np.exp(zeta))**(1./n)
        avepow = avepow + fzeta*powval
        avelog = avelog + zeta*fzeta
        
    return powval,avepow,avelog

def calc_meantipow_nv(off,loglam,chi,n):
    '''
    calculate mean of power transformed topographic index
    needs par-library as input!
    
    See literature:    
    Clark, Martyn P., A. G. Slater, D. E. Rupp, R. A. Woods, Jasper A. Vrugt, H. V. Gupta, Thorsten Wagener, and L. E. Hay. Framework for Understanding Structural Errors (FUSE): A modular framework to diagnose differences between hydrological models. Water Resources Research 44 (2008): 14.
    Original code from Clark, Martyn P.
    '''
    mu=off
    loglambda=loglam                    
    phi= (loglambda-mu)/chi    #Chi -- loglamb is the first parameter (mean) phi eigenlijk
    nn=n

    # loop through the frequency distribution
    LOWERV = 0.
    LOWERP = 0.
    AVELOG = 0.  #testing
    AVEPOW = 0.

    Nbins=20000
    Ti_max=50.
    
    for ibin in range (1,Nbins):
        # get probability for the current bin
        UPPERV = (float(ibin)/Nbins) * Ti_max               # upper value in frequency bin
        GMARG2 = max(0., UPPERV - mu) / chi          # 2nd argument to the Gamma function  Ti_arg = max(0., Ti_log - Ti_off) / chi
        UPPERP = special.gammainc(phi, GMARG2)           # GAMMP is the incomplete Gamma function GAMMP(Ti_shp, GMARG2)
        PROBIN = UPPERP-LOWERP                              # probability of the current bin
        # get the scaled topographic index value
        LOGVAL = 0.5*(LOWERV+UPPERV)                        # log-transformed index for the current bin
        POWVAL = (np.exp(LOGVAL))**(1./nn)        # power-transformed index for the current bin
        AVELOG = AVELOG + LOGVAL*PROBIN                    # ! average log-transformed index (testing)
        AVEPOW += POWVAL*PROBIN                             # average power-transformed index
#        print LOWERV, UPPERV, LOGVAL, POWVAL, AVEPOW
#        !write(*,'(7(f9.3,1x))') lowerv, upperv, logval, powval, avelog, avepow
        # save the lower value and probability
        LOWERV = UPPERV                                     # lower value for the next bin
        LOWERP = UPPERP                                     # cumulative probability for the next bin
   
    return POWVAL,AVEPOW,AVELOG


#calculate value 
nn=10.
mu=3.
   
chi=1.0
phi=2.48
loglambda= chi*phi+mu
powval,avepow,avelog=calc_meantipow_2(mu,loglambda,chi,nn)  

#1 parameter implementation => phi is chi
chi=1.0
phi=2.48
loglambda= chi*phi+mu  
powvala,avepowa,aveloga=calc_meantipow(mu,loglambda,phi,nn)

#3 parameter implementation
chi=1.0
phi=2.48
loglambda= chi*phi+mu  
powvalb,avepowb,avelogb=calc_meantipow_nv(mu,loglambda,chi,nn)
print avepow,avepowa,avepowb
print avelog,aveloga,avelogb

###############################################################################
##FOR OVERLAND FLOW   
###############################################################################
mu=3.
plt.figure() 
plt.subplots_adjust(wspace = 0.2)

plt.subplot(121)   
chi=1.0
loglambda= [5.0,8.0,10.0]
liness=['k-','k--','k-.']

cntt=0
for ll in loglambda:
    phi= (ll-mu)/chi
    zeta=np.arange(.0,18.,0.01)
    fzeta=np.zeros(zeta.size)
    cnt=0
    for i in zeta:
        temp=max(0.,(i-mu)/chi)
        fzeta[cnt]=1./(chi*special.gamma(phi)) * temp**(phi-1) * np.exp(-temp)
        cnt+=1
    
    print cntt,'cnt'
    plt.plot(zeta,1.-fzeta.cumsum()/100.,liness[cntt],label=r'$\lambda$ = '+str(ll))
    cntt+=1
      
#plt.xlabel(r'$\zeta$ ($ln(\alpha / \tan \beta $)')
plt.xlabel(r'$\zeta$')
plt.ylabel(r'$\frac{A_c}{A}$')
plt.legend()

plt.subplot(122)   
chi=[0.1,1.25,3.0]
loglambda= 7.5
liness=['k-','k--','k-.']

cntt=0
for ch in chi:
    phi= (loglambda-mu)/ch
    zeta=np.arange(.0,18.,0.01)
    fzeta=np.zeros(zeta.size)
    cnt=0
    for i in zeta:
        temp=max(0.,(i-mu)/ch)
#        print temp,'temp'
        fzeta[cnt]=1./(ch*special.gamma(phi)) * temp**(phi-1) * np.exp(-temp)
        cnt+=1
    
#    print cntt,'cnt'
    plt.plot(zeta,1.-fzeta.cumsum()/100.,liness[cntt],label=r'$\chi$ = '+str(ch))
    cntt+=1
      
#plt.xlabel(r'$\zeta$ ($ln(\alpha / \tan \beta $)')
plt.xlabel(r'$\zeta$')
#plt.ylabel(r'$\frac{A_c}{A}$')
plt.legend()

#testcase:
S1=np.arange(0.1,499,0.1)
sata=np.zeros(S1.size)

#FOR OVERLAND FLOW: 1 par implementation!
chi=2.48
phi=1.0
loglambda= chi*phi+mu

for i in range(S1.size):
    nozerodivide=1.e-8 #prevent zero dividing 
    Ti_sat = avepow/(S1[i]/(500.+nozerodivide))
    
    if Ti_sat > powval:
        Sat_area = 0.0
    else:
        Ti_log = np.log(Ti_sat**nn)
        Ti_off=mu
        Ti_chi = (loglambda-Ti_off)/chi
        Ti_arg = max(0., Ti_log - Ti_off) / Ti_chi
        sata[i] = 1.0 - special.gammainc(chi, Ti_arg)


#FOR OVERLAND FLOW: 3 par implementation!
chi=1.0
phi=2.48
loglambda= chi*phi+mu
#phi=(loglambda-mu)/chi
#print phi,'phi'
#t1=sp.special.gammainc(phi,(zeta_crit-mu)/chi)
satb=np.zeros(S1.size)
for i in range(S1.size):
    nozerodivide=1.e-8 #prevent zero dividing 
    Ti_sat = avepow/(S1[i]/(500.+nozerodivide))
    
    if Ti_sat > powval:
        Sat_area = 0.0
    else:
        Ti_log = np.log(Ti_sat**nn)
        Ti_off=mu
#        Ti_chi = (loglambda-Ti_off)/chi
        Ti_arg = max(0., Ti_log - Ti_off) / chi
#        satb[i] = 1.0 - special.gammainc(phi, Ti_arg)
        satb[i] = special.gammaincc(phi, Ti_arg)

plt.figure()  
plt.plot(S1,sata)
plt.plot(S1,satb)
plt.title('BEIDE IMPLEMENTATIES EFFECTIEF ANALOOG')

###############################################################################
#CONCLUSION:
    #both ar equal, but chi en phi get interchanged meaning!!
###############################################################################

set_par={}
set_par['mut']=2.

def qtimedelay(set_par,deltim=1.):
    '''
    gamma-function based weight function to control the runoff delay
    '''
    alpha=3.0
    print set_par['mut']
    alamb = alpha/set_par['mut']
    psave=0.0
    set_par['frac_future']=np.zeros(500.)  #Parameter added
    ntdh = set_par['frac_future'].size
    
    deltim=1.
    print 'qtimedelay is calculated with a unit of',deltim,'hours'
    
    for jtim in range(ntdh):
#        print jtim
        tfuture=jtim*deltim
#        print alamb*tfuture
        cumprob= special.gammainc(alpha, alamb*tfuture)# hoeft niet want verschil wordt genomen: /special.gamma(alpha)
#        print cumprob
        set_par['frac_future'][jtim]=max(0.,cumprob-psave)
        if set_par['frac_future'][jtim] > 0.0001:
            print set_par['frac_future'][jtim]
        psave = cumprob
    
    if cumprob < 0.99:
        print 'not enough bins in the frac_future'

    #make sure sum to one
    set_par['frac_future'][:]=set_par['frac_future'][:]/set_par['frac_future'][:].sum()
    
    return set_par

tt=qtimedelay(set_par,deltim=24)
plt.figure()
plt.plot(tt['frac_future'])


###############################################################################
##   plot the distirbution
###############################################################################    
    
mu=3.82
chi=1.0
phi=2.48
loglambda= chi*phi+mu

zeta=np.arange(mu,14.,0.01)
fzeta=np.zeros(zeta.size)

cnt=0
for i in zeta:
   temp=(i-mu)/chi
   fzeta[cnt]=1./(chi*special.gamma(phi)) * temp**(phi-1) * np.exp(-temp)
   cnt+=1
   
plt.plot(zeta,1.-fzeta.cumsum()/100.)
plt.xlabel(r'Topographic Index ($ln(\alpha / \tan \beta $)')
plt.ylabel(r'Ac/A')
  


    
    
    
    
    
    
    