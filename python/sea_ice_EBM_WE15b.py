# Cecilia Bitz's python version:
#
#% This script numerically solves the diffusive energy balance model (EBM)
#% with seasonal variations, sea ice, and stochastic weather noise described
#% in eqn (1) of Wagner and Eisenman (2015b).
#%
#% The original version of this model was introduced and described in Wagner
#% and Eisenman (2015a, hereafter WE15). This version of the model has the
#% addition of stochastic weather noise.
#%
#% Till Wagner and Ian Eisenman, November 2015, wrote matlab version
#% hacked into python by Cecilia Bitz, April 2016
#% minor bug fix Jan 2022 (in eq.A1, S[:,i-1] -> S[:,i])
#%
#% References:
#% T.J.W. Wagner and I. Eisenman (2015a). How climate model complexity
#%   influences sea ice stability. J Climate 28, 3998-4014.
#% T.J.W. Wagner and I. Eisenman (2015b). False alarms: How early warning
#%   signals falsely predict abrupt sea ice loss. Geophys Res Lett 42, 103320341.
import numpy as np
import scipy
# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass
n = 200.    #%grid resolution
dur = 35    #%duration of simulation
sig = 0.5   #%noise amplitude
Fdef = 0.   #%initial Forcing level 
            #%(F = 0 corresponds roughly to pre-industrial levels)
spinup = 15  #%start ramping after 'spinup' years
#%%Model parameters (WE15, Table 1 and Section 2d) -------------------------
D = 0.6      #%diffusivity for heat transport (W m^-2 K^-1)
S1 = 338.    #%insolation seasonal dependence (W m^-2)
A = 193.     #%OLR when T = T_m (W m^-2)
B = 2.1      #%OLR temperature dependence (W m^-2 K^-1)
cw = 9.8     #%ocean mixed layer heat capacity (W yr m^-2 K^-1)
S0 = 420.    #%insolation at equator  (W m^-2)
S2 = 240.    #%insolation spatial dependence (W m^-2)
a0 = 0.7     #%ice-free co-albedo at equator
a2 = 0.1     #%ice=free co-albedo spatial dependence
ai = 0.4     #%co-albedo where there is sea ice
Fb = 4.      #%heat flux from ocean below (W m^-2)
k = 2.       #%sea ice thermal conductivity (W m^-2 K^-1)
Lf = 9.5     #%sea ice latent heat of fusion (W yr m^-3)
cg = 1e-3    #%ghost layer heat capacity(W yr m^-2 K^-1)
tau = 3e-6   #%ghost layer coupling timescale (yr)
#%%time stepping
nt = 1e3           #% time steps per year
dt = 1./nt         #% time step
dF = 1./(20.*nt)   #%%ramping rate
#%%Spatial Grid ------------------------------------------------------------
dx = 1./n
x = np.arange(dx/2., 1.+dx/2., dx)        #%native grid
#%%Diffusion Operator (WE15, Appendix A) -----------------------------------
xb = np.arange(dx, (1.0-dx)+(dx), dx)
llambda = D/dx/dx *(1.-xb*xb)
a = np.insert(-llambda,0,0)
c = np.insert(-llambda,len(llambda),0)
b = -a-c
diffop = -np.diag(b)-np.diag(c[0:int(n)-1], 1)-np.diag(a[1:int(n)], (-1))
#%%Definitions for implicit scheme on Tg
cg_tau = cg/tau
dt_tau = dt/tau
dc = dt_tau*cg_tau
kappa = (1.+dt_tau)*np.eye(n)-dt/cg*diffop

# In[2]:
#%%Seasonal forcing (WE15 eq.3)
ty = np.arange(dt/2., (1.-dt/2.)+(dt), dt)
Spart1=np.tile(S0-S2*x**2,(1,nt))     # makes 1 x 200000
Spart1=Spart1.reshape(int(nt),200)    # makes 1000x200 
Spart1=np.transpose(Spart1)           # finally 200x1000 done right
Spart2=np.tile(S1*np.cos(2*np.pi*ty),(1,int(n)))
Spart2=Spart2.reshape(200,int(nt))
            
Spart3=np.tile(x,(1,nt))
Spart3=Spart3.reshape(int(nt),200)
Spart3=np.transpose(Spart3)
S=Spart1-Spart2*Spart3

# In[3]:
#%%Further definitions
M = B+cg_tau
aw = a0-a2*x**2
kLf = k*Lf
#%%noise timeseries
sig_noise = sig/np.sqrt(dt)
print 'sig_noise'
print sig_noise
noise = sig_noise * plt.randn(1, int(dur*nt))
lp = 1./52.#%1 week 'decorrelation' time
alpha = np.exp(-dt/lp)  # this alpha differs from the one below
nalpha = np.sqrt(1.-alpha**2)
N_red = noise*0.
N_red[0] = noise[0]
for i in np.arange(1, (len(noise))-1):
    N_red[i] = alpha*N_red[i-1] + nalpha*noise[i]
print 'N_red'
print N_red[0:10]
#%%Set up output arrays, saving 100 timesteps/year
E100 = np.zeros((int(n), int(dur*100.)))
T100 = np.zeros((int(n), int(dur*100.)))

# In[4]:
#%%Initial conditions ------------------------------------------------------
T = 10.*np.ones( (int(n), 1))
Tg = T
E = cw*T
p = 0
m = 0
N = 0.  # the noise
F = Fdef
#%%run the model(see WE15_NumericIntegration.pdf)-------------------------
for years in np.arange(1, int(dur+1)):  #%%Loop over Years
    print years
    for i in np.arange(1, int(nt+1)):
        m=m+1
        if years>spinup:
            F=F+dF
            N=N_red[:,m-1]
        if ((p+1)*nt/100)==m:
            p=p+1
            E100[:,(p-1):p]=E
            T100[:,(p-1):p]=T
        alpha=aw.reshape(200,1)*(E>0)+ai*(E<0)
        C =alpha*S[:,i-1].reshape(200,1)+cg_tau*Tg-A+F+N
        #% surface temperature                                              
        T0 =  C/(M-kLf/E)                 #%WE15, eq.A3                  
        T = E/cw*(E>=0)+T0*(E<0)*(T0<0)   #%WE15, eq.9 
        #% Forward Euler on E                                               
        E = E+dt*(C-M*T+Fb)              #%WE15, eq.A2 
        #% Implicit Euler on Tg  %WE15, eq.A1 
        ondiag=(dc/(M-kLf/E)*(T0<0)*(E<0))
        ondiag=ondiag.reshape(int(n))
        thematrix=kappa-np.diag(ondiag)
        rightside=Tg+dt_tau*(E/cw*(E>=0)+
                     (ai*S[:,i].reshape(200,1)
                      -A+F+N)/(M-kLf/E)*(T0<0)*(E<0))
        
        Tg=np.linalg.solve(thematrix,rightside)
print 'Tg'
print Tg[0:4]
print Tg[-4:]

# In[5]:
#%compute ice edge
xi = np.ones([1, p])    # array of the location in x
for i in np.arange(1,p):
    if np.any((E100[:,int(i)-1]<0.)) == 1.:
        iceindexes=np.where(E100[:,int(i)-1]<0.)
        xicy=x[iceindexes]
        xi[0,int(i)-1] = xicy[0]
#%--------------------------------------------------------------------------
#%compute yearly summer and winter ice area and temperature at pole
xialt=xi.reshape(years,100)
T100alt=T100.reshape(int(n),years,100)
SIA_sept = 255.*(1.-xialt[int(spinup):,75]) #%summer ice area in million km^2
SIA_mar = 255.*(1.-xialt[int(spinup):,25]) #%winter ice area in million km^2
Tpole_sept = T100alt[-1:,int(spinup):,75]  #%summer temperature at pole
Tpole_mar = T100alt[-1:,int(spinup):,25]  #%winter temperature at pole
Fv = np.linspace(Fdef, F, int(dur-spinup))  #%yearly forcing without spinup
Tpole_sept=Tpole_sept.reshape(int(years-spinup))
Tpole_mar=Tpole_mar.reshape(int(years-spinup))

# In[10]:
print 'now plot it'
#%plot summer/winter ice areas and temperatures at the pole
#plt.figure(1)
plt.subplot(2, 1, 1)
plt.plot(Fv, SIA_sept, Fv, SIA_mar)
#plt.xlabel('F (W m^{-2)}')
plt.ylabel('Sea Ice Area (10^6 km^2)')
plt.legend(('winter', 'summer'),loc='best')
plt.subplot(2, 1, 2)
plt.plot(Fv, Tpole_sept, Fv, Tpole_mar)
plt.xlabel('F (W m^{-2)}')
plt.ylabel('Pole Temp (deg C)')
plt.show()
