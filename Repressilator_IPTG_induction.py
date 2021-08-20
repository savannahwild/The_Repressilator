#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  6 23:18:21 2021

@author: savannah
"""

# ## simulate the Repressilator using an ODE model  


import numpy as np  
import matplotlib.pyplot as plt  
from scipy.integrate import odeint  
import math

# ### Define repressilator simulation


def sdot_repressilator(s,t,params):  
    
    km, km0, kdm, kp, kdp, K, n, X = params 
    m_tetR, m_lacI, m_cI, p_tetR, p_lacI, p_cI = s
        
    rate_m_tetR_prod = (km*((K**n)/((K**n)+((X*p_lacI)**n))))+km0
    rate_m_lacI_prod = (km*((K**n)/((K**n)+(p_cI**n))))+km0
    rate_m_cI_prod   = (km*((K**n)/((K**n)+(p_tetR**n))))+km0
    
    rate_p_tetR_prod = kp*m_tetR
    rate_p_lacI_prod = kp*m_lacI
    rate_p_cI_prod   = kp*m_cI
    
    rate_m_tetR_loss = kdm*m_tetR
    rate_m_lacI_loss = kdm*m_lacI
    rate_m_cI_loss   = kdm*m_cI
    
    rate_p_tetR_loss = kdp*p_tetR   
    rate_p_lacI_loss = kdp*p_lacI
    rate_p_cI_loss   = kdp*p_cI
    
    dm_tetR = rate_m_tetR_prod - rate_m_tetR_loss
    dm_lacI = rate_m_lacI_prod - rate_m_lacI_loss
    dm_cI   = rate_m_cI_prod - rate_m_cI_loss
    
    dp_tetR = rate_p_tetR_prod - rate_p_tetR_loss
    dp_lacI = rate_p_lacI_prod - rate_p_lacI_loss
    dp_cI   = rate_p_cI_prod - rate_p_cI_loss
    
    ds = [dm_tetR, dm_lacI, dm_cI, dp_tetR, dp_lacI, dp_cI]
    
    return ds  

# ### Parameter values


km = 30
km0 = 0.03
kdm = 0.3466
kp = 6.931
kdp = 0.06931
K = 40
n = 2
X_vals=[0,1]


# ### Initial conditions
p_tetR_vals=[]

m_tetR0 = 5
m_lacI0 = 0
m_cI0   = 0
p_tetR0 = 0
p_lacI0 = 0
p_cI0   = 0
s0 = [ m_tetR0, m_lacI0, m_cI0, p_tetR0, p_lacI0, p_cI0 ]

# ### Time observations


t_max = 1000
t_obs = np.linspace(0,t_max,t_max+1) 

fig = plt.figure()
ax = fig.add_subplot(1,1,1)

### Run simulation
for X in X_vals:  
    params = [ km, km0, kdm, kp, kdp, K, n, X ]
    s_obs = odeint(sdot_repressilator,s0,t_obs,args=(params,))  
    m_tetR_obs = s_obs[:,0]
    m_lacI_obs = s_obs[:,1]
    m_cI_obs = s_obs[:,2]
    
    p_lacI_obs = s_obs[:,4]
    p_cI_obs = s_obs[:,5]
    
    p_tetR_obs = s_obs[:,3]
    ax.plot(t_obs, p_tetR_obs, '-',label='X='+str(X))

ax.set_xlabel('Time (min)')
ax.set_ylabel('Proteins per cell')
ax.set_ylim(0,8000)

ax.set_title('Repressilator simulation with inhibition of tetR where X=0')

#ax.plot(t_obs, p_lacI_obs, '-',label='lacI', color='b')
#ax.plot(t_obs, p_cI_obs, '-',label='cI', color='y')
ax.legend()