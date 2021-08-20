#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  6 22:56:34 2021

@author: savannah
"""

# ## simulate the Repressilator using an ODE model  


import numpy as np  
import matplotlib.pyplot as plt  
from scipy.integrate import odeint  
import math

# ### Define repressilator simulation


def sdot_repressilator(s,t,params):  
    
    km, km0, kdm, kp, kdp, K, n = params 
    m_tetR, m_lacI, m_cI, p_tetR, p_lacI, p_cI = s
        
    rate_m_tetR_prod = (km*((K**n)/((K**n)+(p_lacI**n))))+km0
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
params = [ km, km0, kdm, kp, kdp, K, n ]

# ### Initial conditions


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
# run simulation
s_obs = odeint(sdot_repressilator,s0,t_obs,args=(params,))

p_tetR_obs = s_obs[:,3]
p_lacI_obs = s_obs[:,4]
              
fig = plt.figure(figsize=(10,10))
axP = fig.add_subplot(1,1,1)
axP.plot(p_lacI_obs, p_tetR_obs, 'b-')
axP.set_xlabel('[lacI]')
axP.set_ylabel('[tetR]')

# plot points to show correspondance betwen time evolution and phase plots
max_I_i = np.argmax(p_lacI_obs)
min_I_i = np.argmin(p_lacI_obs)
max_S_i = np.argmax(p_tetR_obs)
min_S_i = np.argmin(p_tetR_obs)

axP.plot(p_lacI_obs[max_I_i], p_tetR_obs[max_I_i],  'ro', label='maximum [lacI]')
axP.plot(p_lacI_obs[min_I_i], p_tetR_obs[min_I_i],  'r^')
axP.plot(p_lacI_obs[max_S_i], p_tetR_obs[max_S_i],  'bo', label='maximum [tetR]')
axP.plot(p_lacI_obs[min_S_i], p_tetR_obs[min_S_i],  'b^')

# plot points to show steady states,
axP.plot(0, 0, 'ko', label='steady states')
axP.plot(p_lacI_obs[min_I_i], p_lacI_obs[min_I_i],'ko')

# plot lines showing null clines,
axP.axvline(p_lacI_obs[min_I_i],ls=':',c='k', label='dp_tetR/dt=0 nullcine') 
axP.axhline(0,ls=':',c='k') 
axP.axhline(p_lacI_obs[min_I_i],ls=':',c='r', label='dp_lacI=0 nullcine') 
axP.axvline(0,ls=':',c='r') 

axP.legend()
axP.set_title('Phase plot of repressilator model')