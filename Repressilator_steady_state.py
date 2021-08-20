#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  6 22:45:06 2021

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

### Parameter values

km = 30
km0 = 0.03
kdm = 0.3466
kp = 6.931
kdp = 0.06931
K = 40
n = 2
params = [km, km0, kdm, kp, kdp, K, n]

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


# Time vs Proteins per cell

vals =[3,4,5]
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
legends = ['tetR', 'lacI', 'cI']
colours = ['r','b','y']

### 
ax.set_title('Repressilator ODE simulation of the steady state')



p_lacI_vals=np.linspace(0,300,100)
p_tetR_vals=[]
for p_lacI0 in p_lacI_vals:
    s0 = [m_tetR0, p_tetR0]
    params = [ km, km0, kdm, kp, kdp, K, n]
    s0 = [ m_tetR0, m_lacI0, m_cI0, p_tetR0, p_lacI0, p_cI0 ]

    #run simulation
    s_obs = odeint(sdot_repressilator,s0,t_obs,args=(params,))

    p_tetR_obs = s_obs[:,1]
    p_tetR_vals.append(p_tetR_obs[-1])

ax.plot(p_lacI_vals, p_tetR_vals)#,label=legends[i-3])
ax.set_xlabel('lacI proteins per cell')
ax.set_ylabel('tetR proteins per cell')