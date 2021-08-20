#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  6 21:41:36 2021

@author: savannah
"""

# ## Model of repression: 
# p_lacI --| p_tetR


import numpy as np  
import matplotlib.pyplot as plt  
from scipy.integrate import odeint  
import math

# at constant p_lacI levels

# ### Define simulation to return steady states reached




def sdot_repression(s,t,params):  
    
    km, km0, kdm, kp, kdp, K, n, p_lacI = params 
    m_tetR, p_tetR = s
    
    rate_m_tetR_prod =  (km*((K**n)/((K**n)+(p_lacI**n))))+km0  
    rate_p_tetR_prod = kp*m_tetR

    rate_m_tetR_loss = kdm*m_tetR
    rate_p_tetR_loss = kdp*p_tetR

    dm_tetR = rate_m_tetR_prod - rate_m_tetR_loss
    dp_tetR = rate_p_tetR_prod - rate_p_tetR_loss

    
    ds = [ dm_tetR, dp_tetR ]
    
    return ds  

# ### Parameter values



km = 30
km0 = 0.03
kdm = 0.3466
kp = 6.931
kdp = 0.06931
K = 40
n = 2
p_lacI_vals = np.linspace(0,300,100)

# ### Initial conditions



m_tetR0 = 0
p_tetR0 = 0

s0 = [m_tetR0, p_tetR0]

# 
# ### time observations

t_max = 1000
t_obs = np.linspace(0,t_max,t_max+1) 

# ### Run simulation

p_tetR_vals=[]
for p_lacI  in p_lacI_vals:
    #set up parameters
    params = [km, km0, kdm, kp, kdp, K, n, p_lacI]

    # run simulation
    s_obs = odeint(sdot_repression,s0,t_obs,args=(params,))  
    m_tetR_obs = s_obs[:,0]
    p_tetR_obs = s_obs[:,1]

    # extract and stroe final protein level
    p_tetR_vals.append(p_tetR_obs[-1])

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.set_title('Repression kinetics plot')
ax.set_xlabel('lacI proteins per cell')
ax.set_ylabel('tetR proteins per cell')
ax.plot(p_lacI_vals, p_tetR_vals)
fig.show()