#!/usr/bin/env python
# coding: utf-8

# In[89]:


import numpy as np  
import matplotlib.pyplot as plt  
from scipy.integrate import odeint  
import math


# In[90]:


# Use ODE's for reaction rates of mRNA and proteins to construct repressilator function
def inducible_repressilator_GFP(s,t,params):  
    
    n, K, kdp, kp, kdm, km0, km, kdp_GFP = params
    m_cI, m_tetR, m_luxR, p_cI, p_tetR, p_luxR, m_GFP, p_GFP = s
    
    if AHL == True:
        rate_m_cI_prod = km * K**n / (K**n + p_luxR**n) + km0
    elif AHL == False:
        rate_m_cI_prod = km + km0
    rate_m_tetR_prod = km * K**n / (K**n + p_cI**n) + km0
    rate_m_luxR_prod = km * K**n / (K**n + p_tetR**n) + km0
    rate_m_GFP_prod  = km * K**n / (K**n + p_tetR**n) + km0
        
    rate_p_cI_prod   = kp * m_cI
    rate_p_tetR_prod = kp * m_tetR
    rate_p_luxR_prod = kp * m_luxR
    rate_p_GFP_prod  = kp * m_GFP
    
    rate_m_cI_loss   = kdm * m_cI
    rate_m_tetR_loss = kdm * m_tetR
    rate_m_luxR_loss = kdm * m_luxR
    rate_m_GFP_loss  = kdm * m_GFP
    
    rate_p_cI_loss   = kdp * p_cI
    rate_p_tetR_loss = kdp * p_tetR
    rate_p_luxR_loss = kdp * p_luxR
    rate_p_GFP_loss  = kdp_GFP * p_GFP
    
    dm_cI   = rate_m_cI_prod - rate_m_cI_loss
    dm_tetR = rate_m_tetR_prod - rate_m_tetR_loss
    dm_luxR = rate_m_luxR_prod - rate_m_luxR_loss
    
    dp_cI   = rate_p_cI_prod - rate_p_cI_loss
    dp_tetR = rate_p_tetR_prod - rate_p_tetR_loss
    dp_luxR = rate_p_luxR_prod - rate_p_luxR_loss
    
    dm_GFP  = rate_m_GFP_prod - rate_m_GFP_loss
    dp_GFP  = rate_p_GFP_prod - rate_p_GFP_loss
    
    ds = [ dm_cI, dm_tetR, dm_luxR, dp_cI, dp_tetR, dp_luxR, dm_GFP, dp_GFP]
    
    return ds  


# In[95]:


# Figure 2: Inducible repressilatory oscillations
n = 2
K = 40
kdp = 0.06931
kp = 6.931
kdm = 0.3466
km0 = 0.03
km = 30
kdp_GFP = 0.01155

params = [ n, K, kdp, kp, kdm, km0, km, kdp_GFP ]

# set time observations
t_max = 1000
t_obs = np.linspace(0,t_max,t_max+1) 

# create figure
fig = plt.figure(figsize=(12,6))
fig.subplots_adjust(hspace=0.3, wspace=0.3)
ax.set_ylim([0.1,1000])

for i in range(1, 3):
    ax = fig.add_subplot(1, 2, i)
    
    if i == 1:
        AHL = True
        # Initial conditions
        m_cI0   = 5
        m_tetR0 = 0
        m_luxR0 = 0
        p_cI0   = 0
        p_tetR0 = 0
        p_luxR0 = 0
        m_GFP0  = 0 
        p_GFP0  = 0 

        s0 = [ m_cI0, m_tetR0, m_luxR0, m_GFP0, p_cI0, p_tetR0, p_luxR0, p_GFP0 ]
        s_obs = odeint(inducible_repressilator_GFP,s0,t_obs,args=(params,))  

        m_cI_obs = s_obs[:,0]
        m_tetR_obs = s_obs[:,1]
        m_luxR_obs = s_obs[:,2]
        p_cI_obs = s_obs[:,3]
        p_tetR_obs = s_obs[:,4]
        p_luxR_obs = s_obs[:,5]
        m_GFP_obs =  s_obs[:,6]
        p_GFP_obs =  s_obs[:,7]

        # plot figure
        ax.plot(t_obs, p_cI_obs, color = 'blue', label=f"= cI")
        ax.plot(t_obs, p_tetR_obs, color = 'yellow', label=f"= tetR")
        ax.plot(t_obs, p_luxR_obs, color = 'red', label=f"= luxR")
        ax.plot(t_obs, p_GFP_obs, color = 'green', label=f"= GFP")
        ax.set(xlabel = 'Time (mins)')
        ax.set(ylabel = 'Protein copies per cell')
        ax.legend()
        plt.grid()
        
    elif i == 2:
        AHL = False
        # Initial conditions
        m_cI0   = 5
        m_tetR0 = 0
        m_luxR0 = 0
        p_cI0   = 0
        p_tetR0 = 0
        p_luxR0 = 0
        m_GFP0  = 0 
        p_GFP0  = 0 

        s0 = [ m_cI0, m_tetR0, m_luxR0, m_GFP0, p_cI0, p_tetR0, p_luxR0, p_GFP0 ]
        s_obs = odeint(inducible_repressilator_GFP,s0,t_obs,args=(params,))  

        m_cI_obs = s_obs[:,0]
        m_tetR_obs = s_obs[:,1]
        m_luxR_obs = s_obs[:,2]
        p_cI_obs = s_obs[:,3]
        p_tetR_obs = s_obs[:,4]
        p_luxR_obs = s_obs[:,5]
        m_GFP_obs =  s_obs[:,6]
        p_GFP_obs =  s_obs[:,7]

        # plot figure
        ax.plot(t_obs, p_cI_obs, color = 'blue', label=f"= cI")
        ax.plot(t_obs, p_tetR_obs, color = 'yellow', label=f"= tetR")
        ax.plot(t_obs, p_luxR_obs, color = 'red', label=f"= luxR")
        ax.plot(t_obs, p_GFP_obs, color = 'green', label=f"= GFP")
        ax.set(xlabel = 'Time (mins)')
        ax.set(ylabel = 'Protein copies per cell')
        ax.legend()
        plt.grid()
plt.show()
plt.savefig('Figure 2')


# In[97]:


# Figure 3 - Hill coefficient variation

AHL = True
K = 40
kdp = 0.06931
kp = 6.931
kdm = 0.3466
km0 = 0.03
km = 30
kdp_GFP = 0.01155

#intitial condtions
m_cI0 = 0
m_tetR0 = 5
m_luxR0   = 0
p_cI0 = 0
p_tetR0 = 0
p_luxR0   = 0
m_GFP0  = 0 
p_GFP0  = 0 

# set time observations
t_max = 1000
t_obs = np.linspace(0,t_max,t_max+1) 

# create figure
fig = plt.figure(figsize=(13,9))
fig.subplots_adjust(hspace=0.2, wspace=0.5)
ax.set_ylim([0,1000])

for i in range(1, 7):
    ax = fig.add_subplot(2, 3, i)

    if i == 1:
        n = 0.5
        params = [ n, K, kdp, kp, kdm, km0, km, kdp_GFP ]
        s0 = [ m_cI0, m_tetR0, m_luxR0, m_GFP0, p_cI0, p_tetR0, p_luxR0, p_GFP0 ]
        s_obs = odeint(inducible_repressilator_GFP,s0,t_obs,args=(params,))  

        m_cI_obs = s_obs[:,0]
        m_tetR_obs = s_obs[:,1]
        m_luxR_obs = s_obs[:,2]
        p_cI_obs = s_obs[:,3]
        p_tetR_obs = s_obs[:,4]
        p_luxR_obs = s_obs[:,5]
        m_GFP_obs =  s_obs[:,6]
        p_GFP_obs =  s_obs[:,7]

        # plot figure
        ax.plot(t_obs, p_cI_obs, color = 'blue', label=f"= cI")
        ax.plot(t_obs, p_tetR_obs, color = 'yellow', label=f"= tetR")
        ax.plot(t_obs, p_luxR_obs, color = 'red', label=f"= luxR")
        ax.plot(t_obs, p_GFP_obs, color = 'green', label=f"= GFP")
        ax.set(xlabel = 'Time (mins)')
        ax.set(ylabel = 'Protein copies per cell')
        ax.legend()
        plt.grid()
        
    elif i == 2:
        n = 1
        params = [ n, K, kdp, kp, kdm, km0, km, kdp_GFP ]
        s0 = [ m_cI0, m_tetR0, m_luxR0, m_GFP0, p_cI0, p_tetR0, p_luxR0, p_GFP0 ]
        s_obs = odeint(inducible_repressilator_GFP,s0,t_obs,args=(params,))  

        m_cI_obs = s_obs[:,0]
        m_tetR_obs = s_obs[:,1]
        m_luxR_obs = s_obs[:,2]
        p_cI_obs = s_obs[:,3]
        p_tetR_obs = s_obs[:,4]
        p_luxR_obs = s_obs[:,5]
        m_GFP_obs =  s_obs[:,6]
        p_GFP_obs =  s_obs[:,7]

        # plot figure
        ax.plot(t_obs, p_cI_obs, color = 'blue', label=f"= cI")
        ax.plot(t_obs, p_tetR_obs, color = 'yellow', label=f"= tetR")
        ax.plot(t_obs, p_luxR_obs, color = 'red', label=f"= luxR")
        ax.plot(t_obs, p_GFP_obs, color = 'green', label=f"= GFP")
        ax.set(xlabel = 'Time (mins)')
        ax.set(ylabel = 'Protein copies per cell')
        ax.legend()
        plt.grid()
        
    elif i == 3:
        n = 1.5
        params = [ n, K, kdp, kp, kdm, km0, km, kdp_GFP ]
        s0 = [ m_cI0, m_tetR0, m_luxR0, m_GFP0, p_cI0, p_tetR0, p_luxR0, p_GFP0 ]
        s_obs = odeint(inducible_repressilator_GFP,s0,t_obs,args=(params,))  

        m_cI_obs = s_obs[:,0]
        m_tetR_obs = s_obs[:,1]
        m_luxR_obs = s_obs[:,2]
        p_cI_obs = s_obs[:,3]
        p_tetR_obs = s_obs[:,4]
        p_luxR_obs = s_obs[:,5]
        m_GFP_obs =  s_obs[:,6]
        p_GFP_obs =  s_obs[:,7]

        # plot figure
        ax.plot(t_obs, p_cI_obs, color = 'blue', label=f"= cI")
        ax.plot(t_obs, p_tetR_obs, color = 'yellow', label=f"= tetR")
        ax.plot(t_obs, p_luxR_obs, color = 'red', label=f"= luxR")
        ax.plot(t_obs, p_GFP_obs, color = 'green', label=f"= GFP")
        ax.set(xlabel = 'Time (mins)')
        ax.set(ylabel = 'Protein copies per cell')
        ax.legend()
        plt.grid()
        
    elif i == 4:
        n = 2
        params = [ n, K, kdp, kp, kdm, km0, km, kdp_GFP ]
        s0 = [ m_cI0, m_tetR0, m_luxR0, m_GFP0, p_cI0, p_tetR0, p_luxR0, p_GFP0 ]
        s_obs = odeint(inducible_repressilator_GFP,s0,t_obs,args=(params,))  

        m_cI_obs = s_obs[:,0]
        m_tetR_obs = s_obs[:,1]
        m_luxR_obs = s_obs[:,2]
        p_cI_obs = s_obs[:,3]
        p_tetR_obs = s_obs[:,4]
        p_luxR_obs = s_obs[:,5]
        m_GFP_obs =  s_obs[:,6]
        p_GFP_obs =  s_obs[:,7]

        # plot figure
        ax.plot(t_obs, p_cI_obs, color = 'blue', label=f"= cI")
        ax.plot(t_obs, p_tetR_obs, color = 'yellow', label=f"= tetR")
        ax.plot(t_obs, p_luxR_obs, color = 'red', label=f"= luxR")
        ax.plot(t_obs, p_GFP_obs, color = 'green', label=f"= GFP")
        ax.set(xlabel = 'Time (mins)')
        ax.set(ylabel = 'Protein copies per cell')
        ax.legend()
        plt.grid()
    
    elif i == 5:
        n = 2.5
        params = [ n, K, kdp, kp, kdm, km0, km, kdp_GFP ]
        s0 = [ m_cI0, m_tetR0, m_luxR0, m_GFP0, p_cI0, p_tetR0, p_luxR0, p_GFP0 ]
        s_obs = odeint(inducible_repressilator_GFP,s0,t_obs,args=(params,))  

        m_cI_obs = s_obs[:,0]
        m_tetR_obs = s_obs[:,1]
        m_luxR_obs = s_obs[:,2]
        p_cI_obs = s_obs[:,3]
        p_tetR_obs = s_obs[:,4]
        p_luxR_obs = s_obs[:,5]
        m_GFP_obs =  s_obs[:,6]
        p_GFP_obs =  s_obs[:,7]

        # plot figure
        ax.plot(t_obs, p_cI_obs, color = 'blue', label=f"= cI")
        ax.plot(t_obs, p_tetR_obs, color = 'yellow', label=f"= tetR")
        ax.plot(t_obs, p_luxR_obs, color = 'red', label=f"= luxR")
        ax.plot(t_obs, p_GFP_obs, color = 'green', label=f"= GFP")
        ax.set(xlabel = 'Time (mins)')
        ax.set(ylabel = 'Protein copies per cell')
        ax.legend()
        plt.grid()
        
    elif i == 6:
        n = 3
        params = [ n, K, kdp, kp, kdm, km0, km, kdp_GFP ]
        s0 = [ m_cI0, m_tetR0, m_luxR0, m_GFP0, p_cI0, p_tetR0, p_luxR0, p_GFP0 ]
        s_obs = odeint(inducible_repressilator_GFP,s0,t_obs,args=(params,))  

        m_cI_obs = s_obs[:,0]
        m_tetR_obs = s_obs[:,1]
        m_luxR_obs = s_obs[:,2]
        p_cI_obs = s_obs[:,3]
        p_tetR_obs = s_obs[:,4]
        p_luxR_obs = s_obs[:,5]
        m_GFP_obs =  s_obs[:,6]
        p_GFP_obs =  s_obs[:,7]

        # plot figure
        ax.plot(t_obs, p_cI_obs, color = 'blue', label=f"= cI")
        ax.plot(t_obs, p_tetR_obs, color = 'yellow', label=f"= tetR")
        ax.plot(t_obs, p_luxR_obs, color = 'red', label=f"= luxR")
        ax.plot(t_obs, p_GFP_obs, color = 'green', label=f"= GFP")
        ax.set(xlabel = 'Time (mins)')
        ax.set(ylabel = 'Protein copies per cell')
        ax.legend()
        plt.grid()
plt.show()
fig.savefig('Figure 3')


# In[98]:


# Figure 4: Bifurcation plot with different Hill coefficient (n) values

AHL = True
km = 30
km0 = 0.03
kdm = 0.3466
kp = 6.931
kdp = 0.06931
K = 40

#intitial condtions
m_cI0 = 0
m_tetR0 = 5
m_luxR0   = 0
p_cI0 = 0
p_tetR0 = 0
p_luxR0   = 0
m_GFP0  = 0 
p_GFP0  = 0

fig3 = plt.figure()
ax = fig3.add_subplot(1,1,1)
ax.set(xlabel = 'Value of n')
ax.set(ylabel = 'Copies of luxR per cell')
ax.set_yscale('log')

s0 = [ m_cI0, m_tetR0, m_luxR0, m_GFP0, p_cI0, p_tetR0, p_luxR0, p_GFP0 ]

# create array of n values
n_vals = np.linspace(0,4,100)

# set time observations
t_max = 10000
t_obs = np.linspace(0,t_max,t_max+1) 
p_luxR_max_vals = []
p_luxR_min_vals = []

for n in n_vals:

    params = [ n, K, kdp, kp, kdm, km0, km, kdp_GFP ]

    # Run simulation
    s_obs = odeint(inducible_repressilator_GFP,s0,t_obs,args=(params,))  

    m_cI_obs = s_obs[:,0]
    m_tetR_obs = s_obs[:,1]
    m_luxR_obs = s_obs[:,2]
    p_cI_obs = s_obs[:,3]
    p_tetR_obs = s_obs[:,4]
    p_luxR_obs = s_obs[:,5]
    m_GFP_obs =  s_obs[:,6]
    p_GFP_obs =  s_obs[:,7]

    # extracting maximum and minimum values at the 600th minute
    osc_max = np.max(s_obs[-600:,3])
    osc_min = np.min(s_obs[-600:,3])
    
    p_luxR_max_vals.append(osc_max)
    p_luxR_min_vals.append(osc_min)
    
ax.plot(n_vals, p_luxR_max_vals, color = 'orange', label=f"= Max values")
ax.plot(n_vals, p_luxR_min_vals, color = 'green', label=f"= Min values")
ax.legend()
plt.grid()
plt.show()            
fig3.savefig('Figure 4')

