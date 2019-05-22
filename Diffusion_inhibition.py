# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 20:22:15 2019

@author: micha
"""

import numpy as np
import matplotlib.pyplot as plt
""" given parameters and introduction of relevant arrays"""
L=10
Da=0.1
Db=10
k0=0.067
delta=1
gamma=1
K=1
N=100
#Ntot=22.6833 #total amount of proteins (active and inactive)
Ntot=80
a_init=Ntot/N#0.2683312 #initial values given in Mori et. al.
b_init=2.0       #units in molecules per length (so multiply by 10 to get total)


#N=100 #number of steps
asteps=L/N
h=asteps/(2*Db) *0.8 #h<=a/(2D) -- h=0.004 here
Rhos=np.zeros((N+1),float)
Rhos2=np.zeros((N+1),float)
"""different options for starting Rho profile"""
#Rhos[0:N+1]=np.linspace(0.1,0.7,101)
#Rhos[0:N+1]=np.linspace(0.7,0.1,101)
#Rhos[0:N+1]=0.4*np.random.rand(101) # all random
#Rhos[0:19]=0.4*np.random.rand(19)  #noise on the left side
Rhos[80:N+1]=0.4*np.random.rand(21) #noise on the right side

Rhos2[0]=a_init
Rhos2[N]=a_init
#plt.plot(Rhos,'.',label='noise')
"""modeling reaction diffusion without inhib to equilibrate"""
def react_ab(a,b):
    #this is the "reaction" part (including cooperativity)
    return b*(k0+gamma*a**2/(K**2+a**2)) - delta*a

t=0
tEnd1=200
t0=1; t1=5; t2=10; t3=30; t4=50; t5=75; t6=100
epsilon=1e-8
const=h*Da/asteps**2
while t<tEnd1:
    Rhos[0]=(Rhos[2]+Rhos[1])/2  
    Rhos[N]=(Rhos[N-2]+Rhos[N-1])/2
    a=sum(Rhos)/L     
    b=(Ntot-a*L)/L

    Rhos2[1:N] = Rhos[1:N] + const*(Rhos[2:N+1]+Rhos[0:N-1]-2*Rhos[1:N]) + h*react_ab(Rhos[1:N],b)
    Rhos,Rhos2 = Rhos2,Rhos
    t += h
  
plt.plot(Rhos,label='t=200s (before inhib)')
"""allow to equibrilate here"""

"""redo the reaction diffusion in the presence of an inhibitor (mimicking 
contact w/ another cell after equilibrium)"""

beta=0.1
j=0.05 
Di=0.1
const_inh=h*Di/asteps**2 
def react_inhib(a,b,inh):
    #this is the RDE, including an inhibitor that propels Rho away
    return b*(k0+gamma*a**2/(K**2+a**2) - beta*inh) - delta*a

"""introduction of inhibitor"""
inh=np.zeros((N+1),float)
inh2=np.zeros((N+1),float)
inh[N]=0.5 #some arbitrary concentration to start with
inh2[N]=inh[N]
tEnd2=400
t0=201;t1=205;t2=210;t3=225;t4=250;t5=275;t6=300
while t<tEnd2: #later I can change this to "while true" to model separation
    Rhos[0]=(Rhos[2]+Rhos[1])/2    
    Rhos[N]=(Rhos[N-2]+Rhos[N-1])/2
    inh[0]=(inh[2]+inh[1])/2
    if inh[N]<j:
        inh[N]=j
    a=sum(Rhos)/L     
    b=(Ntot-a*L)/L

    inh2[1:N] = inh[1:N] + const_inh*(inh[2:N+1]+inh[0:N-1]-2*inh[1:N])
    Rhos2[1:N] = Rhos[1:N] + const*(Rhos[2:N+1]+Rhos[0:N-1]-2*Rhos[1:N]) + h*react_inhib(Rhos[1:N],b,inh[1:N])
    Rhos,Rhos2 = Rhos2,Rhos
    inh,inh2 = inh2,inh
    t += h
    
    """
    if abs(t-t0)<epsilon:
        plt.plot(Rhos,label='t=201s')
    elif abs(t-t1)<epsilon:
        plt.plot(Rhos,label='t=205s')
    elif abs(t-t2)<epsilon:
        plt.plot(Rhos,label='t=210s')
    #elif abs(t-t3)<epsilon:
        #plt.plot(Rhos,label='t=225s')
    elif abs(t-t4)<epsilon:
        plt.plot(Rhos,label='t=250s')
    elif abs(t-t5)<epsilon:
        plt.plot(Rhos,label='t=275s')
    elif abs(t-t6)<epsilon:
        plt.plot(Rhos,label='t=300s')
    """

    if abs(t-t0)<epsilon:
        plt.plot(inh,label='t=201s')
    elif abs(t-t1)<epsilon:
        plt.plot(inh,label='t=205s')
    elif abs(t-t2)<epsilon:
        plt.plot(inh,label='t=210s')
    elif abs(t-t3)<epsilon:
        plt.plot(inh,label='t=225s')
    elif abs(t-t4)<epsilon:
        plt.plot(inh,label='t=250s')
    elif abs(t-t5)<epsilon:
        plt.plot(inh,label='t=275s')
    elif abs(t-t6)<epsilon:
        plt.plot(inh,label='t=300s')
    


plt.plot(Rhos,label='t=400s (Rhos)')
plt.plot(inh,label='final inhib')

plt.xlabel('Intracellular pos. (arb. length)')
plt.ylabel('conc. GTPase (rho)')
plt.title('Inhibitor only')
plt.legend() #(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.savefig('inhibs')
print((a+b)*L)
