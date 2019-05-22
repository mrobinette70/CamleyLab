# -*- coding: utf-8 -*-
"""
this program monitors Rho (active) concentration in the presence 
of a diffusing chemoattractant (chem here)

Things to try:
. Allow the chemoattractant to diffuse over more space; its environment is
    300 long, only look at the middle 100 positions
. find how to determine velocity of a cell
    as the cell moves, the position of the chemoattractant wrt. the cell
    will change!
. find how to monitor collision mechanics between cells

Arbitrary choices here:
    Initial concentration of chemoattractant and starting positions
    Diffusion coefficient of chemoattractant
    Chemoattractant bias toward Rho binding to membrane (in RD equation)
    
4/26:
    Changed the chemoattractant to be nondiffusing and steady
        This still gives the sigmoidal pattern
    Need to:
        Add in inhibitor as another cell comes in contact (start without chemoattractants
        with FirstDiffusionTry.py)
        Then include chemoattractants as in this file
"""

import numpy as np
import matplotlib.pyplot as plt
"""given parameters"""
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
asteps=L/N
h=asteps/(2*Db) *0.8 #h<=a/(2D) -- h=0.004 here

"""initial arrays"""
Rhos=np.zeros((N+1),float)
Rhos2=np.zeros((N+1),float)
chem=np.zeros((N+1),float)
chem2=np.zeros((N+1),float) #it might be good to have a bigger env


#to start with all cytosolic, do not choose any of these Rhos[] arrays

#Rhos[0:N+1]=np.linspace(0.1,0.7,101) #gradient
#Rhos[0:N+1]=np.linspace(0.7,0.1,101) #gradient
Rhos[0:N+1]=0.4*np.random.rand(101)  # all random
#Rhos[0:19]=0.4*np.random.rand(19)    #noise on the left side
#Rhos[80:N+1]=0.4*np.random.rand(21)  #noise on the right side
#chem[0:N+1]=0.1*np.random.rand(101)  #chems mixed
#chem[0:10]=0.6*np.random.rand(10)     #chems on the left side
#chem[0:N+1]=np.linspace(0.5,0,101)
chem[0:19]=np.linspace(0.1,0,19) #steady state, constant flow in
Rhos2[0]=a_init
Rhos2[N]=a_init
plt.plot(Rhos,'.',label='noise')
plt.plot(chem,'o',label='10*chems')

def react_ab(a,b):
    #this is the "reaction" part (including cooperativity)
    return b*(k0+gamma*a**2/(K**2+a**2)) - delta*a

"""reaction-diffusion modeling"""
t=0
tEnd=100
t0=1; t1=5; t2=10; t3=30; t4=50; t5=75
epsilon=1e-8
const=h*Da/asteps**2
while t<tEnd:
    Rhos[0]=(Rhos[2]+Rhos[1])/2     #I should change these two boundaries back when done
    Rhos[N]=(Rhos[N-2]+Rhos[N-1])/2
    #chem[0]=(chem[2]+chem[1])/2
    #chem[N]=(chem[N-2]+chem[N-1])/2
    
    a=sum(Rhos)/L     #number active/L  -- why is this blowing up?
    b=(Ntot-a*L)/L
    
    #chem2[1:N]=chem[1:N] + 0.5*const*(chem[2:N+1]+chem[0:N-1]-2*chem[1:N])
    Rhos2[1:N] = Rhos[1:N] + const*(Rhos[2:N+1]+Rhos[0:N-1]-2*Rhos[1:N]) + 0.02*chem[1:N]+ h*react_ab(Rhos[1:N],b)
    #chem,chem2 = chem2,chem
    Rhos,Rhos2 = Rhos2,Rhos
    t += h
    
    if abs(t-t0)<epsilon:
        plt.plot(Rhos,label='t=1s')
    elif abs(t-t1)<epsilon:
        plt.plot(Rhos,label='t=5s')
    elif abs(t-t2)<epsilon:
        plt.plot(Rhos,label='t=10s')
        plt.plot(chem,'*',label='10*chems (10s)')
    elif abs(t-t3)<epsilon:
        plt.plot(Rhos,label='t=30s')
    elif abs(t-t4)<epsilon:
        plt.plot(Rhos,label='t=50s')
    elif abs(t-t5)<epsilon:
        plt.plot(Rhos,label='t=75s')
       
plt.plot(Rhos,label='t=100s')
#plt.plot(chem,label='10*chem (100s)')
plt.xlabel('Intracellular pos. (arb. length)')
plt.ylabel('conc. GTPase (rho)')
plt.legend()
#plt.ylim([0,1.2])
print((a+b)*L)
 
