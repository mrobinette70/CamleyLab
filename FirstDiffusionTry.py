# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 00:19:03 2019

@author: micha
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
L=10
Da=0.1
Db=10
k0=0.067
delta=1
gamma=1
K=1
N=100
Ntot=80 #Ntot is 22.6833 given, but larger Ntot behaves better
a_init=Ntot/N#0.2683312 #initial values given in Mori et. al.
b_init=2.0       #units in molecules per length (so multiply by 10 to get total)
asteps=L/N
h=(asteps/(0.5*Db)) *0.8 #h<=a/(2D) -- h=0.016 here

Rhos=np.zeros((N+1),float)
Rhos2=np.zeros((N+1),float)

'''potential initial Rho profiles'''
#Rhos[0:N+1]=np.linspace(0.1,0.7,101)
#Rhos[0:N+1]=np.linspace(0.7,0.1,101)
#Rhos[0:N+1]=0.4*np.random.rand(101) # all random
Rhos[0:19]=0.4*np.random.rand(19)  #noise on the left side
#Rhos[90:N+1]=1*np.random.rand(11) #noise on the right side
Rhos2[0]=a_init
Rhos2[N]=a_init

#initial profile plot
plt.figure(1)
plt.plot(Rhos,'.',label='noise')
plt.xlabel('Intracellular pos. (arb. length)')
plt.ylabel('conc. GTPase (rho)')
plt.ylim((0,1.1))
plt.title('Initial Rho profile')
plt.savefig('initialRho')

def react_ab(a,b):
    #this is the "reaction" part (including cooperativity)
    return b*(k0+gamma*a**2/(K**2+a**2)) - delta*a

t=0
count=0
Rhoverall=np.zeros((125,N+1),float)
tEnd=100
t0=1; t1=5; t2=10; t3=30; t4=50; t5=75; t6=100
epsilon=1e-8
const=h*Da/asteps**2
while t<tEnd:
    Rhos[0]=(Rhos[2]+Rhos[1])/2     
    Rhos[N]=(Rhos[N-2]+Rhos[N-1])/2
    a=sum(Rhos)/L     
    b=(Ntot-a*L)/L

    Rhos2[1:N] = Rhos[1:N] + const*(Rhos[2:N+1]+Rhos[0:N-1]-2*Rhos[1:N]) + h*react_ab(Rhos[1:N],b)
    Rhos,Rhos2 = Rhos2,Rhos
    t += h
    """
    #collecting data for heatmap
    if count%50==0:
        Rhoverall[int(count/50)]=Rhos
    count+=1
    """    
#plotting final Rho profile and profiles over time
"""    
plt.figure(2)
plt.imshow(Rhoverall,aspect='auto',cmap='Blues',extent=[0,100,200,0])
#sns.heatmap(Rhoverall,cmap='Blues',label='Rho',cbar_kws={'label':'conc. Rho'})
plt.xlabel('Intracellular pos. (arb. length)')  
plt.ylabel('time')
plt.colorbar(label='conc. Rho')
plt.savefig('Rhoheatmap')
"""
Bs=np.zeros((N+1),float)
Bs+=b
plt.figure(3)  
plt.plot(Rhos,label='active Rho')
plt.plot(0.2*Bs, label='inactive Rho')
plt.xlabel('Intracellular pos. (arb. length)')
plt.ylabel('conc. GTPase (rho)')
plt.title('Rho profile')
plt.legend()
plt.savefig('finalRho')



