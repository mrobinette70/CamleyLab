# -*- coding: utf-8 -*-
"""
Created on Thu May  2 21:16:46 2019

@author: micha
"""

"""
Here I use the functions from RhoFunctions.py to create phase diagrams for
pairs of interesting parameters
"""
import matplotlib.pyplot as plt
import numpy as np
from RhoFunctions import rho_only, rho_CI, turned_around

#parameters to test
js=np.linspace(0,0.15,25)
betas=np.linspace(0,0.04,25)
#Dis=np.linspace(0,0.5,11)
omegas=np.linspace(0,1,25)
chems_max=np.linspace(0,3,25)

Rhos_init=rho_only() #the base cell, right-polarized 
plt.plot(Rhos_init)
""" developing phase diagrams for interesting parameter combinations"""
#zmat=np.zeros([11,11],float)
jandbeta=np.zeros([25,25],float)
chemandomega=np.zeros([25,25],float)
jandchem=np.zeros([25,25],float)
#jandDi=np.zeros([11,11],float)

#rho_CI(Rhos,beta=0.05,j=0.05,Di=0.1,omega=0.7,chem_max=0.5)

for y, Js in enumerate(js):
    for x, Beta in enumerate(betas):
        newRho=rho_CI(Rhos_init,j=Js,beta=Beta)
        #print(newRho)
        yidx=len(js)-y-1
        jandbeta[yidx,x]=turned_around(newRho)
        print('good1')
        print(x)
        print(y)
plt.figure(1)
plt.imshow(jandbeta,extent=[0,0.04,0,0.15],cmap='gray',aspect='auto')
plt.xlabel('beta')
plt.ylabel('j')
plt.title('j vs beta') 
plt.colorbar()
plt.savefig('j vs beta highres')


for y, chm in enumerate(chems_max):
    for x, ome in enumerate(omegas):
        newRho=rho_CI(Rhos_init,chem_max=chm,omega=ome)
        #print(newRho)
        yidx=len(chems_max)-y-1
        chemandomega[yidx,x]=turned_around(newRho)
        print('good2')
        print(x)
        print(y)
plt.figure(2)
plt.imshow(chemandomega,extent=[0,1,0,3],cmap='gray',aspect='auto')
plt.xlabel('omega')
plt.ylabel('max chem conc.')
plt.title('Max chem vs omega')
plt.colorbar()
plt.savefig('Max chem vs omega highres')


for y, Js in enumerate(js):
    for x, chm in enumerate(chems_max):
        newRho=rho_CI(Rhos_init,j=Js,chem_max=chm)
        #print(newRho)
        yidx=len(js)-y-1
        jandchem[yidx,x]=turned_around(newRho)
        print('good3')
        print(x)
        print(y)
plt.figure(3)
plt.imshow(jandchem,extent=[0,3,0,0.15],cmap='gray',aspect='auto')
plt.xlabel('max chem conc.')
plt.ylabel('j')
plt.title('j vs Max chem') 
plt.colorbar()
plt.savefig('j vs Max chem highres')       
    


"""      
for y, Js in enumerate(js):
    for x, DI in enumerate(Dis):
        newRho=rho_CI(Rhos_init,j=Js,Di=DI)
        #print(newRho)
        jandDi[y,x]=turned_around(newRho)
        print('good4')
        print(x)
        print(y)
"""




