# -*- coding: utf-8 -*-
"""
Created on Thu May  2 20:26:08 2019

@author: micha
"""
"""
turning this script into functions, for the purpose of 
making phase diagrams
2.36 min:
    builtin sum 37.43s
"""

def rho_only():
    import numpy as np #lint this
    L=10
    Da=0.1
    Db=10
    k0=0.067
    delta=1
    gamma=1
    K=1
    N=100 #number of steps
    Ntot=80
    a_init=Ntot/N#0.2683312 #initial values given in Mori et. al.
    asteps=L/N
    h=(asteps/(0.5*Db)) *0.8 #h<=a/(2D) -- h=0.004 here
    Rhos=np.zeros((N+1),float)
    Rhos2=np.zeros((N+1),float)
    Rhos[81:N+1]=0.2#*np.random.rand(20) #noise on the right side - causes right-polarization

    Rhos2[0]=a_init
    Rhos2[N]=a_init
    
    """modeling reaction diffusion without inhib or chem to equilibrate"""
    def react_ab(a,b):
            #this is the "reaction" part (including cooperativity)
        return b*(k0+gamma*a**2/(K**2+a**2)) - delta*a

    t=0
    tEnd1=300
    const=h*Da/asteps**2
     
    while t<tEnd1:
        Rhos[0]=(Rhos[2]+Rhos[1])/2  
        Rhos[N]=(Rhos[N-2]+Rhos[N-1])/2
        a=np.sum(Rhos)/L     
        b=(Ntot-a*L)/L

        Rhos2[1:N] = Rhos[1:N] + const*(Rhos[2:N+1]+Rhos[0:N-1]-2*Rhos[1:N]) + h*react_ab(Rhos[1:N],b)
        Rhos,Rhos2 = Rhos2,Rhos
        t += h
        
    return Rhos
        
def rho_CI(Rhos_in,beta=0.05,j=0.05,Di=0.1,omega=0.7,chem_max=0.5): #all of the parameters I'd like to test
    import numpy as np #lint this
    Rhos=Rhos_in.copy()
    L=10    #constants from above made local again
    Da=0.1
    Db=10
    k0=0.067
    delta=1
    gamma=1
    K=1
    N=100 #number of steps
    Ntot=80
    asteps=L/N
    h=(asteps/(0.5*Db)) *0.8 #h<=a/(2D) -- 
    const=h*Da/asteps**2
    const_inh=h*Di/asteps**2 
    
    def react_CI(a,b,inh,chem):
        onrate=b*(k0+gamma*a**2/(K**2+a**2) - beta*inh + omega*chem)
        onrate[onrate<0]=0
        #this is the RDE, including an inhibitor that propels Rho away and a chemoattractant
        return (onrate - delta*a)
    
    """introduction of inhibitor and chemoattractant"""
    Rhos2=np.zeros((N+1),float)
    inh=np.zeros((N+1),float)
    inh2=np.zeros((N+1),float)
    chem=np.zeros((N+1),float)
    
    #chem=np.linspace(chem_max,0,101) 
    chem=np.linspace(0,chem_max,101)
    
    Rhos2[0]=Rhos[0]; Rhos2[N]=Rhos[N]
    inh[N]=j #some arbitrary concentration to start with
    inh2[0]=0; inh2[N]=inh[N]
    t=0
    tEnd2=200 #change this to 200
    
    while t<tEnd2: #later I can change this to "while true" to model separation
        Rhos[0]=(Rhos[2]+Rhos[1])/2    
        Rhos[N]=(Rhos[N-2]+Rhos[N-1])/2
        #Rhos[Rhos<0]=0 #avoids negative values -- all active Rho (a) is gone
        inh[0]=(inh[2]+inh[1])/2
        inh[N]=j/Di + inh[N-1]
        
        a=np.sum(Rhos)/L     
        b=(Ntot-a*L)/L
        
        #diffusion and RD for inhibitor and Rho
        inh2[1:N] = inh[1:N] + const_inh*(inh[2:N+1]+inh[0:N-1]-2*inh[1:N])
        Rhos2[1:N] = Rhos[1:N] + const*(Rhos[2:N+1]+Rhos[0:N-1]-2*Rhos[1:N]) + h*react_CI(Rhos[1:N],b,inh[1:N],chem[1:N])
        #Rhos2[Rhos2<0]=0
        Rhos,Rhos2 = Rhos2,Rhos
        inh,inh2 = inh2,inh
        t += h
    return Rhos#, chem, inh

def turned_around(Rhos):
    #checks if a cell is turned around, starting from right-polarization
    N=len(Rhos)-1
    #print(Rhos)
    left=Rhos[0]
    right=Rhos[N]
    if left>right:
        return 1
    elif right>left:
        return 0
    else:
        return 2
