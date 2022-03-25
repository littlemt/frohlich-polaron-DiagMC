#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 30 19:56:57 2022

@author: Milo Koyama Belmont
"""

import numpy as np

nrandG=np.random.default_rng()
####  for the purposes of calculations hbar=1
#### what is mu? can it be 0? it says tuning paramtere, so i assume constant. 


#H_el = sum_k (hbar k )^2/(2 m ) aDag a
#H_ph = hbar omega_ph Sum( bdag b)
#H _el-ph = Sum_k,q( V(q)(bdag-b)adag a)
#V(q)=i hbar omaga_ph/q (4 pi alpha /V)^.5 (h/(2m pmega_ph)^.25)

#predefine a momentum complex time space

def inBetween(a,b,c,d):
     '''
     takes 2 ordered pairs (a,b),(c,d) and checks if any of the vertex are 
     contained in the other ordered pair. 
 
     Parameters
     ----------
     a : Float
         lower of the ab ordered pair.
     b : Float
         DESCRIPTION.
     c : Float
         DESCRIPTION.
     d : Float
         DESCRIPTION.
 
     Returns
     -------
     bool
         DESCRIPTION.
 
     '''
     w=c<a<d
     x=c<b<d
     y=a<c<b
     z=a<d<b
     if x!=w ^ y!=z:
         return True
     else:
         return False
     
def r(pins,prem,tauR,tauL,tauOne,tauTwo,q,k,mu,alpha,order,m=1,omegaPH=1):
    
    wX=-greensFunction(k, tauR-tauL, mu)
    alphaT=2*2**.5*np.pi*alpha
    wY=-alphaT**2*greensFunction(k, tauOne-tauR, mu)*\
        greensFunction(k,tauR-tauTwo,mu)*greensFunction(k-q,tauTwo-tauOne,mu)\
            *phononProp(tauTwo-tauOne)/q**2/(2*np.pi)**3
    pXY=pins/(tauR-tauL)*omegaPH*np.exp(-omegaPH*(tauTwo-tauOne))\
        *np.exp(-q**2/(2*m)*(tauTwo-tauOne))/2*np.pi*m/(tauTwo-tauOne)**(3/2)
    pYX=prem/(order+1)
    #what are pins prem
    return wY*pYX/(wX*pXY)

def phononProp(tau,omegaPH=1):
    return np.exp(-omegaPH*(tau))

def greensFunction(k,tau,mu,m=1):
    '''
    returns the greens function of the Hamiltonian 
    for this case function is 
    G(k,tau)=-Theta(tau)e^(-(eps_k-mu)tau)

    Parameters
    ----------
    k : TYPE
        DESCRIPTION.
    tau : TYPE
        Complex Time parameter.
    m : TYPE
        mass of the phonon.
    mu : TYPE
        Energy shift used as a tuning parameter.

    Returns
    -------
    int
        DESCRIPTION.

    '''
    epsilonK=k**2/(2*m)
    if tau<=0:
        return 0
    else:
        return -np.exp(-tau*(epsilonK-mu)) 
        #need to define Z_k ask about this, something with partition function

class diagMC:
    def timeSpaceLin(tauMax,Ntau,j):
        '''
        
    
        Parameters
        ----------
        tauMax : TYPE
            DESCRIPTION.
        Ntau : TYPE
            DESCRIPTION.
        j : TYPE
            DESCRIPTION.
    
        Returns
        -------
        TYPE
            DESCRIPTION.
    
        '''
        deltaTau=tauMax/Ntau
        return (j+1/2)*deltaTau
    
    def momSpaceLin(kMax,Nk,j):
        deltaK=kMax/Nk
        return j*deltaK
    
    def Z_k(k):
        #residue of the particle k
        return
        
                
    #should also do for nonlinear spaces but paper says its does not matter
    
            
    
    
    
    def normalise():
        #norm is sum over all momenta of the interal of the 0 order Green func
        return
    
        
    def pChange ():
        #only when expansion is of order 0
        #select new momentum 
        #accept only if metroplis with min of [1,r] where r is G_0(p new,tau)/G_)(p old,tau)
        return
        
    def tauChange(p,tauMax,mu,m=1):
        '''
        only acepts if order=0
        selesct new tau using exponential distrobution
        dispersion eta_p=epsilon_p-mu 
        mu is raondom number [0,1]
        tau = -ln (mu)/abs(eta_p)
        accept only if tau< tau max
    
        Parameters
        ----------
        
        p : Float
            External momentum.
        tauMax : Float
            External time max.
        m : Float, optional
            Mass of the particle. The default is 1.
    
        Returns
        -------
        Float
            new time.
    
        '''
        
    
        
        #tau=nrandG.exponential(scale=beta)
        u=nrandG.uniform(low=0,high=1)
        epsilonP=p**2/(2*m)
        eta_p=epsilonP-mu
        
        tauPrime=-np.log(u)/abs(eta_p)
        if tauPrime<tauMax:
            return tauPrime,1
        else:
            return tauMax,0
        
    def genWorldLine(kMax,tauMax):
        k=1
        #what is list of allowable k?
        qList=[0,tauMax]
        #needs to genrate k, tauMax
        return k,qList
        
    def insertProp(qList,tauMax,k,m,mu,alpha,order,pIns,pRem,omegaPH=1):
        '''
        

        Parameters
        ----------
        qList : TYPE
            DESCRIPTION.
        tauMax : TYPE
            DESCRIPTION.
        k : TYPE
            DESCRIPTION.
        m : TYPE
            DESCRIPTION.
        mu : TYPE
            DESCRIPTION.
        alpha : TYPE
            DESCRIPTION.
        order : TYPE
            DESCRIPTION.
        pIns : TYPE
            DESCRIPTION.
        pRem : TYPE
            DESCRIPTION.
        omegaPH : TYPE, optional
            DESCRIPTION. The default is 1.

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''
        #phononList needs to be list with tuples that are the end points of the electron prop 
        #takes a electron propogator finds nearest neighbors 
        #find its nerest neighbors 
        #select a time in between left and right times uniform 
        #pick a point where there is only k
        tauList=[0]
        for i in qList:
            q,a,b=i
            tauList.extend([a,b])
        
        tauList.sort()
        index=nrandG.integers(len(tauList))
        tauList.append(tauMax)
        tauLeft,tauRight=tauList[index],tauList[index+1]
        tauOne=nrandG.uniform(tauLeft,tauRight)
        u=nrandG.uniform()
        tauTwo=tauOne-np.log(u)/omegaPH
        if tauTwo>tauRight:
            return qList
        mean=np.zeros(3)
        variance=np.diag(m/(tauTwo-tauOne)*np.ones(3))
        q=nrandG.multivariate_normal(mean,variance)
            #does not use box muller but maybe not an issue?
        dummy=q,tauOne,tauTwo
        
        
        x=nrandG.uniform()
        
        
        if x<r(pIns,pRem,tauOne,tauTwo,tauLeft,tauRight,q,k,mu,alpha,order):#put in all the variables
            qList.append(dummy)
        
        return qList,1
    
        
   
        
        
    def removeProp(qList,k,mu,alpha,order,pIns,pRem):
        '''
        

        Parameters
        ----------
        qList : TYPE
            DESCRIPTION.
        k : TYPE
            DESCRIPTION.
        mu : TYPE
            DESCRIPTION.
        alpha : TYPE
            DESCRIPTION.
        order : TYPE
            DESCRIPTION.
        pIns : TYPE
            DESCRIPTION.
        pRem : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''
        #pick random prop
        index=nrandG.integer(len(qList))
        q,tauLeft,tauRight=qList[index]
        
        #currently my list is unordered this would be faster if it was ordered also
        #might be easier if only had qList
        for i in qList:
            dummy,tL,tR=i
            if inBetween(tauLeft,tauRight,tL,tR)==True:
                return qList
            
        if nrandG.uniform(1)<r(pIns,pRem,tR,tL,tauLeft,tauRight,q,k,mu,alpha,order)**-1:
            qList.pop(index)
        return qList,-1
        
        #need to compute inverse r if passes then remove index value from q list
            
            
    
        
    
    def swap(qList,k,tauMax,mu,order):
        '''
        

        Parameters
        ----------
        qList : TYPE
            DESCRIPTION.
        k : TYPE
            DESCRIPTION.
        tauMax : TYPE
            DESCRIPTION.
        mu : TYPE
            DESCRIPTION.
        order : TYPE
            DESCRIPTION.

        Returns
        -------
        qList : TYPE
            DESCRIPTION.

        '''
        
    
        
        index=nrandG.integer(len(qList))
        qOne,tauOne,tauA=qList[index]
        dummy=tauMax
        i=0
        while i<=order:
            if i!=index:
                continue
            qP,tauBP,tauP=qList[i]
            if abs(tauOne-tauP)<dummy:
                index2=i
                dummy=abs(tauOne-tauP)
                qTwo,tauTwo,tauB=qP,tauP,tauBP
        
        kP=k-qOne-qTwo
        
        wX=greensFunction(k, tauOne-tauTwo, mu)*phononProp(abs(tauOne-tauA))/qOne**2\
            *phononProp(abs(tauTwo-tauB))/qTwo**2
        wY=-greensFunction(kP, tauTwo-tauOne, mu)*phononProp(abs(tauOne-tauB))/qTwo**2\
            *phononProp(abs(tauTwo-tauA))/qOne**2
        
        if nrandG.uniform(1)<wY/wX:
            if index<index2:
                qList.pop(index2)
                qList.pop(index)
            else:
                qList.pop(index)
                qList.pop(index2)
            
            phonon1=qOne,tauTwo,tauA
            qList.append(phonon1)
            phonon2=qTwo,tauOne,tauB
            qList.append(phonon2)
            return qList
    
        
    def order(qList):
        return len(qList)
    
        
            
    #current issue this only looks under arcs
    
    #need to create list of nodes
        