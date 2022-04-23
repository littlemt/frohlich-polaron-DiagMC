#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 26 15:23:02 2022

@author: Milo Kamiya Belmont
"""

from dep.FrohlichPoleronDiagMCUpdates import diagMC as FPC
import time
import numpy as np
import numpy.random as nrng
import matplotlib.pyplot as mpl

    

nrg=nrng.default_rng()


def full(tauMax,runTime,p,pExt,mu,k,alpha,m=1):
    
    
    tauList=[]
    qList=[]
    total=sum(p)
    pTau=p[0]/total
    pIns=p[1]/total
    pRem=p[2]/total
    pSwap=p[3]/total
    countT=0
    countI=0
    countR=0
    countS=0
    tau=tauMax
    mcTime=[]
    mcT=0
    
    orderList=[]
    runTime=runTime*3600+time.time()
    
    while time.time()<runTime:
        
        
        order=FPC.order(qList)
        orderList.append(order)
        
        
        x=nrg.uniform()
        
        if x<pTau and order==0:
            
            mcTime.append(mcT)
            tau,i=FPC.tauChange(pExt,tauMax,tau,mu,m)
            tauList.append(tau)
            countT+=i
            mcT=0
        elif pTau<=x<pTau+pIns:
            
            qList,i=FPC.insertProp(qList, tau, k, mu, alpha, order,pIns,pRem)
            countI+=i
            
        elif pTau+pIns<=x<pTau+pIns+pRem and order>=1:
            
            qList,i=FPC.removeProp(qList,tau, k, mu, alpha, order,pIns,pRem)
            countR+=i
            
        elif pTau+pIns+pRem<=x<1 and order>=2:
            
            qList,i=FPC.swap(qList, k, tauMax, mu, order)
            countS+=i
            
        mcT+=1        
            
    count=[countT,countI,countR,countS]
    return tauList,mcTime,qList,count
        
def zero_order(tauMax,runTime,pExt,mu,alpha,m=1):
    tauList=[]
    runTime=runTime*3600+time.time()
    count=0
    tau=tauMax
    while time.time()<runTime:
        tau,i = FPC.tauChange(pExt,tauMax,tau,mu)
        tauList.append(tau)
        count += i
        
    return tauList,count

def first_order(tauMax,runTime,p,pExt,mu,k,alpha,m=1):
    tauList=[]
    qList=[]
    total=sum(p)
    pTau=p[0]/total
    pIns=p[1]/total
    pRem=p[2]/total
    countT=0
    countI=0
    countR=0
    tau=tauMax
    runTime=runTime*3600+time.time()
    while time.time()<runTime:
        x=nrg.uniform()
        order=FPC.order(qList)
        if order==0 and 0<=x<=pTau:
            tau,i = FPC.tauChange(pExt,tauMax,tau,mu,m)
            tauList.append(tau)
            countT += i
        elif pTau<x<=pTau+pIns:
            qList,i=FPC.insertProp(qList, tau, k, mu, alpha, order,pIns,pRem)
            countI+=i
        elif pTau+pIns<x<=1 and order>=1:
            qList,i=FPC.removeProp(qList,tau, k, mu, alpha, order,pIns,pRem)
            countR+=i
    count=[countT,countI,countR]
    return tauList,qList,count



def data_plot(qList):
    length=len(qList)
    
    table=np.ndarray((length,3))
    for i in range(length):
        q,tau1,tau2=qList[i]
        table[i]=[q,tau1,tau2]
        
    return table
        
#add extend 

        
        