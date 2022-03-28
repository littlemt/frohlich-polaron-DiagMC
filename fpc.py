#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 26 15:23:02 2022

@author: Milo Kamiya Belmont
"""

from dep.FredrichPoleronCode import diagMC as FPC
import time
import numpy.random as nrng
import matplotlib.pyplot as mpl

nrg=nrng.default_rng()


def main(tauMax,mu,alpha,P):
    #P is list of probablilitys for each update need to sum to one
    #gen initial conditions 
    for i in range(1000):
        tauList=[]
        qList=[]
        k=1
        order=FPC.order(qList)
        
        x=nrg.uniform()
        
        if order==0:
            tau,i=FPC.changeTau()
            FPC.insert()
            FPC.remove()
        elif order==1:
            FPC.insert()
            FPC.remove()
        else:
            FPC.insert()
            FPC.remove()
            FPC.swap()
    return 
        
def zero_order(tauMax,runTime,p,mu,alpha,m=1):
    tauList=[]
    runTime=runTime*3600+time.time()
    count=0
    while time.time()<runTime:
        tau,i = FPC.tauChange(p,tauMax,mu)
        tauList.append(tau)
        count += i
        
    return tauList,count

def first_order(tauMax,runTime,p,mu,k,alpha,m=1):
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
            tau,i = FPC.tauChange(k,tauMax,mu,m)
            tauList.append(tau)
            countT += i
        elif pTau<x<=pTau+pIns:
            qList,i=FPC.insertProp(qList, tauMax, k, m, mu, alpha, order,pIns,pRem)
            countI+=i
        elif pTau+pIns<x<=1 and FPC.order(qList)>=1:
            qList,i=FPC.removeProp(qList, k, mu, alpha, order,pIns,pRem)
            countR+=i
    count=[countT,countI,countR]
    return tauList,qList,count