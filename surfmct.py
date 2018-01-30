# -*- coding: utf-8 -*-
#!usr/bin/env python2
"""
Created on Sat Jan 13 20:30:18 2018


Surface simulation using perfect matching
min weigth perfect matching

Monte Carlo simulation of many iterations

@author: pedro
"""
import os
import numpy as np
#import networkx as nx
#import matplotlib.pyplot as plt
from surftimehx import *
import time
import sys

L=int(sys.argv[1])
Nit=int(sys.argv[2])

Nsteps=10
x=[]
y=[]
t=[]



print "Work L="+str(L)+", Nit="+str(Nit)
a=.02
b=.04
P=np.linspace(a,b,Nsteps)
for p in P:
    start=time.time()
    per=mcexp(L,p,p,Nit)
    end=time.time()

    x.append(p)
    y.append(per)
    t.append(end-start)
    print p,per,end-start,end-time.time()

'''  #bin search
for i in range(Nsteps):
    p=(a+b)/2.
    start = time.time()
    per=mcexp(L,p,Nit)
    end= time.time()
    x.append(p)
    y.append(per)
    t.append(end-start)
    if per>p:
        b=p
    else:
        a=p
'''

filen="./results/smctL"+str(L)+"Nit"+str(Nit)+".txt"


f=open(filen,"w")
for i in range(len(x)):
    f.write(str(x[i])+" "+str(y[i])+" "+str(t[i])+"\n")
f.close()

