#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 14:59:40 2017

Surface code simulation using blossomV
min weigth perfect matching from Kolmogorov


@author: pedro
"""
import os
import numpy as np
#import networkx as nx
import matplotlib.pyplot as plt

'''
            Write and Read files and Blossom
'''
#saves a graph in a file, preparing it for kolmogorov's code
def savegraph(n,edges,file="graph.txt"):
    f=open(file,"w")
    f.write(str(n)+" "+str(len(edges))+"\n")
    for tup in edges:
        f.write(str(tup[0])+" "+str(tup[1])+" "+str(tup[2])+"\n")
    f.close()
    return
#calls blossomV from Kolmogorov's code
def blossom(graphfile="graph.txt",solfile="solution.txt"):
    order="./blossom5 -e "+graphfile+" -w "+solfile+" -V"
    os.system(order)
    return
#reads the solution given by Kolmogorov's code
def readsol(solfile="solution.txt"):
    f=open(solfile)
    match=[]
    for line in f:
        r=line.split()
        x=int(r[0])
        y=int(r[1])
        match.append((x,y))
    c=match.pop(0)
    f.close()
    return match,c

'''
            Some functions for manipulating indexes in the lattice
            coordinates for the operators
'''
def xy(i,L):    
    x=i%L
    y=i/L
    return x,y
def indexy(x,y,L):
    return x%L+y%L*L       
#coordinates for the qubits 
def xyz(n,Lx,Ly=-1):  ###
    '''
    x=i%L
    y=(i%(L*L))/L
    z=i/(L*L)
    return x,y,z
    '''
    if Ly==-1:
        Ly=Lx
    nk0=(Lx-1)*(Ly-1)
    k=n>=nk0
    i=(1-k)*(n%(Lx-1))+k*((n-nk0)%Lx)
    j=(1-k)*(n/(Lx-1))+k*((n-nk0)/Lx)
    return i,j,k
def indexyz(i,j,z,Lx,Ly=-1):###
    if Ly==-1:
        Ly=Lx
    #return x%L+y%L*L+z*L*L      
    return (1-z)*(i+(Lx-1)*j)+z*(i+Lx*j+(Lx-1)*(Ly-1))     
def mdst(a,b,L):#manhattan distance
    x1=a%L
    x2=b%L
    y1=a/L
    y2=b/L
    dist=np.abs(x1-x2)+np.abs(y1-y2)
    return dist
def distop(a,b,L):#manhattan distance with periodic boundary conditions    
    x1=a%L
    x2=b%L
    y1=a/L
    y2=b/L
    distx=np.abs(x1-x2)
    disty=np.abs(y1-y2)
    dist=min(distx,L-distx)+min(disty,L-disty)
    return dist

class Node:
    def __init__(self,n,Lx,Ly=-1):
        self.n=n
        self.Lx=Lx
        if Ly==-1:
            self.Ly=Lx
        else:
            self.Ly=Ly
        nk0=(Lx-1)*(Ly-1)
        self.k=n/nk0
        self.i=(1-self.k)*(n%(Lx-1))+self.k*((n-nk0)%Lx)
        self.j=(1-self.k)*(n/(Lx-1))+self.k*((n-nk0)/Lx)
        self.op=0
        
class StabX:###
    def __init__(self,n,Lx,Ly=-1):
        self.n=n
        self.Lx=Lx
        if Ly==-1:
            self.Ly=Lx
            Ly=Lx
        else:
            self.Ly=Ly
        self.i=n%(Lx)
        self.j=n/(Lx)
        self.name="X"
        
        #neighbor qubits and neighbor stabilizers for bfs
        
        nop=[]
        nq=[]
        #neighbor qubits for stabilizer measurement
        meas=[]
        nopi=[(1,0),(0,-1),(0,+1),(-1,0)]
        nqi=[(0,0,0),(0,0,1),(0,+1,1),(-1,0,0)]
        
        for i in range(len(nopi)):
            ni=self.i+nopi[i][0]
            nj=self.j+nopi[i][1]
            #print ni,nj, ni>=0 , ni<Lx , nj>=0 , nj<(Ly-1)
            if ni>=0 and ni<Lx and nj>=0 and nj<(Ly-1):
                nn=ni+Lx*nj
                nop.append(nn)
                ni=self.i+nqi[i][0]
                nj=self.j+nqi[i][1]
                nk=nqi[i][2]
                nn=indexyz(ni,nj,nk,Lx,Ly)
                nq.append(nn)
        
        for i in range(len(nopi)):
            ni=self.i+nqi[i][0]
            nj=self.j+nqi[i][1]
            nk=nqi[i][2]
            if nk==0 and ni>-1 and nj>-1:
                if ni<Lx-1 and nj<Ly-1:
                    meas.append(indexyz(ni,nj,nk,Lx,Ly))
            if nk==1 and ni>-1 and nj>-1:
                if ni<Lx and nj<Ly:
                    meas.append(indexyz(ni,nj,nk,Lx,Ly))
        
        #neighbor qubits and neighbor stabilizers for bfs
        #print meas
        self.nop=nop
        self.nq=nq
        #neighbor qubits for stabilizer measurement
        self.meas=meas
class StabZ:###
    def __init__(self,n,Lx,Ly=-1):
        self.n=n
        self.Lx=Lx
        if Ly==-1:
            self.Ly=Lx
            Ly=Lx
        else:
            self.Ly=Ly
        self.i=n%(Lx-1)
        self.j=n/(Lx-1)
        self.name="Z"
        #neighbor qubits and neighbor stabilizers for bfs
        
        self.nop=[]
        self.nq=[]
        #neighbor qubits for stabilizer measurement
        self.meas=[]
        nopi=[(0,1),(-1,0),(+1,0),(0,-1)]
        nqi=[(0,0,0),(0,0,1),(+1,0,1),(0,-1,0)]
        for i in range(len(nopi)):
            ni=self.i+nopi[i][0]
            nj=self.j+nopi[i][1]
            if ni>=0 and ni<(Lx-1) and nj>=0 and nj<Ly:
                nn=ni+(Lx-1)*nj
                self.nop.append(nn)
                ni=self.i+nqi[i][0]
                nj=self.j+nqi[i][1]
                nk=nqi[i][2]
                nn=indexyz(ni,nj,nk,Lx,Ly)
                self.nq.append(nn)
        
        for i in range(len(nopi)):
            ni=self.i+nqi[i][0]
            nj=self.j+nqi[i][1]
            nk=nqi[i][2]
            if nk==0 and ni>-1 and nj>-1:
                if ni<Lx-1 and nj<Ly-1:
                    self.meas.append(indexyz(ni,nj,nk,Lx,Ly))
            if nk==1 and ni>-1 and nj>-1:
                if ni<Lx and nj<Ly:
                    self.meas.append(indexyz(ni,nj,nk,Lx,Ly))
        

'''
            simulation functions
'''


#applies noise to a state with probability p for each qbit
def noisestate(L,p,Ly=-1):###
    if Ly==-1:
        Ly=L
    nq=(L-1)*(Ly-1)+Ly*L
    
    err=[0]*nq
    for i in range(nq):
        r=np.random.random()
        if r<p:
            r=np.random.random()
            if (r<1./3.):
                pi=1
            elif(r<2./3):
                pi=2
            else:
                pi=3
            err[i]=pi
    return err
def plotnoisestate(err,L,col='r',mize=8):###
    Ly=L
    nq=(L-1)*(Ly-1)+Ly*L
    for i in range(nq):
        x,y,z=xyz(i,L)
        
        if err[i]==1:
            plt.plot(x+.5*(1-z),y+0.5*(1-z),marker='$X$',color=col,markersize=mize)
        if err[i]==2:
            plt.plot(x+.5*(1-z),y+0.5*(1-z),marker='$Y$',color=col,markersize=mize)
        if err[i]==3:
            plt.plot(x+.5*(1-z),y+0.5*(1-z),marker='$Z$',color=col,markersize=mize)
    return
    

#evaluates the syndrome for a certain noise state:
    
def syndrome(err,Zc,Xc,L):### can be optimized by putting Xc,Zc as input
    #Ly=L
    nstab=L*(L-1)
    Z=[0]*nstab
    X=[0]*nstab
    '''
    Zc=[]
    Xc=[]
    #list of stabilizers to be measured
    for i in range(nstab):
        Zc.append(StabZ(i,L))
        Xc.append(StabX(i,L))
    '''
        
    for k in range(nstab):
        neigx=Xc[k].meas #list of qubits involved
        for ind in neigx:
            if err[ind]>1:#if err==Y or Z
                X[k]+=1
        X[k]=X[k]%2
        neigz=Zc[k].meas #list of qubits involved
        for ind in neigz:
            if err[ind]>0 and err[ind]<3:#if err==Y or Z
                Z[k]+=1
        Z[k]=Z[k]%2
    return X,Z
def plotsyndrome(X,Z,L,colx='y',colz='g',mize=8):###
    for k in range(L*(L-1)):
        i,j=xy(k,L)
        if X[k]>0:
            plt.plot(i,j+.5,marker='d',color=colx,markersize=mize)
        i,j=xy(k,L-1)
        if Z[k]>0:
            plt.plot(i+.5,j,marker='d',color=colz,markersize=mize)
    return


def distToBorderX(n,L):###
    #X stabilizers are linked to rough edges (top and bottom)
    j=n/L
    return min(j+1,L-1-j)

def distToBorderZ(n,L):###
    #X stabilizers are linked to smooth edges (laterals)
    i=n%(L-1)
    return min(i+1,L-i-1)


    
#evaluates min weight perfect matching for the syndrome

def blossomopX(X,L):###
    nx=sum(X)
    xl=[]
    sol=[]
    edges=[]    
    ner=nx
    
    if nx>0:
        #a node for each syndrome
        for i in range(len(X)):
            if X[i]==1:
                xl.append(i)
        #edges between nodes and closest vertical border
        for n in range(nx):
            dist=distToBorderX(xl[n],L)
            edges.append((n,n+ner,dist))
            
        #edges between nodes
        for e1 in range(nx-1):
            for e2 in range(e1+1,nx):
                dist=mdst(xl[e1],xl[e2],L)
                edges.append((e1,e2,dist))
                edges.append((e1+ner,e2+ner,0))
                
        savegraph(nx*2,edges)
        blossom()
        match,d=readsol()
            
        sol=[]
        for edg in match:
            if edg[0]<ner or edg[1]<ner:#only record interesting edges
                e1=xl[edg[0]%ner]#+L*L*(edg[0]/ner)
                e2=xl[edg[1]%ner]#+L*L*(edg[1]/ner)
                sol.append((e1,e2))
                
    return sol
          
def blossomopZ(X,L):###
    nx=sum(X)
    xl=[]
    sol=[]
    edges=[]    
    ner=nx
    
    if nx>0:
        #a node for each syndrome
        for i in range(len(X)):
            if X[i]==1:
                xl.append(i)
        #edges between nodes and closest vertical border
        for n in range(nx):
            dist=distToBorderZ(xl[n],L)
            edges.append((n,n+ner,dist))
            
        #edges between nodes
        for e1 in range(nx-1):
            for e2 in range(e1+1,nx):
                dist=mdst(xl[e1],xl[e2],L-1)
                edges.append((e1,e2,dist))
                edges.append((e1+ner,e2+ner,0))
                
        savegraph(nx*2,edges)
        blossom()
        match,d=readsol()
            
        sol=[]
        for edg in match:
            if edg[0]<ner or edg[1]<ner:
                e1=xl[edg[0]%ner]#+L*L*(edg[0]/ner)
                e2=xl[edg[1]%ner]#+L*L*(edg[1]/ner)
                sol.append((e1,e2))
                
    return sol
              

def bfsx(a,b,S,L):#
    if a==b:
        name=S[a].name
        path=[]           
        i=S[a].i
        j=S[a].j
        if name=="Z":
            #trace the shortest horizontal path to the border 
            if (i+1)<(L-i-1):
                #path from left to 'a'                
                for n in range(i+1):
                    path.append(indexyz(n,j,1,L))
            else:
                #path from right to 'a'                
                for n in range(L-i-1):
                    path.append(indexyz(L-1-n,j,1,L))
 
        if name=="X":
            #trace the shortest vertical path to the border
            if (j+1)<(L-j-1):
                #path bottom to 'a'
                for n in range(j+1):
                    path.append(indexyz(i,n,1,L))
            else:
                #path from top to 'a'
                for n in range(L-j-1):
                    path.append(indexyz(i,L-1-n,1,L))
                    
        return path
    front=[a]
    paths=[[]]
    visited=set()
    while front:
        node=front.pop(0)
        path=paths.pop(0)
        neigs=S[node].nop
        qneigs=S[node].nq
        visited.add(node)
        for n in range(len(neigs)):
            if neigs[n] == b:
                path+=[qneigs[n]]
                return path
            if neigs[n] not in visited:
                front.append(neigs[n])
                paths.append(path+[qneigs[n]])
    return []     
def plotpathx(a,b,path,L,col='y',marker='$Z$', mize=7):###
    for n in path:
        i,j,z=xyz(n,L)
        x=i+.5*(1-z)
        y=j+.5*(1-z)
        #xa,ya=xy(a,L)
        #xb,yb=xy(b,L)
        #plt.plot([xa,xb],[ya+.5,yb+.5],'r-')
        plt.plot(x,y,color=col,marker=marker,markersize=mize)
    
    
    
    
def plotpathz(a,b,path,L,col='g',marker='$X$', mize=4):###
    for n in path:
        i,j,z=xyz(n,L)
        x=i+.5*(1-z)
        y=j+.5*(1-z)
        #xa,ya=xy(a,L-1)
        #xb,yb=xy(b,L-1)
        #plt.plot([xa+.5,xb+.5],[ya,yb],'b-')
        plt.plot(x,y,color=col,marker=marker,markersize=mize)
 
    
    
def plotmatch(solX,solZ,L):###
    nstab=L*(L-1)
    Sx=[]
    Sz=[]
    #list of stabilizers to be measured
    for i in range(nstab):
        Sz.append(StabZ(i,L))
        Sx.append(StabX(i,L))
    for pair in solX:
        a=pair[0]
        b=pair[1]
        path=bfsx(a,b,Sx,L)
        plotpathx(a,b,path,L)
    for pair in solZ:
        a=pair[0]
        b=pair[1]
        path=bfsx(a,b,Sz,L)
        plotpathz(a,b,path,L)
        
        
        
'''
            final result error+corrections
'''
        
def corrx(solX,S,L):###
    corx=[0]*((L-1)*(L-1)+L*L)
    for pair in solX:
        a=pair[0]
        b=pair[1]
        path=bfsx(a,b,S,L)
        for q in path:
            corx[q]=3
    return corx
def corrz(solX,S,L):###
    corx=[0]*((L-1)*(L-1)+L*L)
    for pair in solX:
        a=pair[0]
        b=pair[1]
        path=bfsx(a,b,S,L)
        for q in path:
            corx[q]=1
    return corx
#product of pauli matrices
def pprod(a,b):###
    if a*b==0:
        return a+b
    if a==b:
        return 0
    prod=[3,2,1]
    return prod[(a+b)%3]
    
def ercor(err,corx,corz):###
    sol=[0]*len(err)
    for i in range(len(err)):
        sol[i]=pprod(pprod(err[i],corx[i]),corz[i])
    return sol
        
        
'''
            checking logical errors
'''        
def conmz(a):###
    if a==1 or a==2:
        return 1
    else:
        return 0
def conmx(a):###
    if a>1:
        return 1
    else:
        return 0
def checkl(erc,L):###
    vx=0
    vz=0
    hx=0
    hz=0    
    for i in range(L):
        hx+=conmx(erc[indexyz(i,0,1,L)])
        vz+=conmz(erc[indexyz(0,i,1,L)])
    #for i in range(L-1):
        #hx+=conmz(erc[indexyz(i,0,1,L)])
        #vx+=conmx(erc[indexyz(0,i,0,L)])
    hx=hx%2
    #hz=hz%2
    #vx=vx%2
    vz=vz%2
    return min(hx+vz,1),hz,vx


'''
            plotting functions
'''


#plots the entire lattice
def latplot(Lx,Ly=-1,col='b'):###
    if Ly==-1:
        Ly=Lx
    L=Lx
    #square lattice plot
    
    for n in range(L):
        plt.plot([n,n],[-.5,Ly-.5],col)
    for n in range(Ly-1):
        plt.plot([0,Lx-1],[n+.5,n+.5],col)
    return


def experimentplot(L,p,Ly=-1):###
    nstab=L*(L-1)
    Sx=[]
    Sz=[]
    #list of stabilizers to be measured
    for i in range(nstab):
        Sz.append(StabZ(i,L))
        Sx.append(StabX(i,L))
    
    
    plt.figure(0,figsize=(2,2))
    plt.title("Surface code")
    plt.clf()
    plt.subplot(2,2,1)
    plt.title('Errors')
    latplot(L)
    
    err=noisestate(L,p)
    plotnoisestate(err,L)
    
    
    
    plt.subplot(2,2,2)
    plt.title("Syndrome")
    latplot(L)
    plotnoisestate(err,L)
    X,Z=syndrome(err,Sz,Sx,L)
    plotsyndrome(X,Z,L)
    #print sum(X), sum(Z)
    
    
    
    plt.subplot(2,2,3)
    plt.title("Matching")
    latplot(L)
    solX=blossomopX(X,L)
    solZ=blossomopZ(Z,L)
    plotsyndrome(X,Z,L)
    #plotmatchX(solX,L)
    #plotmatchZ(solZ,L)
    
    
    plotmatch(solX,solZ,L)
    plt.subplot(2,2,4)
    
    latplot(L)
    corx=corrx(solX,Sx,L)
    corz=corrz(solZ,Sz,L)
    erc=ercor(err,corx,corz)
    plotnoisestate(erc,L)
    X,Z=syndrome(erc,Sz,Sx,L)
    plotsyndrome(X,Z,L)
    print erc
    logic=checkl(erc,L)
    
    if logic[0]>0:
        plt.title("Final result:\n logical error")
    else:
        plt.title("Final result:\n No logical error")
    
    plt.show()
    return logic[0]

def experiment(L,p):
        
    nstab=L*(L-1)
    Sx=[]
    Sz=[]
    #list of stabilizers to be measured
    for i in range(nstab):
        Sz.append(StabZ(i,L))
        Sx.append(StabX(i,L))
    
    
    
    err=noisestate(L,p) 
    
    X,Z=syndrome(err,Sz,Sx,L)
    
    solX=blossomopX(X,L)
    solZ=blossomopZ(Z,L)    
    
    corx=corrx(solX,Sx,L)
    corz=corrz(solZ,Sz,L)
    erc=ercor(err,corx,corz)
    X,Z=syndrome(erc,Sz,Sx,L)
    logic=checkl(erc,L)
    
    return logic[0]


#THRESHOLD SEARCH

def errorfind(err,Sz,Sx,L):
    X,Z=syndrome(err,Sz,Sx,L)
    
    solX=blossomopX(X,L)
    solZ=blossomopZ(Z,L)    
    
    corx=corrx(solX,Sx,L)
    corz=corrz(solZ,Sz,L)
    erc=ercor(err,corx,corz)
    X,Z=syndrome(erc,Sz,Sx,L)
    logic=checkl(erc,L)
    
    return logic[0]
    
#creates a state with c errors from nqor list of nodes
def errorc(nqor,c):    
    nq=len(nqor)
    err=[0]*nq
    for i in range(c):
        err[nqor[i]]=1
    return err
    
def threshold(L):
     
    nstab=L*(L-1)
    nqs=(L-1)*(L-1)+L*L
    Sx=[]
    Sz=[]
    #list of stabilizers to be measured
    for i in range(nstab):
        Sz.append(StabZ(i,L))
        Sx.append(StabX(i,L))
    
    #we create a random permutation of the indexes of the qubits
    qor=[]
    qbag=range(nqs)
    while qbag:
        r=np.random.randint(0,len(qbag))
        qor.append(qbag.pop(r))
    
    #we start a binary search between 0 and nq/5
    b=nqs/3
    a=0    
    c=(b+a)/2
    logic=1
    while (b-a)>1:     
        c=(b+a)/2
        #print a,c,b,logic
        err=errorc(qor,c)
        logic=errorfind(err,Sz,Sx,L)
        if logic>0:
            b=c
        else:
            a=c
    return 1.*c/nqs
    
def threshold2(L):
     
    nstab=L*(L-1)
    nqs=(L-1)*(L-1)+L*L
    Sx=[]
    Sz=[]
    #list of stabilizers to be measured
    for i in range(nstab):
        Sz.append(StabZ(i,L))
        Sx.append(StabX(i,L))
    
    #we create a random permutation of the indexes of the qubits
    qor=[]
    qbag=range(nqs)
    while qbag:
        r=np.random.randint(0,len(qbag))
        qor.append(qbag.pop(r))
    
    #we start the search with a small c, until there is any logic error
    c=nqs/20
    logic=1
    while c<nqs:     
        c+=1
        err=errorc(qor,c)
        logic=errorfind(err,Sz,Sx,L)
        if logic>0:
            return c
    return c
    
    
def phreshold(L,plotres=1):
     
    nstab=L*(L-1)
    nqs=(L-1)*(L-1)+L*L
    Sx=[]
    Sz=[]
    #list of stabilizers to be measured
    for i in range(nstab):
        Sz.append(StabZ(i,L))
        Sx.append(StabX(i,L))
    
    #we create a random permutation of the indexes of the qubits
    qor=[]
    qbag=range(nqs)
    while qbag:
        r=np.random.randint(0,len(qbag))
        qor.append(qbag.pop(r))
    
    #we start a binary search between 0 and nq/5
    N=range(nqs)
    sol=[0]*nqs    
    for c in N:
        
        err=errorc(qor,c)
        logic=errorfind(err,Sz,Sx,L)
        sol[c]=logic
    if plotres==1:
        plt.plot(N,sol,'o')
    return sol
    
    
def ploterr(err,L):###
    nstab=L*(L-1)
    Sx=[]
    Sz=[]
    #list of stabilizers to be measured
    for i in range(nstab):
        Sz.append(StabZ(i,L))
        Sx.append(StabX(i,L))
    
    
    plt.figure(2,figsize=(2,2))
    plt.title("Surface code")
    plt.clf()
    plt.subplot(2,2,1)
    plt.title('Errors')
    latplot(L)
    
    plotnoisestate(err,L)
    
    
    
    plt.subplot(2,2,2)
    plt.title("Syndrome")
    latplot(L)
    plotnoisestate(err,L)
    X,Z=syndrome(err,Sz,Sx,L)
    plotsyndrome(X,Z,L)
    #print sum(X), sum(Z)
    
    
    
    plt.subplot(2,2,3)
    plt.title("Matching")
    latplot(L)
    solX=blossomopX(X,L)
    solZ=blossomopZ(Z,L)
    plotsyndrome(X,Z,L)
    #plotmatchX(solX,L)
    #plotmatchZ(solZ,L)
    
    
    plotmatch(solX,solZ,L)
    plt.subplot(2,2,4)
    
    latplot(L)
    corx=corrx(solX,Sx,L)
    corz=corrz(solZ,Sz,L)
    erc=ercor(err,corx,corz)
    plotnoisestate(erc,L)
    X,Z=syndrome(erc,Sz,Sx,L)
    plotsyndrome(X,Z,L)
    #print erc
    logic=checkl(erc,L)
    
    if logic[0]>0:
        plt.title("Final result:\n logical error")
    else:
        plt.title("Final result:\n No logical error")
    
    plt.show()
    return logic[0]

def nploterr(err,L):###same as plotter, does not plot anything
    nstab=L*(L-1)
    Sx=[]
    Sz=[]
    #list of stabilizers to be measured
    for i in range(nstab):
        Sz.append(StabZ(i,L))
        Sx.append(StabX(i,L))
    
    X,Z=syndrome(err,Sz,Sx,L)
    
    
    solX=blossomopX(X,L)
    solZ=blossomopZ(Z,L)
    
    corx=corrx(solX,Sx,L)
    corz=corrz(solZ,Sz,L)
    erc=ercor(err,corx,corz)
    X,Z=syndrome(erc,Sz,Sx,L)
    logic=checkl(erc,L)
    
    return logic[0]



'''
x=[]
for i in range(1000):
    x.append(threshold2(17))
plt.clf()
plt.hist(x,100)
'''


'''
plt.clf()
print threshold2(7)
'''


'''
x=np.array(phreshold(7,0))
for i in range(1000):
    x+=np.array(phreshold(7,0))
plt.plot(x*1.0/sum(x))
'''


'''

L=7
nstab=L*(L-1)
nqs=(L-1)*(L-1)+L*L
Sx=[]
Sz=[]
#list of stabilizers to be measured
for i in range(nstab):
    Sz.append(StabZ(i,L))
    Sx.append(StabX(i,L))

#we create a random permutation of the indexes of the qubits
qor=[]
qbag=range(nqs)
while qbag:
    r=np.random.randint(0,len(qbag))
    qor.append(qbag.pop(r))

#we start a binary search between 0 and 1/2
b=nqs/2
a=0    
c=(b+a)/2
logic=1


err=errorc(qor,c)
ploterr(err,L)




c=(b+a)/2
#print a,c,b,logic
err=errorc(qor,c)
logic=errorfind(err,Sz,Sx,L)
ploterr(err,L)
if logic>0:
    b=c
else:
    a=c
'''