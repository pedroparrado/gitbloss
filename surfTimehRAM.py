#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 17:34:54 2018

@author: pedro
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 16:21:22 2017


Surface code simulation using blossomV
min weigth perfect matching from Kolmogorov
with measurement errors


@author: pedro
"""
import os
import numpy as np
import surfCodeh as scd
import surfCodehRAM as scdR
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
    order="./blossom5 -e "+graphfile+" -w "+solfile+ " -V"
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

def ijtX(n,L):    #]
    i=n%L
    j=(n%(L*(L-1)))/L
    t=n/(L*(L-1))
    return i,j,t
def ijtZ(n,L):    #]
    i=n%(L-1)
    j=(n%(L*(L-1)))/(L-1)
    t=n/(L*(L-1))
    return i,j,t
def indexX(i,j,t,L):#]
    return i+j*L+t*L*(L-1)    
def indexZ(i,j,t,L):#]
    return i+j*(L-1)+t*L*(L-1)    




   
#coordinates for the qubits 
def ijzt(n,L):  #]
    
    nq=(L-1)*(L-1)+L*L
    t=n/nq
    z=(n%nq)/((L-1)*(L-1))
    if z>1:
        z=1
    np=(n%nq)%((L-1)*(L-1))
    npp=n%nq-(L-1)*(L-1)
    if z==0:
        i=np%(L-1)
        j=np/(L-1)
    else:
        i=npp%L
        j=npp/L
        
    return i,j,z,t

def xyq(n,L):#]
    phi=60*np.pi/180.
    vi=np.array((1,0))
    vj=np.array((np.sin(phi),np.cos(phi)))*.75
    vt=np.array((0,1))
    i,j,z,t=ijzt(n,L)
    if z==0:
        x,y=vi*(i+.5)+vj*(j+1)+vt*t
    else:
        x,y=vi*i+vj*(j+.5)+vt*t
    return x,y
def xyX(n,L):#]
    phi=60*np.pi/180.
    vi=np.array((1,0))
    vj=np.array((np.sin(phi),np.cos(phi)))*.75
    vt=np.array((0,1))
    i,j,t=ijtX(n,L)
    x,y=vi*i+vj*(j+1)+vt*t
    return x,y
def xyZ(n,L):#]
    phi=60*np.pi/180.
    vi=np.array((1,0))
    vj=np.array((np.sin(phi),np.cos(phi)))*.75
    vt=np.array((0,1))
    i,j,t=ijtZ(n,L)
    x,y=vi*(i+.5)+vj*(j+.5)+vt*t
    return x,y
def indexq(i,j,z,t,L):#]
    nq=(L-1)*(L-1)+L*L
    if z==0:
        return i+j*(L-1)+nq*t
    else:
        return i+j*L+(L-1)*(L-1)+nq*t 
 
def mdst(a,b,L):#manhattan distance #]
    global wq,wp
    x1,y1,t1=ijtX(a,L)
    x2,y2,t2=ijtX(b,L)
    dist=np.abs(x1-x2)*wp+wp*np.abs(y1-y2)+wq*np.abs(t1-t2)
    return int(dist*10000)
def mdstZ(a,b,L):#manhattan distance #]
    global wq,wp
    x1,y1,t1=ijtZ(a,L)
    x2,y2,t2=ijtZ(b,L)
    dist=np.abs(x1-x2)*wp+wp*np.abs(y1-y2)+wq*np.abs(t1-t2)
    return int(dist*10000)


class Node:#]?
    def __init__(self,n,L):
        self.n=n
        self.L=L
        i,j,z,t=ijzt(n,L)
        self.i=i
        self.j=j
        self.z=z
        self.t=t
        self.op=0
        
class StabX:#]
    def __init__(self,n,L):
        self.n=n
        self.L=L
        i,j,t=ijtX(n,L)
        self.i=i
        self.j=j
        self.t=t
        self.name="X"
        
        #neighbor qubits and neighbor stabilizers for bfs
        
        nop=[]
        nq=[]
        #neighbor qubits for stabilizer measurement
        meas=[]
        
        nopi=[(1,0,0),(0,-1,0),(0,+1,0),(-1,0,0),(0,0,1),(0,0,-1)]
        nqi=[(0,0,0),(0,0,1),(0,+1,1),(-1,0,0)]
        
        #adding neighbors for the pathfinding
        for i in range(len(nopi)):
            ni=self.i+nopi[i][0]
            nj=self.j+nopi[i][1]
            nt=self.t+nopi[i][2]
            #print ni,nj, ni>=0 , ni<Lx , nj>=0 , nj<(Ly-1)
            if ni>=0 and ni<L and nj>=0 and nj<(L-1) and nt>=0 and nt<L:
                nn=indexX(ni,nj,nt,L)
                nop.append(nn)
                if i<4:#meaning, not moving in a temporal direction
                    ni=self.i+nqi[i][0]
                    nj=self.j+nqi[i][1]
                    nz=nqi[i][2]
                    nt=self.t
                    nn=indexq(ni,nj,nz,nt,L)
                    nq.append(nn)
        #adding neighbors for stabilizer measurement    
        for i in range(len(nopi)-2):#last 2 are temporal neighbors
            ni=self.i+nqi[i][0]
            nj=self.j+nqi[i][1]
            nz=nqi[i][2]
            nt=self.t
            if nz==0 and ni>-1 and nj>-1:
                if ni<L-1 and nj<L-1:
                    meas.append(indexq(ni,nj,nz,nt,L))
            if nz==1 and ni>-1 and nj>-1:
                if ni<L and nj<L:
                    meas.append(indexq(ni,nj,nz,nt,L))
        
        #neighbor qubits and neighbor stabilizers for bfs
        #print meas
        self.nop=nop
        self.nq=nq
        #neighbor qubits for stabilizer measurement
        self.meas=meas
class StabZ:#]
    def __init__(self,n,L):
        self.n=n
        self.L=L
        i,j,t=ijtZ(n,L)
        self.i=i
        self.j=j
        self.t=t
        self.name="Z"
        #neighbor qubits and neighbor stabilizers for bfs
        
        nop=[]
        nq=[]
        #neighbor qubits for stabilizer measurement
        meas=[]
        nopi=[(0,1,0),(-1,0,0),(+1,0,0),(0,-1,0),(0,0,1),(0,0,-1)]
        nqi=[(0,0,0),(0,0,1),(+1,0,1),(0,-1,0)]

        #adding neighbors for the pathfinding
        for i in range(len(nopi)):
            ni=self.i+nopi[i][0]
            nj=self.j+nopi[i][1]
            nt=self.t+nopi[i][2]
            #print ni,nj, ni>=0 , ni<Lx , nj>=0 , nj<(Ly-1)
            if ni>=0 and ni<L and nj>=0 and nj<(L-1) and nt>=0 and nt<L:
                nn=indexZ(ni,nj,nt,L)
                nop.append(nn)
                if i<4:#meaning, not moving in a temporal direction
                    ni=self.i+nqi[i][0]
                    nj=self.j+nqi[i][1]
                    nz=nqi[i][2]
                    nt=self.t
                    nn=indexq(ni,nj,nz,nt,L)
                    nq.append(nn)
        #adding neighbors for stabilizer measurement    
        for i in range(len(nopi)-2):#last 2 are temporal neighbors
            ni=self.i+nqi[i][0]
            nj=self.j+nqi[i][1]
            nz=nqi[i][2]
            nt=self.t
            if nz==0 and ni>-1 and nj>-1:
                if ni<L-1 and nj<L-1:
                    meas.append(indexq(ni,nj,nz,nt,L))
            if nz==1 and ni>-1 and nj>-1:
                if ni<L and nj<L:
                    meas.append(indexq(ni,nj,nz,nt,L))
        
        #neighbor qubits and neighbor stabilizers for bfs
        #print meas
        self.nop=nop
        self.nq=nq
        #neighbor qubits for stabilizer measurement
        self.meas=meas
        
        
'''
            simulation functions
'''


#applies noise to a state with probability p for each qbit
def noisestate(L,p,q):#]3
    
    nq=((L-1)*(L-1)+L*L)*L
    
    #err contains the index of the pauli matrix applied
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
            
    nop=L*(L-1)*L
    
    #1 if measurement error
    merX=[0]*nop
    merZ=[0]*nop
    for i in range(nop):
        r=np.random.random()
        if r<q:
            merX[i]=1
        r=np.random.random()
        if r<q:
            merZ[i]=1
    return err,merX,merZ



def plotnoisestateq(err,L,col='r',mize=8):#]3
    nq=len(err)
    for i in range(nq):
        x,y=xyq(i,L)
        
        if err[i]==1:
            plt.plot(x,y,marker='$X$',color=col,markersize=mize)
        if err[i]==2:
            plt.plot(x,y,marker='$Y$',color=col,markersize=mize)
        if err[i]==3:
            plt.plot(x,y,marker='$Z$',color=col,markersize=mize)
    return

def plotnoisestateX(err,L,col='r',colne='b',mize=8):#]
    nq=len(err)
    for i in range(nq):
        x,y=xyX(i,L)
        
        #if err[i]==0:
            #plt.plot(x,y,marker='$o$',color=colne,markersize=mize)
        if err[i]==1:
            plt.plot(x,y,marker='o',color=col,markersize=mize)
    return
def plotnoisestateZ(err,L,col='r',colne='b',mize=8):#]
    nq=len(err)
    for i in range(nq):
        x,y=xyZ(i,L)
        
        #if err[i]==0:
            #plt.plot(x,y,marker='$o$',color=colne,markersize=mize)
        if err[i]==1:
            plt.plot(x,y,marker='o',color=col,markersize=mize)
    return


#evaluates the syndrome for a certain noise state:
    
#evaluates the syndrome at t=0
def syndrome0(erc,Zc,Xc,L):#]3 
    #Ly=L
    nstab=L*(L-1)*L
    Z=[0]*nstab
    X=[0]*nstab
    
    #errc are the places where an error occur
    #in err, we compute the error state at each time
    #(an X in t0 and an Z in t1 will mean err=1 in t0 and err=2 in t1)  
    err=[0]*len(erc)
    nq=((L-1)*(L-1)+L*L)
    for i in range(len(erc)):
        #first errors are just erc
        if i<nq:
            err[i]=erc[i]
        else:
            err[i]=pprod(erc[i],err[i-nq])    
    

    
        
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

#generates the nodes of the graphs for the MinWeightPerfectMatching
def syndromet(meX,meZ,X,Z,L):#]3
    #mX: measures of the stabilizers at each timestep
    #meX: errors in the measures
    #sX: final nodes in the graph
    mX=[0]*len(meX)
    mZ=[0]*len(meZ)
    sX=[0]*len(meX)
    sZ=[0]*len(meZ)
    
    nst=L*(L-1)
    for n in range(len(meX)):
        mX[n]=(X[n]+meX[n])%2
        mZ[n]=(Z[n]+meZ[n])%2    
        
        #now we compute the nodes
        #if the measures are different from the previous step, 
        #then we have a node (sX=1)
        if n<nst:
            sX[n]=mX[n]
            sZ[n]=mZ[n]
        else:
            if(mX[n]==mX[n-nst]):
                sX[n]=0
            else:
                sX[n]=1
            if(mZ[n]==mZ[n-nst]):
                sZ[n]=0
            else:
                sZ[n]=1
            
    return sX,sZ
            


def plotsyndromet(X,Z,L,only='x',colx='y',colz='g',mize=8):#]
    for k in range(len(X)):
        if only=='x':
            x,y=xyX(k,L)
            if X[k]>0:
                plt.plot(x,y,marker='d',color=colx,markersize=mize)
        else:
                
            x,y=xyZ(k,L)
            if Z[k]>0:
                plt.plot(x,y,marker='d',color=colz,markersize=mize)
    return


def distToBorderX(n,L):#]
    #X stabilizers are linked to rough edges (top and bottom)
    global wp,wq
    i,j,t=ijtX(n,L)
    return int(min(min(j+1,L-1-j)*wp,(L-t)*wq)*10000)

def distToBorderZ(n,L):#]
    global wp,wq
    #X stabilizers are linked to smooth edges (laterals)
    i,j,t=ijtZ(n,L)
    return int(min(min(i+1,L-i-1)*wp,(L-t)*wq)*10000)


    
#evaluates min weight perfect matching for the syndrome

def blossomopX(X,L):#]
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
                
        name=str(np.random.randint(0,10**9))
        gfile="/dev/shm/"+name+"g.txt"
        sfile="/dev/shm/"+name+"s.txt"
        savegraph(nx*2,edges,gfile)
        blossom(gfile,sfile)
        match,d=readsol(sfile)
        os.remove(gfile)
        os.remove(sfile)
            
        sol=[]
        for edg in match:
            if edg[0]<ner or edg[1]<ner:#only record interesting edges
                e1=xl[edg[0]%ner]#+L*L*(edg[0]/ner)
                e2=xl[edg[1]%ner]#+L*L*(edg[1]/ner)
                sol.append((e1,e2))
                
    return sol
          
def blossomopZ(X,L):#]
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
                dist=mdstZ(xl[e1],xl[e2],L)
                edges.append((e1,e2,dist))
                edges.append((e1+ner,e2+ner,0))
                
        name=str(np.random.randint(0,10**9))
        gfile="/dev/shm/"+name+"g.txt"
        sfile="/dev/shm/"+name+"s.txt"
        savegraph(nx*2,edges,gfile)
        blossom(gfile,sfile)
        match,d=readsol(sfile)
        os.remove(gfile)
        os.remove(sfile)
            
        sol=[]
        for edg in match:
            if edg[0]<ner or edg[1]<ner:
                e1=xl[edg[0]%ner]#+L*L*(edg[0]/ner)
                e2=xl[edg[1]%ner]#+L*L*(edg[1]/ner)
                sol.append((e1,e2))
                
    return sol
              

def bfsx(a,b,S,L):#]
    #if a==b, means that the match is to the border
    if a==b:
        name=S[a].name
        path=[]           
        i=S[a].i
        j=S[a].j
        t=S[a].t
        if name=="Z":
            #trace the shortest horizontal path to the border 
            
            #in this case the shortest is to the top
            if (L-t)<min((i+1),(L-i-1)):
                return path
            
            if (i+1)<(L-i-1):
                #path from left to 'a'                
                for n in range(i+1):
                    path.append(indexq(n,j,1,t,L))
            else:
                #path from right to 'a'                
                for n in range(L-i-1):
                    path.append(indexq(L-1-n,j,1,t,L))
 
        if name=="X":
            #trace the shortest vertical path to the border
            #in this case the shortest is to the top
            if (L-t)<min((j+1),(L-j-1)):
                return path
            if (j+1)<(L-j-1):
                #path bottom to 'a'
                for n in range(j+1):
                    path.append(indexq(i,n,1,t,L))
            else:
                #path from top to 'a'
                for n in range(L-j-1):
                    path.append(indexq(i,L-1-n,1,t,L))
                    
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
        nqn=len(qneigs)
        for n in range(len(neigs)):
            if neigs[n] == b:
                if n<nqn:
                    path+=[qneigs[n]]
                return path
            if neigs[n] not in visited:
                front.append(neigs[n])
                
                if n<nqn:
                    paths.append(path+[qneigs[n]])
                else:
                    paths.append(path)
    return []     
def plotpathx(a,b,path,L,col='y',marker='$Z$', mize=7):#]
    for n in path:
        #i,j,z,t=ijzt(n,L)
        x,y=xyq(n,L)
        #xa,ya=xy(a,L)
        #xb,yb=xy(b,L)
        #plt.plot([xa,xb],[ya+.5,yb+.5],'r-')
        plt.plot(x,y,color=col,marker=marker,markersize=mize)
    
    
    
    
def plotpathz(a,b,path,L,col='g',marker='$X$', mize=4):#]
    for n in path:
        #i,j,z=xyz(n,L)
        #x=i+.5*(1-z)
        #y=j+.5*(1-z)
        x,y=xyq(n,L)
        #xa,ya=xy(a,L-1)
        #xb,yb=xy(b,L-1)
        #plt.plot([xa+.5,xb+.5],[ya,yb],'b-')
        plt.plot(x,y,color=col,marker=marker,markersize=mize)
 
    
    
def plotmatch(solX,solZ,L):#]
    nstab=L*(L-1)*L
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
    
def plotmatchX(solX,L):#]
    nstab=L*(L-1)*L
    Sx=[]
    #list of stabilizers to be measured
    for i in range(nstab):
        Sx.append(StabX(i,L))
    for pair in solX:
        a=pair[0]
        b=pair[1]
        xa,ya=xyX(a,L)
        xb,yb=xyX(b,L)
        plt.plot([xa,xb],[ya,yb],'rp-')
        path=bfsx(a,b,Sx,L)
        plotpathx(a,b,path,L)

    
def plotmatchZ(solZ,L):#]
    nstab=L*(L-1)*L
    Sz=[]
    #list of stabilizers to be measured
    for i in range(nstab):
        Sz.append(StabZ(i,L))
    for pair in solZ:
        a=pair[0]
        b=pair[1]
        xa,ya=xyZ(a,L)
        xb,yb=xyZ(b,L)
        plt.plot([xa,xb],[ya,yb],'rp-')
        path=bfsx(a,b,Sz,L)
        plotpathz(a,b,path,L)      
        
        
'''
            final result error+corrections
'''
        
def corrx(solX,S,L):#]
    corx=[0]*((L-1)*(L-1)+L*L)
    for pair in solX:
        a=pair[0]
        b=pair[1]
        path=bfsx(a,b,S,L)
        for q in path:
            #i,j,z,t=ijzt(q,L)
            #n=indexq(i,j,z,0,L)
            n=q%len(corx)
            corx[n]=pprod(corx[n],3)
    return corx
def corrz(solX,S,L):#]
    corx=[0]*((L-1)*(L-1)+L*L)
    for pair in solX:
        a=pair[0]
        b=pair[1]
        path=bfsx(a,b,S,L)
        for q in path:
            #i,j,z,t=ijzt(q,L)
            #n=indexq(i,j,z,0,L)
            n=q%len(corx)
            corx[n]=pprod(corx[n],1)
    return corx
#product of pauli matrices
def pprod(a,b):#]
    if a*b==0:
        return a+b
    if a==b:
        return 0
    prod=[3,2,1]
    return prod[(a+b)%3]
    
def ercor(ers,corx,corz):#]3
    nq=len(corx)
    err=[0]*nq
    sol=[0]*nq
    
    for i in range(len(ers)):
        err[i%nq]=pprod(err[i%nq],ers[i])
    
    
    for i in range(nq):
        sol[i]=pprod(pprod(err[i],corx[i]),corz[i])
    return sol
        
        
'''
            checking logical errors
'''        
def conmz(a):#]
    if a==1 or a==2:
        return 1
    else:
        return 0
def conmx(a):#]
    if a>1:
        return 1
    else:
        return 0
def checkl(erc,L):#]
    vx=0
    vz=0
    hx=0
    hz=0    
    for i in range(L):
        hx+=conmx(erc[indexq(i,0,0,L)])
        vz+=conmz(erc[indexq(0,i,0,L)])
    #for i in range(L-1):
        #hx+=conmz(erc[indexyz(i,0,1,L)])
        #vx+=conmx(erc[indexyz(0,i,0,L)])
    hx=hx%2
    #hz=hz%2
    #vx=vx%2
    vz=vz%2
    return max(hx,vz),hz,vx


'''
            plotting functions
'''


#plots the entire lattice
def latplot(L,col='b'):#]
    #square lattice plot
    
    for n in range(L):
        plt.plot([n,n],[-.5,L-.5],col)
    for n in range(L-1):
        plt.plot([0,L-1],[n+.5,n+.5],col)
    return
def plotnoisestate2D(err,L,col='r',mize=8):#]!
    nq=(L-1)*(L-1)+L*L
    for i in range(nq):
        x,y,z,t=ijzt(i,L)
        
        if err[i]==1:
            plt.plot(x+.5*(1-z),y+0.5*(1-z),marker='$X$',color=col,markersize=mize)
        if err[i]==2:
            plt.plot(x+.5*(1-z),y+0.5*(1-z),marker='$Y$',color=col,markersize=mize)
        if err[i]==3:
            plt.plot(x+.5*(1-z),y+0.5*(1-z),marker='$Z$',color=col,markersize=mize)
    return
    

def latplotX(L,col='b',mize=7):###

    #square lattice plot
    phi=60*np.pi/180.
    vi=np.array((1,0))
    vj=np.array((np.sin(phi),np.cos(phi)))*.75
    #vt=np.array((0,1))
    
    #vertical qbits
    for n in range(L):
        plt.plot([vi[0]*n,(vi[0]*n+vj[0]*(L))],[0.,L*vj[1]],col)
    #horizontal qbits
    for n in range(1,L):
        plt.plot([vj[0]*n,(vi[0]*(L-1)+vj[0]*n)],[vj[1]*(n),n*vj[1]],col)
    
    for n in range(L*(L-1)):
        x1,y1=xyX(n,L)
        i,j,t=ijtX(n,L)
        nn=indexX(i,j,(L-1),L)
        x2,y2=xyX(nn,L)
        plt.plot([x1,x2],[y1,y2],'b-.')   
    nq=L*L*(L-1)
    for i in range(nq):
        x,y=xyX(i,L)
        plt.plot(x,y,marker='$o$',color=col,markersize=mize) 
    
    return

def latplotZ(L,col='b',mize=7):###

    #square lattice plot
    phi=60*np.pi/180.
    vi=np.array((1,0))
    vj=np.array((np.sin(phi),np.cos(phi)))*.75
    #vt=np.array((0,1))
    
    #vertical qbits
    for n in range(L):
        plt.plot([vi[0]*n,(vi[0]*n+vj[0]*(L))],[0.,L*vj[1]],col)
    #horizontal qbits
    for n in range(1,L):
        plt.plot([vj[0]*n,(vi[0]*(L-1)+vj[0]*n)],[vj[1]*(n),n*vj[1]],col)
    
    for n in range(L*(L-1)):
        x1,y1=xyZ(n,L)
        i,j,t=ijtZ(n,L)
        nn=indexZ(i,j,(L-1),L)
        x2,y2=xyZ(nn,L)
        plt.plot([x1,x2],[y1,y2],'b-.')  
    nq=L*L*(L-1)
    for i in range(nq):
        x,y=xyZ(i,L)
        plt.plot(x,y,marker='$o$',color=col,markersize=mize)   
    
    return

def experimentplot(L,p,q):###
    nstab=L*(L-1)*L
    Sx=[]
    Sz=[]
    global wp
    wp=np.log((1-p)/p)
    global wq
    wq=np.log((1-q)/q)
    #list of stabilizers to be measured
    for i in range(nstab):
        Sz.append(StabZ(i,L))
        Sx.append(StabX(i,L))
    
    
    plt.figure(0,figsize=(2,3))
    plt.title("Surface code")
    plt.clf()
    plt.subplot(2,3,1)
    plt.title('Errors, X Stabilizers')
    latplotX(L)
    
    err,meX,meZ=noisestate(L,p,q)
    plotnoisestateq(err,L)
    plotnoisestateX(meX,L)
    plt.subplot(2,3,4)
    plt.title('Errors, Z Stabilizers')
    
    
    latplotZ(L)
    plotnoisestateq(err,L)
    plotnoisestateZ(meZ,L)
    
    
    plt.subplot(2,3,2)
    plt.title('Syndrome measurements X')
    X,Z=syndrome0(err,Sz,Sx,L)
    synX,synZ=syndromet(meX,meZ,X,Z,L)
    latplotX(L)
    plotnoisestateq(err,L)
    plotsyndromet(synX,synZ,L)


    
    plt.subplot(2,3,5)
    plt.title('Syndrome measurements Z')
    latplotZ(L)
    plotnoisestateq(err,L)
    plotsyndromet(synX,synZ,L,'z')
    
    
    plt.subplot(2,3,3)
    plt.title("Matching X")
    latplotX(L)
    solX=blossomopX(synX,L)
    solZ=blossomopZ(synZ,L)
    #plotsyndrome(X,Z,L)
    plotmatchX(solX,L)
    #plotmatchZ(solZ,L)
    
    plt.subplot(2,3,6)
    plt.title("Matching Z")
    latplotZ(L)
    #plotsyndrome(X,Z,L)
    plotmatchZ(solZ,L)
    #plotmatchZ(solZ,L)
    
    
    #plt.figure(1)
    #plt.clf()
    corx=corrx(solX,Sx,L)
    corz=corrz(solZ,Sz,L)
    erc=ercor(err,corx,corz)
    #latplot(L)
    #plotnoisestate2D(erc,L)
    return scd.ploterr(erc,L)
  
def experiment(L,p,q):###
    nstab=L*(L-1)*L
    Sx=[]
    Sz=[]
    global wp
    wp=np.log((1-p)/p)
    global wq
    wq=np.log((1-q)/q)
    #list of stabilizers to be measured
    for i in range(nstab):
        Sz.append(StabZ(i,L))
        Sx.append(StabX(i,L))
    
    
    err,meX,meZ=noisestate(L,p,q)
    
    
    X,Z=syndrome0(err,Sz,Sx,L)
    synX,synZ=syndromet(meX,meZ,X,Z,L)

    
    solX=blossomopX(synX,L)
    solZ=blossomopZ(synZ,L)
    
    corx=corrx(solX,Sx,L)
    corz=corrz(solZ,Sz,L)
    erc=ercor(err,corx,corz)
    
    return scdR.nploterr(erc,L)

def mcexp(L,p,q,Nit):
    r=0.0
    for i in range(Nit):
        r+=experiment(L,p,q)
    return r/Nit
#experimentplot(4,.04,.04)


'''
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
            return 1.*c/nqs
    return 1.*c/nqs
    
    
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
    
    
    plt.figure(0,figsize=(2,2))
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

'''

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